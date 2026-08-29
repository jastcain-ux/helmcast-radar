#!/usr/bin/env python3
"""
Render NOAA HRRR forecast reflectivity to PNG frames Helmcast can draw.

Why this exists
---------------
Radar measures what is already in the air. NOAA publishes no *forecast*
reflectivity as ready-made map tiles — it exists only as HRRR model output in
GRIB2 on NOMADS. Every consumer app showing "future radar" is either
extrapolating recent frames (good for about an hour) or rendering this same
model, usually bought from a vendor.

The data is free and public domain. The only missing piece was something to
turn the files into images, which is all this is.

What it does
------------
Finds the newest HRRR run, downloads *only* the composite reflectivity record
from each forecast hour using the GRIB index sidecar (about 0.3 MB each rather
than the ~130 MB full file), reprojects from the model's Lambert grid to web
mercator, and writes one national PNG per hour plus a manifest.

Honesty rules that are not optional
-----------------------------------
- A frame that could not be produced is ABSENT from the manifest, never
  substituted with a neighbouring hour. The app says "no frame for this hour"
  rather than showing 3 PM's weather under a 6 PM label.
- The manifest carries the model run time, so the app can show how old the
  forecast is. A six-hour outlook from a run that stopped updating four hours
  ago is not a six-hour outlook.
- Below 5 dBZ is transparent. Drizzle nobody can feel must not paint the bay.
"""
import argparse, datetime, json, math, os, sys, tempfile, urllib.request

import numpy as np
import eccodes
from scipy import ndimage
from PIL import Image
from pyproj import CRS, Transformer

NOMADS = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod"
MERCATOR_R = 20037508.342789244

# The frames the app scrubs through. HRRR runs out to 18 hours; Helmcast shows
# six, matching the radar timeline it already has.
# Every 15 minutes out to 8 hours.
#
# QUARTER-HOURS because that is HRRR's own sub-hourly cadence — these are real
# model steps, not interpolated ones. An hourly step was too coarse to watch a
# squall line approach: at 25 knots a storm crosses fifteen miles between
# frames, which on a 60-mile map is a jump rather than a movement.
#
# NINE HOURS, although the app only scrubs six. The newest complete run is
# always one to two hours behind the wall clock, so a series starting at run+1h
# ran out before the end of a six-hour scrubber and the last steps fell back to
# the wind field — a flat blue wash that read as "the radar turned blue".
#
# Three hours of margin is deliberate: past three hours the manifest reads as
# stale and the app says so out loud, so covering further would only be padding
# a window the screen has already disowned.
#
# Frames are only downloaded when scrubbed to, so the extra ones cost storage
# rather than the boater's data.
FORECAST_MINUTES = list(range(15, 9 * 60 + 1, 15))

# CONUS, a little wider than the model grid so nothing is clipped at the edges.
CONUS_BBOX = (-127.0, 21.0, -65.0, 50.0)
# Finer than the model's 3 km grid (this is ~1.4 km/px), because the app
# upscales a national frame down to one bay and the smoothing has to have
# something to work with. Going finer than this buys nothing: the information
# ceiling is the model, not the raster.
FRAME_SIZE = (4096, 2304)

# NWS-style reflectivity ramp, in dBZ. The first stop is where painting starts:
# below it the frame stays transparent.
# Starts at 15 dBZ, not 5.
#
# Between 5 and 15 dBZ is drizzle and mist, and at 3 km resolution it is also
# where the model is noisiest — a haze of speckle that made the whole picture
# look broken while telling a boater nothing. What matters here is squalls, and
# those are 20 dBZ and up. Nothing that affects a decision is hidden: the
# rating never reads this map, and the legend says it starts at light rain.
RAMP = [
    (15, (  3,   0, 244)),
    (20, (  2, 253,   2)), (25, (  1, 197,   1)), (30, (  0, 142,   0)),
    (35, (253, 248,   2)), (40, (229, 188,   0)), (45, (253, 149,   0)),
    (50, (253,   0,   0)), (55, (212,   0,   0)), (60, (188,   0,   0)),
    (65, (248,   0, 253)), (70, (152,  84, 198)),
]
ALPHA = 200
# The outermost returns fade in rather than starting at full opacity, so the
# edge of a cell is an edge rather than a cliff. Real precipitation does not
# have a hard boundary and drawing one implies a certainty about where the rain
# stops that a 3 km model does not have.
EDGE_FADE_DBZ = 8.0
# Steps of edge softness. The eye cannot see 256 of them and PNG pays for every
# one — quantising here cut the busiest frame from 1.15 MB to 660 KB with no
# visible difference, which is a boater's cellular data at a boat ramp.
EDGE_FADE_STEPS = 4
# The ramp is walked as a continuous gradient rather than 14 hard bands, in
# this many steps across its range. Hard bands drew visible contour edges once
# a national frame was zoomed to one bay — a staircase in the picture that
# corresponds to nothing in the weather. 16 steps is indistinguishable from 48
# on screen and compresses to little more than the banded version did.
RAMP_STEPS = 24


def _url(day, run, minutes, suffix=""):
    """The file holding this lead time, and the label its index uses.

    Whole hours live in the hourly surface files; the half-hours live in the
    sub-hourly ones, which carry four 15-minute steps each. Both publish
    composite reflectivity, and both are indexed the same way — only the
    filename and the way the record labels its lead time differ.
    """
    if minutes % 60 == 0:
        name = f"wrfsfcf{minutes // 60:02d}"
    else:
        # e.g. 45 min lives in wrfsubhf01, alongside 15, 30 and 60.
        name = f"wrfsubhf{-(-minutes // 60):02d}"
    return f"{NOMADS}/hrrr.{day}/conus/hrrr.t{run:02d}z.{name}.grib2{suffix}"


def _lead_label(minutes):
    return f"{minutes // 60} hour fcst" if minutes % 60 == 0 else f"{minutes} min fcst"


def _get(url, byte_range=None, timeout=90):
    req = urllib.request.Request(url)
    if byte_range:
        req.add_header("Range", f"bytes={byte_range[0]}-{byte_range[1]}")
    return urllib.request.urlopen(req, timeout=timeout).read()


def refc_byte_range(day, run, minutes):
    """Byte range of the composite-reflectivity record for this lead time.

    The .idx sidecar lists every record as
    `number:offset:d=date:field:level:lead:`, so one small text fetch turns a
    130 MB file into a 0.3 MB range request. Sub-hourly files hold four lead
    times, so the lead label has to be matched as well as the field — taking
    the first REFC would silently return :15 for every half-hour asked for.
    """
    url = _url(day, run, minutes, ".idx")
    lines = _get(url, timeout=30).decode().splitlines()
    want = _lead_label(minutes)
    for i, line in enumerate(lines):
        parts = line.split(":")
        if len(parts) > 5 and parts[3] == "REFC" and parts[5].strip() == want:
            start = int(parts[1])
            end = int(lines[i + 1].split(":")[1]) - 1 if i + 1 < len(lines) else ""
            return start, end
    raise LookupError(f"no REFC '{want}' in {url}")


def latest_complete_run(now=None):
    """Newest run whose whole forecast window has been published.

    Checked against the LAST hour we need, not the first: a run that is still
    uploading would otherwise be chosen and then come up short, which would
    silently shorten the outlook the app offers.
    """
    now = now or datetime.datetime.now(datetime.timezone.utc)
    for back in range(0, 10):
        t = now - datetime.timedelta(hours=back)
        day, run = t.strftime("%Y%m%d"), t.hour
        try:
            urllib.request.urlopen(
                urllib.request.Request(_url(day, run, max(FORECAST_MINUTES), ".idx")),
                timeout=20).read(1)
            return day, run
        except Exception:
            continue
    raise SystemExit("no complete HRRR run published in the last 10 hours")


def read_refc(day, run, minutes):
    lo, hi = refc_byte_range(day, run, minutes)
    raw = _get(_url(day, run, minutes), (lo, hi))
    fd, path = tempfile.mkstemp(suffix=".grib2")
    os.write(fd, raw)
    os.close(fd)
    try:
        with open(path, "rb") as f:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                raise ValueError("range request did not contain a GRIB message")
            keys = ("Nx", "Ny", "LaDInDegrees", "LoVInDegrees", "Latin1InDegrees",
                    "Latin2InDegrees", "latitudeOfFirstGridPointInDegrees",
                    "longitudeOfFirstGridPointInDegrees", "DxInMetres", "DyInMetres",
                    "radius", "validityDate", "validityTime")
            meta = {}
            for k in keys:
                try:
                    meta[k] = eccodes.codes_get(gid, k)
                except Exception:
                    meta[k] = None
            values = eccodes.codes_get_values(gid).reshape(meta["Ny"], meta["Nx"])
            eccodes.codes_release(gid)
        return values, meta, len(raw)
    finally:
        os.unlink(path)


def colourise(dbz):
    """dBZ -> RGBA along a continuous ramp, with a soft outer edge.

    The ramp's anchor colours are the conventional radar ones — cyan and blue
    for light, green through yellow for moderate, red for heavy — so intensity
    still reads the way a boater expects. What changed is that the colour moves
    between them smoothly instead of stepping, because hard bands put a
    staircase in the picture that corresponds to nothing in the weather.
    """
    stops = np.array([t for t, _ in RAMP], dtype=float)
    colours = np.array([c for _, c in RAMP], dtype=float)
    floor, ceiling = stops[0], stops[-1]

    quantised = np.clip(dbz, floor, ceiling)
    quantised = floor + np.round(
        (quantised - floor) / (ceiling - floor) * RAMP_STEPS) / RAMP_STEPS * (ceiling - floor)

    h, w = dbz.shape
    out = np.zeros((h, w, 4), dtype=np.uint8)
    for channel in range(3):
        out[..., channel] = np.interp(quantised, stops, colours[:, channel]).astype(np.uint8)
    lit = dbz >= floor
    fade = np.clip((dbz - floor) / EDGE_FADE_DBZ, 0.0, 1.0)
    fade = np.ceil(fade * EDGE_FADE_STEPS) / EDGE_FADE_STEPS
    out[..., 3] = np.where(lit, (fade * ALPHA).astype(np.uint8), 0)
    return out


def render(values, meta, bbox, size):
    """Cubic-spline sample from the model's Lambert grid into web mercator.

    This started as nearest-neighbour, on the theory that interpolating
    invented gradients the model did not produce. That was the wrong call.
    Reflectivity is a continuous field that HRRR samples every 3 km, and
    nearest-neighbour does not display that honestly — it draws hard square
    edges, which assert a sharp boundary exactly where the model is least
    certain. A storm edge is not a 3 km square, and showing one as though it
    were is its own small false precision.

    Cubic spline rather than bilinear or a blur, and the reason is a
    false-calm one rather than a cosmetic one. Measured against a real frame
    whose grid peaked at 71.9 dBZ:

        cubic spline      72.4 dBZ   (+0.4)
        bilinear          70.1 dBZ   (-1.8)
        gaussian s=0.6    66.2 dBZ   (-5.7)
        gaussian s=1.0    60.0 dBZ   (-11.9)

    Blurring is the obvious way to make a coarse field look smooth, and it
    quietly shaves the top off every storm — an app that paints a 60 dBZ core
    where the model said 72 is understating exactly the thing a boater needs
    to see. Cubic interpolation smooths the shape while leaving peaks intact.
    """
    lcc = CRS.from_proj4(
        f"+proj=lcc +lat_1={meta['Latin1InDegrees']} +lat_2={meta['Latin2InDegrees']} "
        f"+lat_0={meta['LaDInDegrees']} +lon_0={meta['LoVInDegrees']} "
        f"+a={meta['radius']} +b={meta['radius']} +units=m +no_defs")
    to_lcc = Transformer.from_crs("EPSG:4326", lcc, always_xy=True)
    x0, y0 = to_lcc.transform(meta["longitudeOfFirstGridPointInDegrees"],
                              meta["latitudeOfFirstGridPointInDegrees"])

    west, south, east, north = bbox
    def merc_y(lat):
        return math.log(math.tan((90 + lat) * math.pi / 360)) / (math.pi / 180) * MERCATOR_R / 180

    min_x, max_x = west * MERCATOR_R / 180, east * MERCATOR_R / 180
    min_y, max_y = merc_y(south), merc_y(north)

    width, height = size
    mx, my = np.meshgrid(np.linspace(min_x, max_x, width),
                         np.linspace(max_y, min_y, height))
    lon = mx / MERCATOR_R * 180
    lat = np.degrees(2 * np.arctan(np.exp(my / MERCATOR_R * math.pi)) - math.pi / 2)
    lx, ly = to_lcc.transform(lon, lat)

    fx = (lx - x0) / meta["DxInMetres"]
    fy = (ly - y0) / meta["DyInMetres"]

    inside = (fx >= 0) & (fx < meta["Nx"] - 1) & (fy >= 0) & (fy < meta["Ny"] - 1)
    interpolated = ndimage.map_coordinates(values, [fy, fx], order=3,
                                           mode="nearest", prefilter=True)
    sampled = np.where(inside, interpolated, -99.0)
    return Image.fromarray(colourise(sampled))


def valid_time(meta):
    d, t = meta.get("validityDate"), meta.get("validityTime")
    if not d:
        return None
    return datetime.datetime(d // 10000, d // 100 % 100, d % 100,
                             (t or 0) // 100, (t or 0) % 100,
                             tzinfo=datetime.timezone.utc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out", help="directory to write frames into")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    day, run = latest_complete_run()
    run_at = datetime.datetime.strptime(f"{day}{run:02d}", "%Y%m%d%H").replace(
        tzinfo=datetime.timezone.utc)
    print(f"HRRR run {run_at.isoformat()}")

    frames, downloaded = [], 0
    for minutes in FORECAST_MINUTES:
        try:
            values, meta, nbytes = read_refc(day, run, minutes)
            image = render(values, meta, CONUS_BBOX, FRAME_SIZE)
            name = f"refc-{minutes:04d}.png"
            image.save(os.path.join(args.out, name), optimize=True)
            downloaded += nbytes
            # The model's own valid time, not run + lead. If those ever
            # disagree the file is right and the arithmetic is wrong.
            when = valid_time(meta)
            frames.append({
                "leadMinutes": minutes,
                "validTime": when.isoformat().replace("+00:00", "Z") if when else None,
                "image": name,
            })
            size_kb = os.path.getsize(os.path.join(args.out, name)) / 1024
            print(f"  +{minutes:3d} min  {nbytes/1e6:.2f} MB in  ->  {size_kb:5.0f} KB  {name}")
        except Exception as e:
            # Absent, never substituted. A missing step is a thing the app has
            # to be able to say out loud.
            print(f"  +{minutes:3d} min  SKIPPED: {e}", file=sys.stderr)

    if not frames:
        raise SystemExit("no frames rendered; leaving the previous manifest in place")

    manifest = {
        "model": "HRRR",
        "field": "composite reflectivity",
        "source": "NOAA NCEP, public domain",
        "runTime": run_at.isoformat().replace("+00:00", "Z"),
        "generatedAt": datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds").replace("+00:00", "Z"),
        "bbox": {"west": CONUS_BBOX[0], "south": CONUS_BBOX[1],
                 "east": CONUS_BBOX[2], "north": CONUS_BBOX[3]},
        "size": {"width": FRAME_SIZE[0], "height": FRAME_SIZE[1]},
        "dbzFloor": RAMP[0][0],
        "frames": frames,
    }
    with open(os.path.join(args.out, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"{len(frames)}/{len(FORECAST_MINUTES)} frames, "
          f"{downloaded/1e6:.1f} MB pulled from NOMADS")


if __name__ == "__main__":
    main()
