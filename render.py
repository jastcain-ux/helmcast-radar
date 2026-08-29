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
from PIL import Image
from pyproj import CRS, Transformer

NOMADS = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod"
MERCATOR_R = 20037508.342789244

# The frames the app scrubs through. HRRR runs out to 18 hours; Helmcast shows
# six, matching the radar timeline it already has.
FORECAST_HOURS = range(1, 7)

# CONUS, a little wider than the model grid so nothing is clipped at the edges.
CONUS_BBOX = (-127.0, 21.0, -65.0, 50.0)
FRAME_SIZE = (3200, 1800)

# NWS-style reflectivity ramp, in dBZ. The first stop is where painting starts:
# below it the frame stays transparent.
RAMP = [
    ( 5, (  4, 233, 231)), (10, (  1, 159, 244)), (15, (  3,   0, 244)),
    (20, (  2, 253,   2)), (25, (  1, 197,   1)), (30, (  0, 142,   0)),
    (35, (253, 248,   2)), (40, (229, 188,   0)), (45, (253, 149,   0)),
    (50, (253,   0,   0)), (55, (212,   0,   0)), (60, (188,   0,   0)),
    (65, (248,   0, 253)), (70, (152,  84, 198)),
]
ALPHA = 200


def _url(day, run, fh, suffix=""):
    return f"{NOMADS}/hrrr.{day}/conus/hrrr.t{run:02d}z.wrfsfcf{fh:02d}.grib2{suffix}"


def _get(url, byte_range=None, timeout=90):
    req = urllib.request.Request(url)
    if byte_range:
        req.add_header("Range", f"bytes={byte_range[0]}-{byte_range[1]}")
    return urllib.request.urlopen(req, timeout=timeout).read()


def refc_byte_range(day, run, fh):
    """Byte range of the composite-reflectivity record.

    The .idx sidecar lists every record as `number:offset:d=date:field:...`,
    so one small text fetch turns a 130 MB file into a 0.3 MB range request.
    """
    lines = _get(_url(day, run, fh, ".idx"), timeout=30).decode().splitlines()
    for i, line in enumerate(lines):
        parts = line.split(":")
        if len(parts) > 3 and parts[3] == "REFC":
            start = int(parts[1])
            end = int(lines[i + 1].split(":")[1]) - 1 if i + 1 < len(lines) else ""
            return start, end
    raise LookupError(f"no REFC record in {_url(day, run, fh, '.idx')}")


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
                urllib.request.Request(_url(day, run, max(FORECAST_HOURS), ".idx")),
                timeout=20).read(1)
            return day, run
        except Exception:
            continue
    raise SystemExit("no complete HRRR run published in the last 10 hours")


def read_refc(day, run, fh):
    lo, hi = refc_byte_range(day, run, fh)
    raw = _get(_url(day, run, fh), (lo, hi))
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
    h, w = dbz.shape
    out = np.zeros((h, w, 4), dtype=np.uint8)
    for threshold, rgb in RAMP:
        m = dbz >= threshold
        out[m, 0], out[m, 1], out[m, 2] = rgb
        out[m, 3] = ALPHA
    return out


def render(values, meta, bbox, size):
    """Nearest-neighbour sample from the Lambert grid into web mercator.

    Nearest neighbour rather than interpolation on purpose: reflectivity is a
    measurement of intensity, and smoothing it invents gradients between cells
    that the model did not produce.
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

    ix = np.rint((lx - x0) / meta["DxInMetres"]).astype(int)
    iy = np.rint((ly - y0) / meta["DyInMetres"]).astype(int)
    inside = (ix >= 0) & (ix < meta["Nx"]) & (iy >= 0) & (iy < meta["Ny"])

    sampled = np.full((height, width), -99.0)
    sampled[inside] = values[iy[inside], ix[inside]]
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
    for fh in FORECAST_HOURS:
        try:
            values, meta, nbytes = read_refc(day, run, fh)
            image = render(values, meta, CONUS_BBOX, FRAME_SIZE)
            name = f"refc-f{fh:02d}.png"
            image.save(os.path.join(args.out, name), optimize=True)
            downloaded += nbytes
            when = valid_time(meta)
            frames.append({
                "forecastHour": fh,
                "validTime": when.isoformat().replace("+00:00", "Z") if when else None,
                "image": name,
            })
            size_kb = os.path.getsize(os.path.join(args.out, name)) / 1024
            print(f"  f{fh:02d}  {nbytes/1e6:.2f} MB in  ->  {size_kb:5.0f} KB  {name}")
        except Exception as e:
            # Absent, never substituted. A missing hour is a thing the app has
            # to be able to say out loud.
            print(f"  f{fh:02d}  SKIPPED: {e}", file=sys.stderr)

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
    print(f"{len(frames)}/{len(list(FORECAST_HOURS))} frames, "
          f"{downloaded/1e6:.1f} MB pulled from NOMADS")


if __name__ == "__main__":
    main()
