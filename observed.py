#!/usr/bin/env python3
"""
Render NOAA MRMS *measured* reflectivity to PNG frames SeaWise can draw.

Why this exists
---------------
The forecast half of SeaWise's radar timeline has always been rendered here,
and it looks like weather. The measured half was NOAA's nowCOAST WMS — someone
else's pre-drawn tiles, painted nearest-neighbour with hard colour bands, which
at the zoom a boater actually uses reads as a grid of coloured squares.

That is not a resolution problem and buying a sharper source does not fix it.
It is a *rendering* problem, and this file already contains the fix: the same
cubic sampling, continuous ramp and soft edge that made the forecast frames
stop looking broken.

So the measured half is rendered here too, from the same kind of input.

What it does
------------
Takes NOAA MRMS MergedReflectivityQCComposite — the national radar mosaic, one
frame every two minutes, public domain on NOAA's open bucket — and writes one
PNG per timeline step plus a manifest, using `render.colourise` unchanged.

Why that matters beyond looks
-----------------------------
The two halves of the timeline now come out of one renderer with one palette
and one floor. They used to be matched by hand, and when they drifted apart
rain appeared to *thin out* exactly at "now" as the scrubber crossed over —
false calm at the instant of most interest. That class of bug is now
impossible rather than merely fixed.

Honesty rules carried over unchanged
------------------------------------
- A frame that could not be produced is ABSENT from the manifest, never
  substituted with a neighbouring time.
- Each frame carries the time it was *observed*, not when it was fetched.
- Below 5 dBZ is transparent, the same floor the forecast frames use.
"""
import argparse, multiprocessing, datetime, gzip, io, json, math, os, re, sys, urllib.request

import numpy as np
import eccodes
from scipy import ndimage
from PIL import Image

from render import ALPHA, CONUS_BBOX, MERCATOR_R, colourise, RAMP
from smoothing import smooth

BUCKET = "https://noaa-mrms-pds.s3.amazonaws.com"
PRODUCT = "MergedReflectivityQCComposite_00.50"

# How far back the scrubber can reach, and how finely.
#
# Ten minutes matches what the app already thinned NOAA's four-minute mosaic
# to: closer spacing than that is more frames than anyone drags through, and
# the difference between two frames ten minutes apart is a movement you can
# see. Two hours is the window a boater actually uses — watching a line
# approach — rather than the seven hours nowCOAST offered and nobody scrubbed.
STEP_MINUTES = 10
HISTORY_MINUTES = 120

# Regional cells, not one national frame.
#
# This was a national 5120 px frame first, and it was the wrong shape of
# answer. A boater looks at about 35 miles of water; out of a frame spanning
# the whole country that view is *46 pixels wide*, and blown up to a phone
# screen it is mush. The renderer was right and the delivery was wrong.
#
# So each step is rendered once per cell instead. A cell is 5 x 4 degrees —
# roughly 480 x 440 km — drawn at 2400 x 1920, about 0.2 km/px. That is finer
# than MRMS's own 1 km grid, so nothing the mosaic knows is lost, and the cubic
# smoothing is baked in before the file reaches the phone: the app then
# magnifies a picture with no blocks in it rather than a grid.
#
# 2400 rather than 1200 was measured, not guessed. At 1200 a 35-mile view is
# 134 px and still soft; at 2400 it is 269 px and reads clean. The cost of the
# difference is 107 KB against 42 — and it is still *lighter* than the 400 KB
# national frame it replaces, which could not serve this zoom at all.
from cells import (CELL_ORIGINS, CELL_SIZE, CELL_SPAN,  # noqa: F401
                   MID_ORIGINS, MID_SIZE, MID_SPAN,
                   NATIONAL_BBOX, NATIONAL_ID, NATIONAL_SIZE)


# Set before the pool forks, read inside the workers. Module-level because a
# forked child inherits them for free; passing the decoded grid as an argument
# would pickle tens of megabytes per cell.
_VALUES = _META = _SIZE = _OUT = _STAMP = None


def _render_one(cell):
    """Render one cell of the frame the pool was forked for."""
    name, bbox = cell
    try:
        # The national frame is the model's own resolution rather than a
        # cell's: at continental zoom the extra pixels carry nothing and cost
        # a boater bytes.
        # Each tier at its own size: close cells sharp, the middle tier at the
        # density the forecast half publishes, the national frame at the
        # model's own resolution.
        if name == NATIONAL_ID:
            px = NATIONAL_SIZE
        elif name.startswith("wide-"):
            px = MID_SIZE
        else:
            px = _SIZE
        image = render(_VALUES, _META, bbox, px)
        image.save(os.path.join(_OUT, f"{name}-{_STAMP}.png"), optimize=True)
        return name, True, None
    except Exception as e:            # noqa: BLE001 - reported, never substituted
        return name, False, str(e)


def cells():
    """Every cell as (id, bbox) with bbox west/south/east/north.

    The national frame comes last, and is the same whole-domain picture the
    forecast half publishes. NOAA's own national tiles served this zoom for one
    release and could not animate: every step is a different time, so each one
    refetched twenty-odd tiles over the network and the map flashed once per
    step. One image does not.
    """
    out = [(name, (west, south, west + CELL_SPAN[0], south + CELL_SPAN[1]))
           for name, west, south in CELL_ORIGINS]
    # The middle tier, prefixed so its ids cannot collide with a close cell's.
    # A view too wide for a close cell lands here instead of falling all the way
    # to the national frame, which is blocky over a few hundred miles.
    out += [(f"wide-{name}", (west, south, west + MID_SPAN[0], south + MID_SPAN[1]))
            for name, west, south in MID_ORIGINS]
    out.append((NATIONAL_ID, NATIONAL_BBOX))
    return out

# Written as a 256-colour PNG rather than RGBA.
#
# Measured against a real frame it is 543 KB instead of 1,593 KB — a third of
# the bytes — and the two are indistinguishable side by side, because the
# picture only ever contains ramp colours at a fixed set of alphas. The
# boater's data allowance is part of the product; a picture that costs three
# times as much to look at is not a better picture.
PALETTE_COLOURS = 256   # no longer used to write frames; kept for the nowcast import
# No smoothing on measured frames. A 1 km core is often one pixel, and on real
# frames a smooth of a quarter cell cut small cores by 4 to 8 dBZ — the table is
# in `smoothing.py`. The measured half gets its smoothness from dither-free
# RGBA output and the cubic sample alone (D-65, C-15).
SMOOTH_CELLS = 0.0

# Values below this are "no echo" or "no coverage" rather than weak returns.
# MRMS uses -99 for the first and -999 for the second; both must be pulled up
# before interpolation or cubic ringing throws huge negative overshoot into the
# edge of every cell.
NO_ECHO = -30.0


def _get(url, timeout=90):
    return urllib.request.urlopen(urllib.request.Request(url), timeout=timeout).read()


def list_frames(day):
    """Every published frame for a UTC day, oldest first."""
    url = f"{BUCKET}/?list-type=2&prefix=CONUS/{PRODUCT}/{day}/&max-keys=1000"
    body = _get(url, timeout=60).decode()
    return re.findall(r"<Key>([^<]+)</Key>", body)


def observed_time(key):
    """The moment the mosaic describes, read off the filename."""
    m = re.search(r"_(\d{8})-(\d{6})\.grib2", key)
    if not m:
        return None
    return datetime.datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S").replace(
        tzinfo=datetime.timezone.utc)


def wanted_keys(now=None):
    """One key per timeline step, newest last.

    Frames arrive about every two minutes, so each step takes the closest
    frame at or before it. A step with nothing within half a step is left out
    rather than filled from further away — a gap the app can state is better
    than a frame under the wrong label.
    """
    now = now or datetime.datetime.now(datetime.timezone.utc)
    keys = []
    for day in {(now - datetime.timedelta(minutes=HISTORY_MINUTES)).strftime("%Y%m%d"),
                now.strftime("%Y%m%d")}:
        try:
            keys += list_frames(day)
        except Exception as e:
            print(f"  listing {day} failed: {e}", file=sys.stderr)
    stamped = sorted((observed_time(k), k) for k in keys if observed_time(k))
    if not stamped:
        return []

    newest = stamped[-1][0]
    tolerance = datetime.timedelta(minutes=STEP_MINUTES / 2)
    out, seen = [], set()
    steps = range(HISTORY_MINUTES, -1, -STEP_MINUTES)
    for back in steps:
        target = newest - datetime.timedelta(minutes=back)
        candidates = [(abs(t - target), t, k) for t, k in stamped if t <= target + tolerance]
        if not candidates:
            continue
        gap, t, k = min(candidates)
        if gap <= tolerance and k not in seen:
            seen.add(k)
            out.append((t, k))
    return out


def read_frame(key):
    """dBZ on the mosaic's own lat/lon grid, plus what it takes to place it."""
    raw = _get(f"{BUCKET}/{key}", timeout=120)
    data = gzip.decompress(raw) if key.endswith(".gz") else raw
    path = os.path.join(os.environ.get("TMPDIR", "/tmp"), "mrms-frame.grib2")
    with open(path, "wb") as f:
        f.write(data)
    try:
        with open(path, "rb") as f:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                raise ValueError("not a GRIB message")
            meta = {k: eccodes.codes_get(gid, k) for k in (
                "Ni", "Nj",
                "latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees",
                "iDirectionIncrementInDegrees", "jDirectionIncrementInDegrees",
                "jScansPositively")}
            values = eccodes.codes_get_values(gid).reshape(meta["Nj"], meta["Ni"])
            eccodes.codes_release(gid)
        return values, meta, len(raw)
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def render(values, meta, bbox, size):
    """Cubic-spline sample from the mosaic's lat/lon grid into web mercator.

    The same sampling the forecast frames use, and for the same reason: a
    storm edge is not a square, and nearest-neighbour draws one anyway. See
    `render.render` for the measurements that settled cubic over bilinear and
    over a blur — blurring shaves the top off every core, which understates
    exactly what a boater needs to see.
    """
    # No-coverage and no-echo are not weak returns. Left at -999 they ring
    # violently through a cubic filter and paint a dark halo round every cell.
    field = np.where(values < NO_ECHO, NO_ECHO, values)

    lon0 = meta["longitudeOfFirstGridPointInDegrees"]
    if lon0 > 180:
        lon0 -= 360
    lat0 = meta["latitudeOfFirstGridPointInDegrees"]
    di = meta["iDirectionIncrementInDegrees"]
    dj = meta["jDirectionIncrementInDegrees"]
    # MRMS scans north to south, so a row index counts *down* from lat0.
    south_up = bool(meta.get("jScansPositively"))

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

    fx = (lon - lon0) / di
    fy = (lat - lat0) / dj if south_up else (lat0 - lat) / dj

    inside = (fx >= 0) & (fx < meta["Ni"] - 1) & (fy >= 0) & (fy < meta["Nj"] - 1)
    interpolated = ndimage.map_coordinates(field, [fy, fx], order=3,
                                           mode="nearest", prefilter=True)
    sampled = np.where(inside, interpolated, -99.0)
    sampled, before, after = smooth(sampled, di, bbox, size, SMOOTH_CELLS)
    if before - after > 0.3:
        print(f"smooth: peak {before:.1f} -> {after:.1f} dBZ on {size[0]}x{size[1]}", flush=True)
    return Image.fromarray(colourise(sampled))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/observed", help="directory to write frames into")
    ap.add_argument("--width", type=int, default=CELL_SIZE[0], help="cell width in pixels")
    ap.add_argument("--limit", type=int, default=0, help="render at most this many steps")
    ap.add_argument("--only", default="", help="render just this cell, for checking one place")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    size = (args.width, args.width * CELL_SIZE[1] // CELL_SIZE[0])

    steps = wanted_keys()
    if args.limit:
        steps = steps[-args.limit:]
    if not steps:
        raise SystemExit("no MRMS frames listed; leaving the previous manifest in place")

    wanted_cells = [c for c in cells() if not args.only or c[0] == args.only]

    frames, downloaded, reused, rendered = [], 0, 0, 0
    for when, key in steps:
        stamp = when.strftime("%Y%m%dT%H%M") + "Z"
        # Which cells still need this step. One GRIB read serves all of them,
        # so the download only happens if at least one is missing.
        todo = [(name, bbox) for name, bbox in wanted_cells
                if not os.path.exists(os.path.join(args.out, f"{name}-{stamp}.png"))]
        values = meta = None
        if todo:
            try:
                values, meta, nbytes = read_frame(key)
                downloaded += nbytes
            except Exception as e:
                # Absent, never substituted.
                print(f"  {when:%H:%M}Z  SKIPPED: {e}", file=sys.stderr)
                continue

        # Cells are rendered in parallel across the runner's cores.
        #
        # Measured 2026-09-02 on a real run: this step took **689 seconds**,
        # eleven and a half minutes, and it is the reason the app kept showing
        # "NOAA's live radar isn't responding". The frames were correct; they
        # were simply older than the 30-minute staleness limit by the time the
        # run finished, so the app refused to draw them — which is the right
        # behaviour and must not be loosened.
        #
        # The work is 43 independent reprojections of one grid onto one cell
        # each, run one after another on a four-core machine. Nothing about it
        # was serial except the loop.
        #
        # `fork` is the default on Linux, so each worker inherits the decoded
        # grid rather than having it pickled across — the array is tens of
        # megabytes and sending it per cell would cost more than the drawing.
        for name, bbox in wanted_cells:
            image_name = f"{name}-{stamp}.png"
            if os.path.exists(os.path.join(args.out, image_name)):
                frames.append({"observedTime": when.isoformat().replace("+00:00", "Z"),
                               "cell": name, "image": image_name})
                reused += 1

        missing = [(n, b) for n, b in wanted_cells
                   if not os.path.exists(os.path.join(args.out, f"{n}-{stamp}.png"))]
        if missing:
            global _VALUES, _META, _SIZE, _OUT, _STAMP
            _VALUES, _META, _SIZE, _OUT, _STAMP = values, meta, size, args.out, stamp
            with multiprocessing.Pool(processes=min(len(missing), os.cpu_count() or 2)) as pool:
                for name, ok, err in pool.imap_unordered(_render_one, missing):
                    if ok:
                        frames.append({"observedTime": when.isoformat().replace("+00:00", "Z"),
                                       "cell": name,
                                       "image": f"{name}-{stamp}.png"})
                        rendered += 1
                    else:
                        print(f"  {when:%H:%M}Z {name}: SKIPPED: {err}", file=sys.stderr)

        if todo:
            kb = sum(os.path.getsize(os.path.join(args.out, f"{n}-{stamp}.png"))
                     for n, _ in wanted_cells
                     if os.path.exists(os.path.join(args.out, f"{n}-{stamp}.png"))) / 1024
            print(f"  {when:%H:%M}Z  {len(todo)} cells  ->  {kb:6.0f} KB total")

    if not frames:
        raise SystemExit("no frames rendered; leaving the previous manifest in place")

    manifest = {
        "product": "MRMS MergedReflectivityQCComposite",
        "field": "composite reflectivity",
        "source": "NOAA MRMS, public domain",
        "generatedAt": datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds").replace("+00:00", "Z"),
        "stepMinutes": STEP_MINUTES,
        "dbzFloor": RAMP[0][0],
        "smoothing": SMOOTH_CELLS,
        "alpha": ALPHA,
        "size": {"width": size[0], "height": size[1]},
        "cells": [{"id": name,
                   "bbox": {"west": b[0], "south": b[1], "east": b[2], "north": b[3]}}
                  for name, b in wanted_cells],
        "frames": frames,
    }
    with open(os.path.join(args.out, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"{len(steps)} steps x {len(wanted_cells)} cells: "
          f"{rendered} rendered, {reused} reused, {downloaded/1e6:.1f} MB pulled from MRMS")


if __name__ == "__main__":
    main()
