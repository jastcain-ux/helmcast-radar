#!/usr/bin/env python3
"""Extrapolated measured radar: the last two MRMS frames, motion-tracked forward.

The reference apps' forecast stays crisp for the first hour or two because it
is not a model at all — it is the newest measured frame moved along the motion
seen between it and the one before. Ours went to HRRR's 3 km grid at
+15 minutes, which is why our forecast half could never look as sharp as our
measured half however well it was drawn. D-64.

**What this is, and is not.** It moves rain; it does not grow or kill it. A
steady system tracks for an hour or two; developing or collapsing convection
will not, and every app's nowcast shares that limit. The app captions these
frames "projected from the last hour" and never "forecast" — D-57. The manifest
says `"source": "extrapolation"` so nothing downstream can mistake it.

**Motion.** Block-wise phase correlation between the two newest frames on a
pooled grid: one vector per ~250 km block, smoothed, with blocks that have too
little echo to correlate inheriting their neighbours'. A single national vector
was tried first in thought and rejected — a front over Texas and a sea-breeze
line over Florida do not move together.

**Cost.** Wide tier at half resolution plus the national frame, eight leads.
The close cells are deliberately absent: an advected field is already a smear
of the source, and 480 px/degree of a smear costs render time for nothing. The
app's tier chooser then picks the wide frame at close zoom, which is the right
picture. Stamp-keyed on the newest MRMS frame, so a ten-minute tick that brings
no new frame renders nothing.
"""
import argparse, datetime, json, os, sys
from multiprocessing import Pool

import numpy as np
from PIL import Image
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import observed                                            # noqa: E402
import render as forecast                                   # noqa: E402
from cells import (MID_ORIGINS, MID_SPAN, MID_SIZE, NATIONAL_BBOX,  # noqa: E402
                   NATIONAL_ID, NATIONAL_SIZE)

LEAD_MINUTES = [15, 30, 45, 60, 75, 90, 105, 120]
# Pooling before correlation: the mosaic is 0.01 degrees and storm motion is a
# few kilometres a minute, so tracking at full resolution costs memory for no
# precision. Sixteen cells is ~16 km per pooled pixel.
POOL = 16
# Correlation blocks, in pooled pixels: ~16 of them ~250 km across CONUS.
BLOCK = 16
# Below this fraction of echo a block cannot be correlated and borrows motion.
MIN_ECHO_FRACTION = 0.02
# Wide tier at half resolution; see the module docstring.
NOWCAST_MID_SIZE = (MID_SIZE[0] // 2, MID_SIZE[1] // 2)
NOWCAST_NATIONAL_SIZE = (NATIONAL_SIZE[0] // 2, NATIONAL_SIZE[1] // 2)


def _pool(field):
    h, w = field.shape
    h2, w2 = h - h % POOL, w - w % POOL
    lit = np.where(field > observed.NO_ECHO, field, 0.0)[:h2, :w2]
    return lit.reshape(h2 // POOL, POOL, w2 // POOL, POOL).mean(axis=(1, 3))


def _phase_shift(a, b):
    """Integer shift (dy, dx) that best moves `a` onto `b`, by phase correlation."""
    fa, fb = np.fft.fft2(a), np.fft.fft2(b)
    cross = fb * np.conj(fa)
    denom = np.abs(cross)
    denom[denom == 0] = 1.0
    corr = np.fft.ifft2(cross / denom).real
    peak = np.unravel_index(np.argmax(corr), corr.shape)
    dy, dx = peak
    if dy > a.shape[0] // 2: dy -= a.shape[0]
    if dx > a.shape[1] // 2: dx -= a.shape[1]
    return dy, dx, float(corr.max())


def motion_field(prev, curr, minutes_apart):
    """Per-pixel (dy, dx) in grid cells per minute, on the full mosaic grid."""
    pa, pb = _pool(prev), _pool(curr)
    H, W = pa.shape
    vy = np.zeros((H // BLOCK + 1, W // BLOCK + 1))
    vx = np.zeros_like(vy)
    ok = np.zeros_like(vy, dtype=bool)
    for by in range(0, H, BLOCK):
        for bx in range(0, W, BLOCK):
            a = pa[by:by + BLOCK, bx:bx + BLOCK]
            b = pb[by:by + BLOCK, bx:bx + BLOCK]
            if a.size == 0 or (a > 0).mean() < MIN_ECHO_FRACTION:
                continue
            dy, dx, _ = _phase_shift(a, b)
            vy[by // BLOCK, bx // BLOCK] = dy
            vx[by // BLOCK, bx // BLOCK] = dx
            ok[by // BLOCK, bx // BLOCK] = True
    if not ok.any():
        return None
    # Blocks with no echo borrow the nearest tracked motion, then the whole
    # field is smoothed so neighbouring blocks do not tear the picture apart.
    idx = ndimage.distance_transform_edt(~ok, return_distances=False, return_indices=True)
    vy, vx = vy[tuple(idx)], vx[tuple(idx)]
    vy, vx = ndimage.gaussian_filter(vy, 1.0), ndimage.gaussian_filter(vx, 1.0)
    # Back to the full grid, in cells per minute.
    scale = POOL / max(minutes_apart, 1.0)
    zoom = (prev.shape[0] / vy.shape[0], prev.shape[1] / vy.shape[1])
    vy_full = ndimage.zoom(vy, zoom, order=1) * scale
    vx_full = ndimage.zoom(vx, zoom, order=1) * scale
    # `zoom` rounds its output size; the advection indexes this against the
    # full grid and a one-row mismatch would raise on the runner where it
    # cannot be seen. Trim or pad to the grid exactly.
    def fit(a):
        out = np.zeros(prev.shape, dtype=a.dtype)
        h, w = min(a.shape[0], prev.shape[0]), min(a.shape[1], prev.shape[1])
        out[:h, :w] = a[:h, :w]
        return out
    return fit(vy_full), fit(vx_full)


def advect(field, vy, vx, minutes):
    """The field carried along its motion for `minutes`. Semi-Lagrangian: each
    output cell looks back along the flow for its value."""
    H, W = field.shape
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    src_y = yy - vy * minutes
    src_x = xx - vx * minutes
    lit = np.where(field > observed.NO_ECHO, field, observed.NO_ECHO)
    return ndimage.map_coordinates(lit, [src_y, src_x], order=1,
                                   mode="constant", cval=observed.NO_ECHO)


_FIELD = _META = _OUT = _STAMP = None


def _render_one(cell):
    name, bbox = cell
    try:
        px = NOWCAST_NATIONAL_SIZE if name == NATIONAL_ID else NOWCAST_MID_SIZE
        image = observed.render(_FIELD, _META, bbox, px)
        image.quantize(colors=observed.PALETTE_COLOURS, method=Image.FASTOCTREE,
                       dither=Image.FLOYDSTEINBERG) \
             .save(os.path.join(_OUT, f"{name}-{_STAMP}.png"), optimize=True)
        return name, True, None
    except Exception as e:                                  # noqa: BLE001
        return name, False, str(e)


def tiers():
    out = [(f"wide-{n}", (w, s, w + MID_SPAN[0], s + MID_SPAN[1]))
           for n, w, s in MID_ORIGINS]
    out.append((NATIONAL_ID, NATIONAL_BBOX))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/nowcast")
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    steps = observed.wanted_keys()
    if len(steps) < 2:
        raise SystemExit("fewer than two MRMS frames; nothing to extrapolate")
    (t_prev, k_prev), (t_curr, k_curr) = steps[-2], steps[-1]
    stamp = t_curr.strftime("%Y%m%dT%H%M") + "Z"
    manifest_path = os.path.join(args.out, "manifest.json")

    # Nothing new since last run: the manifest already describes this stamp.
    if os.path.exists(manifest_path):
        try:
            if json.load(open(manifest_path)).get("runTime") == \
               t_curr.isoformat().replace("+00:00", "Z"):
                print(f"nowcast for {stamp} already published; nothing to do")
                return
        except Exception:                                   # noqa: BLE001
            pass

    prev, meta_prev, _ = observed.read_frame(k_prev)
    curr, meta, _ = observed.read_frame(k_curr)
    if prev.shape != curr.shape:
        raise SystemExit("MRMS grid changed between frames; refusing to track")
    apart = (t_curr - t_prev).total_seconds() / 60.0
    motion = motion_field(prev, curr, apart)
    if motion is None:
        raise SystemExit("no echo anywhere to track; nothing to extrapolate")
    vy, vx = motion
    print(f"tracked {stamp}: frames {apart:.0f} min apart, "
          f"median motion {np.median(np.hypot(vy, vx)) * 60 * 1.1:.0f} km/h")

    wanted = [c for c in tiers() if not args.only or c[0] == args.only]
    frames = []
    global _FIELD, _META, _OUT, _STAMP
    _META, _OUT = meta, args.out
    for minutes in LEAD_MINUTES:
        _FIELD = advect(curr, vy, vx, minutes)
        _STAMP = f"{stamp}-{minutes:03d}"
        when = t_curr + datetime.timedelta(minutes=minutes)
        with Pool() as pool:
            results = pool.map(_render_one, wanted)
        for name, ok, err in results:
            if ok:
                frames.append({"cell": name, "leadMinutes": minutes,
                               "validTime": when.isoformat().replace("+00:00", "Z"),
                               "image": f"{name}-{_STAMP}.png"})
            else:
                print(f"  +{minutes:3d} {name}: SKIPPED: {err}", file=sys.stderr)
        print(f"  +{minutes:3d} min -> {sum(1 for _, ok, _ in results if ok)} frames")

    if not frames:
        raise SystemExit("no nowcast frames rendered; leaving the previous manifest in place")

    manifest = {
        "model": "MRMS extrapolation",
        "field": "composite reflectivity",
        "source": "extrapolation",
        "runTime": t_curr.isoformat().replace("+00:00", "Z"),
        "generatedAt": datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds").replace("+00:00", "Z"),
        "size": {"width": NOWCAST_MID_SIZE[0], "height": NOWCAST_MID_SIZE[1]},
        "dbzFloor": forecast.RAMP[0][0],
        "rampSteps": forecast.RAMP_STEPS,
        "quantiser": forecast.QUANTISER,
        "palette": forecast.PALETTE,
        "cells": [{"id": n, "bbox": {"west": b[0], "south": b[1], "east": b[2], "north": b[3]}}
                  for n, b in wanted],
        "frames": frames,
    }
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"{len(LEAD_MINUTES)} leads x {len(wanted)} tiers = {len(frames)} nowcast frames")


if __name__ == "__main__":
    main()
