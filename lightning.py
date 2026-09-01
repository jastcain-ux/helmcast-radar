#!/usr/bin/env python3
"""Observed lightning for SeaWise, from NOAA's GOES Geostationary Lightning Mapper.

Public domain, no key, no licence, no per-user cost — unlike every commercial
lightning network, which is quote-priced and would cost more than the app earns.

**Points, not a field.** The radar halves are rendered as images because
reflectivity is a continuous field. Lightning is a list of coordinates, so it
travels as JSON: far smaller, and it lets the app fade a strike as it ages and
size it by zoom rather than baking those decisions into a raster.

What this is and is not: GLM is a satellite looking down, detecting total
lightning — in-cloud as well as cloud-to-ground, so it sees *more* than a
ground network. What it gives up is position. Roughly 8-10 km, against about
150 m for NLDN. Good enough for "there is lightning on your water"; not good
enough to say which dock was hit, and the app must not imply otherwise.
"""
import argparse, io, json, os, re, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np

# GOES-19 became GOES-East in 2025, replacing GOES-16; GOES-18 is West. Both
# publish GLM. East alone covers every coast SeaWise supports except the
# Pacific, so West is fetched too.
SATELLITES = ["noaa-goes19", "noaa-goes18"]

# Everything SeaWise forecasts for: CONUS coasts, the Gulf, the Great Lakes,
# and the inland water in between. Trimming here is what keeps the payload
# small — the raw files cover a hemisphere.
BOUNDS = dict(south=23.0, north=50.0, west=-126.0, east=-64.0)

S3 = "https://{bucket}.s3.amazonaws.com/"

# The same cells the measured radar already publishes in, and for the same
# reason at a different scale: a national file is 332 KB and only twelve of its
# 7,648 strikes were anywhere near Galveston Bay. A boater downloads their own
# cell — typically a few KB, about 40 KB in the busiest weather in the country.
#
# Imported from the radar renderer rather than copied, so a cell can never
# exist for radar and not for lightning. `test_cells.py` already guards the
# layout against gaps.
# Flat beside the radar scripts in the renderer repo; one directory over in the
# app repo. Try the neighbour first so the deployed layout needs no path fixing.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "hrrr-radar"))
from cells import CELL_ORIGINS, CELL_SPAN  # noqa: E402


def list_keys(bucket, prefix):
    url = S3.format(bucket=bucket) + "?list-type=2&max-keys=1000&prefix=" + prefix
    with urllib.request.urlopen(url, timeout=60) as r:
        return re.findall(r"<Key>([^<]+)</Key>", r.read().decode())


def started_at(key):
    """The s-timestamp in the filename: sYYYYDDDHHMMSSt, tenths of a second."""
    m = re.search(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})", key)
    if not m:
        return None
    y, doy, hh, mm, ss = (int(g) for g in m.groups())
    return (datetime(y, 1, 1, tzinfo=timezone.utc)
            + timedelta(days=doy - 1, hours=hh, minutes=mm, seconds=ss))


def recent_keys(bucket, minutes, now):
    """Keys for the last `minutes`, crossing the hour boundary when needed."""
    keys = []
    for back in range(0, minutes // 60 + 2):
        t = now - timedelta(hours=back)
        prefix = f"GLM-L2-LCFA/{t.year}/{t.timetuple().tm_yday:03d}/{t.hour:02d}/"
        try:
            keys += list_keys(bucket, prefix)
        except Exception as e:
            print(f"  {bucket} {prefix}: {e}", file=sys.stderr)
    cutoff = now - timedelta(minutes=minutes)
    return [k for k in keys if (started_at(k) or now) >= cutoff]


def flashes_in(bucket, key):
    """Flashes inside BOUNDS from one 20-second file, as (lat, lon, unix_ts)."""
    try:
        with urllib.request.urlopen(S3.format(bucket=bucket) + key, timeout=60) as r:
            raw = r.read()
        with h5py.File(io.BytesIO(raw), "r") as f:
            lat = f["flash_lat"][:].astype("float64")
            lon = f["flash_lon"][:].astype("float64")
            # Seconds since the file's own epoch attribute; the filename start
            # time is accurate to the 20-second window, which is finer than the
            # app displays.
            t0 = started_at(key)
    except Exception as e:
        print(f"  skip {key.split('/')[-1]}: {e}", file=sys.stderr)
        return []

    keep = ((lat >= BOUNDS["south"]) & (lat <= BOUNDS["north"])
            & (lon >= BOUNDS["west"]) & (lon <= BOUNDS["east"]))
    ts = int(t0.timestamp()) if t0 else 0
    # Rounded to 3 decimals — about 100 m, well inside GLM's own 8-10 km
    # accuracy, and it nearly halves the JSON.
    return [(round(float(a), 3), round(float(o), 3), ts)
            for a, o in zip(lat[keep], lon[keep])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--minutes", type=int, default=30,
                    help="how far back to publish")
    ap.add_argument("--fetch-minutes", type=int, default=15,
                    help="how far back to download when carrying strikes "
                         "forward; overlaps the last run on purpose")
    ap.add_argument("--out", default="public/lightning")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    now = datetime.now(timezone.utc)

    # Strikes already published, so a ten-minute run fetches ten minutes of
    # files rather than the whole thirty-minute window three times over. Same
    # principle as the radar reusing frames already on disk.
    #
    # Deduplicated as a set: the fetch window deliberately overlaps the last
    # run, because a file that was still being written when the previous run
    # listed the bucket would otherwise be lost for good.
    kept = set()
    cutoff = (now - timedelta(minutes=args.minutes)).timestamp()
    for name, _, _ in CELL_ORIGINS:
        path = os.path.join(args.out, f"{name}.json")
        if not os.path.exists(path):
            continue
        try:
            with open(path) as f:
                for s in json.load(f).get("strikes", []):
                    if s["t"] >= cutoff:
                        kept.add((s["lat"], s["lon"], s["t"]))
        except Exception as e:
            print(f"  ignoring unreadable {name}.json: {e}", file=sys.stderr)
    if kept:
        print(f"carried forward: {len(kept)} strikes still inside the window")

    strikes = list(kept)
    # Overlap the previous run rather than fetching exactly since it: a file
    # still being written when the last run listed the bucket would otherwise
    # never be picked up.
    fetch_minutes = args.minutes if not kept else args.fetch_minutes
    for bucket in SATELLITES:
        keys = recent_keys(bucket, fetch_minutes, now)
        print(f"{bucket}: {len(keys)} files in the last {fetch_minutes} min")
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for got in pool.map(lambda k: flashes_in(bucket, k), keys):
                strikes.extend(got)

    # Newest last, so the app can fade by age without sorting.
    # Deduplicate across the overlap, then oldest first so the app can fade by
    # age without sorting.
    strikes = sorted(set(strikes), key=lambda s: s[2])
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    span_lon, span_lat = CELL_SPAN
    stamp = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    index, total_bytes = [], 0

    for name, west, south in CELL_ORIGINS:
        east, north = west + span_lon, south + span_lat
        inside = [(a, o, t) for a, o, t in strikes
                  if south <= a <= north and west <= o <= east]
        cell = {
            "source": "NOAA GOES GLM",
            "generatedAt": stamp,
            "windowMinutes": args.minutes,
            # Stated rather than assumed: the app draws an area, not a pinpoint.
            # GLM locates a flash to roughly 10 km, and a map that implies a
            # struck dock would be claiming precision the satellite never had.
            "positionAccuracyKm": 10,
            "bounds": {"south": south, "north": north, "west": west, "east": east},
            "strikes": [{"lat": a, "lon": o, "t": t} for a, o, t in inside],
        }
        path = os.path.join(out_dir, f"{name}.json")
        with open(path, "w") as f:
            json.dump(cell, f, separators=(",", ":"))
        size = os.path.getsize(path)
        total_bytes += size
        index.append({"cell": name, "south": south, "north": north,
                      "west": west, "east": east, "strikes": len(inside)})

    # One small index so the app can pick its cell without guessing at the
    # layout, and so a dead run is visible as a stale generatedAt rather than
    # as an empty map.
    with open(os.path.join(out_dir, "index.json"), "w") as f:
        json.dump({"generatedAt": stamp, "windowMinutes": args.minutes,
                   "cells": index}, f, separators=(",", ":"))

    busiest = max(index, key=lambda c: c["strikes"]) if index else None
    # Published, not fetched. Most of what comes back is over the interior,
    # Mexico or open ocean beyond the cells and is never written — reporting
    # the national count made a healthy run look like it was losing strikes.
    published = sum(c["strikes"] for c in index)
    print(f"\n{published} strikes published across {len(CELL_ORIGINS)} cells "
          f"({total_bytes/1024:.0f} KB total); {len(strikes)} fetched nationally")
    if busiest:
        b = os.path.getsize(os.path.join(out_dir, busiest["cell"] + ".json"))
        print(f"busiest cell: {busiest['cell']} "
              f"{busiest['strikes']} strikes ({b/1024:.0f} KB) — "
              f"what one boater actually downloads")


if __name__ == "__main__":
    main()
