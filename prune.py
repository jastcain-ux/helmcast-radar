#!/usr/bin/env python3
"""Drop frames the timeline can no longer reach.

The published directory is carried between runs, so nothing removes old
frames on its own. Left alone it would grow by ~400 KB every ten minutes for
ever — about 58 MB a day for the measured half alone.

Anything not named in the manifest goes. The manifest is written last by
`observed.py` and `render.py` and lists exactly the frames the app can ask
for, so it is the only thing that needs to be right.

Run it over the forecast directory too. Forecast frame names are derived from
the cell and the lead time, so they normally overwrite in place — but the
moment the cell list changes, every frame belonging to a retired cell is
orphaned and would be carried for ever.
"""
import argparse, json, os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    args = ap.parse_args()

    manifest_path = os.path.join(args.dir, "manifest.json")
    if not os.path.exists(manifest_path):
        print("no manifest; pruning nothing")
        return
    with open(manifest_path) as f:
        keep = {frame["image"] for frame in json.load(f).get("frames", [])}
    keep.add("manifest.json")

    removed = 0
    for name in os.listdir(args.dir):
        path = os.path.join(args.dir, name)
        # The forecast frames sit in the root of `public/` and the measured
        # ones in `public/observed`, so pruning the root must step over that
        # directory rather than trying to unlink it.
        if name in keep or os.path.isdir(path):
            continue
        os.unlink(path)
        removed += 1
    print(f"kept {len(keep) - 1} frames, removed {removed}")


if __name__ == "__main__":
    main()
