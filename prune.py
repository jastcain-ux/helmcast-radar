#!/usr/bin/env python3
"""Drop measured frames the timeline can no longer reach.

The published directory is committed, so nothing removes old frames on its
own. Left alone the repo and the CDN would grow by ~400 KB every ten minutes
for ever — about 58 MB a day.

Anything not named in the manifest goes. The manifest is written last by
`observed.py` and lists exactly the frames the app can ask for, so it is the
only thing that needs to be right.
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
        if name in keep:
            continue
        os.unlink(os.path.join(args.dir, name))
        removed += 1
    print(f"kept {len(keep) - 1} frames, removed {removed}")


if __name__ == "__main__":
    main()
