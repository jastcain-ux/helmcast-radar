# Forecast radar frames

Renders NOAA HRRR composite reflectivity into PNG frames the app draws over
the map for the next six hours.

## Why the app needs this at all

Radar measures what is already in the air. NOAA publishes observed
reflectivity as map tiles, but **no forecast reflectivity** — that exists only
as HRRR model output in GRIB2 on NOMADS. Every weather app showing "future
radar" is either extrapolating recent frames (useful for about an hour, then
it falls apart) or rendering this same model, usually bought from a vendor
whose pricing is quote-only.

The data is free and public domain. The missing piece was something to turn
the files into images. That is all this is.

## What it costs

About 16 MB pulled from NOMADS per run, because the GRIB index sidecar lets
us range-request only the reflectivity record — 0.3 MB instead of the 130 MB
full file. Frames every 15 minutes out to 9 hours, although the app scrubs six. Quarter
hours are HRRR's own sub-hourly cadence, so each step is a real model frame
rather than an interpolated one. The extra hours cover run lag: the newest
complete run is always one to two hours behind the clock, so a six-hour series
ran out before the end of a six-hour scrubber.

Note that sub-hourly GRIB files hold four lead times each, so the index has to
be matched on the lead label as well as the field name — taking the first REFC
record returns :15 for every quarter hour asked for.

Output is 36 national PNGs at 4096x2304 — about 1.4 km per pixel, finer than
the model's own 3 km grid because the app upscales a national frame down to a
single bay. Overwritten hourly. Rendering takes a few seconds.

## Running it by hand

```bash
pip install -r requirements.txt
python render.py --out out
```

Writes `refc-f01.png` … `refc-f06.png` and `manifest.json` into `out/`.

## Setting up the hourly job

The workflow in `workflow/hrrr-radar.yml` is meant to live in a **separate
public repository**, because the app repo is private and private repos get
neither unlimited Actions minutes nor free Pages hosting. Nothing sensitive is
published by doing this — it is a script reading public NOAA data.

1. Create a public repo, e.g. `helmcast-radar`
2. Copy `render.py`, `requirements.txt`, and `workflow/hrrr-radar.yml`
   (as `.github/workflows/hrrr-radar.yml`) into it
3. In that repo: **Settings → Pages → Source → GitHub Actions**
4. Run the workflow once by hand from the Actions tab to check it

Frames then appear at `https://<user>.github.io/helmcast-radar/manifest.json`.

## Before launch

GitHub Pages is fine for testing and early users — its bandwidth allowance
covers tens of thousands of sessions a month — but it is not intended as a
production CDN for a commercial app. Move the publish step to object storage
(Cloudflare R2 has no egress fees) when there are enough subscribers for it to
matter. Only the last step of the workflow changes; the renderer does not.

## Rules that are not negotiable

- **A frame that could not be rendered is absent from the manifest**, never
  replaced by a neighbouring hour. The app has to be able to say "no frame for
  this hour" rather than show 3 PM's weather under a 6 PM label.
- **The manifest carries the model run time**, so the app can show how old the
  forecast is. Six hours of outlook from a run that stopped updating four
  hours ago is not six hours of outlook.
- **Below 5 dBZ is transparent.** Drizzle nobody can feel must not paint the
  bay green.
- **Bilinear sampling, and a soft outer edge.** This started as
  nearest-neighbour on the theory that smoothing invents gradients the model
  did not produce. That was wrong. Reflectivity is a continuous field sampled
  every 3 km, and hard square cells assert a sharp boundary exactly where the
  model is least certain — a storm edge is not a 3 km square. The colour
  *bands* stay discrete, so nothing is invented about intensity; only the shape
  is smoothed.

## It is a model, not a measurement

HRRR is 3 km. It will show a squall line crossing the bay, and it will
sometimes put that line in the wrong place or an hour off. It gets the same
treatment as every other forecast in this app: its own palette, labelled as
forecast, never blended into the measured half of the timeline.


## observed.py — the measured half

`render.py` draws the forecast frames. `observed.py` draws the **measured**
ones, from NOAA MRMS composite reflectivity, and exists because of a look
rather than a gap in the data.

The measured radar used to be NOAA's nowCOAST WMS: someone else's pre-drawn
tiles, painted nearest-neighbour with hard colour bands. At the zoom a boater
uses it read as a grid of coloured squares. That is not a resolution problem —
a sharper source would only have given sharper squares — it is a rendering
problem, and the fix was already in this directory. `observed.py` imports
`render.colourise` unchanged and gets the same cubic sampling, continuous ramp
and soft edge that made the forecast frames stop looking broken.

Two things fall out of that beyond the picture:

- **The palettes cannot drift.** They used to be matched by hand, and when they
  drifted rain appeared to thin out exactly at "now" as the scrubber crossed
  over. One renderer means one ramp and one floor.
- **Frames are regional cells, not one national picture.** That was the first
  version and it failed at the only zoom that matters: a 35-mile view out of a
  country-wide frame is 46 pixels. Cells are 5 x 4 degrees at 2400 px, 25 of
  them over the coasts, the Great Lakes and the Gulf bays, so the same view is
  269 px and reads clean — and a cell is ~75 KB against the national frame's
  400 KB.

Frames are 256-colour PNGs: 398 KB against 1,593 KB for RGBA, and
indistinguishable side by side, because the picture only ever contains ramp
colours at a fixed set of alphas.

```bash
python observed.py --out public/observed      # 13 steps x 25 cells, 2 hours
python observed.py --only florida-ne --limit 1   # one cell, for checking a place
python prune.py --dir public/observed         # drop anything the manifest dropped
```

Frames already on disk are reused, so a ten-minute run pulls one new frame
rather than re-rendering the lot.
