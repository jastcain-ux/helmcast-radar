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

About 2.2 MB pulled from NOMADS per run, because the GRIB index sidecar lets
us range-request only the reflectivity record — 0.3 MB instead of the 130 MB
full file. Output is six national PNGs totalling roughly 1.7 MB, overwritten
hourly. Rendering takes a few seconds.

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
- **Nearest-neighbour sampling, never interpolation.** Smoothing reflectivity
  invents gradients between cells that the model did not produce.

## It is a model, not a measurement

HRRR is 3 km. It will show a squall line crossing the bay, and it will
sometimes put that line in the wrong place or an hour off. It gets the same
treatment as every other forecast in this app: its own palette, labelled as
forecast, never blended into the measured half of the timeline.
