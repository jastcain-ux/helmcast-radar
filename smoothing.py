"""A light smooth of the reprojected reflectivity field, after the cubic upsample.

D-65 and C-15. The cubic sample keeps every model value, which at 200 px per
degree shows a 3 km model's cells as soft squares. Smoothing *after* the
upsample, at a fraction of one model cell, rounds those off. It also shaves
peaks, and the fraction is chosen from what that costs on real frames
(2026-09-02, whole published frames, every local maximum above 45 dBZ):

    HRRR 3 km, 200 px/deg   0.25 cell: max 0.1  p99 1.0    0.3: 0.1 / 1.3    0.4: 0.3 / 2.0
    MRMS 1 km, 480 px/deg   0.25 cell: max 3.9  p99 7.5    0.3: 5.4 / 9.9    0.4: 8.2 / 14.2

A measured 1 km core is often a single pixel, and any smooth at all understates
it. So the smooth is followed by `restore_peaks`, which puts back what the
smooth took off the tops: the loss of the local maximum, spread over half a
cell, is added back around every peak. Measured on the same frames at 0.7 of a
cell (D-66, Jason's choice against four candidates):

    HRRR 3 km forecast   peak 0.0   p99 small-core loss 1.6 dBZ
    MRMS 1 km measured   peak 0.0   p99 small-core loss 0.9 dBZ

Every product now uses 0.7 with peaks restored. A frame that still loses more
than a third of a ramp step at its peak is named in the run log.

Pure numpy and scipy, so `test_render.py` can import it without eccodes.
"""
import numpy as np
from scipy import ndimage

OUTSIDE = -99.0


def sigma_px(cell_deg, bbox, size, cells):
    """Gaussian sigma in output pixels for `cells` of one model cell."""
    west, _, east, _ = bbox
    return cells * cell_deg * (size[0] / (east - west))


def smooth(sampled, cell_deg, bbox, size, cells):
    """Smooth the field where it is inside the model grid; outside stays OUTSIDE.

    Returns (field, peak_before, peak_after) so the caller can log the cost.
    Outside pixels are filled with the field's own minimum before filtering so
    the sentinel never bleeds in as a dark rim, then restored.
    """
    inside = sampled > OUTSIDE + 0.5
    if not inside.any():
        return sampled, float("nan"), float("nan")
    before = float(sampled[inside].max())
    sigma = sigma_px(cell_deg, bbox, size, cells)
    if sigma <= 0.05:
        return sampled, before, before
    field = np.where(inside, sampled, sampled[inside].min())
    field = ndimage.gaussian_filter(field, sigma, mode="nearest")
    after = float(field[inside].max())
    return np.where(inside, field, OUTSIDE), before, after


def restore_peaks(original, smoothed, cell_px):
    """Add back what the smooth took off the tops.

    The local maximum of the smoothed field sits below the original's near
    every core; that difference, spread over half a cell so it does not draw
    a dot, is added around each peak. Values only go up, never down, so no
    core is painted weaker than the field said.
    """
    k = max(int(round(2 * cell_px)) | 1, 3)
    gain = np.clip(ndimage.maximum_filter(original, k) - ndimage.maximum_filter(smoothed, k),
                   0.0, None)
    return smoothed + ndimage.gaussian_filter(gain, max(cell_px * 0.5, 0.5), mode="nearest")


def shape(sampled, cell_deg, bbox, size, cells):
    """Smooth at `cells` of a model cell, then restore the peaks.

    Same contract as `smooth`: outside stays OUTSIDE, returns the field and
    the peak before and after so the caller can log the cost.
    """
    inside = sampled > OUTSIDE + 0.5
    if not inside.any():
        return sampled, float("nan"), float("nan")
    before = float(sampled[inside].max())
    sigma = sigma_px(cell_deg, bbox, size, cells)
    if sigma <= 0.05:
        return sampled, before, before
    field = np.where(inside, sampled, sampled[inside].min())
    smoothed = ndimage.gaussian_filter(field, sigma, mode="nearest")
    west, _, east, _ = bbox
    out = restore_peaks(field, smoothed, cell_deg * size[0] / (east - west))
    after = float(out[inside].max())
    return np.where(inside, out, OUTSIDE), before, after
