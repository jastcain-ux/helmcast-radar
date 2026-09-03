"""A light smooth of the reprojected reflectivity field, after the cubic upsample.

D-65 and C-15. The cubic sample keeps every model value, which at 200 px per
degree shows a 3 km model's cells as soft squares. Smoothing *after* the
upsample, at a fraction of one model cell, rounds those off. It also shaves
peaks, and the fraction is chosen from what that costs on real frames
(2026-09-02, whole published frames, every local maximum above 45 dBZ):

    HRRR 3 km, 200 px/deg   0.25 cell: max 0.1  p99 1.0    0.3: 0.1 / 1.3    0.4: 0.3 / 2.0
    MRMS 1 km, 480 px/deg   0.25 cell: max 3.9  p99 7.5    0.3: 5.4 / 9.9    0.4: 8.2 / 14.2

A measured 1 km core is often a single pixel, and any smooth at all understates
it — so the measured half and the nowcast use a strength of 0 and rely on the
dither-free RGBA output for their smoothness. The forecast uses 0.3: a tenth of
a dBZ at the peak, about one ramp step on the smallest cores.

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
