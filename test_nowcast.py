#!/usr/bin/env python3
"""The nowcast's motion tracking, on synthetic storms whose motion is known.

This file exists because two real bugs were caught by a dry run before the
runner ever executed the code, and would have been invisible there: the
renderer publishes with a cheerful 200 whatever the vectors say.

- Pooling 16x meant a phase-correlation peak could only land on multiples of
  16 km; a true 6-row, 10-column shift came back as 0 and 16.
- Empty windows borrowed the nearest tracked vector before smoothing, so a
  storm whose neighbours were nearer another storm had that storm's opposite
  motion averaged into its own — down to a fifth.

Needs numpy and scipy. The GRIB-reading modules are stubbed so this runs
anywhere those two install: `python -m unittest test_nowcast`.
"""
import sys, types, unittest

import numpy as np

# Stub what nowcast.py imports at module load; only the names it touches.
_obs = types.ModuleType("observed")
_obs.NO_ECHO, _obs.PALETTE_COLOURS = -30.0, 256
_obs.render = _obs.wanted_keys = _obs.read_frame = lambda *a, **k: None
sys.modules.setdefault("observed", _obs)
_ren = types.ModuleType("render")
_ren.RAMP, _ren.RAMP_STEPS, _ren.QUANTISER, _ren.PALETTE = [(15, (0, 0, 0))], 96, "x", "x"
sys.modules.setdefault("render", _ren)

import nowcast  # noqa: E402

H, W = 512, 768
YY, XX = np.mgrid[0:H, 0:W]
NE = -30.0


def storm(cy, cx, peak=45.0):
    return peak * np.exp(-(((YY - cy) / 18.0) ** 2 + ((XX - cx) / 28.0) ** 2))


def field(*storms):
    f = np.maximum.reduce([storm(*s) for s in storms])
    return np.where(f > 5.0, f, NE)


def peak_at(a, rows=slice(None)):
    sub = a[rows]
    y, x = np.unravel_index(np.argmax(sub), sub.shape)
    return y + (rows.start or 0), x


class MotionTracking(unittest.TestCase):

    def test_one_storm_is_tracked_and_carried_forward(self):
        prev, curr = field((240, 300)), field((246, 310))
        vy, vx = nowcast.motion_field(prev, curr, minutes_apart=10.0)
        self.assertEqual(vy.shape, prev.shape)
        self.assertAlmostEqual(vy[246, 310] * 10, 6.0, delta=1.5)
        self.assertAlmostEqual(vx[246, 310] * 10, 10.0, delta=1.5)
        adv = nowcast.advect(curr, vy, vx, 30.0)
        y, x = peak_at(adv)
        self.assertLessEqual(abs(y - 264), 3)
        self.assertLessEqual(abs(x - 340), 4)
        self.assertGreater(adv.max(), 40.0, "advection must not shave the core")

    def test_two_storms_moving_oppositely_keep_their_own_motion(self):
        prev = field((120, 150), (400, 600))
        curr = field((124, 158), (396, 592))
        vy, vx = nowcast.motion_field(prev, curr, 10.0)
        self.assertAlmostEqual(vy[124, 158] * 10, 4.0, delta=1.5)
        self.assertAlmostEqual(vx[124, 158] * 10, 8.0, delta=1.5)
        self.assertAlmostEqual(vy[396, 592] * 10, -4.0, delta=1.5)
        self.assertAlmostEqual(vx[396, 592] * 10, -8.0, delta=1.5)
        adv = nowcast.advect(curr, vy, vx, 30.0)
        ay, ax = peak_at(adv, slice(0, 256))
        by, bx = peak_at(adv, slice(256, H))
        self.assertLessEqual(abs(ay - 136), 4); self.assertLessEqual(abs(ax - 182), 5)
        self.assertLessEqual(abs(by - 384), 4); self.assertLessEqual(abs(bx - 568), 5)

    def test_a_small_shift_does_not_round_to_nothing(self):
        prev, curr = field((240, 300)), field((241, 302))
        vy, vx = nowcast.motion_field(prev, curr, 10.0)
        self.assertAlmostEqual(vy[241, 302] * 10, 1.0, delta=1.0)
        self.assertAlmostEqual(vx[241, 302] * 10, 2.0, delta=1.0)

    def test_no_echo_anywhere_means_no_field(self):
        empty = np.full((H, W), NE)
        self.assertIsNone(nowcast.motion_field(empty, empty, 10.0))


if __name__ == "__main__":
    unittest.main()
