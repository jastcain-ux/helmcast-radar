"""The smoothing helper, pinned. Pure numpy; runs without eccodes."""
import unittest
import numpy as np
from smoothing import smooth, sigma_px, OUTSIDE

BBOX = (-100.0, 25.0, -90.0, 32.0)
SIZE = (2000, 1400)          # 200 px per degree, the forecast close tier
HRRR = 0.027                 # degrees per 3 km cell


class SmoothingTests(unittest.TestCase):

    def test_sigma_is_a_fraction_of_one_model_cell_in_output_pixels(self):
        self.assertAlmostEqual(sigma_px(HRRR, BBOX, SIZE, 0.3), 0.3 * HRRR * 200, places=6)

    def test_zero_strength_changes_nothing(self):
        f = np.random.default_rng(1).uniform(-10, 60, (80, 80))
        out, before, after = smooth(f, HRRR, BBOX, SIZE, 0.0)
        np.testing.assert_array_equal(out, f)
        self.assertEqual(before, after)

    def test_a_broad_core_keeps_its_peak(self):
        f = np.full((200, 200), 20.0); f[80:120, 80:120] = 60.0
        out, before, after = smooth(f, HRRR, BBOX, SIZE, 0.3)
        self.assertEqual(before, 60.0)
        self.assertAlmostEqual(after, 60.0, places=5)
        self.assertAlmostEqual(out.max(), 60.0, places=5)

    def test_outside_the_grid_stays_outside_and_never_bleeds_in(self):
        f = np.full((100, 100), 40.0); f[:, 50:] = OUTSIDE
        out, _, _ = smooth(f, HRRR, BBOX, SIZE, 0.3)
        np.testing.assert_array_equal(out[:, 50:], OUTSIDE)
        self.assertGreaterEqual(out[:, :50].min(), 40.0 - 1e-6,
                                "a -99 sentinel next to real data used to ring a dark rim round every cell")

    def test_the_cost_on_a_single_model_cell_is_reported_not_hidden(self):
        # One 3 km cell at 70 dBZ over 20: the case the measured half refuses
        # to smooth (strength 0). The helper must report the loss so the run
        # log shows it rather than a frame quietly understating a core.
        f = np.full((60, 60), 20.0); f[27:33, 27:33] = 70.0
        _, before, after = smooth(f, HRRR, BBOX, SIZE, 0.3)
        self.assertEqual(before, 70.0)
        self.assertLess(after, 70.0)
        self.assertGreater(after, 55.0)


if __name__ == "__main__":
    unittest.main()
