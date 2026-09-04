"""The measured timeline's frame choice, pinned. Stdlib only."""
import datetime
import unittest
from steps import select

UTC = datetime.timezone.utc


def frames(first, count, every=2, seconds=40):
    """MRMS-shaped keys: one every two minutes, forty seconds past the minute."""
    out = []
    for i in range(count):
        t = first + datetime.timedelta(minutes=every * i, seconds=seconds)
        out.append((t, f"MRMS_{t:%Y%m%d-%H%M%S}.grib2.gz"))
    return out


class StepSelectionTests(unittest.TestCase):
    start = datetime.datetime(2026, 9, 3, 20, 0, tzinfo=UTC)

    def test_one_frame_per_clock_boundary_plus_the_newest(self):
        available = frames(self.start, 90)                      # 20:00:40 … 22:58:40
        chosen = select(available, 10, 120)
        times = [t for t, _ in chosen]
        self.assertEqual(times[-1], available[-1][0], "the newest frame is the last step")
        boundaries = [t for t in times[:-1]]
        for t in boundaries:
            self.assertEqual(t.second, 40)
            self.assertEqual(t.minute % 10, 0, f"{t:%H:%M:%S} is the frame nearest a clock boundary")
        self.assertEqual(len(chosen), 14)                       # 13 boundaries, 20:50 … 22:50, and 22:58:40

    def test_stamps_do_not_drift_when_the_newest_frame_moves_two_minutes(self):
        earlier = select(frames(self.start, 90), 10, 120)       # newest 22:58:40
        later = select(frames(self.start, 91), 10, 120)         # newest 23:00:40
        reused = {k for _, k in earlier} & {k for _, k in later}
        # Every boundary frame still inside the two-hour window is chosen
        # again. Two of the earlier run's frames go: its newest, 22:58:40,
        # because 23:00:40 now sits on the boundary; and 20:50:40, which the
        # window has moved past. One frame to draw, not fourteen.
        self.assertEqual(len(reused), len(earlier) - 2)
        self.assertEqual(len(later) - len(reused), 1, "one new frame to draw, not fourteen")

    def test_the_old_anchoring_drifted(self):
        # The fault, kept as a fact. Under the old rule the grid started at the
        # newest frame's own minute, so a :58 run and a :00 run named their
        # frames :x8 and :x0 and shared nothing — every cell redrawn, every run.
        def old_rule(stamped):
            newest = stamped[-1][0]
            tolerance = datetime.timedelta(minutes=5)
            chosen = set()
            for back in range(120, -1, -10):
                target = newest - datetime.timedelta(minutes=back)
                candidates = [(abs(t - target), t, k) for t, k in stamped if t <= target + tolerance]
                if candidates:
                    gap, _, k = min(candidates)
                    if gap <= tolerance:
                        chosen.add(k)
            return chosen
        a, b = old_rule(frames(self.start, 90)), old_rule(frames(self.start, 91))
        self.assertEqual(len(a), 13)
        self.assertEqual(a & b, set())

    def test_a_boundary_with_no_frame_within_half_a_step_is_left_out(self):
        available = frames(self.start, 90)
        # Remove every frame between 21:26 and 21:34: the 21:30 boundary has
        # nothing within five minutes.
        thinned = [(t, k) for t, k in available
                   if not (datetime.datetime(2026, 9, 3, 21, 26, tzinfo=UTC) <= t
                           <= datetime.datetime(2026, 9, 3, 21, 34, tzinfo=UTC))]
        chosen = select(thinned, 10, 120)
        minutes = {(t.hour, t.minute) for t, _ in chosen}
        self.assertNotIn((21, 30), minutes)
        self.assertIn((21, 20), minutes)
        self.assertIn((21, 40), minutes)

    def test_a_newest_frame_on_the_boundary_is_not_listed_twice(self):
        available = frames(self.start, 91)                      # newest 23:00:40, on the boundary
        chosen = select(available, 10, 120)
        keys = [k for _, k in chosen]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual(chosen[-1][0], available[-1][0])

    def test_a_newest_frame_just_past_a_boundary_is_not_a_second_step(self):
        available = frames(self.start, 92)                      # newest 23:02:40, 2 min past 23:00:40
        chosen = select(available, 10, 120)
        self.assertEqual(chosen[-1][0], available[-2][0], "23:00:40 stands for now; 23:02:40 adds a stutter")

    def test_nothing_in_gives_nothing_out(self):
        self.assertEqual(select([], 10, 120), [])


if __name__ == "__main__":
    unittest.main()
