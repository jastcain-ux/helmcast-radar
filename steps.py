"""Which MRMS frames make the measured timeline — and why their names must not drift.

Pure Python, no numpy, so it is testable on any machine and importable by the
renderer without its heavy dependencies.

The renderer reuses any frame whose PNG already exists from the last run and
draws only the missing ones. That reuse silently failed for most runs: the
ten-minute step grid was anchored on the newest MRMS frame's own minute, and
MRMS publishes every two minutes while runs land at irregular intervals, so
the grid — and every frame's file stamp with it — drifted between runs. A run
found no file under the new names and redrew all thirteen steps of every cell
(832 frames, 522–689 s), which is what kept the measured half arriving at the
app close to the 30-minute limit (HelmCast B-20, diagnosis of 2026-09-03).

Anchoring the grid on the clock — :00, :10, :20 — makes a step's frame, and so
its stamp, the same in every run once the frames around that boundary exist.
Only the newest step is ever new.
"""
import datetime


def select(stamped, step_minutes, history_minutes):
    """The timeline's frames, oldest first, as (time, key) pairs.

    `stamped` is [(observed time, key)] for every frame available. One frame per
    ten-minute clock boundary back through the history — the frame nearest that
    boundary within half a step, or nothing when none is that close: a gap the
    app can state beats a frame under the wrong label. Then the newest frame
    itself, so the timeline is as fresh as MRMS is, unless it is within half a
    step of the last boundary's frame, where a second near-identical picture
    would only make the loop stutter and give the projection a two-minute
    baseline to measure motion from.
    """
    if not stamped:
        return []
    stamped = sorted(stamped)
    newest_time, newest_key = stamped[-1]
    step = datetime.timedelta(minutes=step_minutes)
    tolerance = step / 2
    anchor = newest_time.replace(second=0, microsecond=0)
    anchor -= datetime.timedelta(minutes=anchor.minute % step_minutes)
    out, seen = [], set()
    for back in range(history_minutes, -1, -step_minutes):
        target = anchor - datetime.timedelta(minutes=back)
        candidates = [(abs(t - target), t, k) for t, k in stamped if abs(t - target) <= tolerance]
        if not candidates:
            continue
        _, t, k = min(candidates)
        if k not in seen:
            seen.add(k)
            out.append((t, k))
    if newest_key not in seen and (not out or newest_time - out[-1][0] > tolerance):
        out.append((newest_time, newest_key))
    return out
