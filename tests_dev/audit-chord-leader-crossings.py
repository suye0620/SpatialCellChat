"""Measure strict segment crossings in both recorded chord-label prototypes.

Only proper crossings are counted; touching, near-coincident strokes, width, text,
and visual readability require separate checks.
"""
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / 'tests_dev/figures/P2/visium'
OLD = FIG / 'chord-external-prototype-2026-11-11'
NEW = FIG / 'chord-external-unfold-2026-11-11'


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def old_paths():
    paths = {}
    for row in rows(OLD / 'positions.csv'):
        x, y, r, label_y = (float(row[k]) for k in ('x', 'y', 'radius', 'label.y'))
        side = 1 if row['side'] == 'right' else -1
        paths[row['sector']] = [(x*r, y*r), (x*1.01, y*1.01),
                                (side*1.12, label_y), (side*1.2, label_y)]
    return paths


def new_paths():
    paths = defaultdict(list)
    for row in rows(NEW / 'paths.csv'):
        path = paths[row['sector']]
        assert int(row['vertex']) == len(path) + 1
        path.append((float(row['x']), float(row['y'])))
    positions = {r['sector']: r for r in rows(NEW / 'positions.csv')}
    assert set(paths) == set(positions) and len(paths) == 111
    for sector, path in paths.items():
        row = positions[sector]
        x, y, radius, label_y = (float(row[k]) for k in ('x', 'y', 'radius', 'label.y'))
        assert len(path) == 34
        assert abs(path[0][0]-x*radius) < 1e-12 and abs(path[0][1]-y*radius) < 1e-12
        assert abs(path[-1][1]-label_y) < 1e-12
    return paths


def orientation(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


def crosses(a, b, c, d):
    if (max(a[0], b[0]) < min(c[0], d[0]) or max(c[0], d[0]) < min(a[0], b[0]) or
            max(a[1], b[1]) < min(c[1], d[1]) or max(c[1], d[1]) < min(a[1], b[1])):
        return False
    return orientation(a, b, c)*orientation(a, b, d) < 0 and orientation(c, d, a)*orientation(c, d, b) < 0


def count(paths):
    entries = list(paths.items())
    intersections = 0
    pairs = 0
    for i, (_, left) in enumerate(entries):
        for _, right in entries[i+1:]:
            hits = sum(crosses(a, b, c, d)
                       for a, b in zip(left, left[1:])
                       for c, d in zip(right, right[1:]))
            intersections += hits
            pairs += bool(hits)
    return intersections, pairs


def main():
    old = old_paths()
    new = new_paths()
    assert set(old) == set(new)
    old_count, old_pairs = count(old)
    new_count, new_pairs = count(new)
    result = [{'layout': 'straight', 'strict_crossings': old_count, 'pairs_with_crossings': old_pairs},
              {'layout': 'unfolded', 'strict_crossings': new_count, 'pairs_with_crossings': new_pairs}]
    with (NEW / 'crossing-summary.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=result[0])
        writer.writeheader()
        writer.writerows(result)
    print(result)
    if old_count == 0 or new_count != 0:
        raise SystemExit('Crossing comparison did not satisfy its focused checks')


if __name__ == '__main__':
    main()
