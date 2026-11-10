"""Check strict intersections in paths drawn by the production outside helper."""
import csv
import importlib.util
from collections import defaultdict
from pathlib import Path

spec = importlib.util.spec_from_file_location('leader_crossings', Path(__file__).with_name('audit-chord-leader-crossings.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
base = Path(__file__).resolve().parent / 'figures/P2/chord-production-outside-2026-11-11'
results = []
for dataset in ('simulated', 'visium', 'slideseq'):
    for case in ('outside-large', 'outside-small', 'outside-small-legend-overlap'):
        path = base / dataset / (case + '-paths.csv')
        lines = defaultdict(list)
        with path.open(encoding='utf-8-sig', newline='') as f:
            for row in csv.DictReader(f):
                vertices = lines[row['sector']]
                if int(row['vertex']) != len(vertices) + 1:
                    raise ValueError(f'Out-of-order path vertices: {path}')
                vertices.append((float(row['x']), float(row['y'])))
        if not lines or any(len(vertices) != 34 for vertices in lines.values()):
            raise ValueError(f'Unexpected path count or length: {path}')
        n, pairs = module.count(lines)
        results.append({'dataset': dataset, 'case': case, 'sectors': len(lines),
                        'strict_crossings': n, 'pairs_with_crossings': pairs})
        print(results[-1])
with (base / 'crossing-summary.csv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=results[0]); writer.writeheader(); writer.writerows(results)
if any(row['strict_crossings'] for row in results):
    raise SystemExit('At least one production outside layout has strict crossings')
