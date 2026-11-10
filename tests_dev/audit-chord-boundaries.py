"""Count proper crossings for the saved boundary-study polylines."""
import csv
from collections import defaultdict
from pathlib import Path
import importlib.util
spec = importlib.util.spec_from_file_location('leader_crossings', Path(__file__).with_name('audit-chord-leader-crossings.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
count = module.count

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / 'tests_dev/figures/P2/chord-external-boundaries-2026-11-11'
result = []
for dataset in ('simulated', 'visium', 'slideseq'):
    with (BASE / dataset / 'summary.csv').open(encoding='utf-8-sig', newline='') as f:
        summary = list(csv.DictReader(f))
    for row in summary:
        if row['status'] != 'DRAWN-REVIEW':
            raise RuntimeError(f"Case failed: {dataset}, {row['case']}: {row['error']}")
        case = BASE / dataset / row['case']
        paths = defaultdict(list)
        with (case / 'paths.csv').open(encoding='utf-8-sig', newline='') as f:
            for point in csv.DictReader(f):
                path = paths[point['sector']]
                if int(point['vertex']) != len(path) + 1:
                    raise ValueError('Inconsistent path vertices')
                path.append((float(point['x']), float(point['y'])))
        if len(paths) != int(row['sectors']):
            raise ValueError('Path and sector counts do not agree')
        crossings, pairs = count(paths)
        result.append({'dataset': dataset, 'case': row['case'], 'sectors': len(paths),
                       'strict_crossings': crossings, 'pairs_with_crossings': pairs})
        print(result[-1], flush=True)
with (BASE / 'crossing-summary.csv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=result[0]); writer.writeheader(); writer.writerows(result)
if any(r['strict_crossings'] for r in result):
    raise SystemExit('At least one case has strict crossings')
