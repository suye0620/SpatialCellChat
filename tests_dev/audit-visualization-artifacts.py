"""Inventory historical visualization captures referenced by migration summaries.

This validates files and local widget assets, not scientific or visual correctness.
"""
import argparse
import csv
import hashlib
from html.parser import HTMLParser
from pathlib import Path
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent.parent
FIGURES = ROOT / 'tests_dev' / 'figures'
OUTPUT = FIGURES / 'artifact-inventory-2026-11-10.csv'


class AssetParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.assets = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'script' and attrs.get('src'):
            self.assets.append(attrs['src'])
        if tag == 'link' and attrs.get('href'):
            self.assets.append(attrs['href'])


def rel(path):
    return path.relative_to(ROOT).as_posix()


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check-tracked', action='store_true')
    args = parser.parse_args()
    items = {}
    errors = []

    def add(path, source, status, kind):
        path = path.resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            errors.append(f'{source}: missing or external path {path}')
            return
        key = rel(path)
        if key not in items:
            items[key] = {'path': key, 'kind': kind, 'status': status, 'sources': set()}
        items[key]['sources'].add(source)

    for dataset in ('simulated', 'visium', 'slideseq'):
        summaries = [FIGURES / f'summary-{dataset}-P{phase}.csv' for phase in range(4)]
        summaries += [FIGURES / f'summary-{dataset}-P0-spatialC2CDistPlot-fresh.csv']
        for summary in summaries:
            if not summary.is_file():
                errors.append(f'missing summary: {rel(summary)}')
                continue
            for row in rows(summary):
                status = row.get('status', '')
                if status not in ('PASS', 'VISUAL-REVIEW'):
                    continue
                if row.get('file'):
                    paths = [ROOT / row['file']]
                else:
                    base = FIGURES / row['figure']
                    paths = sorted(base.parent.glob(base.name + '*.pdf')) if 'netVisual-export-' in base.name else [base.with_suffix('.pdf')]
                    if 'netVisual-export-' in base.name and len(paths) != 2:
                        errors.append(f'{rel(summary)}: expected two exports for {row["figure"]}, got {len(paths)}')
                for path in paths:
                    add(path, rel(summary), status, 'pdf')

    for summary in sorted(FIGURES.glob('summary-all-P2-*.csv')):
        for row in rows(summary):
            if row.get('status') == 'PASS':
                add((FIGURES / row['figure']).with_suffix('.pdf'), rel(summary), 'PASS', 'pdf')

    for dataset in ('simulated', 'visium', 'slideseq'):
        summaries = [FIGURES / 'P4' / dataset / 'summary.csv',
                     FIGURES / 'P3' / dataset / 'netVisual_CommunFlowWidget-summary.csv']
        for summary in summaries:
            for row in rows(summary):
                if row.get('status') != 'HTML-WRITTEN':
                    continue
                path = ROOT / row['file']
                add(path, rel(summary), 'HTML-WRITTEN', 'html')
                if not path.is_file():
                    continue
                html = AssetParser()
                html.feed(path.read_text(encoding='utf-8'))
                for url in html.assets:
                    parsed = urlsplit(url)
                    if parsed.scheme or parsed.netloc or url.startswith('#'):
                        continue
                    asset = path.parent / unquote(parsed.path)
                    add(asset, rel(path), 'LOCAL-ASSET', 'asset')
                    if asset.name in ('plotly-latest.min.js', 'typedarray.min.js'):
                        add(asset.parent / 'LICENSE', rel(asset), 'LICENSE', 'license')

    diagnostics = [
        FIGURES / 'P3' / 'slideseq' / 'display-diagnostics-2026-11-11' / 'summary.csv',
        FIGURES / 'P3' / 'streamline-cancellation-2026-11-11' / 'summary.csv',
        FIGURES / 'P3' / 'streamline-patch-experiment-2026-11-11' / 'summary.csv',
        FIGURES / 'P2' / 'visium' / 'chord-label-experiment-2026-11-11' / 'summary.csv',
        FIGURES / 'P2' / 'visium' / 'chord-external-prototype-2026-11-11' / 'summary.csv',
    ]
    for summary in diagnostics:
        if not summary.is_file():
            errors.append(f'missing diagnostic summary: {rel(summary)}')
            continue
        for row in rows(summary):
            add(ROOT / row['file'], rel(summary), 'DIAGNOSTIC', 'pdf')

    # Producer-derived flow comparisons use the same field for both stats.
    production = FIGURES / 'P3' / 'streamline-production-comparison-2026-11-11'
    for dataset in ('simulated', 'visium', 'slideseq'):
        summary = production / f'{dataset}-summary.csv'
        if not summary.is_file():
            errors.append(f'missing production comparison summary: {rel(summary)}')
            continue
        for row in rows(summary):
            tag = 'default' if float(row['res']) == 0.5 else 'res8'
            for kind in ('original', 'adapted'):
                add(production / f'{dataset}-{tag}-{kind}.pdf', rel(summary), 'DIAGNOSTIC', 'pdf')

    # Explicit arrow scale comparisons on producer-derived fields.
    arrow_scale = FIGURES / 'P3' / 'arrow-display-scale-2026-11-11'
    for dataset in ('simulated', 'visium', 'slideseq'):
        summary = arrow_scale / f'{dataset}-summary.csv'
        if not summary.is_file():
            errors.append(f'missing arrow scale comparison summary: {rel(summary)}')
            continue
        for scale in (1, 4):
            add(arrow_scale / f'{dataset}-scale-{scale}.pdf', rel(summary), 'DIAGNOSTIC', 'pdf')

    # Browser screenshots and explicitly documented visual-review previews are
    # included as supporting evidence, separate from PASS or visual acceptance.
    for phase in ('P0', 'P2', 'P3', 'P4'):
        for path in (FIGURES / phase).rglob('*.png'):
            if 'preview' in path.name.lower() or 'browser' in path.name.lower():
                add(path, 'visual-review capture', 'REVIEW', 'image')

    tracked = {}
    if args.check_tracked:
        index = subprocess.check_output(['git', 'ls-files', '--stage', '-z'], cwd=ROOT).decode('utf-8')
        for entry in index.split('\0'):
            if '\t' in entry:
                info, name = entry.split('\t', 1)
                mode, oid, stage = info.split()
                if stage == '0':
                    tracked[name] = oid
    result = []
    for key, item in sorted(items.items()):
        path = ROOT / key
        if item['kind'] == 'pdf':
            with path.open('rb') as f:
                start = f.read(5)
                f.seek(max(0, path.stat().st_size - 2048))
                end = f.read()
            if start != b'%PDF-' or b'%%EOF' not in end:
                errors.append(f'invalid PDF header or trailer: {key}')
        content = path.read_bytes()
        if args.check_tracked:
            if key not in tracked:
                errors.append(f'not tracked or staged: {key}')
            else:
                blob = b'blob ' + str(len(content)).encode('ascii') + b'\0' + content
                if hashlib.sha1(blob).hexdigest() != tracked[key]:
                    errors.append(f'index content differs from inventory bytes: {key}')
        digest = hashlib.sha256(content).hexdigest()
        result.append({'path': key, 'kind': item['kind'], 'status': item['status'],
                       'size_bytes': path.stat().st_size, 'sha256': digest,
                       'sources': ' | '.join(sorted(item['sources']))})
    with OUTPUT.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['path', 'kind', 'status', 'size_bytes', 'sha256', 'sources'])
        writer.writeheader()
        writer.writerows(result)
    from collections import Counter
    print('Artifacts:', len(result), dict(Counter(r['kind'] for r in result)))
    print('Inventory:', rel(OUTPUT))
    for error in errors:
        print('ERROR:', error, file=sys.stderr)
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
