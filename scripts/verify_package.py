#!/usr/bin/env python3
"""Check packaging and recorded numerical evidence without extra packages."""
from pathlib import Path
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT / 'provenance/source-manifest.json').read_text())
checked = 0
for record in manifest['files']:
    file = ROOT / record['path']
    assert file.is_file(), f'Missing copied source: {file}'
    if file.suffix == '.m':
        assert hashlib.sha256(file.read_bytes()).hexdigest() == record['sha256'], file
        checked += 1

for file in ROOT.rglob('*.md'):
    if any(part in {'output', '.git', '.local'} for part in file.parts):
        continue
    for target in re.findall(r'\]\(([^)]+)\)', file.read_text()):
        target = target.split('#')[0]
        if not target or '://' in target or target.startswith('mailto:'):
            continue
        assert (file.parent / target).exists(), f'Broken link in {file}: {target}'

for stem in ['starfish-overview', 'starfish-hierarchy', 'starfish-validation']:
    for suffix in ['.png', '.pdf']:
        file = ROOT / 'figures' / (stem + suffix)
        assert 1000 < file.stat().st_size < 10_000_000, f'Unexpected figure size: {file}'
        if suffix == '.pdf':
            assert file.read_bytes().startswith(b'%PDF-'), file

data = json.loads((ROOT / 'verification/starfish-full.json').read_text())
assert data['levels'] == 5 and data['finestEquivalent'] == 1024
assert 0 < data['visibleCells'] < data['uniformCells']
assert data['uncoveredBandCells'] == 0
assert data['relativeResidual'] < 1e-10
assert max(data['uniformRelativeResiduals']) < 1e-10
assert data['relativeL2VsUniform1024'] < data['acceptanceRelativeL2']
assert isinstance(data['linfVsUniform1024'], (int, float))
assert math.isfinite(data['linfVsUniform1024'])
assert data['balanceDefect'] < 1e-10
suite = (ROOT / 'verification/suite-output.txt').read_text()
assert 'ALL BSAM PORTFOLIO VERIFICATION SUITES PASSED' in suite
assert 'README EXAMPLES PASSED' in (ROOT / 'verification/readme-examples.txt').read_text()
print(f'PASS: {checked} copied MATLAB files unchanged; links, figures, and numerical evidence checked.')
