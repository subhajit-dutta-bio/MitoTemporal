"""Create local folders and fixed analysis-configuration files, never observations.

Run from a fresh clone before obtaining external inputs. Does not download data.
"""
from pathlib import Path
import json
from gene_sets import CANONICAL_MODULES, CONTROL_GENE_SETS, COMPARATOR_GENE_SETS, LEGACY_MODULES

ROOT = Path(__file__).resolve().parent

def main():
    folders = ['replication', 'source/hummel2', 'revision2', 'revision3', 'revision4',
               'revision6', 'revision7', 'revision8', 'revision9', 'revision10',
               'revision11_polished', 'revision12', 'revision13', 'revision14']
    for folder in folders:
        for sub in ['analysis', 'metadata', 'raw', 'deliverables/figures', 'deliverables/source_data']:
            (ROOT / folder / sub).mkdir(parents=True, exist_ok=True)
    for name, values in {
        'revision6/analysis/canonical_modules.json': CANONICAL_MODULES,
        'revision4/analysis/control_gene_sets.json': CONTROL_GENE_SETS,
        'revision10/fixed_comparator_gene_sets.json': COMPARATOR_GENE_SETS,
        'source/hummel2/gene_sets/mito_modules.json': LEGACY_MODULES,
        'replication/scripts/mito_modules.json': LEGACY_MODULES,
    }.items():
        path = ROOT / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            if json.loads(path.read_text(encoding='utf-8')) != values:
                raise RuntimeError(f'Existing configuration differs: {name}; preserve it and resolve manually.')
        else:
            path.write_text(json.dumps(values, indent=2) + '\n', encoding='utf-8')
    print('Local folders and fixed gene-set configuration ready. No datasets downloaded.')

if __name__ == '__main__':
    main()
