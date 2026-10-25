"""Re-render inherited context with its new supplementary numbering."""
from pathlib import Path
import importlib.util,sys
r=Path(__file__).resolve().parent;root=r.parent
spec=importlib.util.spec_from_file_location('legacy_common',root/'revision7/common.py');legacy=importlib.util.module_from_spec(spec);sys.modules['legacy_common']=legacy;spec.loader.exec_module(legacy)
source=(root/'revision7/build_figures.py').read_text()
source=source.replace('from common import *','from legacy_common import *\nr=Path(__file__).resolve().parent\ndef legacy_name(n):return {3:"Supplementary_Figure_S8",4:"Supplementary_Figure_S9"}.get(n,f"Figure_{n}")')
source=source.replace("f'FIGURE {n}'","legacy_name(n).replace('_',' ').upper()")
source=source.replace("f'Figure_{n}_{letter}.csv.gz'","f'{legacy_name(n)}_{letter}.csv.gz'")
source=source.replace("figure=f'Figure_{n}'","figure=legacy_name(n)")
source=source.replace("f'Figure_{n}.{ext}'","f'{legacy_name(n)}.{ext}'")
source=source.replace('main_panel_source_manifest.csv','inherited_context_panel_manifest.csv')
source=source.replace('Created four main figures, eight panels each','Updated inherited plots and supplementary numbering')
exec(compile(source,'revision7/build_figures.py preserved plot logic','exec'),dict(__file__=__file__))
