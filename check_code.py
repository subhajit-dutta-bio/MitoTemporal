"""Syntax and data-free synthetic checks of the published numerical helpers.

AST extraction avoids running data-dependent top-level analysis statements.
All numerical arrays below are synthetic test fixtures, not study observations.
"""
from pathlib import Path
import ast
import unittest
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from gene_sets import CANONICAL_MODULES, CONTROL_GENE_SETS, COMPARATOR_GENE_SETS

ROOT = Path(__file__).resolve().parent

def functions(path, names, **namespace):
    tree = ast.parse((ROOT / path).read_text(encoding='utf-8-sig'))
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    if {n.name for n in nodes} != set(names):
        raise AssertionError(f'Missing helper in {path}')
    scope = dict(np=np, pd=pd, stats=stats, sm=sm, **namespace)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), scope)
    return scope

class CodeChecks(unittest.TestCase):
    def test_all_python_sources_parse(self):
        paths = [p for p in ROOT.rglob('*.py') if '.venv' not in p.parts and 'python_packages' not in p.parts]
        for path in paths:
            with self.subTest(path=str(path.relative_to(ROOT))):
                ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path))

    def test_canonical_import_and_controls(self):
        genes = CANONICAL_MODULES['Import_TIM_TOM']
        self.assertEqual(len(genes), 24)
        self.assertEqual(len(set(genes)), 24)
        for name in ['ER_programme', 'Proliferation']:
            self.assertFalse(set(genes) & set(CONTROL_GENE_SETS[name]))
        self.assertFalse(set(genes) & set(COMPARATOR_GENE_SETS['Nonimport_mito']))

    def test_bh_matches_independent_library(self):
        p = np.array([.7, .004, .09, .03, 1., .004])
        expected = multipletests(p, method='fdr_bh')[1]
        for path in ['revision7/common.py', 'revision8/common.py', 'revision9/common.py',
                     'revision10/common.py', 'revision14/v14_common.py']:
            with self.subTest(path=path):
                np.testing.assert_allclose(functions(path, ['bh'])['bh'](p), expected)

    def test_extension_bh_preserves_untested_slots(self):
        p = np.array([.04, np.nan, .001, np.inf])
        result = functions('revision14/v14_common.py', ['bh'])['bh'](p)
        np.testing.assert_allclose(result[[0, 2]], multipletests(p[[0, 2]], method='fdr_bh')[1])
        self.assertTrue(np.isnan(result[[1, 3]]).all())

    def test_source_score_coverage_gate(self):
        score = functions('revision8/common.py', ['scores'], sets={})['scores']
        a = pd.DataFrame({'complete': [1., 2., 3., 4.], 'missing': [1., np.nan, np.nan, 4.]},
                         index=['A', 'B', 'C', 'D'])
        result = score(a, {'module': ['A', 'B', 'C', 'D'], 'single': ['A']})
        self.assertAlmostEqual(result.loc['complete', 'module'], 2.5)
        self.assertTrue(np.isnan(result.loc['missing', 'module']))
        self.assertEqual(result.loc['missing', 'single'], 1.)

    def test_extension_score_measured_panel_gate(self):
        score = functions('revision14/v14_common.py', ['module_score'])['module_score']
        a = pd.DataFrame({'complete': np.arange(6.), 'undercovered': [1., 2., 3., 4., np.nan, np.nan]},
                         index=list('ABCDEF'))
        result = score(a, list('ABCDEF'))
        self.assertAlmostEqual(result['complete'], 2.5)
        self.assertTrue(np.isnan(result['undercovered']))

    def test_contrast_helpers_against_scipy(self):
        ns = functions('revision10/analysis_utils.py', ['welch', 'onesample'])
        x = np.array([.2, .6, .1, .5, .9]); y = np.array([.1, .3, .2, .4])
        w = ns['welch'](x, y)
        self.assertAlmostEqual(w['p'], stats.ttest_ind(x, y, equal_var=False).pvalue)
        self.assertAlmostEqual(w['estimate'], x.mean()-y.mean())
        t = ns['onesample'](x)
        self.assertAlmostEqual(t['p'], stats.ttest_1samp(x, 0).pvalue)
        self.assertTrue(t['low'] < t['estimate'] < t['high'])

    def test_regression_eligibility_and_rank(self):
        fit = functions('revision14/v14_common.py', ['fit_ols'])['fit_ols']
        small = pd.DataFrame({'y': [1., 2.], 'x': [0., 1.]})
        self.assertEqual(fit(small, 'y', 'x')['status'], 'unavailable')
        collinear = pd.DataFrame({'y': np.arange(20.)**2, 'x': np.arange(20.), 'z': np.arange(20.)})
        result = fit(collinear, 'y', 'x', ['z'])
        self.assertEqual(result['reason'], 'rank-deficient fixed model')

if __name__ == '__main__':
    unittest.main(verbosity=2)
