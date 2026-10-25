#!/usr/bin/env python3
"""
02_qc_pipeline.py — Quality Control & Mitochondrial Scoring
============================================================
Reads raw count matrices, applies:
  1. Adaptive QC (NO hard MT% cutoff — key for this project)
  2. Doublet detection (scrublet)
  3. Per-cell mitochondrial module scoring
  4. QC figures
Outputs: processed/{dataset}_qc.h5ad
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from scipy.io import mmread
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
sc.settings.verbosity = 2

# ── Paths ──
BASE = Path(os.environ.get("MITOTEMPORAL_SINGLECELL_ROOT", str(Path(__file__).resolve().parents[1])))
RAW = BASE / "raw"
PROCESSED = BASE / "processed"
GENE_SETS = BASE / "gene_sets"
FIGURES = BASE / "figures"

for d in [PROCESSED, FIGURES]:
    d.mkdir(parents=True, exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# ============================================================================
# LOAD MITOCHONDRIAL GENE MODULES
# ============================================================================
def load_mito_modules():
    """Load curated mitochondrial gene modules."""
    modules_file = GENE_SETS / "mito_modules.json"
    with open(modules_file) as f:
        modules = json.load(f)
    log(f"Loaded {len(modules)} mito modules ({sum(len(v) for v in modules.values())} genes)")
    return modules


# ============================================================================
# DATASET LOADING FUNCTIONS
# ============================================================================

def load_gse167036():
    """Load GSE167036 — Liu et al. HR+ BC longitudinal NAC."""
    log("Loading GSE167036 (mtx + csv format)...")
    d = RAW / "GSE167036"

    # Read sparse matrix
    mtx = mmread(str(d / "GSE167036_cell_counts_matrix.mtx.gz"))
    mtx = csr_matrix(mtx.T)  # genes × cells → cells × genes

    # Read barcodes and features
    barcodes = pd.read_csv(d / "GSE167036_barcodes.csv.gz", index_col=0)
    features = pd.read_csv(d / "GSE167036_features.csv.gz", index_col=0)

    # Build AnnData — matrix is genes(25301) x cells(118845), need NO transpose
    # because mmread gives (rows=25301, cols=118845) and we need cells x genes
    adata = ad.AnnData(X=mtx)
    adata.obs_names = barcodes["cell_id"].values.astype(str)
    adata.var_names = features["x"].values.astype(str)
    adata.var_names_make_unique()

    # Load metadata — set cell_id as index for alignment
    meta = pd.read_csv(d / "GSE167036_meta_all.csv.gz", index_col=0)
    meta = meta.set_index("cell_id")
    common = adata.obs_names.intersection(meta.index)
    log(f"  Barcode-metadata overlap: {len(common)}/{adata.shape[0]}")
    adata = adata[common].copy()
    meta_aligned = meta.loc[common]
    for col in meta_aligned.columns:
        adata.obs[col] = meta_aligned[col].values

    adata.obs["dataset"] = "GSE167036"
    log(f"  GSE167036: {adata.shape[0]} cells × {adata.shape[1]} genes")
    log(f"  Obs columns: {list(adata.obs.columns[:10])}")
    return adata


def load_gse169246():
    """Load GSE169246 — Zhang et al. TNBC (RNA only)."""
    log("Loading GSE169246 (10x mtx format)...")
    d = RAW / "GSE169246"

    mtx = mmread(str(d / "GSE169246_TNBC_RNA.counts.mtx.gz"))
    mtx = csr_matrix(mtx.T)  # transpose if genes × cells

    barcodes = pd.read_csv(d / "GSE169246_TNBC_RNA.barcode.tsv.gz",
                           header=None, sep="\t")
    features = pd.read_csv(d / "GSE169246_TNBC_RNA.feature.tsv.gz",
                           header=None, sep="\t")

    adata = ad.AnnData(X=mtx)
    adata.obs_names = barcodes[0].values.astype(str)

    # Features may have gene_id and gene_name columns
    if features.shape[1] >= 2:
        adata.var["gene_id"] = features[0].values.astype(str)
        adata.var_names = features[1].values.astype(str)
    else:
        adata.var_names = features[0].values.astype(str)
    adata.var_names_make_unique()

    adata.obs["dataset"] = "GSE169246"
    log(f"  GSE169246: {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata


def load_gse176078():
    """Load GSE176078 — Wu et al. Breast Cancer Atlas."""
    log("Loading GSE176078 (sparse mtx + barcodes + genes)...")
    d = RAW / "GSE176078" / "Wu_etal_2021_BRCA_scRNASeq"

    mtx_file = d / "count_matrix_sparse.mtx"
    barcodes_file = d / "count_matrix_barcodes.tsv"
    genes_file = d / "count_matrix_genes.tsv"
    meta_file = d / "metadata.csv"

    # Read sparse matrix
    log("  Reading sparse matrix (this may take a few minutes for 2.3 GB)...")
    mtx = mmread(str(mtx_file))
    mtx = csr_matrix(mtx.T)  # genes x cells -> cells x genes
    log(f"  Matrix shape after transpose: {mtx.shape}")

    # Read barcodes and genes
    barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")
    genes = pd.read_csv(genes_file, header=None, sep="\t")

    # Build AnnData
    adata = ad.AnnData(X=mtx)
    adata.obs_names = barcodes[0].values.astype(str)

    # Genes file may have gene_id and gene_name
    if genes.shape[1] >= 2:
        adata.var["gene_id"] = genes[0].values.astype(str)
        adata.var_names = genes[1].values.astype(str)
    else:
        adata.var_names = genes[0].values.astype(str)
    adata.var_names_make_unique()

    # Load metadata
    log("  Loading metadata...")
    meta = pd.read_csv(meta_file, index_col=0)
    common = adata.obs_names.intersection(meta.index)
    if len(common) > 0:
        adata = adata[common].copy()
        meta_aligned = meta.loc[common]
        for col in meta_aligned.columns:
            adata.obs[col] = meta_aligned[col].values
        log(f"  Metadata attached: {len(meta_aligned.columns)} columns")
    else:
        log("  ⚠️ No barcode overlap with metadata — checking format...")
        # Try without index alignment
        if len(meta) == adata.shape[0]:
            for col in meta.columns:
                adata.obs[col] = meta[col].values
            log(f"  Metadata attached by position: {len(meta.columns)} columns")

    adata.obs["dataset"] = "GSE176078"
    log(f"  GSE176078: {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata


# ============================================================================
# QC PIPELINE
# ============================================================================

def run_qc(adata, dataset_name, mito_modules):
    """
    Adaptive QC pipeline.
    CRITICAL: No hard MT% cutoff. MT% is kept as a continuous covariate.
    """
    log(f"\n{'='*60}")
    log(f"QC: {dataset_name}")
    log(f"{'='*60}")
    log(f"Input: {adata.shape[0]} cells × {adata.shape[1]} genes")

    # ── 1. Basic gene/cell filtering ──
    log("Step 1: Basic filtering...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    log(f"  After min_genes=200, min_cells=3: {adata.shape}")

    # ── 2. QC metrics (MT%, ribo%, hemoglobin) ──
    log("Step 2: QC metrics...")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.match("^RP[SL]\\d")
    adata.var["hb"] = adata.var_names.str.match("^HB[ABDEGM]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=False, inplace=True
    )

    log(f"  MT% range:   {adata.obs['pct_counts_mt'].min():.1f} – {adata.obs['pct_counts_mt'].max():.1f}")
    log(f"  MT% median:  {adata.obs['pct_counts_mt'].median():.1f}")
    log(f"  Genes/cell:  {adata.obs['n_genes_by_counts'].median():.0f} median")
    log(f"  UMIs/cell:   {adata.obs['total_counts'].median():.0f} median")

    # ── 3. Adaptive filtering (MAD-based, NO hard MT% cutoff) ──
    log("Step 3: Adaptive MAD-based filtering...")

    def mad_outlier(series, nmads=5):
        """Flag outliers using median absolute deviation."""
        median = series.median()
        mad = np.median(np.abs(series - median))
        if mad == 0:
            mad = series.std()
        lower = median - nmads * 1.4826 * mad
        upper = median + nmads * 1.4826 * mad
        return lower, upper

    # n_genes: remove extremely low/high
    g_lo, g_hi = mad_outlier(np.log1p(adata.obs["n_genes_by_counts"]), nmads=5)
    gene_mask = (
        (np.log1p(adata.obs["n_genes_by_counts"]) > g_lo) &
        (np.log1p(adata.obs["n_genes_by_counts"]) < g_hi)
    )

    # total_counts: remove extremely low/high
    c_lo, c_hi = mad_outlier(np.log1p(adata.obs["total_counts"]), nmads=5)
    count_mask = (
        (np.log1p(adata.obs["total_counts"]) > c_lo) &
        (np.log1p(adata.obs["total_counts"]) < c_hi)
    )

    # MT%: ADAPTIVE upper bound only — use generous 5 MAD
    # This is INTENTIONALLY permissive to retain metabolically active cancer cells
    mt_median = adata.obs["pct_counts_mt"].median()
    mt_mad = np.median(np.abs(adata.obs["pct_counts_mt"] - mt_median))
    mt_upper = mt_median + 5 * 1.4826 * mt_mad
    mt_upper = max(mt_upper, 50)  # Never filter below 50% — cancer cells can be >40%
    mt_mask = adata.obs["pct_counts_mt"] < mt_upper

    log(f"  Gene filter (5 MAD): {gene_mask.sum()}/{len(gene_mask)} pass")
    log(f"  Count filter (5 MAD): {count_mask.sum()}/{len(count_mask)} pass")
    log(f"  MT% filter (upper={mt_upper:.1f}%): {mt_mask.sum()}/{len(mt_mask)} pass")

    combined_mask = gene_mask & count_mask & mt_mask
    n_removed = (~combined_mask).sum()
    adata = adata[combined_mask].copy()
    log(f"  Removed {n_removed} cells → {adata.shape[0]} remaining")

    # ── 4. Doublet detection (scrublet) ──
    log("Step 4: Doublet detection (scrublet)...")
    try:
        sc.pp.scrublet(adata, batch_key="dataset" if "dataset" in adata.obs else None)
        n_doublets = adata.obs["predicted_doublet"].sum()
        log(f"  Detected {n_doublets} doublets ({100*n_doublets/adata.shape[0]:.1f}%)")
        # Flag but do NOT remove — keep for downstream inspection
        adata.obs["qc_doublet"] = adata.obs["predicted_doublet"].astype(bool)
    except Exception as e:
        log(f"  ⚠️ Scrublet failed: {e}. Skipping doublet detection.")
        adata.obs["qc_doublet"] = False

    # ── 5. Store raw counts before normalization ──
    log("Step 5: Storing raw counts layer...")
    adata.layers["raw_counts"] = adata.X.copy()

    # ── 6. Normalize + log1p ──
    log("Step 6: Normalize + log1p...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ── 7. Score mitochondrial modules ──
    log("Step 7: Scoring mitochondrial modules...")
    for module_name, genes in mito_modules.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) >= 2:
            score_key = f"mito_{module_name}"
            try:
                sc.tl.score_genes(adata, gene_list=present, score_name=score_key)
                log(f"  ✅ {module_name}: {len(present)}/{len(genes)} genes, "
                    f"range=[{adata.obs[score_key].min():.3f}, {adata.obs[score_key].max():.3f}]")
            except Exception as e:
                log(f"  ⚠️ {module_name}: scoring failed ({e})")
        else:
            log(f"  ⚠️ {module_name}: only {len(present)}/{len(genes)} genes found, skipping")

    # ── 8. HVG selection ──
    log("Step 8: HVG selection...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3",
                                 layer="raw_counts")
    log(f"  {adata.var['highly_variable'].sum()} HVGs selected")

    # ── 9. PCA + neighbors + UMAP (for QC visualization) ──
    log("Step 9: PCA + UMAP...")
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)

    # ── 10. QC Figures ──
    log("Step 10: Generating QC figures...")
    try:
        fig_dir = FIGURES / dataset_name
        fig_dir.mkdir(exist_ok=True)

        # Violin plots
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))
        for ax, metric in zip(axes, ["n_genes_by_counts", "total_counts",
                                       "pct_counts_mt", "pct_counts_ribo"]):
            ax.hist(adata.obs[metric], bins=50, alpha=0.7, edgecolor="black", linewidth=0.3)
            ax.set_title(metric)
            ax.axvline(adata.obs[metric].median(), color="red", ls="--", lw=1)
        plt.tight_layout()
        plt.savefig(fig_dir / "qc_distributions.png", dpi=150, bbox_inches="tight")
        plt.close()

        # UMAP colored by MT%
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        sc.pl.umap(adata, color="pct_counts_mt", ax=axes[0], show=False,
                    title="MT%", frameon=False)
        sc.pl.umap(adata, color="n_genes_by_counts", ax=axes[1], show=False,
                    title="Genes/cell", frameon=False)
        if "qc_doublet" in adata.obs:
            sc.pl.umap(adata, color="qc_doublet", ax=axes[2], show=False,
                        title="Doublets", frameon=False)
        plt.tight_layout()
        plt.savefig(fig_dir / "umap_qc.png", dpi=150, bbox_inches="tight")
        plt.close()

        # Mito module score UMAPs
        mito_cols = [c for c in adata.obs.columns if c.startswith("mito_")]
        if mito_cols:
            n_cols = min(len(mito_cols), 6)
            fig, axes = plt.subplots(1, n_cols, figsize=(5*n_cols, 4))
            if n_cols == 1:
                axes = [axes]
            for ax, col in zip(axes, mito_cols[:n_cols]):
                sc.pl.umap(adata, color=col, ax=ax, show=False,
                            title=col.replace("mito_", ""), frameon=False)
            plt.tight_layout()
            plt.savefig(fig_dir / "umap_mito_modules.png", dpi=150, bbox_inches="tight")
            plt.close()

        log(f"  Figures saved to: {fig_dir}")
    except Exception as e:
        log(f"  ⚠️ Figure generation failed: {e}")

    # ── Summary ──
    log(f"\nQC Summary for {dataset_name}:")
    log(f"  Final: {adata.shape[0]} cells × {adata.shape[1]} genes")
    log(f"  MT% median: {adata.obs['pct_counts_mt'].median():.1f}%")
    log(f"  Doublets flagged: {adata.obs['qc_doublet'].sum()}")
    log(f"  Mito modules scored: {len([c for c in adata.obs.columns if c.startswith('mito_')])}")

    return adata


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    log("MitoTemporal — QC Pipeline")
    log(f"Scanpy: {sc.__version__}, AnnData: {ad.__version__}")
    log(f"Threads: {sc.settings.n_jobs}")
    log("")

    mito_modules = load_mito_modules()

    # ── Process each dataset ──
    datasets = {
        "GSE167036": load_gse167036,
        "GSE169246": load_gse169246,
        "GSE176078": load_gse176078,
    }

    qc_summary = {}

    for name, loader in datasets.items():
        try:
            adata = loader()
            if adata is None:
                log(f"⚠️ Skipping {name} — loader returned None")
                continue

            adata = run_qc(adata, name, mito_modules)

            # Save
            out_path = PROCESSED / f"{name}_qc.h5ad"
            adata.write_h5ad(str(out_path))
            log(f"✅ Saved: {out_path} ({os.path.getsize(out_path)/1e9:.2f} GB)")

            qc_summary[name] = {
                "n_cells": adata.shape[0],
                "n_genes": adata.shape[1],
                "mt_pct_median": float(adata.obs["pct_counts_mt"].median()),
                "doublets_flagged": int(adata.obs["qc_doublet"].sum()),
                "mito_modules_scored": len([c for c in adata.obs.columns
                                            if c.startswith("mito_")]),
            }
        except Exception as e:
            log(f"❌ FAILED processing {name}: {e}")
            import traceback
            traceback.print_exc()
            qc_summary[name] = {"error": str(e)}

    # Write QC summary
    summary_file = PROCESSED / "qc_summary.json"
    with open(summary_file, "w") as f:
        json.dump(qc_summary, f, indent=2)
    log(f"\nQC summary saved: {summary_file}")

    log("\n" + "=" * 60)
    log("QC PIPELINE COMPLETE")
    log("Next: sbatch scripts/03_integrate.sh")
    log("=" * 60)
