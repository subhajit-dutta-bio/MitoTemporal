#!/usr/bin/env python3
"""
03_scvi_integrate.py — scVI Batch Integration (GPU)
====================================================
Integrates QC'd datasets using scVI (VAE-based, CUDA-accelerated).
Batch keys: dataset + patient_id
Outputs: processed/integrated_scvi.h5ad
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

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ──
BASE = Path(os.environ.get("MITOTEMPORAL_SINGLECELL_ROOT", str(Path(__file__).resolve().parents[1])))
PROCESSED = BASE / "processed"
FIGURES = BASE / "figures" / "integration"
FIGURES.mkdir(parents=True, exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# ============================================================================
# 1. CHECK GPU
# ============================================================================
def check_gpu():
    import torch
    log(f"PyTorch: {torch.__version__}")
    log(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        log(f"GPU: {torch.cuda.get_device_name(0)}")
        log(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    else:
        log("⚠️ No GPU — scVI will run on CPU (slower but functional)")


# ============================================================================
# 2. LOAD QC'D DATASETS
# ============================================================================
def load_qc_datasets():
    """Load all QC'd h5ad files and concatenate."""
    datasets = ["GSE167036", "GSE169246", "GSE176078"]
    adatas = []

    for name in datasets:
        path = PROCESSED / f"{name}_qc.h5ad"
        if not path.exists():
            log(f"⚠️ {name}_qc.h5ad not found — skipping")
            continue

        log(f"Loading {name}...")
        a = sc.read_h5ad(str(path))
        a.obs["dataset"] = name
        log(f"  {a.shape[0]} cells × {a.shape[1]} genes")

        # Ensure raw counts are in .X for scVI
        if "raw_counts" in a.layers:
            a.X = a.layers["raw_counts"].copy()
            log("  Restored raw counts to .X")

        adatas.append(a)

    if len(adatas) == 0:
        raise RuntimeError("No QC'd datasets found!")

    # Concatenate
    log(f"\nConcatenating {len(adatas)} datasets...")
    adata = ad.concat(adatas, join="inner", label="dataset_source",
                      keys=[a.obs["dataset"].iloc[0] for a in adatas])
    log(f"Combined: {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Ensure batch key exists
    if "dataset" not in adata.obs.columns:
        adata.obs["dataset"] = adata.obs["dataset_source"]

    return adata


# ============================================================================
# 3. PREPARE FOR SCVI
# ============================================================================
def prepare_for_scvi(adata):
    """Prepare AnnData for scVI: HVG selection on raw counts."""
    log("\nPreparing for scVI...")

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Normalize for HVG selection only
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Select HVGs per batch to avoid batch-driven HVG bias
    log("Selecting HVGs (batch-aware)...")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=4000, flavor="seurat_v3",
        layer="counts", batch_key="dataset"
    )
    n_hvg = adata.var["highly_variable"].sum()
    log(f"  {n_hvg} HVGs selected")

    # Restore raw counts to .X for scVI
    adata.X = adata.layers["counts"].copy()

    return adata


# ============================================================================
# 4. RUN SCVI
# ============================================================================
def run_scvi(adata):
    """Train scVI model for batch integration."""
    import scvi

    log("\nSetting up scVI model...")
    log(f"scvi-tools: {scvi.__version__}")

    # Register AnnData
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="dataset",
    )

    # Create model
    model = scvi.model.SCVI(
        adata,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb",  # negative binomial for scRNA-seq
    )

    log(f"Model parameters: {sum(p.numel() for p in model.module.parameters()):,}")

    # Train
    log("Training scVI (this takes 20-60 min on H100)...")
    model.train(
        max_epochs=400,
        early_stopping=True,
        early_stopping_patience=20,
        batch_size=256,
        plan_kwargs={"lr": 1e-3},
    )

    # Extract latent representation
    log("Extracting latent representation...")
    adata.obsm["X_scVI"] = model.get_latent_representation()

    # Get normalized expression
    log("Computing denoised expression...")
    adata.layers["scvi_normalized"] = model.get_normalized_expression(
        library_size=1e4
    )

    # Save model
    model_dir = BASE / "models" / "scvi"
    model_dir.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=True)
    log(f"  Model saved: {model_dir}")

    return adata, model


# ============================================================================
# 5. POST-INTEGRATION ANALYSIS
# ============================================================================
def post_integration(adata):
    """Compute neighbors, UMAP, clustering on scVI latent space."""
    log("\nPost-integration analysis...")

    # Neighbors on scVI latent space
    sc.pp.neighbors(adata, use_rep="X_scVI", n_pcs=30)

    # UMAP
    log("Computing UMAP...")
    sc.tl.umap(adata)

    # Leiden clustering at multiple resolutions
    for res in [0.5, 1.0, 1.5]:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()
        log(f"  Leiden (res={res}): {n_clusters} clusters")

    # Re-normalize .X for downstream scoring (log-normalized)
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Re-score mito modules on integrated data
    log("Re-scoring mitochondrial modules on integrated data...")
    mito_modules_file = BASE / "gene_sets" / "mito_modules.json"
    if mito_modules_file.exists():
        with open(mito_modules_file) as f:
            mito_modules = json.load(f)
        for module_name, genes in mito_modules.items():
            present = [g for g in genes if g in adata.var_names]
            if len(present) >= 2:
                score_key = f"mito_{module_name}"
                try:
                    sc.tl.score_genes(adata, gene_list=present, score_name=score_key)
                except Exception:
                    pass

    return adata


# ============================================================================
# 6. INTEGRATION FIGURES
# ============================================================================
def plot_integration(adata):
    """Generate integration QC figures."""
    log("\nGenerating integration figures...")

    try:
        # UMAP by dataset
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        sc.pl.umap(adata, color="dataset", ax=axes[0], show=False,
                    title="Dataset", frameon=False)
        sc.pl.umap(adata, color="leiden_1.0", ax=axes[1], show=False,
                    title="Leiden (1.0)", frameon=False, legend_loc="on data",
                    legend_fontsize=6)
        sc.pl.umap(adata, color="pct_counts_mt", ax=axes[2], show=False,
                    title="MT%", frameon=False)
        plt.tight_layout()
        plt.savefig(FIGURES / "umap_integration.png", dpi=200, bbox_inches="tight")
        plt.close()

        # Mito module UMAPs
        mito_cols = sorted([c for c in adata.obs.columns if c.startswith("mito_")])
        if mito_cols:
            n = min(len(mito_cols), 8)
            ncols = 4
            nrows = (n + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
            axes_flat = axes.flatten() if hasattr(axes, 'flatten') else [axes]
            for i, col in enumerate(mito_cols[:n]):
                sc.pl.umap(adata, color=col, ax=axes_flat[i], show=False,
                            title=col.replace("mito_",""), frameon=False, vmin=-1, vmax=2)
            for j in range(i+1, len(axes_flat)):
                axes_flat[j].set_visible(False)
            plt.tight_layout()
            plt.savefig(FIGURES / "umap_mito_modules_integrated.png",
                        dpi=200, bbox_inches="tight")
            plt.close()

        log(f"  Figures saved to: {FIGURES}")

    except Exception as e:
        log(f"  ⚠️ Figure error: {e}")


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    log("MitoTemporal — scVI Integration Pipeline")
    log("")
    check_gpu()
    log("")

    # Load
    adata = load_qc_datasets()

    # Prepare
    adata = prepare_for_scvi(adata)

    # Integrate
    adata, model = run_scvi(adata)

    # Post-process
    adata = post_integration(adata)

    # Plots
    plot_integration(adata)

    # Save
    out_path = PROCESSED / "integrated_scvi.h5ad"
    log(f"\nSaving integrated object: {out_path}")
    adata.write_h5ad(str(out_path))
    log(f"  Size: {os.path.getsize(out_path) / 1e9:.2f} GB")

    # Summary
    log("\n" + "=" * 60)
    log("INTEGRATION COMPLETE")
    log(f"  Cells: {adata.shape[0]}")
    log(f"  Genes: {adata.shape[1]}")
    log(f"  Datasets: {adata.obs['dataset'].nunique()}")
    log(f"  Latent dim: {adata.obsm['X_scVI'].shape[1]}")
    log(f"  Clusters (res=1.0): {adata.obs['leiden_1.0'].nunique()}")
    log(f"  Mito modules: {len([c for c in adata.obs.columns if c.startswith('mito_')])}")
    log(f"  Output: {out_path}")
    log("=" * 60)
    log("\nModule 1 COMPLETE. Ready for Module 2 (Longitudinal Analysis).")
