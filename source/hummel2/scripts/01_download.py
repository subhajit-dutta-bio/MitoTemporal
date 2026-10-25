#!/usr/bin/env python3
"""
01_download.py — Download GEO datasets + mitochondrial gene sets
Run via: sbatch scripts/01_download.sh
"""

import os
import sys
import subprocess
import hashlib
import json
import gzip
import shutil
from datetime import datetime
from pathlib import Path

# ── Paths ──
BASE = Path(os.environ.get("MITOTEMPORAL_SINGLECELL_ROOT", str(Path(__file__).resolve().parents[1])))
RAW = BASE / "raw"
GENE_SETS = BASE / "gene_sets"
LOG = BASE / "logs"

for d in [RAW, GENE_SETS, LOG]:
    d.mkdir(parents=True, exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def wget(url, dest_dir, filename=None):
    """Download a file with wget. Returns path to downloaded file."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["wget", "-q", "--show-progress", "-c",
           "--timeout=120", "--tries=3",
           "-P", str(dest_dir)]
    if filename:
        cmd += ["-O", str(dest_dir / filename)]
    cmd.append(url)
    log(f"  wget: {url.split('/')[-1]}")
    subprocess.run(cmd, check=True)
    if filename:
        return dest_dir / filename
    return dest_dir / url.split("/")[-1]

def md5sum(filepath):
    """Compute MD5 checksum."""
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

def file_size_human(filepath):
    size = os.path.getsize(filepath)
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


# ============================================================================
# DATASET 1: GSE167036 — Liu et al. (PRIMARY LONGITUDINAL)
# ============================================================================
def download_gse167036():
    log("=" * 60)
    log("DATASET 1: GSE167036 — Liu et al. HR+ BC Longitudinal NAC")
    log("=" * 60)

    dest = RAW / "GSE167036"
    dest.mkdir(exist_ok=True)

    ftp = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167036/suppl/"
    files = [
        "GSE167036_barcodes.csv.gz",
        "GSE167036_cell_counts_matrix.mtx.gz",
        "GSE167036_features.csv.gz",
        "GSE167036_meta_all.csv.gz",
        "GSE167036_meta_tcr.csv.gz",
    ]

    for f in files:
        target = dest / f
        if target.exists():
            log(f"  Already exists: {f} ({file_size_human(target)})")
        else:
            wget(ftp + f, dest)

    # Verify
    for f in files:
        fp = dest / f
        if fp.exists():
            log(f"  ✅ {f}: {file_size_human(fp)}")
        else:
            log(f"  ❌ MISSING: {f}")
            sys.exit(1)

    log("GSE167036 download complete.\n")


# ============================================================================
# DATASET 2: GSE169246 — Zhang et al. (TNBC chemo ± anti-PD-L1)
# ============================================================================
def download_gse169246():
    log("=" * 60)
    log("DATASET 2: GSE169246 — Zhang et al. TNBC Anti-PD-L1")
    log("=" * 60)

    dest = RAW / "GSE169246"
    dest.mkdir(exist_ok=True)

    ftp = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/suppl/"

    # RNA files (we need these; ATAC files are bonus)
    rna_files = [
        "GSE169246_TNBC_RNA.barcode.tsv.gz",
        "GSE169246_TNBC_RNA.counts.mtx.gz",
        "GSE169246_TNBC_RNA.feature.tsv.gz",
    ]
    atac_files = [
        "GSE169246_TNBC_ATAC.barcode.tsv.gz",
        "GSE169246_TNBC_ATAC.counts.mtx.gz",
        "GSE169246_TNBC_ATAC.feature.tsv.gz",
    ]

    for f in rna_files + atac_files:
        target = dest / f
        if target.exists():
            log(f"  Already exists: {f} ({file_size_human(target)})")
        else:
            wget(ftp + f, dest)

    for f in rna_files:
        fp = dest / f
        if fp.exists():
            log(f"  ✅ {f}: {file_size_human(fp)}")
        else:
            log(f"  ❌ MISSING (critical): {f}")
            sys.exit(1)

    log("GSE169246 download complete.\n")


# ============================================================================
# DATASET 3: GSE176078 — Wu et al. (REFERENCE atlas)
# ============================================================================
def download_gse176078():
    log("=" * 60)
    log("DATASET 3: GSE176078 — Wu et al. Breast Cancer Atlas")
    log("=" * 60)

    dest = RAW / "GSE176078"
    dest.mkdir(exist_ok=True)

    ftp = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/"

    # The processed tar.gz has the count matrices we need
    files = [
        "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz",
    ]

    for f in files:
        target = dest / f
        if target.exists():
            log(f"  Already exists: {f} ({file_size_human(target)})")
        else:
            wget(ftp + f, dest)

    # Extract the tar.gz
    tar_file = dest / files[0]
    if tar_file.exists():
        log(f"  ✅ {files[0]}: {file_size_human(tar_file)}")
        # Check if already extracted
        extracted_marker = dest / ".extracted"
        if not extracted_marker.exists():
            log("  Extracting tar.gz...")
            subprocess.run(["tar", "xzf", str(tar_file), "-C", str(dest)],
                           check=True)
            extracted_marker.touch()
            log("  Extraction complete.")
        else:
            log("  Already extracted.")
    else:
        log(f"  ❌ MISSING: {files[0]}")
        sys.exit(1)

    log("GSE176078 download complete.\n")


# ============================================================================
# MITOCHONDRIAL GENE SETS
# ============================================================================
def download_gene_sets():
    log("=" * 60)
    log("MITOCHONDRIAL GENE SETS")
    log("=" * 60)

    # ── MitoCarta 3.0 ──
    log("Downloading MitoCarta 3.0...")
    mc_dest = GENE_SETS / "mitocarta3"
    mc_dest.mkdir(exist_ok=True)

    # Primary URL (Broad Institute)
    mc_urls = [
        "https://personal.broadinstitute.org/scalMDiRecto/MitoCarta3.0/Human.MitoCarta3.0.xls",
        "https://www.broadinstitute.org/files/shared/metabolism/mitocarta/human.mitocarta3.0.html",
    ]

    mc_file = mc_dest / "Human.MitoCarta3.0.xls"
    mc_html = mc_dest / "human.mitocarta3.0.html"

    # Try XLS first, fall back to HTML
    if not mc_file.exists():
        try:
            wget(mc_urls[0], mc_dest)
            log(f"  ✅ MitoCarta XLS: {file_size_human(mc_file)}")
        except subprocess.CalledProcessError:
            log("  ⚠️ XLS download failed (URL may have moved)")
            log("  Trying HTML table fallback...")
            try:
                wget(mc_urls[1], mc_dest, "human.mitocarta3.0.html")
                log(f"  ✅ MitoCarta HTML: {file_size_human(mc_html)}")
            except subprocess.CalledProcessError:
                log("  ⚠️ HTML also failed — will use embedded gene lists")
    else:
        log(f"  Already exists: MitoCarta XLS ({file_size_human(mc_file)})")

    # ── KEGG OXPHOS + TCA pathways ──
    log("\nFetching KEGG pathway gene lists...")
    kegg_pathways = {
        "hsa00190": "OXPHOS",
        "hsa00020": "TCA_cycle",
        "hsa04141": "protein_processing_ER",  # for mtUPR context
    }

    import urllib.request

    for pathway_id, pathway_name in kegg_pathways.items():
        outfile = GENE_SETS / f"kegg_{pathway_name}_{pathway_id}.tsv"
        if outfile.exists():
            log(f"  Already exists: {outfile.name}")
            continue
        try:
            url = f"https://rest.kegg.jp/link/hsa/{pathway_id}"
            urllib.request.urlretrieve(url, str(outfile))
            n_lines = sum(1 for _ in open(outfile))
            log(f"  ✅ KEGG {pathway_name}: {n_lines} gene entries")
        except Exception as e:
            log(f"  ⚠️ KEGG {pathway_name} failed: {e}")

    # ── Curated mito module gene sets (embedded — reliable) ──
    log("\nWriting curated mitochondrial module gene sets...")
    mito_modules = {
        "OXPHOS_Complex_I": [
            "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
            "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8",
            "NDUFV1", "NDUFV2", "NDUFV3", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4",
            "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11",
            "NDUFA12", "NDUFA13", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5",
            "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11", "NDUFC1", "NDUFC2",
            "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4",
        ],
        "OXPHOS_Complex_II": [
            "SDHA", "SDHB", "SDHC", "SDHD", "SDHAF1", "SDHAF2",
        ],
        "OXPHOS_Complex_III": [
            "MT-CYB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRB", "UQCRH",
            "UQCRQ", "UQCR10", "UQCR11",
        ],
        "OXPHOS_Complex_IV": [
            "MT-CO1", "MT-CO2", "MT-CO3", "COX4I1", "COX4I2", "COX5A", "COX5B",
            "COX6A1", "COX6A2", "COX6B1", "COX6B2", "COX6C", "COX7A1", "COX7A2",
            "COX7B", "COX7C", "COX8A", "NDUFA4",
            "COX10", "COX11", "COX14", "COX15", "COX16", "COX17", "COX18", "COX19", "COX20",
            "SCO1", "SCO2", "SURF1",
        ],
        "OXPHOS_Complex_V": [
            "MT-ATP6", "MT-ATP8", "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
            "ATP5PB", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5PD", "ATP5PF", "ATP5PO",
            "ATP5IF1", "ATP5MF", "ATP5MG",
        ],
        "TCA_cycle": [
            "CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G",
            "OGDH", "OGDHL", "DLST", "DLD", "SUCLG1", "SUCLG2", "SUCLA2",
            "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH1", "MDH2", "PC",
        ],
        "Fission": [
            "DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2",
        ],
        "Fusion": [
            "MFN1", "MFN2", "OPA1",
        ],
        "Mitophagy": [
            "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "PHB2",
            "OPTN", "NDP52", "SQSTM1", "TBK1",
        ],
        "Biogenesis": [
            "TFAM", "TFB1M", "TFB2M", "POLG", "POLG2", "TWNK",
            "PPARGC1A", "PPARGC1B", "NRF1", "GABPA", "GABPB1",
            "PPRC1", "ESRRA",
        ],
        "mtUPR": [
            "HSPD1", "HSPE1", "CLPP", "CLPX", "LONP1", "AFG3L2",
            "SPG7", "YME1L1", "HSPA9", "DNAJA3",
            "ATF4", "ATF5", "DDIT3", "CHOP",
        ],
        "mtDNA_encoded": [
            "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
            "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB",
        ],
        "Import_TIM_TOM": [
            "TOMM5", "TOMM6", "TOMM7", "TOMM20", "TOMM22", "TOMM34", "TOMM40",
            "TOMM40L", "TOMM70", "TIMM8A", "TIMM8B", "TIMM9", "TIMM10", "TIMM10B",
            "TIMM13", "TIMM17A", "TIMM17B", "TIMM21", "TIMM22", "TIMM23",
            "TIMM44", "TIMM50", "PAM16", "GRPEL1",
        ],
        "Membrane_potential_proxy": [
            "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5IF1",
            "SLC25A4", "SLC25A5", "SLC25A6",  # ANT1/2/3
            "VDAC1", "VDAC2", "VDAC3",
            "UCP2", "UCP3",
        ],
    }

    modules_file = GENE_SETS / "mito_modules.json"
    with open(modules_file, "w") as f:
        json.dump(mito_modules, f, indent=2)
    log(f"  ✅ Wrote {len(mito_modules)} modules ({sum(len(v) for v in mito_modules.values())} total genes)")
    log(f"     Saved to: {modules_file}")

    # Summary
    log("\nModule summary:")
    for name, genes in mito_modules.items():
        log(f"  {name:30s} → {len(genes):3d} genes")

    log("\nGene set download complete.\n")


# ============================================================================
# MANIFEST
# ============================================================================
def write_manifest():
    """Write a manifest of all downloaded files with checksums."""
    log("Writing download manifest...")
    manifest = {}
    for root, dirs, files in os.walk(str(RAW)):
        for f in files:
            fp = os.path.join(root, f)
            rel = os.path.relpath(fp, str(BASE))
            manifest[rel] = {
                "size_bytes": os.path.getsize(fp),
                "size_human": file_size_human(fp),
                "md5": md5sum(fp),
            }
    for root, dirs, files in os.walk(str(GENE_SETS)):
        for f in files:
            fp = os.path.join(root, f)
            rel = os.path.relpath(fp, str(BASE))
            manifest[rel] = {
                "size_bytes": os.path.getsize(fp),
                "size_human": file_size_human(fp),
                "md5": md5sum(fp),
            }

    manifest_file = BASE / "download_manifest.json"
    with open(manifest_file, "w") as f:
        json.dump(manifest, f, indent=2)
    log(f"  ✅ Manifest: {len(manifest)} files → {manifest_file}")

    # Total size
    total = sum(v["size_bytes"] for v in manifest.values())
    log(f"  Total downloaded: {total / 1e9:.2f} GB")


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    log("MitoTemporal — Data Download Pipeline")
    log(f"Base directory: {BASE}")
    log("")

    download_gse167036()
    download_gse169246()
    download_gse176078()
    download_gene_sets()
    write_manifest()

    log("=" * 60)
    log("ALL DOWNLOADS COMPLETE")
    log("Next step: sbatch scripts/02_qc.sh")
    log("=" * 60)
