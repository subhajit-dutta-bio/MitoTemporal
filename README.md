# MitoTemporal

Analysis code accompanying **Mitochondrial import transcription tracks estrogen receptor activity in breast cancer**.

Author: Subhajit Dutta. Correspondence: subhajit@duttasubhajit.com.

## Repository scope

This is a **code-only repository**. It contains analysis and plotting scripts, software requirements, fixed gene-set constants defining the scoring algorithms, and methods documentation. It contains **no raw or processed datasets, expression matrices, sample records, clinical crosswalks, numerical result tables, figure files, manuscript files, or reviewer archives**. Public source inputs must be obtained separately from their original repositories. Download scripts retrieve inputs to the user's own computer; downloaded files are not part of this repository.

The original analysis-stage directory names are retained because scripts use relative paths and share helpers across stages. They identify code dependencies, not separate manuscript versions to cite. Main analysis entry points and external-input requirements are mapped in [WORKFLOW.md](WORKFLOW.md).

## Setup

Use a separate Python environment. The pinned requirements record the manuscript analysis environment; some packages may require a compatible Python/platform combination.

```bash
python -m venv .venv
# Activate .venv for your operating system.
python -m pip install -r requirements.txt
python prepare_configuration.py
python check_code.py
```

The configuration step creates empty local input/output folders and serializes the fixed gene-set definitions from `gene_sets.py`. It does not obtain any observations. Single-cell processing additionally requires the packages listed in `requirements-singlecell.txt`; GPU/scVI dependencies should be installed in a separate suitable environment.

## Running analyses

1. Obtain the original public source files and annotation resources identified in the relevant `fetch_*.py` / `download_*.py` scripts and in [WORKFLOW.md](WORKFLOW.md). Inspect the selected script before running it: downloads can be large, and source resources remain subject to their providers' terms.
2. Place inputs at the relative paths used by that stage. HGNC annotation is expected at `revision6/metadata/hgnc_complete_set.txt`. Earlier stages create matrices, pairing manifests and scores consumed by later stages. These intermediate files are generated locally and are not distributed here.
3. Execute one stage in its own Python process, for example `python revision8/analyze_e2.py`. Stages reuse a module named `common`; avoid importing several stages into the same persistent interpreter.
4. Run downstream specificity, sensitivity and plotting scripts only after their required local inputs have been generated. The code does not supply a one-command pipeline or an offline data bundle. Plotting scripts retain manuscript annotation strings, but load the plotted measurements from external/local analysis files.

Use `MITOTEMPORAL_DEPMAP_DIR` and `MITOTEMPORAL_CCLE_PROTEOME` to point the supplementary analyses at separately obtained files. Single-cell scripts accept `MITOTEMPORAL_SINGLECELL_ROOT` / `MITOTEMPORAL_REPLICATION_ROOT`. Defaults and expected filenames are documented in the code. Poppler's `pdftoppm` must be on `PATH` for the vector-figure relabelling step.

## Analysis conventions and provenance

Scripts retain the original scoring, sample-selection, paired-patient/culture contrasts, missingness rules, comparator definitions and multiple-testing families. The fixed gene-symbol lists are algorithm configuration rather than measured data. The three extension plans are included as methods documentation; they describe local prespecification after earlier project results were known, not prospective registration of the whole study.

Publication preparation changed machine-specific paths to configurable local paths and renamed plan-document extensions from `.txt` to `.md` without changing their bytes. It did not recompute or select biological results. Historical plotting scripts are retained to document figure provenance; the current four-main-figure assembly uses `revision12/build_figures.py` and `revision13/build_rebalanced_figure3.py`.

## Verification

`python check_code.py` parses all Python files and runs synthetic checks of the published scoring and statistical helpers without downloading or loading study data. Passing these checks is not a claim that every data-dependent pipeline has been rerun from a fresh clone. Those runs require the separately obtained inputs and the stage order described above.

## Citation and archival status

Use [CITATION.cff](CITATION.cff) to cite this code repository and identify the exact commit used. A Zenodo DOI has not yet been assigned; GitHub–Zenodo archiving will be configured separately. No preprint DOI is claimed here.
