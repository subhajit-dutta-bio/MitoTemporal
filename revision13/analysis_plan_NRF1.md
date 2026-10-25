NRF1 loss/rescue extension: fixed before expression download or scoring
25 October 2026. Local prespecification after revision 11 and source-design review;
not prospective registration of the project. The publication's general NRF1/exon
findings are known; canonical import-module outcomes have not been inspected.

Question: Does nuclear respiratory factor 1 (HGNC NRF1, not NFE2L1) regulate the
unchanged 24-gene TIM/TOM transcriptional module? Can siRNA-resistant full-length
NRF1 rescue loss, and does its exon-7 isoform matter? These experiments lack an
ER/E2 factorial perturbation and are not a test of ER-to-NRF1 mediation in breast
cancer or a measurement of protein-import activity.

Selection: GSE292988 provides six HEK293 Flp-In groups with three deposited RNA
replicates each: siNT, siNRF1, siNT+FL, siNRF1+FL, siNT+DE7, siNRF1+DE7.
Kumari et al., Nat Commun 2026;17:1218, doi:10.1038/s41467-026-68774-w,
documents siRNA-resistant inducible full-length or exon-7-deleted constructs.
The article specifies dox 6 h after transfection and harvest 48 h after induction;
GEO instead says dox after 24 h. Record this unresolved timing discrepancy.
Replicate numbers are not sufficient evidence of cross-condition pairing. The
unit is a deposited RNA sample; independent cell-line engineering is not assumed.

Orthogonal validation GSE292986: endogenous exon-7 splice editing in HEK293T,
three NRF1-targeted and three intergenic-control RNA samples per platform:
SaCas9-ABE, SpCas9-ABE, SpCas9-CBE. Same publication, not an independent study.
All three platforms retained irrespective of editing efficiency or module effect.

Processing: raw integer gene counts, no DE-gene filtering. Strip Ensembl version
suffixes, map via supplied symbols if present, otherwise inherited HGNC approved
Ensembl IDs. Sum duplicate unique approved genes; exclude ambiguous identifiers.
For each dataset, retain genes with CPM >=1 in at least 3 samples. DESeq-style
median-of-ratios size factors over genes positive in every sample; use
log2(normalized count+1) as the source scale. Sample rank sensitivity is gene
percentile rank over the same expressed-gene universe. Require >=12 of 24 import
genes; same measured genes across all groups. Report module coverage, library
sizes, size factors and all sample scores. NRF1 total-gene RNA and TFAM are
descriptive target/context measures; transgene rescue makes total NRF1 RNA an
imperfect measure of endogenous depletion. No external spike-ins are available.

Inference: independent group-means OLS with HC3 covariance and residual-t
intervals, df=12 in rescue and df=12 in editing. Two-sided tests, unadjusted 95%
confidence intervals. Genes are never treated as independent biological samples.
Three primary rescue contrasts, one BH family per module and scale:
  A. siNRF1 - siNT (hypothesised negative).
  B. (siNRF1+FL - siNT+FL) - (siNRF1 - siNT), rescue interaction
     (hypothesised positive; controls general FL-expression effects).
  C. siNRF1+DE7 - siNRF1+FL (isoform comparison, no direction assumed).
Eight secondary contrasts, separate eight-test BH families per module/scale:
  FL rescue vs KD; DE7 rescue vs KD; FL rescue vs siNT; DE7 rescue vs siNT;
  siNT+FL vs siNT; siNT+DE7 vs siNT;
  DE7 rescue interaction; DE7-minus-FL rescue-interaction difference.
Orthogonal editing: targeted-minus-intergenic per platform, three-test BH family
per module and scale. No pooled platform effect or meta-analysis.

Context modules unchanged: non-import nuclear mitochondrial genes, cytosolic
ribosome, translation factors, proliferation and ER programme. Use exactly the
revision-10 gene definitions and disjointness conventions. Three primary rescue
contrasts and three editing contrasts on direct import-minus-four-comparator
scores form a single 24-test specificity family per scale. Context-module tests
use the same within-dataset primary and secondary families as import. No
post-scoring selection of genes or comparator sets.

Fixed diagnostics: classic OLS; leave-one-sample-out effect ranges; replicate-label
fixed-block OLS sensitivity with HC3 (df=10), explicitly unverified cross-arm
matching, never replacing primary inference. Gene-wise fold changes and
concordance heatmaps descriptive only. Present every planned contrast and scale.

Display/interpretation: full rescue groups, primary effects and editing contrasts
in a supplementary mechanistic figure; update a main regulatory figure only if
the evidence is material. A positive rescue establishes context-specific NRF1
control of this RNA score, not ER mediation, import flux or cancer phenotype.
No new ER-mediated import claim from import-engagement maps or non-ER stress
screens. Keep searched direct-import datasets in the candidate-screen record.
