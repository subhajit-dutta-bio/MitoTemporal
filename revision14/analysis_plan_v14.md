PUBLIC FUNCTION / CISTROME / PROTEOME / SPECIFICITY EXTENSION
25 October 2026. Fixed before inspecting any new import-module outcomes.
Prior MitoTemporal results and general source-publication findings are known.
Cached input matrices are present but their import/ER/dependency outcomes have
not been read in this extension. This is local prespecification, not project-wide
prospective registration. Retain all planned results, regardless of direction.

PRIORITY AND CLAIMS
Execute Tier 1 (dependency and cistrome), then Tier 2 protein validation, then the
genome-wide ranking. Johansson is conditional as specified below. Public CRISPR
knockout scores measure essentiality in cultured models. A context association
does not show that ER causes the dependency or that import flux changed. ChIP
occupancy/co-occurrence is not mediation; sparse ER binding plus NRF1 binding
cannot by itself establish an indirect NRF1 pathway. Do not replace an unavailable
primary endpoint with a favourable proxy. Main-figure promotion depends on
relevance and evidential strength, not significance alone.

COMMON FIXED DEFINITIONS
Canonical 24-gene TIM/TOM module and revision10 fixed comparator sets unchanged.
Nonimport_mito contains 167 nuclear genes; use the full measured fixed set,
including the prior disjointness rules. Canonical ER programme is the existing
seven-gene set, not ESR1 alone. HGNC mapping inherited from revision6. For a module,
require at least half of the fixed members and >=3 genes (import >=12/24).
Unique approved gene symbols; ambiguous protein groups excluded, duplicate
gene-level protein quantities median-collapsed. Missing abundance is not zero.
Expression/protein values already on log scale are not re-logged.
Statistical unit: one independent cell model/patient/biological trajectory, never
one gene or technical replicate. Two-sided inference with directional hypotheses
specified below. Unadjusted 95% intervals; BH families explicitly listed. Separate
families are not project-wide FDR control. Source/rank sensitivities both retained.

TIER 1A: DEPMAP / CCLE
Pin coherent cached DepMap Public 24Q2 (Figshare article 25880521) CRISPR Chronos
gene effects, log2(TPM+1) expression and Model metadata; verify provenance/hashes.
This is a versioned release, not a claim to use the newest release. No switching
release on the basis of results. Breast = Model OncotreeLineage Breast, cancer
cell lines only; exclude normal, engineered or treatment-derived variants where
explicitly labelled. Retain one model per PatientID (lexicographically first
eligible ModelID), with missing PatientID treated as unique. One basal/default
expression profile per model. Explicit model/alias crosswalk and exclusions.

ER labels must come from an independent published receptor/protein annotation
(priority Neve 2006 primary cell-line table; explicit positive/negative only).
Ambiguous/unknown labels excluded from binary test, retained in continuous tests.
Luminal subtype or ESR1 expression is not substituted for measured receptor status.
Use existing ER-programme expression score continuously; no outcome-selected cut.
Quartiles may be displayed descriptively. Require >=5 models per ER group for
binary inference, >=20 eligible breast models for continuous breast inference.
Use >=12 measured import dependency genes, at least 80% of the available common
module per model; same gene panel across groups. Negative gene effect means more
essential. Median genome-wide gene effect is an overall screening-severity control.

Eight primary tests in one BH family:
D1. Breast ER+ minus ER- difference in (import minus other mitochondrial gene
    effect), hypothesis negative; independent-group HC3 OLS.
D2. Within breast, standardized ER-programme slope for the same differential
    dependency, adjusted for global gene-effect median and proliferation RNA;
    hypothesis negative, HC3 OLS residual-t inference.
D3. Pan-cancer import mean gene effect versus ESR1 gene effect, adjusted for global
    gene-effect median and lineage fixed effects; hypothesis positive co-dependency.
D4-D6. Same models for TOMM20, TOMM40 and TIMM23, each versus ESR1 gene effect;
    hypothesis positive. These genes fixed by user request, no selection.
D7. Pan-cancer import-minus-other-mitochondrial dependency versus ESR1 expression,
    adjusted for lineage, global gene-effect median and proliferation RNA;
    hypothesis negative.
D8. Same as D7 with the ER-programme score instead of ESR1, hypothesis negative.
Lineages with <10 eligible models are combined into Other, based only on n.
Continuous variables standardized in the eligible model set; intercept and
lineage indicators not standardized. Global scores exclude ESR1 and fixed import
genes to avoid mechanical self-correlation. Record rank/design/condition numbers.
If rank-deficient or below n gate, mark unavailable rather than drop covariates.

Context/sensitivities: raw import and comparator group differences, unadjusted
Spearman co-dependencies and ESR1-expression correlations descriptive; all 24
gene-specific ER differential-dependency tests form a separate 24-test BH family.
Leave-one-breast-model-out coefficient ranges; source-rank ER programme sensitivity;
excluding Breast in D3-D8; adjust D3-D6 for nonimport mitochondrial dependency as
an explicit sensitivity. No pooled gene-as-sample t test. Broad essentials or
lineage-confounded correlations do not qualify as ER-selective import dependence.

TIER 1B: ER CISTROME
Design-selected MCF7 GSE99680 wild-type vehicle/E2 45-minute ChIP-seq, two deposited
replicates/condition (GSM3595410-13). Check source paper and processing; if the
E2 peak sets are already vehicle-subtracted, do not treat condition peak counts
as a calibrated differential-binding assay. Require reproducible overlap across
both condition replicates for primary occupancy; union is a fixed sensitivity.
Map GRCh37 peaks to GRCh38 using UCSC chain with unique, same-chromosome,
orientation-consistent endpoint mapping; log failed/ambiguous intervals.
Use inherited GENCODE v47 protein-coding TSS annotation and NRF1 ENCSR135ANT
binding. Gene-level co-occurrence means both factors have promoter support; it
does not imply simultaneous binding in the same cells or direct protein contact.
Promoter = +/-1 kb of any annotated TSS; +/-2 kb sensitivity.
Independent enhancer/promoter links: MCF7 POLR2A ChIA-PET ENCSR059HDE, or a
verified ENCODE processed loop file from that experiment. Replicate-supported
intrachromosomal links <=1 Mb; one anchor overlaps a promoter and the other is
distal (>1 kb from that gene's TSS). Call these linked distal regions, not proven
causal enhancers. If no usable matched-build loop data exist, mark distal tests
unavailable rather than substitute arbitrary nearest-gene assignments.
Background = GSE78167 baseline mean TPM>=1, valid TSS and relevant assay coverage.
Match 10,000 random nonimport gene panels on baseline-expression decile and TSS
count bins (1,2,3-5,>=6); no outcome-based matching. Primary enrichment permutation
P=(1+number null counts>=observed)/(10001), effect fold over matched mean.
Fisher odds ratios/context provided separately. Freeze seed 20261025.

Seven primary tests, one BH family:
C1-C4. Import promoter and promoter-or-linked-distal ER occupancy in vehicle and
       E2, versus expression/TSS-matched background, hypothesised enrichment.
C5-C6. NRF1-plus-ER promoter co-occurrence in import genes, vehicle/E2, same matched
       background, hypothesised enrichment.
C7. Among fixed expressed import genes, E2-bound (promoter or linked distal) minus
    unbound one-hour induction in GSE78167: calculate the gene-set difference
    separately in each of three biological trajectories and test their mean with
    t2 inference. Require >=4 genes per group; otherwise unavailable. Genes are
    not treated as biological replicates. Hypothesis bound genes respond earlier.
Secondary: 1/2-hour mean contrast; onset times (first >=0.2 log2 increase followed
by another >=0.2 timepoint) descriptive; peak-distance and each replicate summaries.
Absence of detectable ER binding does not establish an indirect response.

TIER 2: INDEPENDENT PROTEOMICS
Nusinow 2020 CCLE normalized protein matrix (375 source lines) and metadata,
cached from Gygi laboratory's public resource; verify source/hash and aliases.
One cell line/model after averaging repeat TMT-plex measurements, not n=375
breast samples. Use independent ER labels as above. Require import >=12 proteins
and a constant panel quantified in >=80% of each comparison group; each sample
must have >=12 and >=80% of this panel. Mean deposited log quantities, no imputation.
Matched RNA/protein comparison uses exactly the same import genes. Do not infer
ER negativity from undetected ESR1 protein. Report actual breast and ER group n.
Three CCLE primary tests, one BH family:
P1. Breast import protein ER+ minus ER-, adjusted for nonimport mito protein;
    hypothesis positive (HC3).
P2. Within breast import protein versus ESR1 protein, adjusted for nonimport mito
    protein; hypothesis positive (HC3), complete detected pairs only.
P3. Within breast matched import RNA/protein Spearman correlation, hypothesis
    positive; model bootstrap 95% interval (5,000, seed20261025).
Raw group contrast and protein/dependency coupling in the exact model intersection
are fixed contextual analyses, clearly distinct from independent replication.

Asleh 2022, Nat Commun doi:10.1038/s41467-022-28524-0, 300 source specimens in
two clinical-era cohorts (178 and122), selected 75 per PAM50 subtype overall.
Use author-processed quantitative proteins and patient metadata, verify overlap
and cohort labels; no pooling raw sample numbers as independent replication.
Primary protein gate as above; if <12 canonical import proteins, report descriptive
coverage/observed-subset summaries, no substitute whole-module P or survival claim.
Three Asleh primary tests, one BH family if eligible:
A1. Luminal A/B versus basal import protein, adjusted for nonimport mitochondrial
    protein and cohort-era indicator; hypothesis positive, HC3.
A2. ER status versus import protein after same adjustment, if deposited receptor
    labels are available and both groups >=10; hypothesis positive.
A3. Continuous import-protein score (perSD) versus breast-cancer-specific/disease-
    specific survival among luminal cases, Cox with cohort-era strata, age and
    LumB indicator. No optimal cutpoints. Outcome field must be unambiguously
    breast-cancer-specific; otherwise mark unavailable. Direction exploratory,
    report two-sided interval. Require >=20 events and >=10 events per fitted
    coefficient; otherwise unadjusted stratified model only as exploratory.
Report unadjusted source/rank subtype profiles and separate era estimates as
sensitivities. Different eras/treatment regimens limit endocrine-specific claims.

Johansson 2019 Oslo2/PXD008841/PXD011385 conditional third proteome cohort:
only proceed if CCLE and Asleh both meet coverage and have the same sign of the
unadjusted luminal/ER protein contrast (regardless of significance). This is a
predeclared concordance gate, not choosing among multiple datasets for a positive
result. If triggered, test luminal versus basal import protein adjusted for other
mitochondrial protein, one primary test (two-sided). Verify independent donors;
if overlap cannot be resolved, no claim of a third independent replication.

TIER 3: GENOME-WIDE MODULE RANKING
MSigDB2025.1.Hs: Hallmark, GO biological process and Reactome; all approved sets
of 15-500 genes after HGNC canonicalization, deduplicate exact gene memberships
with a full alias map. Do not force a predetermined number such as 5,000.
Use unchanged existing expression matrices and sample/patient contrasts:
GSE78167 fixed6/8/12/24h mean vs matched baseline;
GSE89888 four vehicle mutant-minus-WT contexts (MCF7/T47D Y537S/D538G);
GSE231397 E2+fulvestrant vs E2;
six AI cohorts with their original paired populations; additionally POETIC treated
minus untreated change as a controlled sensitivity. No donor pooling/meta-P.
Score all eligible modules, retain >=half fixed members and >=10 genes, require
>=80% of original eligible sample/patient n for a ranking context. Import retains
its existing >=12/24 gate. No DE-gene filtering or outcome-selected pathways.
Primary rank = signed mean within-sample gene-percentile-rank change, with
sign + for E2/activating-mutant retention and - for receptor/AI blockade.
Rank=(1+number comparison modules with strictly greater directional effect),
top-tail percent=100*rank/(N+1); lower is stronger. Ties midrank; retain N and
effective module/sample coverage in every context. Secondary ranking by unit-
level t statistic (paired trajectory/patient or independent culture as appropriate)
and mean source-scale change, shown separately. Module tests get BH across all
eligible modules plus import within each context; dependence/redundancy among
gene sets precludes treating ranks as independent replications. Size-matched
15-40 gene sets and excluding >=50% import-overlap sets are fixed sensitivities.
Report all contexts, not just top ranks. A large induction shared with translation
pathways is not import exclusivity. Ranking cannot replace the fixed direct
comparator test or establish a functional mechanism.

RELEASE
Integrate results honestly into a new revision preserving revision13. Keep four
main figures; move full screens, all families and secondary diagnostics to the
supplement. Prior numerical results unchanged. Update abstract/cover only for
material, adequately supported advances; no invented causal or clinical claims.
Provide data/source/identifier manifests, plan hashes, code and compact reviewer
ZIP under50MB; raw bulk inputs may remain locally with extraction recipes.
