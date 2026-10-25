GSE8597 targeted mechanistic extension - fixed before expression download/scoring
25 October 2026 Europe/Berlin; 24 October 2026 23:28 UTC.

This plan follows known revision-11 results and browser inspection of GEO and the
original publication design. It is a local extension plan, not project-wide
prospective registration. No GSE8597 expression values/module effects have been
inspected, apart from the GEO sample viewer's default first 20 probes, none used
to choose the module or endpoint.

Question: Does the 24-hour estradiol-associated canonical TIM/TOM RNA response
change during protein-synthesis blockade? This tests the later steady-state RNA
response, not the 3-5-hour peak, direct ER DNA binding, NRF1 mediation or import flux.

Source: GSE8597 / Bourdeau et al., Nucleic Acids Res 2008;36:76-93,
doi:10.1093/nar/gkm945. MCF7 vehicle/E2 crossed with vehicle/CHX, four deposited RNA
samples per condition. CHX 10 micrograms/ml begins one hour before 25 nM E2 for
24 hours. The paper reports quadruplicate samples, not verified cross-arm pairing.

Processing: deposited RMA log2 intensities, no additional logarithm or DE-gene
filter. GPL570 annotation retained locally; drop control/unmapped/ambiguous probe
symbols, convert unambiguous HGNC aliases using the inherited mapping and median
collapse to unique approved genes. Canonical 24-gene import module unchanged;
at least 12 measured genes are required, complete across all 16 samples. Use
the same measured genes in all groups. Within-array gene percentile rank is a
separate sensitivity using the same complete-gene universe.

Primary contrasts per scale, all two-sided, in one three-test BH family:
1. E2 - vehicle without CHX (hypothesised positive).
2. E2+CHX - CHX (remaining E2 effect during blockade; direction not constrained).
3. (E2+CHX - CHX) - (E2 - vehicle) interaction (hypothesised negative if the later
   response requires new protein synthesis).
Use factorial OLS with HC3 residual-t inference (df=12 for 16 complete arrays),
unadjusted 95% intervals and all three contrasts. Treat deposited samples as
independent; do not infer pairing solely from matching rep-number labels.
Report the replicate unit as deposited RNA samples, without claiming independent
engineering derivations or culture batches. Classic homoskedastic OLS and leave-
one-sample-out estimate ranges are fixed diagnostic sensitivities, not replacement
primary P values. Do not infer absence from a nonsignificant interaction.

Context: unchanged non-import mitochondrial, cytosolic ribosome, translation-
factor, proliferation and ER-programme sets. Three contrasts per programme/scale
in separate three-test families. Direct import-minus-the-four-inherited-comparator
outcomes form a twelve-test BH family per scale. No new gene sets are selected
after scoring. NRF1/TFAM gene responses may be shown descriptively without a new
gene-as-replicate test.

Display: all four groups and all primary contrasts; full fixed-gene heatmap and
source/rank sensitivities. Preserve complete tables regardless of direction.
Place this targeted extension in Supplementary Figure S17; any main-text statement
must distinguish the observed result from mechanistic interpretation.

Interpretation: CHX is pleiotropic, can change RNA stability and ER abundance, and
was applied for a prolonged interval. An interaction supports translation-sensitive
regulation under these conditions, not a unique regulator or a proven secondary
transcriptional cascade. Persistence under CHX is compatible with pre-existing
regulators, not proof of direct ER occupancy. Broad control responses must be
reported alongside import to avoid claiming import specificity from CHX alone.
