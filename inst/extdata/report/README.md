# Reference files of pai_deg_report()

## chem_common_response_reference_v2.csv

How each of the 5,425 compounds of the L1000 Chem library usually ranks in unrelated perturbation
experiments. Source: 12 contrasts from 9 public GEO series, scored with the Chem library at the
package settings. For each compound, its within-contrast percentile (rank of the activity score /
5,425; 1 = most Concordant) is averaged within a group of contrasts that compare the same thing;
the median over the 9 groups is `reference_percentile`.

| contrast (input column) | GEO series | cells | arms compared (treated vs control) | group |
|---|---|---|---|---|
| GSE278231_Carboplatin | GSE278231 | MCF-7 | carboplatin vs control | carboplatin |
| GSE181517_Doxorubicin_g1 | GSE181517 | hiPSC-derived cardiomyocytes, RARG-WT/S427L | 1 uM doxorubicin vs DMSO | doxorubicin |
| GSE181517_Doxorubicin_g2 | GSE181517 | hiPSC-derived cardiomyocytes, RARG-WT/WT | 1 uM doxorubicin vs DMSO | doxorubicin |
| GSE164637_Doxorubicin_g1 | GSE164637 | hiPSC-derived cardiomyocytes, no doxorubicin | RARG S427L carriers vs controls | RARG S427L genotype |
| GSE164637_Doxorubicin_g2 | GSE164637 | hiPSC-derived cardiomyocytes, 1 uM doxorubicin | RARG S427L carriers vs controls | RARG S427L genotype |
| GSE268568_Indisulam | GSE268568 | CAL27, wild type and LATS1/2 knockout (both in both arms) | indisulam vs DMSO (1) | indisulam |
| GSE276457_Palbociclib | GSE276457 | HL-60 | palbociclib vs untreated | palbociclib |
| GSE209538_Spermidine | GSE209538 | primary melanocytes | 1 uM spermidine vs PBS | spermidine |
| GSE247920_Tamoxifen_g1 | GSE247920 | MDA-MB-231 | medium of tamoxifen-treated vs DMSO-treated MCF7 cells | tamoxifen-conditioned medium |
| GSE247920_Tamoxifen_g2 | GSE247920 | BT-549 | medium of tamoxifen-treated vs DMSO-treated MCF7 cells | tamoxifen-conditioned medium |
| GSE134127_Cadmium | GSE134127 | MCF7 and the MCF7-derived lines Cd7 and Cd12 (all three in both arms) | 100 nM ICI-182,780 (fulvestrant) vs mock | ICI-182,780 (fulvestrant) |
| GSE245491_Berbamine | GSE245491 | H9 | 1 uM PA4 vs DMSO | PA4 |

The arms were read from the GEO sample records. The input column names keep the labels the
contrasts had when they were collected. For five contrasts these labels do not describe what the
arms compare: GSE164637 (two genotype contrasts at equal doxorubicin exposure), GSE134127
(fulvestrant, not cadmium) and GSE247920 (two conditioned-medium contrasts). For GSE245491 the
compound is recorded as PA4; the records do not say whether PA4 is berbamine. (1) Genotypes as given by the sample titles; the genotype
characteristic in GEO disagrees with the titles for 6 of the 12 samples.

| column | meaning |
|---|---|
| compound | Chem library compound name |
| reference_percentile | median percentile over the 9 groups |
| n_groups | groups with a percentile |
| n_groups_top10 / n_groups_bottom10 | groups in which the compound was in the top / bottom 10% |
| common_response | `usually_top` (reference >= 0.80, 193 compounds), `usually_bottom` (<= 0.20, 138), or empty |

Input: the percentile matrix of the 12 contrasts (`PERCENTILE_MATRIX.csv`), SHA-256
e95bf708fa4399acf49f35e1127535cb7a64e5095f70e52851262f2a055d8fb1; built with
`data-raw/report/build_common_response_reference.py`.
File SHA-256: f881417955cb84b0f43cbcc7de6e9393d28c5ee87899f03960c28f50033a6899.

## hallmark_proliferation_sets.csv

HALLMARK_E2F_TARGETS and HALLMARK_G2M_CHECKPOINT (200 genes each), human gene symbols, MSigDB
2026.1.Hs (read with msigdbr 26.1.1). MSigDB gene sets are distributed under the Creative Commons
Attribution 4.0 International License (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp);
copyright (c) 2004-2025 Broad Institute.
Cite: Liberzon A, et al. The Molecular Signatures Database Hallmark gene set collection.
Cell Systems 2015;1:417-425.
File SHA-256: 302e8d80b48b7ae3b247f46b64750f8da087993f9c483b6afee7a7f1abcf6405.
