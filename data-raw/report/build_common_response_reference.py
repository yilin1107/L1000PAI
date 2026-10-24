"""Chem common-response reference v2 for the DEG report.

How each Chem library compound usually ranks in unrelated perturbation experiments. Input: the
percentile matrix of 12 public GEO contrasts scored with the Chem library (one row per library term,
one column per contrast; percentile = rank of the activity score / 5,425, 1 = most Concordant).
Each column name keeps the label the contrast had when it was collected; GROUPS below records what
the two arms of each contrast actually compare, read from the GEO sample records (see
inst/extdata/report/README.md). Percentiles are averaged within a group, and the median over the
9 groups is the reference percentile. A compound whose reference is >= 0.80 (<= 0.20) usually sits
near the top (bottom) of unrelated perturbation experiments.

Usage, from the package root:
    python data-raw/report/build_common_response_reference.py PERCENTILE_MATRIX.csv
The input's SHA-256 is checked against the value in inst/extdata/report/README.md.
"""
import argparse, hashlib
import numpy as np, pandas as pd

INPUT_SHA256 = "e95bf708fa4399acf49f35e1127535cb7a64e5095f70e52851262f2a055d8fb1"
OUT = "inst/extdata/report/chem_common_response_reference_v2.csv"
GROUPS = {  # matrix column -> what its two arms compare
    "GSE278231_Carboplatin": "carboplatin",
    "GSE181517_Doxorubicin_g1": "doxorubicin",
    "GSE181517_Doxorubicin_g2": "doxorubicin",
    "GSE164637_Doxorubicin_g1": "RARG S427L genotype",          # same doxorubicin exposure in both arms
    "GSE164637_Doxorubicin_g2": "RARG S427L genotype",
    "GSE268568_Indisulam": "indisulam",
    "GSE276457_Palbociclib": "palbociclib",
    "GSE209538_Spermidine": "spermidine",
    "GSE247920_Tamoxifen_g1": "tamoxifen-conditioned medium",   # medium of tamoxifen-treated MCF7 cells
    "GSE247920_Tamoxifen_g2": "tamoxifen-conditioned medium",
    "GSE134127_Cadmium": "ICI-182,780 (fulvestrant)",           # MCF7 and its Cd7 / Cd12 lines in both arms
    "GSE245491_Berbamine": "PA4",
}

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("percentile_matrix", help="percentile matrix of the 12 contrasts (PERCENTILE_MATRIX.csv)")
a = ap.parse_args()
sha = hashlib.sha256(open(a.percentile_matrix, "rb").read()).hexdigest()
if sha != INPUT_SHA256:
    raise SystemExit(f"{a.percentile_matrix}: SHA-256 {sha} differs from the documented input {INPUT_SHA256}")

M = pd.read_csv(a.percentile_matrix).set_index("term")
if sorted(M.columns) != sorted(GROUPS):
    raise SystemExit("The matrix columns are not the 12 documented contrasts.")
E = M.T.groupby(M.columns.map(GROUPS)).mean().T
assert E.shape[1] == 9
out = pd.DataFrame({"compound": E.index, "reference_percentile": E.median(axis=1).round(4).values,
                    "n_groups": E.notna().sum(axis=1).values,
                    "n_groups_top10": (E >= 0.9).sum(axis=1).values, "n_groups_bottom10": (E <= 0.1).sum(axis=1).values})
out["common_response"] = np.where(out.reference_percentile >= 0.8, "usually_top", np.where(out.reference_percentile <= 0.2, "usually_bottom", ""))
out.to_csv(OUT, index=False)
print(out.common_response.value_counts().to_dict(), "sha256", hashlib.sha256(open(OUT, "rb").read()).hexdigest())
