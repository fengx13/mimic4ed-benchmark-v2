"""Stage 3: 5-year diagnosis lookback, Charlson/Elixhauser flags, 5y diagnosis history.

The 48 comorbidity flags are packed into a single int64 bitmask so that rolling
them up (diagnosis -> admission -> ED stay) is a bitwise-OR reduction over
sorted arrays rather than a 48-column groupby over millions of rows.
"""
import sys, time, pathlib
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from codemaps import CCI_9, CCI_10, ECI_9, ECI_10, CCI_CATS, ECI_CATS

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "derived"
LOOKBACK = pd.Timedelta(days=365.25 * 5)
_t0 = time.time()
def log(msg): print(f"[{time.time()-_t0:7.1f}s] {msg}", flush=True)


def or_reduce(group_ids, masks):
    """Bitwise-OR `masks` within each run of equal `group_ids` (must be sorted)."""
    starts = np.flatnonzero(np.r_[True, group_ids[1:] != group_ids[:-1]])
    return group_ids[starts], np.bitwise_or.reduceat(masks, starts)


m = pd.read_pickle(OUT / "_stage2.pkl")

log("loading hospital diagnoses")
dx = pd.read_csv(ROOT / "data/hosp/diagnoses_icd.csv.gz")
dx["icd_code"] = dx["icd_code"].astype("string").str.strip().str.upper()

log("building comorbidity bitmask on distinct codes")
uniq = dx[["icd_code", "icd_version"]].drop_duplicates().reset_index(drop=True)
is9 = uniq["icd_version"].to_numpy() == 9
code = uniq["icd_code"]
FLAG_COLS, mask = [], np.zeros(len(uniq), dtype="int64")
for prefix, m9, m10, cats in [("cci", CCI_9, CCI_10, CCI_CATS),
                              ("eci", ECI_9, ECI_10, ECI_CATS)]:
    for cat in cats:
        hit = np.where(is9, code.str.startswith(tuple(m9[cat]), na=False),
                       code.str.startswith(tuple(m10[cat]), na=False))
        mask |= hit.astype("int64") << len(FLAG_COLS)
        FLAG_COLS.append(f"{prefix}_{cat}")
uniq["mask"] = mask
log(f"{len(uniq):,} distinct codes -> {len(FLAG_COLS)} flags packed into int64")

dx = dx.merge(uniq[["icd_code", "icd_version", "mask"]], on=["icd_code", "icd_version"], how="left")
dx["mask"] = dx["mask"].fillna(0).astype("int64")
dxs = dx.sort_values("hadm_id")
h_id, h_mask = or_reduce(dxs["hadm_id"].to_numpy(), dxs["mask"].to_numpy())
hadm_mask = pd.Series(h_mask, index=h_id, name="hadm_mask")

log("building ED-stay x prior-admission pairs")
adm = pd.read_csv(ROOT / "data/hosp/admissions.csv.gz",
                  usecols=["subject_id", "hadm_id", "admittime", "dischtime"],
                  parse_dates=["admittime", "dischtime"])
pairs = m[["subject_id", "stay_id", "intime"]].merge(adm, on="subject_id", how="inner")
# Strictly prior: the admission must have ENDED before this ED visit began, so
# the index visit's own admission can never leak into its own history.
keep = (pairs["dischtime"] < pairs["intime"]) & (pairs["dischtime"] >= pairs["intime"] - LOOKBACK)
pairs = pairs.loc[keep, ["stay_id", "hadm_id", "dischtime"]].copy()
log(f"{len(pairs):,} (stay, prior admission) pairs in the {LOOKBACK.days}-day window")

pairs["mask"] = pairs["hadm_id"].map(hadm_mask).fillna(0).astype("int64")
ps = pairs.sort_values("stay_id")
s_id, s_mask = or_reduce(ps["stay_id"].to_numpy(), ps["mask"].to_numpy())
stay_mask = m["stay_id"].map(pd.Series(s_mask, index=s_id)).fillna(0).astype("int64").to_numpy()
for i, c in enumerate(FLAG_COLS):
    m[c] = ((stay_mask >> i) & 1).astype("int8")
log("comorbidity flags attached")

# ---------------------------------------------------------- 5-year dx history
log("exploding 5-year diagnosis history")
d_icd = pd.read_csv(ROOT / "data/hosp/d_icd_diagnoses.csv.gz")
d_icd["icd_code"] = d_icd["icd_code"].astype("string").str.strip().str.upper()
hist = (pairs[["stay_id", "hadm_id", "dischtime"]]
        .merge(dx[["hadm_id", "seq_num", "icd_code", "icd_version"]], on="hadm_id")
        .merge(d_icd, on=["icd_code", "icd_version"], how="left")
        .sort_values(["stay_id", "dischtime", "seq_num"], kind="stable"))
log(f"{len(hist):,} historical diagnosis rows")
hist.rename(columns={"dischtime": "dx_time", "long_title": "dx_name",
                     "icd_code": "dx_code", "icd_version": "dx_code_type"}
            )[["stay_id", "hadm_id", "dx_time", "seq_num", "dx_code",
               "dx_code_type", "dx_name"]].to_csv(
                   OUT / "dx_history_5y.csv.gz", index=False, compression="gzip")

log("joining history into per-stay strings")
sid = hist["stay_id"].to_numpy()
bounds = np.flatnonzero(np.r_[True, sid[1:] != sid[:-1]])
groups = np.r_[bounds, len(sid)]
uniq_sid = sid[bounds]


def join_col(values):
    v = np.asarray(values, dtype=object)
    return ["|".join(v[groups[i]:groups[i + 1]]) for i in range(len(uniq_sid))]


agg = pd.DataFrame({
    "stay_id": uniq_sid,
    "dx_name_5y": join_col(hist["long_title"].fillna("").astype(str).to_numpy()),
    "dx_code_5y": join_col(hist["icd_code"].fillna("").astype(str).to_numpy()),
    "dx_code_type_5y": join_col(hist["icd_version"].astype(str).to_numpy()),
    "dx_time_5y": join_col(hist["dischtime"].dt.strftime("%Y-%m-%d %H:%M:%S").to_numpy()),
    "n_dx_5y": np.diff(groups),
})
m = m.merge(agg, on="stay_id", how="left")
m["n_dx_5y"] = m["n_dx_5y"].fillna(0).astype("int32")

m.to_pickle(OUT / "_stage3.pkl")
log(f"stage 3 complete -> {m.shape}")
