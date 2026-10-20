#!/usr/bin/env python3
"""Build the MIMIC-IV-ED master analysis dataset (one row per ED stay).

Usage:  python3 scripts/build_master.py
Writes: derived/master_ed.parquet, derived/master_ed.csv.gz,
        derived/dx_history_5y.parquet, derived/column_dictionary.csv
"""
import sys, time, pathlib
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from codemaps import CCI_9, CCI_10, ECI_9, ECI_10, CCI_CATS, ECI_CATS
import scores as sc

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "derived"
OUT.mkdir(exist_ok=True)
LOOKBACK_YEARS = 5
SEC = 10_000_000_000  # subject-block stride for the integer sort key

_t0 = time.time()
def log(msg):
    print(f"[{time.time()-_t0:7.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- utilities
def _key(subject_codes, times):
    """Exact int64 key that sorts by (subject, time); enables global searchsorted."""
    secs = (times.astype("datetime64[s]").astype("int64")
            if hasattr(times, "astype") else times)
    return subject_codes.astype("int64") * SEC + secs


def prior_counts(base_code, base_t, ev_code, ev_t, windows):
    """For each base row, count events in [t - window, t) for the same subject."""
    order = np.lexsort((ev_t, ev_code))
    ev_keys = _key(ev_code[order], ev_t[order])
    now = _key(base_code, base_t)
    upto = np.searchsorted(ev_keys, now, "left")
    out = {}
    for w in windows:
        back = _key(base_code, base_t - np.int64(w * 86400))
        out[w] = upto - np.searchsorted(ev_keys, back, "left")
    return out


def next_event(base_code, base_t, ev_code, ev_t):
    """Time of each subject's next event strictly after base_t (NaT if none)."""
    order = np.lexsort((ev_t, ev_code))
    ev_code_s, ev_t_s = ev_code[order], ev_t[order]
    ev_keys = _key(ev_code_s, ev_t_s)
    pos = np.searchsorted(ev_keys, _key(base_code, base_t), "right")
    ok = (pos < len(ev_keys)) & (ev_code_s[np.clip(pos, 0, len(ev_keys) - 1)] == base_code)
    res = np.full(len(base_t), np.iinfo("int64").min, dtype="int64")
    res[ok] = ev_t_s[pos[ok]]
    return pd.to_datetime(np.where(ok, res, np.nan), unit="s")


def flags_from_codes(dx, mapping9, mapping10, prefix, group_key):
    """Boolean comorbidity flags per `group_key` from a long diagnosis frame."""
    code = dx["icd_code"].astype("string").str.strip().str.upper()
    v9 = dx["icd_version"].to_numpy() == 9
    out = pd.DataFrame({group_key: dx[group_key].to_numpy()})
    for cat in mapping9:
        pat9 = tuple(mapping9[cat])
        pat10 = tuple(mapping10[cat])
        hit = np.where(v9, code.str.startswith(pat9, na=False),
                       code.str.startswith(pat10, na=False))
        out[f"{prefix}_{cat}"] = hit
    return out.groupby(group_key).max()


# ---------------------------------------------------------------- 1. cohort
log("loading ED stays + patients + admissions")
ed = pd.read_csv(ROOT / "data/ed/edstays.csv.gz",
                 parse_dates=["intime", "outtime"])
pat = pd.read_csv(ROOT / "data/hosp/patients.csv.gz", parse_dates=["dod"])
adm = pd.read_csv(ROOT / "data/hosp/admissions.csv.gz",
                  parse_dates=["admittime", "dischtime", "deathtime",
                               "edregtime", "edouttime"])

m = ed.rename(columns={"race": "ethnicity"}).merge(
        pat[["subject_id", "gender", "anchor_age", "anchor_year",
             "anchor_year_group", "dod"]],
        on="subject_id", how="left", suffixes=("_ed", ""))
m = m.drop(columns=["gender_ed"])
m = m.merge(adm[["hadm_id", "admittime", "dischtime", "deathtime", "edregtime",
                 "edouttime", "insurance", "hospital_expire_flag",
                 "admission_type", "discharge_location"]],
            on="hadm_id", how="left")

m["in_year"] = m["intime"].dt.year
m["age"] = m["anchor_age"] + (m["in_year"] - m["anchor_year"])
log(f"cohort: {len(m):,} ED stays, {m.subject_id.nunique():,} patients")

# ---------------------------------------------------------------- 2. outcomes
log("deriving mortality / disposition outcomes")
m["outcome_inhospital_mortality"] = (m["hospital_expire_flag"] == 1).astype("int8")
m["ed_death"] = (m["disposition"] == "EXPIRED").astype("int8")
m["before_ed_mortality"] = (m["dod"].notna() & (m["dod"] < m["intime"].dt.normalize())).astype("int8")
m["time_death_discharge"] = m["dod"] - m["dischtime"]
m["time_death_ed_depart"] = m["dod"] - m["outtime"]
m["time_death_ed_arrival"] = m["dod"] - m["intime"]
m["ed_los"] = m["outtime"] - m["intime"]
m["ed_los_hours"] = m["ed_los"].dt.total_seconds() / 3600
m["outcome_hospitalization"] = m["hadm_id"].notna().astype("int8")

log("linking ICU transfers")
icu = pd.read_csv(ROOT / "data/icu/icustays.csv.gz", parse_dates=["intime", "outtime"])
icu_s = icu.sort_values(["subject_id", "intime"])
nxt_icu = next_event(m["subject_id"].to_numpy(),
                     m["outtime"].to_numpy().astype("datetime64[s]").astype("int64"),
                     icu_s["subject_id"].to_numpy(),
                     icu_s["intime"].to_numpy().astype("datetime64[s]").astype("int64"))
m["intime_icu_after_ed"] = nxt_icu
m["time_to_icu_transfer"] = m["intime_icu_after_ed"] - m["outtime"]
m["time_to_icu_transfer_hours"] = m["time_to_icu_transfer"].dt.total_seconds() / 3600
m["outcome_icu_transfer_12h"] = ((m["time_to_icu_transfer_hours"] >= 0)
                                 & (m["time_to_icu_transfer_hours"] <= 12)).astype("int8")
m["outcome_critical"] = ((m["outcome_icu_transfer_12h"] == 1)
                         | (m["outcome_inhospital_mortality"] == 1)).astype("int8")

# ---------------------------------------------------------------- 3. history
log("counting prior ED / hospital / ICU utilisation")
b_code = m["subject_id"].to_numpy()
b_t = m["intime"].to_numpy().astype("datetime64[s]").astype("int64")
assert b_t.min() > 0, "negative epoch seconds would break the integer sort key"

for label, ev in [("ed", (ed["subject_id"].to_numpy(),
                          ed["intime"].to_numpy().astype("datetime64[s]").astype("int64"))),
                  ("hosp", (adm["subject_id"].to_numpy(),
                            adm["admittime"].to_numpy().astype("datetime64[s]").astype("int64"))),
                  ("icu", (icu["subject_id"].to_numpy(),
                           icu["intime"].to_numpy().astype("datetime64[s]").astype("int64")))]:
    cnt = prior_counts(b_code, b_t, ev[0], ev[1], [30, 90, 365])
    for w, v in cnt.items():
        m[f"n_{label}_{w}d"] = v.astype("int32")

nxt_ed = next_event(b_code, b_t, ed["subject_id"].to_numpy(),
                    ed["intime"].to_numpy().astype("datetime64[s]").astype("int64"))
m["next_ed_visit_time"] = nxt_ed
m["next_ed_visit_time_diff"] = m["next_ed_visit_time"] - m["outtime"]
m["next_ed_visit_time_diff_days"] = m["next_ed_visit_time_diff"].dt.total_seconds() / 86400
m["outcome_ed_revisit_3d"] = ((m["next_ed_visit_time_diff_days"] >= 0)
                              & (m["next_ed_visit_time_diff_days"] <= 3)).astype("int8")

m.to_pickle(OUT / "_stage1.pkl")
log(f"stage 1 complete -> {m.shape}")
