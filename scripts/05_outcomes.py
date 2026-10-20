"""Stage 5: the 5 benchmark outcomes (ACS, sepsis, AHF, stroke, PE), severity
scores, chest-pain cohort flags, and ventilation timing.

This is the source project's stage4_outcomes.py + stage8_cp_outcomes.py,
trimmed to only the 5 target outcomes. Every PREDICTOR column either script
produced (scores, dxes/dxes_code, cohort_chest_pain, dx_chest_pain,
first_vent_time, time_to_invasive_vent_hours) is kept unchanged; only the
OUTCOME columns are reduced — from 23 icd_codes.csv phenotypes + 14 chest-pain
differential conditions + 4 schema-alias duplicates down to exactly:

    outcome_sepsis    (exact match, icd_codes.csv dx_sepsis)
    outcome_hf_acute  (exact match, icd_codes.csv dx_hf_acute -- "AHF")
    outcome_stroke    (exact match, icd_codes.csv dx_stroke)
    outcome_cp_acs    (prefix match, cp_codes.py "acs" -- the broader,
                       code-audited ACS composite: STEMI + NSTEMI +
                       unstable angina + other acute ischemic heart disease)
    outcome_cp_pe     (prefix match, cp_codes.py "pe")

outcome_cp_acs / outcome_cp_pe keep their ORIGINAL column names from the full
master dataset (rather than being renamed to outcome_acs / outcome_pe) so a
column-by-column diff against master_ed.csv lines up exactly.
"""
import sys, time, pathlib
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import scores as sc
from cp_codes import CP_9, CP_10
from outcome_codes_exact import EXACT_OUTCOMES

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "derived"
_t0 = time.time()
def log(msg): print(f"[{time.time()-_t0:7.1f}s] {msg}", flush=True)

m = pd.read_pickle(OUT / "_stage3.pkl")

# ------------------------------------------------- current-visit diagnoses
log("loading current-visit diagnoses (ED + linked admission)")
eddx = pd.read_csv(DATA / "ed/diagnosis.csv.gz")
eddx["icd_code"] = eddx["icd_code"].astype("string").str.strip().str.upper()
eddx = eddx.sort_values(["stay_id", "seq_num"])
g = eddx.groupby("stay_id")
m = m.merge(pd.DataFrame({
    "dxes": g["icd_title"].agg(lambda s: "|".join(s.fillna("").astype(str))),
    "dxes_code": g["icd_code"].agg(lambda s: "|".join(s.fillna("").astype(str))),
}), on="stay_id", how="left")

hdx = pd.read_csv(DATA / "hosp/diagnoses_icd.csv.gz",
                  usecols=["hadm_id", "icd_code", "icd_version"])
hdx["icd_code"] = hdx["icd_code"].astype("string").str.strip().str.upper()

# ------------------------------------------------- sepsis / AHF / stroke (exact)
log("flagging sepsis / AHF / stroke (exact ICD match)")
ed_key = eddx["icd_code"] + "_" + eddx["icd_version"].astype(str)
h_key = hdx["icd_code"] + "_" + hdx["icd_version"].astype(str)
for name, (label, c9, c10) in EXACT_OUTCOMES.items():
    keys = {f"{c}_9" for c in c9} | {f"{c}_10" for c in c10}
    ed_hit = eddx.loc[ed_key.isin(keys), "stay_id"].unique()
    h_hit = hdx.loc[h_key.isin(keys), "hadm_id"].unique()
    m[f"outcome_{name}"] = (m["stay_id"].isin(ed_hit) | m["hadm_id"].isin(h_hit)).astype("int8")

# ------------------------------------------------- scores (unchanged from source)
log("computing severity scores")
ALERT_GCS = pd.Series(15.0, index=m.index)
ALERT_AVPU = pd.Series(0.0, index=m.index)
NO_SUPP_O2 = pd.Series(0.0, index=m.index)
temp_c = sc.f_to_c(m["triage_temperature"])

m["score_CCI"] = sc.score_cci(m, m["age"])
m["score_CART"] = sc.score_cart(m["triage_resprate"], m["triage_heartrate"],
                                m["triage_dbp"], m["age"])
m["score_REMS"] = sc.score_rems(m["age"], m["triage_MAP"], m["triage_heartrate"],
                                m["triage_resprate"], m["triage_o2sat"], ALERT_GCS)
m["score_NEWS"] = sc.score_news(m["triage_resprate"], m["triage_o2sat"], temp_c,
                                m["triage_sbp"], m["triage_heartrate"],
                                pd.Series(1.0, index=m.index), NO_SUPP_O2)
m["score_NEWS2"] = sc.score_news2(m["triage_resprate"], m["triage_o2sat"], temp_c,
                                  m["triage_sbp"], m["triage_heartrate"],
                                  pd.Series(1.0, index=m.index), NO_SUPP_O2)
m["score_MEWS"] = sc.score_mews(m["triage_sbp"], m["triage_heartrate"],
                                m["triage_resprate"], temp_c, ALERT_AVPU)
m["triage_sofa"] = sc.score_sofa_vitals(m["triage_MAP"], ALERT_GCS, m["triage_o2sat"])
m["score_ESRP1"] = np.nan  # placeholder, no definition was ever supplied (source project)

# ------------------------------------------------- ACS / PE (prefix, from cp_codes.py)
log("flagging ACS / PE (prefix match, cp_codes.py -- post pneumothorax/hemothorax audit fix)")


def hit(df, key, cond):
    v9 = df["icd_version"].to_numpy() == 9
    c = df["icd_code"]
    h = np.where(v9, c.str.startswith(tuple(CP_9[cond]), na=False),
                 c.str.startswith(tuple(CP_10[cond]), na=False))
    return df.loc[h, key].unique()


for cond in ("acs", "pe"):
    col = f"outcome_cp_{cond}"
    m[col] = (m["stay_id"].isin(hit(eddx, "stay_id", cond))
              | m["hadm_id"].isin(hit(hdx, "hadm_id", cond))).astype("int8")

# ------------------------------------------------- chest-pain cohort predictor flags
log("computing cohort_chest_pain / dx_chest_pain predictors")
CP9, CP10 = ["7865"], ["R071", "R072", "R078", "R079"]
v9 = eddx["icd_version"].to_numpy() == 9
edcp = eddx.loc[np.where(v9, eddx["icd_code"].str.startswith(tuple(CP9), na=False),
                         eddx["icd_code"].str.startswith(tuple(CP10), na=False)), "stay_id"].unique()
m["dx_chest_pain"] = m["stay_id"].isin(edcp).astype("int8")
m["cohort_chest_pain"] = ((m["chiefcom_chest_pain"] == 1) | (m["dx_chest_pain"] == 1)).astype("int8")

# ------------------------------------------------- ventilation timing (predictor only)
log("timing invasive ventilation (kept as a predictor; no ARDS-timed outcome in this benchmark)")
pe = pd.read_csv(DATA / "icu/procedureevents.csv.gz",
                 usecols=["hadm_id", "itemid", "starttime"], parse_dates=["starttime"])
vent = pe[pe["itemid"].isin([225792, 224385])]  # Invasive Ventilation, Intubation
first_vent = vent.groupby("hadm_id")["starttime"].min().rename("first_vent_time")
m = m.merge(first_vent, on="hadm_id", how="left")
m["time_to_invasive_vent_hours"] = (m["first_vent_time"] - m["outtime"]).dt.total_seconds() / 3600

m.to_pickle(OUT / "_stage5.pkl")
log(f"stage 5 complete -> {m.shape}")
for name in list(EXACT_OUTCOMES) + []:
    c = f"outcome_{name}"
    print(f"  {c:16s} n={m[c].sum():7,d}  rate={m[c].mean()*100:5.2f}%")
for c in ["outcome_cp_acs", "outcome_cp_pe"]:
    print(f"  {c:16s} n={m[c].sum():7,d}  rate={m[c].mean()*100:5.2f}%")
