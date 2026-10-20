"""Stage 2: fix ICU linkage, add triage / vital-sign / medication / chief-complaint features."""
import sys, time, pathlib, re
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import scores as sc

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "derived"
_t0 = time.time()
def log(msg): print(f"[{time.time()-_t0:7.1f}s] {msg}", flush=True)

m = pd.read_pickle(OUT / "_stage1.pkl")

# -------------------------------------------------- ICU transfer, same admission
log("re-linking ICU transfer to the SAME hospital admission")
icu = pd.read_csv(ROOT / "data/icu/icustays.csv.gz", parse_dates=["intime"])
first_icu = (icu.sort_values("intime").groupby("hadm_id")["intime"].first()
             .rename("intime_icu_after_ed"))
m = m.drop(columns=["intime_icu_after_ed"]).merge(first_icu, on="hadm_id", how="left")
# An ICU intime before ED departure means the ICU stay is not a transfer *from*
# this ED visit; null it rather than emit a negative interval.
m.loc[m["intime_icu_after_ed"] < m["outtime"], "intime_icu_after_ed"] = pd.NaT
m["time_to_icu_transfer"] = m["intime_icu_after_ed"] - m["outtime"]
m["time_to_icu_transfer_hours"] = m["time_to_icu_transfer"].dt.total_seconds() / 3600
m["outcome_icu_transfer_12h"] = ((m["time_to_icu_transfer_hours"] >= 0)
                                 & (m["time_to_icu_transfer_hours"] <= 12)).astype("int8")
m["outcome_critical"] = ((m["outcome_icu_transfer_12h"] == 1)
                         | (m["outcome_inhospital_mortality"] == 1)).astype("int8")
m["has_hosp_link"] = m["anchor_age"].notna().astype("int8")

# -------------------------------------------------- triage
log("merging triage vitals")
tri = pd.read_csv(ROOT / "data/ed/triage.csv.gz")
VITALS = ["temperature", "heartrate", "resprate", "o2sat", "sbp", "dbp"]
tri = sc.clean_vitals(tri, VITALS)
tri["pain"] = sc.parse_pain(tri["pain"])
tri = tri.rename(columns={c: f"triage_{c}" for c in VITALS + ["pain", "acuity"]})
m = m.merge(tri.drop(columns=["subject_id"]), on="stay_id", how="left")
m["triage_MAP"] = (m["triage_sbp"] + 2 * m["triage_dbp"]) / 3

# -------------------------------------------------- serial ED vitals
log("first / last recorded ED vital signs")
vs = pd.read_csv(ROOT / "data/ed/vitalsign.csv.gz", parse_dates=["charttime"])
vs = sc.clean_vitals(vs, VITALS)
vs["pain"] = sc.parse_pain(vs["pain"])
vs = vs.sort_values(["stay_id", "charttime"])
# first/last NON-NULL observation per vital, not first/last row - a row can be
# charted with only some vitals filled in.
agg = vs.groupby("stay_id")[VITALS + ["pain"]]
first = agg.first().add_suffix("_first")
last = agg.last().add_suffix("_last")
obs = first.join(last)
ren = {}
for c in VITALS:
    ren[f"{c}_first"] = f"triage_{c}_obs_first"
    ren[f"{c}_last"] = f"triage_{c}_obs_last"
ren["pain_first"], ren["pain_last"] = "triage_pain_first", "triage_pain_last"
m = m.merge(obs.rename(columns=ren), on="stay_id", how="left")

# Glasgow Coma Scale is absent from every MIMIC-IV-ED table (verified: no
# gcs/glasgow/coma/avpu column in edstays, triage, vitalsign, pyxis, medrecon,
# diagnosis). Columns are emitted as all-NaN to preserve the requested schema.
m["triage_gcs_first"] = np.nan
m["triage_gcs_last"] = np.nan

# -------------------------------------------------- medications
log("counting ED medications")
py = pd.read_csv(ROOT / "data/ed/pyxis.csv.gz", usecols=["stay_id", "name"])
mr = pd.read_csv(ROOT / "data/ed/medrecon.csv.gz", usecols=["stay_id", "name"])
m["n_med"] = m["stay_id"].map(py.groupby("stay_id")["name"].nunique()).fillna(0).astype("int16")
m["n_medrecon"] = m["stay_id"].map(mr.groupby("stay_id")["name"].nunique()).fillna(0).astype("int16")

# -------------------------------------------------- chief complaint
log("flagging chief complaints")
CC = {
    "chest_pain":          r"chest\s*pain|cp\b",
    "abdominal_pain":      r"abd(ominal)?\s*pain",
    "headache":            r"headache|h/?a\b",
    "shortness_of_breath": r"\bsob\b|short(ness)?\s*of\s*breath|dyspnea|difficulty breathing",
    "back_pain":           r"back\s*pain",
    "cough":               r"cough",
    "nausea_vomiting":     r"nausea|vomit|\bn/?v\b|emesis",
    "fever_chills":        r"fever|chills|febrile",
    "syncope":             r"syncope|fainting|passed out|loc\b",
    "dizziness":           r"dizz|vertigo|light\s*head",
}
cctext = m["chiefcomplaint"].astype("string").fillna("").str.lower()
for name, pat in CC.items():
    m[f"chiefcom_{name}"] = cctext.str.contains(pat, regex=True, na=False).astype("int8")

m.to_pickle(OUT / "_stage2.pkl")
log(f"stage 2 complete -> {m.shape}")
