"""Stage 7: order columns, write the benchmark dataset and column dictionary.

Column order is IDENTICAL to the source project's master_ed.csv for every
predictor column; the outcome block is reduced to exactly 5 columns:
outcome_sepsis, outcome_hf_acute, outcome_stroke, outcome_cp_acs, outcome_cp_pe.
"""
import sys, time, pathlib
import pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from codemaps import CCI_CATS, ECI_CATS

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "derived"
_t0 = time.time()
def log(m): print(f"[{time.time()-_t0:7.1f}s] {m}", flush=True)

m = pd.read_pickle(OUT / "_stage6.pkl").copy()
VIT = ["temperature", "heartrate", "resprate", "o2sat", "sbp", "dbp"]
cols = set(m.columns)

BENCHMARK_OUTCOMES = ["outcome_sepsis", "outcome_hf_acute", "outcome_stroke",
                      "outcome_cp_acs", "outcome_cp_pe"]

ORDER = [c for c in (
    ["subject_id", "hadm_id", "stay_id", "intime", "outtime",
     "anchor_age", "gender", "anchor_year", "anchor_year_group", "dod",
     "admittime", "dischtime", "deathtime", "ethnicity", "edregtime", "edouttime",
     "insurance", "in_year", "age",
     "cohort_chest_pain", "dx_chest_pain",
     "outcome_inhospital_mortality", "time_death_discharge",
     "ed_death", "time_death_ed_depart", "time_death_ed_arrival",
     "before_ed_mortality", "ed_los", "ed_los_hours",
     "intime_icu_after_ed", "time_to_icu_transfer", "time_to_icu_transfer_hours",
     "outcome_icu_transfer_12h", "outcome_hospitalization", "outcome_critical",
     "first_vent_time", "time_to_invasive_vent_hours",
     "n_ed_30d", "n_ed_90d", "n_ed_365d",
     "next_ed_visit_time", "next_ed_visit_time_diff", "next_ed_visit_time_diff_days",
     "outcome_ed_revisit_3d",
     "n_hosp_30d", "n_hosp_90d", "n_hosp_365d", "n_icu_30d", "n_icu_90d", "n_icu_365d"]
    + [f"triage_{v}" for v in VIT] + [f"triage_{v}_missing" for v in VIT]
    + [f"triage_{v}_obs_first" for v in VIT] + [f"triage_{v}_obs_first_missing" for v in VIT]
    + [f"triage_{v}_obs_last" for v in VIT] + [f"triage_{v}_obs_last_missing" for v in VIT]
    + ["triage_pain", "triage_pain_missing", "triage_pain_first", "triage_pain_first_missing",
       "triage_pain_last", "triage_pain_last_missing",
       "triage_gcs_first", "triage_gcs_last",
       "triage_acuity", "triage_MAP", "triage_MAP_missing", "n_vitals_imputed", "triage_sofa",
       "chiefcomplaint"]
    + [f"chiefcom_{c}" for c in ["chest_pain", "abdominal_pain", "headache",
       "shortness_of_breath", "back_pain", "cough", "nausea_vomiting", "fever_chills",
       "syncope", "dizziness"]]
    + [f"cci_{c}" for c in CCI_CATS] + [f"eci_{c}" for c in ECI_CATS]
    + ["n_med", "n_medrecon",
       "score_CCI", "score_CART", "score_REMS", "score_NEWS", "score_NEWS2",
       "score_MEWS", "score_ESRP1",
       "dxes", "dxes_code"]
    + BENCHMARK_OUTCOMES
    + ["has_labs_ed", "n_labs_ed", "n_lab_assays_ed", "n_labs_before_ed_departure",
       "time_to_first_lab_min"]
    + sorted(c for c in cols if c.startswith("lab_"))
    + ["disposition", "arrival_transport", "admission_type", "discharge_location",
       "hospital_expire_flag", "has_hosp_link", "n_dx_5y",
       "dx_name_5y", "dx_code_5y", "dx_code_type_5y", "dx_time_5y"]
) if c in cols]

extra = [c for c in m.columns if c not in ORDER]
if extra:
    log(f"appending {len(extra)} unlisted columns: {extra}")
    ORDER += extra
m = m[ORDER]
log(f"final shape {m.shape}")

m.to_pickle(OUT / "master_benchmark5.pkl")
log("writing master_benchmark5.csv.gz")
m.to_csv(OUT / "master_benchmark5.csv.gz", index=False, compression="gzip")
log("writing master_benchmark5_slim.csv.gz")
m.drop(columns=[c for c in ["dx_name_5y", "dx_code_5y", "dx_code_type_5y", "dx_time_5y"]
                if c in m.columns]).to_csv(OUT / "master_benchmark5_slim.csv.gz", index=False,
                                           compression="gzip")

log("writing column dictionary")
rows = []
for c in m.columns:
    s = m[c]; k = s.dtype.kind
    rows.append({"column": c, "dtype": str(s.dtype), "n_missing": int(s.isna().sum()),
                "pct_missing": round(float(s.isna().mean()) * 100, 2),
                "min": s.min() if k in "iuf" else "",
                "median": round(float(s.median()), 3) if k in "iuf" and s.notna().any() else "",
                "max": s.max() if k in "iuf" else ""})
pd.DataFrame(rows).to_csv(OUT / "column_dictionary.csv", index=False)
log("done")
print(f"\nFinal dataset: {m.shape[0]:,} rows x {m.shape[1]} columns")
print("\nOutcome prevalence:")
for c in BENCHMARK_OUTCOMES:
    print(f"  {c:16s} n={m[c].sum():7,d}  rate={m[c].mean()*100:5.2f}%")
