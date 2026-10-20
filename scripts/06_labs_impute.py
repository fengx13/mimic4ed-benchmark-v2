"""Stage 9: pivot ED-window labs to wide columns, median-impute vitals, rescore, assemble."""
import sys, re, time, pathlib, json
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import scores as sc

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "derived"
_t0 = time.time()
def log(m): print(f"[{time.time()-_t0:7.1f}s] {m}", flush=True)

m = pd.read_pickle(OUT / "_stage5.pkl").copy()
N = len(m)
labs = pd.read_pickle(OUT / "_ed_labs.pkl")

# ------------------------------------------------------------------ lab panel
log("selecting lab panel")
labs = labs[labs["valuenum"].notna()]
cov = labs.groupby("itemid")["stay_id"].nunique()
# Chest-pain-critical assays are kept regardless of how rarely they are ordered:
# a rarely-drawn troponin or D-dimer is exactly the signal this project needs.
CP_ASSAYS = [51003, 50915, 51196, 50963, 50910, 50911, 50908, 50813, 50820, 50818,
             50821, 50802, 50817, 50805, 50814, 50889, 50954, 51214, 50808, 50852,
             51000, 50907, 50904, 50905, 50903, 50964, 50956, 50867]
panel = sorted(set(cov[cov >= 0.05 * N].index) | (set(CP_ASSAYS) & set(cov.index)))
labs = labs[labs["itemid"].isin(panel)]
log(f"panel: {len(panel)} assays, {len(labs):,} numeric results")

meta = (labs.groupby("itemid")
        .agg(label=("label", "first"), fluid=("fluid", "first"),
             category=("category", "first"), n_results=("valuenum", "size"),
             n_stays=("stay_id", "nunique")).reset_index())


def slug(lbl, fluid):
    s = re.sub(r"[^a-z0-9]+", "_", f"{lbl}_{fluid}".lower()).strip("_")
    return re.sub(r"_+", "_", s)


meta["slug"] = [slug(l, f) for l, f in zip(meta["label"], meta["fluid"])]
# Some assays share a label across itemids (e.g. D-Dimer 50915 and 51196 are run by
# different methods and may not share units). Disambiguate rather than merge them.
dup = meta["slug"].duplicated(keep=False)
meta.loc[dup, "slug"] = meta.loc[dup, "slug"] + "_" + meta.loc[dup, "itemid"].astype(str)
assert meta["slug"].is_unique
meta.to_csv(OUT / "lab_panel_dictionary.csv", index=False)
name = dict(zip(meta["itemid"], meta["slug"]))

log("pivoting labs to one static value per assay")
labs = labs.sort_values(["stay_id", "itemid", "charttime"])
# One static result per assay per ED stay: if the assay was resulted at all during
# the capture window, its earliest result is that stay's value. No first/min/max
# series and no per-assay timing - the result is a fixed attribute of the visit.
agg = labs.groupby(["stay_id", "itemid"])["valuenum"].first().reset_index()
agg["col"] = agg["itemid"].map(name)
W = agg.pivot(index="stay_id", columns="col", values="valuenum")
W.columns = [f"lab_{c}" for c in W.columns]
W = W.astype("float32")
log(f"wide lab block: {W.shape}")
m = m.merge(W, on="stay_id", how="left")

allb = pd.read_pickle(OUT / "_ed_labs.pkl")
m["n_labs_ed"] = m["stay_id"].map(allb.groupby("stay_id").size()).fillna(0).astype("int32")
m["n_lab_assays_ed"] = m["stay_id"].map(allb.groupby("stay_id")["itemid"].nunique()).fillna(0).astype("int16")
m["time_to_first_lab_min"] = m["stay_id"].map(allb.groupby("stay_id")["offset_min"].min()).astype("float32")
m["has_labs_ed"] = (m["n_labs_ed"] > 0).astype("int8")
ed_only = allb[allb["in_ed_window"] == 1]
m["n_labs_before_ed_departure"] = m["stay_id"].map(ed_only.groupby("stay_id").size()).fillna(0).astype("int32")

# ------------------------------------------------------------- vital imputation
log("median-imputing vitals")
VIT = ["temperature", "heartrate", "resprate", "o2sat", "sbp", "dbp"]
VCOLS = ([f"triage_{v}" for v in VIT]
         + [f"triage_{v}_obs_first" for v in VIT]
         + [f"triage_{v}_obs_last" for v in VIT]
         + ["triage_pain", "triage_pain_first", "triage_pain_last"])
medians = {}
for c in VCOLS:
    med = float(m[c].median())
    medians[c] = med
    m[f"{c}_missing"] = m[c].isna().astype("int8")
    m[c] = m[c].fillna(med)
(OUT / "imputation_medians.json").write_text(json.dumps(medians, indent=1))
m["n_vitals_imputed"] = m[[f"{c}_missing" for c in VCOLS]].sum(axis=1).astype("int8")
# MAP is recomputed from the imputed pressures rather than imputed on its own,
# so MAP stays consistent with the SBP/DBP in the same row.
m["triage_MAP_missing"] = ((m["triage_sbp_missing"] == 1) | (m["triage_dbp_missing"] == 1)).astype("int8")
m["triage_MAP"] = (m["triage_sbp"] + 2 * m["triage_dbp"]) / 3
# triage_acuity is deliberately NOT imputed: its missingness is informative
# (46.8% critical rate, higher than ESI 1).

log("recomputing scores on imputed vitals")
ALERT_GCS = pd.Series(15.0, index=m.index)
temp_c = sc.f_to_c(m["triage_temperature"])
m["score_CART"] = sc.score_cart(m["triage_resprate"], m["triage_heartrate"], m["triage_dbp"], m["age"])
m["score_REMS"] = sc.score_rems(m["age"], m["triage_MAP"], m["triage_heartrate"],
                                m["triage_resprate"], m["triage_o2sat"], ALERT_GCS)
m["score_NEWS"] = sc.score_news(m["triage_resprate"], m["triage_o2sat"], temp_c, m["triage_sbp"],
                                m["triage_heartrate"], pd.Series(1.0, index=m.index),
                                pd.Series(0.0, index=m.index))
m["score_NEWS2"] = m["score_NEWS"]
m["score_MEWS"] = sc.score_mews(m["triage_sbp"], m["triage_heartrate"], m["triage_resprate"],
                                temp_c, pd.Series(0.0, index=m.index))
m["triage_sofa"] = sc.score_sofa_vitals(m["triage_MAP"], ALERT_GCS, m["triage_o2sat"])

m = m.copy()
m.to_pickle(OUT / "_stage6.pkl")
log(f"stage 9 complete -> {m.shape}")
print("\nlab coverage of the chest-pain panel (% of all stays):")
for c in ["lab_troponin_t_blood", "lab_d_dimer_blood_50915", "lab_ntprobnp_blood",
          "lab_creatine_kinase_mb_isoenzyme_blood", "lab_lactate_blood",
          "lab_creatinine_blood", "lab_hematocrit_blood"]:
    if c in m.columns:
        print(f"  {c:48s} {m[c].notna().mean()*100:5.1f}%")
print(f"\nvitals imputed: {(m.n_vitals_imputed>0).mean()*100:.1f}% of stays had >=1 vital imputed")
print("score missingness after imputation:",
      {c: round(float(m[c].isna().mean()*100), 2) for c in
       ["score_CART", "score_REMS", "score_NEWS", "score_MEWS", "triage_sofa"]})
