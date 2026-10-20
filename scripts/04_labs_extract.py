"""Stage 7: extract every lab resulted during an ED stay window -> long format.

Scans hosp/labevents.csv.gz (2.4 GB gz, ~158M rows) in chunks and keeps only rows
whose subject has an ED stay whose [intime, outtime] window contains the charttime.
Window membership is resolved with one searchsorted per chunk using the same exact
int64 (subject, time) key used elsewhere in the pipeline.
"""
import time, pathlib
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "derived"
SEC = 10_000_000_000
CHUNK = 4_000_000
_t0 = time.time()
def log(m): print(f"[{time.time()-_t0:7.1f}s] {m}", flush=True)

ed = pd.read_csv(ROOT / "data/ed/edstays.csv.gz", usecols=["subject_id", "stay_id", "intime", "outtime"],
                 parse_dates=["intime", "outtime"])
ed = ed.sort_values(["subject_id", "intime"]).reset_index(drop=True)
ed_subj = ed["subject_id"].to_numpy()
ed_in = ed["intime"].to_numpy().astype("datetime64[s]").astype("int64")
ed_out = ed["outtime"].to_numpy().astype("datetime64[s]").astype("int64")
# Capture window runs to ED departure OR 24 h from arrival, whichever is later: a
# troponin drawn in the ED frequently results after the patient moves to a ward bed,
# and clipping at outtime would discard it. Every row keeps offset_min and an
# in_ed_window flag so triage-available labs can still be isolated exactly.
ed_cap = np.maximum(ed_out, ed_in + 24 * 3600)
ed_stay = ed["stay_id"].to_numpy()
ed_key = ed_subj.astype("int64") * SEC + ed_in
subj_set = np.unique(ed_subj)
log(f"{len(ed):,} ED windows over {len(subj_set):,} subjects")

USE = ["subject_id", "itemid", "charttime", "value", "valuenum", "valueuom", "flag",
       "ref_range_lower", "ref_range_upper"]
kept, n_seen, n_chunk = [], 0, 0
reader = pd.read_csv(ROOT / "data/hosp/labevents.csv.gz", usecols=USE, chunksize=CHUNK,
                     parse_dates=["charttime"])
for ch in reader:
    n_chunk += 1
    n_seen += len(ch)
    ch = ch[ch["subject_id"].isin(subj_set)]
    if ch.empty:
        continue
    t = ch["charttime"].to_numpy().astype("datetime64[s]").astype("int64")
    k = ch["subject_id"].to_numpy().astype("int64") * SEC + t
    # index of the latest ED stay whose intime is <= this charttime
    pos = np.searchsorted(ed_key, k, "right") - 1
    ok = pos >= 0
    pos_c = np.clip(pos, 0, len(ed_key) - 1)
    ok &= (ed_subj[pos_c] == ch["subject_id"].to_numpy()) & (t <= ed_cap[pos_c])
    if not ok.any():
        continue
    sub = ch.loc[ok].copy()
    p = pos_c[ok]
    sub["stay_id"] = ed_stay[p]
    sub["offset_min"] = (t[ok] - ed_in[p]) / 60.0
    sub["in_ed_window"] = (t[ok] <= ed_out[p]).astype("int8")
    kept.append(sub)
    if n_chunk % 5 == 0:
        log(f"chunk {n_chunk}: {n_seen:,} rows scanned, {sum(len(k) for k in kept):,} kept")

labs = pd.concat(kept, ignore_index=True)
log(f"DONE scanning {n_seen:,} rows -> {len(labs):,} ED-window lab results")
d = pd.read_csv(ROOT / "data/hosp/d_labitems.csv.gz")
labs = labs.merge(d, on="itemid", how="left")
labs = labs.sort_values(["stay_id", "itemid", "charttime"])
labs.to_csv(OUT / "ed_labs_long.csv.gz", index=False, compression="gzip")
labs.to_pickle(OUT / "_ed_labs.pkl")
log(f"wrote ed_labs_long.csv.gz  ({labs.stay_id.nunique():,} stays, "
    f"{labs.itemid.nunique():,} distinct assays)")
