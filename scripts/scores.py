"""Vital-sign cleaning and bedside clinical severity scores.

All score functions are vectorised over pandas Series and return float Series
with NaN wherever a required input is missing (no silent zero-filling, so a
missing vital never masquerades as a normal one).
"""
import numpy as np
import pandas as pd

# Plausible physiologic ranges. Values outside are set to NaN, not clipped:
# MIMIC-IV-ED triage contains transcription artefacts (dbp up to 661672,
# temperature up to 986) that would otherwise dominate every score.
RANGES = {
    "temperature": (70.0, 115.0),   # Fahrenheit, as recorded
    "heartrate":   (20.0, 300.0),
    "resprate":    (2.0, 80.0),
    "o2sat":       (30.0, 100.0),
    "sbp":         (30.0, 350.0),
    "dbp":         (10.0, 250.0),
}


def clean_vitals(df, cols=None):
    """NaN-out physiologically impossible values in place-safe fashion."""
    out = df.copy()
    for c in (cols or RANGES):
        if c in out.columns:
            lo, hi = RANGES[c]
            v = pd.to_numeric(out[c], errors="coerce")
            out[c] = v.where((v >= lo) & (v <= hi))
    return out


def parse_pain(s):
    """ED pain is free text ('7', '10', 'critical', 'uta', '13'). Keep 0-10 only."""
    v = pd.to_numeric(pd.Series(s).astype("string").str.strip(), errors="coerce")
    return v.where((v >= 0) & (v <= 10))


def f_to_c(temp_f):
    return (pd.to_numeric(temp_f, errors="coerce") - 32.0) * 5.0 / 9.0


def _band(x, cuts, points):
    """Right-open banding: points[i] applies where cuts[i] <= x < cuts[i+1]."""
    x = pd.to_numeric(x, errors="coerce")
    out = pd.Series(np.nan, index=x.index, dtype="float64")
    for i, p in enumerate(points):
        lo, hi = cuts[i], cuts[i + 1]
        out = out.mask((x >= lo) & (x < hi), float(p))
    return out.where(x.notna())


INF = np.inf


def score_cart(resprate, heartrate, dbp, age):
    """Cardiac Arrest Risk Triage (Churpek 2012). Range 0-57."""
    rr = _band(resprate, [-INF, 21, 24, 26, 30, INF], [0, 8, 12, 15, 22])
    hr = _band(heartrate, [-INF, 110, 140, INF], [0, 4, 13])
    bp = _band(dbp, [-INF, 35, 40, 50, INF], [13, 6, 4, 0])
    ag = _band(age, [-INF, 55, 70, INF], [0, 4, 9])
    return rr + hr + bp + ag


def score_rems(age, mapp, heartrate, resprate, o2sat, gcs):
    """Rapid Emergency Medicine Score (Olsson 2004). Range 0-26."""
    ag = _band(age, [-INF, 45, 55, 65, 75, INF], [0, 2, 3, 5, 6])
    mp = _band(mapp, [-INF, 50, 70, 110, 130, 160, INF], [4, 2, 0, 2, 3, 4])
    hr = _band(heartrate, [-INF, 40, 55, 70, 110, 140, 180, INF], [4, 3, 2, 0, 2, 3, 4])
    rr = _band(resprate, [-INF, 6, 10, 12, 25, 35, 50, INF], [4, 2, 1, 0, 1, 3, 4])
    ox = _band(o2sat, [-INF, 75, 86, 90, INF], [4, 3, 1, 0])
    gc = _band(gcs, [-INF, 5, 8, 11, 15, INF], [4, 3, 2, 1, 0])
    return ag + mp + hr + rr + ox + gc


def score_news(resprate, o2sat, temp_c, sbp, heartrate, alert, supp_o2):
    """National Early Warning Score (RCP 2012). Range 0-20."""
    rr = _band(resprate, [-INF, 9, 12, 21, 25, INF], [3, 1, 0, 2, 3])
    ox = _band(o2sat, [-INF, 92, 94, 96, INF], [3, 2, 1, 0])
    tp = _band(temp_c, [-INF, 35.1, 36.1, 38.1, 39.1, INF], [3, 1, 0, 1, 2])
    bp = _band(sbp, [-INF, 91, 101, 111, 220, INF], [3, 2, 1, 0, 3])
    hr = _band(heartrate, [-INF, 41, 51, 91, 111, 131, INF], [3, 1, 0, 1, 2, 3])
    return rr + ox + tp + bp + hr + (1 - alert) * 3 + supp_o2 * 2


def score_news2(resprate, o2sat, temp_c, sbp, heartrate, alert, supp_o2):
    """NEWS2 (RCP 2017), SpO2 Scale 1.

    Scale 2 applies only to patients with known hypercapnic respiratory failure
    and a prescribed 88-92% target; that flag does not exist in MIMIC-IV-ED, so
    Scale 1 is used for everyone. With supplemental-oxygen and consciousness
    data also absent, this is numerically identical to NEWS here - the column is
    retained for schema compatibility, not because it adds information.
    """
    return score_news(resprate, o2sat, temp_c, sbp, heartrate, alert, supp_o2)


def score_mews(sbp, heartrate, resprate, temp_c, avpu):
    """Modified Early Warning Score (Subbe 2001). Range 0-14."""
    bp = _band(sbp, [-INF, 71, 81, 101, 200, INF], [3, 2, 1, 0, 2])
    hr = _band(heartrate, [-INF, 40, 51, 101, 111, 130, INF], [2, 1, 0, 1, 2, 3])
    rr = _band(resprate, [-INF, 9, 15, 21, 30, INF], [2, 0, 1, 2, 3])
    tp = _band(temp_c, [-INF, 35, 38.5, INF], [2, 0, 2])
    return bp + hr + rr + tp + avpu


def score_sofa_vitals(mapp, gcs, o2sat):
    """Vitals-only SOFA surrogate - NOT a true SOFA score.

    True SOFA needs PaO2/FiO2, platelets, bilirubin, creatinine and vasopressor
    dose, none of which exist in the ED module. This computes only the three
    domains derivable from triage vitals:
      - cardiovascular: MAP >= 70 -> 0, MAP < 70 -> 1 (vasopressor tiers 2-4
        are unreachable without infusion data)
      - CNS: from GCS
      - respiratory: SpO2/FiO2 surrogate for PaO2/FiO2, assuming room air
        (FiO2 = 0.21) because O2 delivery is not recorded. S/F is mapped onto
        the SOFA P/F cutpoints via Rice (2007): S/F = 64 + 0.84 * P/F, so
        P/F 400/300/200/100 correspond to S/F 400/316/232/148.
    Maximum attainable is 9, versus 24 for real SOFA. Interpret as an ordinal
    triage-severity proxy only.
    """
    cardio = pd.Series(np.where(pd.to_numeric(mapp, errors="coerce") < 70, 1.0, 0.0),
                       index=pd.Series(mapp).index).where(pd.notna(mapp))
    cns = _band(gcs, [-INF, 6, 10, 13, 15, INF], [4, 3, 2, 1, 0])
    sf = pd.to_numeric(o2sat, errors="coerce") / 0.21
    resp = _band(sf, [-INF, 148, 232, 316, 400, INF], [4, 3, 2, 1, 0])
    return cardio + cns + resp


def score_cci(flags, age):
    """Charlson Comorbidity Index with age points (Charlson 1987/1994).

    `flags` is a DataFrame of 0/1 columns named cci_<Category>. Severity
    hierarchy is applied: DM2 cancels DM1, Liver2 cancels Liver1, Cancer2
    cancels Cancer1.
    """
    from codemaps import CCI_WEIGHTS, CCI_HIERARCHY
    f = {c: flags[f"cci_{c}"].astype("float64") for c in CCI_WEIGHTS}
    for severe, mild in CCI_HIERARCHY:
        f[mild] = f[mild] * (1 - f[severe])
    total = sum(f[c] * w for c, w in CCI_WEIGHTS.items())
    age_pts = _band(age, [-INF, 50, 60, 70, 80, INF], [0, 1, 2, 3, 4]).fillna(0)
    return total + age_pts
