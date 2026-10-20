"""Chest-pain differential: ICD-9-CM and ICD-10-CM code sets (dot-less prefixes).

Each entry is matched by prefix against MIMIC's dot-less icd_code. Every prefix was
verified against hosp/d_icd_diagnoses to confirm it exists and that its long_title
denotes the intended condition (see scripts/verify_cp_codes.py).
"""
CP_9 = {
  "stemi":            ["4100","4101","4102","4103","4104","4105","4106","4108"],
  "nstemi":           ["4107"],
  "unstable_angina":  ["4111"],
  "mi_any":           ["410"],
  "acs":              ["410","4111","4118"],
  "pe":               ["4151","6732"],
  "aortic_dissection":["4410"],
  "pneumonia":        ["480","481","482","483","484","485","486","4870","5070"],
  "ards":             ["51882"],   # matches icd_codes.csv; ICD-9 has no dedicated ARDS code
  # 8602/8603 ("traumatic hemothorax ... without/with open wound") deliberately
  # excluded: ICD-9-CM groups hemothorax into the same category-860 block as
  # pneumothorax for historical reasons, but it denotes blood in the pleural
  # space with no air at all. Verified against hosp/diagnoses_icd: 136 of 139
  # admissions coded 8602/8603 carry no air-pneumothorax code (8600/01/04/05)
  # whatsoever, and would otherwise be a false-positive pneumothorax outcome.
  "pneumothorax":     ["5120","5121","5128","8600","8601","8604","8605"],
  "pneumomediastinum":["5181"],    # 9587 excluded: subcutaneous, not mediastinal, emphysema
  "pericardial":      ["4200","4209","4230","4231","4232","4233","4238","4239","393","4110"],
  "tamponade":        ["4233"],
  "esoph_rupture":    ["5304","86222","86232"],
}
CP_10 = {
  "stemi":            ["I210","I211","I212","I213"],
  "nstemi":           ["I214"],
  "unstable_angina":  ["I200"],
  "mi_any":           ["I21","I22"],
  "acs":              ["I21","I22","I200","I248","I249"],
  "pe":               ["I26","O882"],
  "aortic_dissection":["I710"],
  "pneumonia":        ["J12","J13","J14","J15","J16","J17","J18","J690",
                       "J09X1","J1000","J1008","J1100","J1108"],
  "ards":             ["J80"],     # J95x is post-surgical insufficiency, not ARDS
  "pneumothorax":     ["J930","J9311","J9312","J9381","J9382","J9383","J939",
                       "J95811","S270"],
  "pneumomediastinum":["J982"],
  "pericardial":      ["I30","I31","I32","I092","I241"],
  "tamponade":        ["I314"],
  "esoph_rupture":    ["K223","S27813"],  # perforation + laceration only, not contusion
}
LABELS = {
  "stemi": "ST-elevation myocardial infarction",
  "nstemi": "Non-ST-elevation myocardial infarction",
  "unstable_angina": "Unstable angina",
  "mi_any": "Myocardial infarction, any type",
  "acs": "Acute coronary syndrome (composite)",
  "pe": "Pulmonary embolism",
  "aortic_dissection": "Aortic dissection",
  "pneumonia": "Pneumonia",
  "ards": "Acute respiratory distress syndrome",
  "pneumothorax": "Pneumothorax",
  "pneumomediastinum": "Pneumomediastinum",
  "pericardial": "Pericardial disease (pericarditis, effusion, tamponade)",
  "tamponade": "Cardiac tamponade",
  "esoph_rupture": "Esophageal rupture / perforation",
}
assert set(CP_9) == set(CP_10) == set(LABELS)
