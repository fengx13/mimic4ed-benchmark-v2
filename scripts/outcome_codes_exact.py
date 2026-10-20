"""Sepsis, acute heart failure (AHF) and stroke: EXACT-match ICD-9-CM /
ICD-10-CM code lists, copied verbatim from the source project's icd_codes.csv
(disease rows dx_sepsis, dx_hf_acute, dx_stroke). ACS and PE are NOT here —
they use cp_codes.py's PREFIX-match definitions (outcome_cp_acs / outcome_cp_pe),
which were specifically audited code-by-code against hosp/d_icd_diagnoses.csv.gz
(see README) and are the broader, validated composite definitions.

Membership is EXACT, not prefix: these lists enumerate every code intended.
"""

SEPSIS_9 = ["0380", "03810", "03811", "03812", "03819", "0382", "0383", "03840",
            "03841", "03842", "03843", "03844", "03849", "0388", "0389",
            "67020", "67022", "67024", "78552", "99591", "99592"]
SEPSIS_10 = ["A40", "A400", "A401", "A403", "A408", "A409", "A41", "A410",
             "A4101", "A4102", "A411", "A412", "A413", "A414", "A415", "A4150",
             "A4151", "A4152", "A4153", "A4159", "A418", "A4181", "A4189",
             "A419", "O0337", "O0387", "O0487", "O85", "R6520", "R6521"]
# Septicaemia, sepsis, severe sepsis, septic shock, and obstetric/abortion-related sepsis.

AHF_9 = ["42821", "42823", "42831", "42833", "42841", "42843"]
AHF_10 = ["I5021", "I5023", "I5031", "I5033", "I5041", "I5043", "I50811", "I50813"]
# Codes explicitly labelled acute or acute-on-chronic heart failure only
# (chronic/unspecified heart failure codes are deliberately excluded).

STROKE_9 = ["430", "431", "4320", "4321", "4329", "43301", "43311", "43321",
            "43331", "43381", "43391", "43401", "43411", "43491", "436"]
STROKE_10 = [
    "I60", "I600", "I6000", "I6001", "I6002", "I601", "I6010", "I6011", "I6012",
    "I602", "I6020", "I6021", "I6022", "I603", "I6030", "I6031", "I6032", "I604",
    "I605", "I6050", "I6051", "I6052", "I606", "I607", "I608", "I609", "I61",
    "I610", "I611", "I612", "I613", "I614", "I615", "I616", "I618", "I619", "I62",
    "I620", "I6200", "I6201", "I6202", "I6203", "I621", "I629", "I63", "I630",
    "I6300", "I6301", "I63011", "I63012", "I63013", "I63019", "I6302", "I6303",
    "I63031", "I63032", "I63033", "I63039", "I6309", "I631", "I6310", "I6311",
    "I63111", "I63112", "I63113", "I63119", "I6312", "I6313", "I63131", "I63132",
    "I63133", "I63139", "I6319", "I632", "I6320", "I6321", "I63211", "I63212",
    "I63213", "I63219", "I6322", "I6323", "I63231", "I63232", "I63233", "I63239",
    "I6329", "I633", "I6330", "I6331", "I63311", "I63312", "I63313", "I63319",
    "I6332", "I63321", "I63322", "I63323", "I63329", "I6333", "I63331", "I63332",
    "I63333", "I63339", "I6334", "I63341", "I63342", "I63343", "I63349", "I6339",
    "I634", "I6340", "I6341", "I63411", "I63412", "I63413", "I63419", "I6342",
    "I63421", "I63422", "I63423", "I63429", "I6343", "I63431", "I63432", "I63433",
    "I63439", "I6344", "I63441", "I63442", "I63443", "I63449", "I6349", "I635",
    "I6350", "I6351", "I63511", "I63512", "I63513", "I63519", "I6352", "I63521",
    "I63522", "I63523", "I63529", "I6353", "I63531", "I63532", "I63533", "I63539",
    "I6354", "I63541", "I63542", "I63543", "I63549", "I6359", "I636", "I638",
    "I6381", "I6389", "I639",
]
# Non-traumatic cerebrovascular events: subarachnoid haemorrhage, intracerebral
# haemorrhage, and cerebral infarction. Verified: all 161 ICD-10 codes matched a
# keyword check against "infarct|hemorrhage|occlusion|cerebrovascular|stenosis".

EXACT_OUTCOMES = {
    "sepsis":   ("Sepsis", SEPSIS_9, SEPSIS_10),
    "hf_acute": ("Acute heart failure (AHF)", AHF_9, AHF_10),
    "stroke":   ("Stroke (non-traumatic)", STROKE_9, STROKE_10),
}
