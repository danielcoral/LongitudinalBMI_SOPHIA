# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Project folder
projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

#-------------------------------------------------------------------------------
# Weight data

# Weight data at recruitment (instance 0)
wrec <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Weight_Recruitment_participant.csv"
  )
)
# Rename columns
renamer <- c(
  weight = "p21002_i0",
  height = "p50_i0",
  sex = "p31",
  date_entry = "p53_i0",
  age = "p21003_i0"
)
wrec <- wrec |>
  rename(all_of(renamer)) |>
  drop_na() |>
  filter(
    age >= 40,
    age <= 70,
    height >= 150,
    weight >= 50,
    weight <= 200,
    weight / ((height / 100) ^ 2) >= 20,
    weight / ((height / 100) ^ 2) <= 50
  )

# Weight data at instance 1
w1i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Weight_Instance1_participant.csv"
  )
)
# Rename columns
renamer <- c(
  height = "p50_i0",
  weight = "p21002_i1",
  age = "p21003_i1",
  date_visit = "p53_i1"
)
w1i <- w1i |>
  select(-p31) |>
  rename(all_of(renamer)) |>
  drop_na() |>
  filter(
    age >= 40,
    age <= 70,
    height >= 150,
    weight >= 50,
    weight <= 200,
    weight / ((height / 100) ^ 2) >= 20,
    weight / ((height / 100) ^ 2) <= 50
  )

head(w1i)

nrow(w1i)

# Weight data at instance 2
w2i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Weight_Instance2_participant.csv"
  )
)
# Rename columns
renamer <- c(
  height = "p50_i2",
  weight = "p21002_i2",
  age = "p21003_i2",
  date_visit = "p53_i2"
)
w2i <- w2i |>
  select(-p31) |>
  rename(all_of(renamer)) |>
  drop_na() |>
  filter(
    age >= 40,
    age <= 70,
    height >= 150,
    weight >= 50,
    weight <= 200,
    weight / ((height / 100) ^ 2) >= 20,
    weight / ((height / 100) ^ 2) <= 50
  )

# Weight data from GP records
wgpc <- read_tsv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "OnlyWeights_gp.tsv"
  ),
  col_types = "ncDnnn"
) |>
  select(
    eid = `participant$eid`,
    event_dt,
    weight = value1
  ) |>
  drop_na() |>
  filter(
    weight >= 50,
    weight <= 200
  ) |>
  inner_join(wrec[, c("eid", "height")]) |>
  filter(
    weight / ((height / 100) ^ 2) >= 20,
    weight / ((height / 100) ^ 2) <= 50
  ) |>
  select(-height)

#-------------------------------------------------------------------------------
# Longitudinal BMI datasets

# Until recruitment
bmilongrec <- wrec |>
  # Joining longitudinal data with data at recruitment
  transmute(
    eid,
    event_dt = date_entry,
    weight = weight
  ) |>
  drop_na |>
  bind_rows(wgpc) |>
  # Only including data up to recruitment
  inner_join(wrec, by = "eid", suffix = c("", ".ref")) |>
  transmute(
    eid, sex,
    time = as.numeric(event_dt - date_entry) / 365.25,
    bmi0 = weight.ref / ((height / 100)^2),
    bmi = weight / ((height / 100)^2),
    bmi = ifelse(time == 0 & bmi != bmi0, NA, bmi)
  ) |>
  filter(time >= -5, time <= 0) |>
  drop_na() |>
  mutate(
    percbase = 100 * bmi / bmi0,
    over6m = time <= -.5
  ) |>
  filter(percbase > 50, percbase < 150) |>
  group_by(eid) |>
  filter(any(over6m)) |>
  filter(n() > 1) |>
  ungroup() |>
  select(-c(percbase, over6m))
# Saving
write_tsv(
  bmilongrec,
  file.path(projfld, "data", "UKB", "processed_data", "bmilongrec.tsv")
)

# Until instance 1
bmilong1i <- bind_rows(
  w1i |>
    transmute(
      eid,
      event_dt = date_visit,
      weight = weight
    ),
  wrec |>
    transmute(
      eid,
      event_dt = date_entry,
      weight = weight
    ),
  wgpc
) |>
  drop_na |>
  # Only including data up to recruitment
  inner_join(w1i, by = "eid", suffix = c("", ".ref")) |>
  inner_join(wrec[, c("eid", "sex")], by = join_by(eid)) |>
  transmute(
    eid, sex,
    time = as.numeric(event_dt - date_visit) / 365.25,
    bmi0 = weight.ref / ((height / 100)^2),
    bmi = weight / ((height / 100)^2),
    bmi = ifelse(time == 0 & bmi != bmi0, NA, bmi)
  ) |>
  filter(time >= -5, time <= 0) |>
  drop_na() |>
  mutate(
    percbase = 100 * bmi / bmi0,
    over6m = time <= -.5
  ) |>
  filter(percbase > 50, percbase < 150) |>
  group_by(eid) |>
  filter(any(over6m)) |>
  filter(n() > 1) |>
  ungroup() |>
  select(-c(percbase, over6m))
# Saving
write_tsv(
  bmilong1i,
  file.path(projfld, "data", "UKB", "processed_data", "bmilong1i.tsv")
)

# Until instance 2
bmilong2i <- bind_rows(
  w2i |>
    transmute(
      eid,
      event_dt = date_visit,
      weight = weight
    ),
  w1i |>
    transmute(
      eid,
      event_dt = date_visit,
      weight = weight
    ),
  wrec |>
    transmute(
      eid,
      event_dt = date_entry,
      weight = weight
    ),
  wgpc
) |>
  drop_na |>
  # Only including data up to recruitment
  inner_join(w2i, by = "eid", suffix = c("", ".ref")) |>
  inner_join(wrec[, c("eid", "sex")], by = join_by(eid)) |>
  transmute(
    eid, sex,
    time = as.numeric(event_dt - date_visit) / 365.25,
    bmi0 = weight.ref / ((height / 100)^2),
    bmi = weight / ((height / 100)^2),
    bmi = ifelse(time == 0 & bmi != bmi0, NA, bmi)
  ) |>
  filter(time >= -5, time <= 0) |>
  drop_na() |>
  mutate(
    percbase = 100 * bmi / bmi0,
    over6m = time <= -.5
  ) |>
  filter(percbase > 50, percbase < 150) |>
  group_by(eid) |>
  filter(any(over6m)) |>
  filter(n() > 1) |>
  ungroup() |>
  select(-c(percbase, over6m))
# Saving
write_tsv(
  bmilong2i,
  file.path(projfld, "data", "UKB", "processed_data", "bmilong2i.tsv")
)

#-------------------------------------------------------------------------------
# Biomarker data

# Biomarker data at recruitment
bmrec <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "BMOut_Recruitment_participant.csv"
  )
)
# Rename columns
renamer <- c(
  dbp_a1 = "p4079_i0_a0",
  dbp_a2 = "p4079_i0_a1",
  dbp_m1 = "p94_i0_a0",
  dbp_m2 = "p94_i0_a1",
  sbp_a1 = "p4080_i0_a0",
  sbp_a2 = "p4080_i0_a1",
  sbp_m1 = "p93_i0_a0",
  sbp_m2 = "p93_i0_a1",
  wstc = "p48_i0",
  hipc = "p49_i0",
  ALT = "p30620_i0",
  AST = "p30650_i0",
  CRP = "p30710_i0",
  GGT = "p30730_i0",
  GLU = "p30740_i0",
  HDL = "p30760_i0",
  LDL = "p30780_i0",
  TG = "p30870_i0",
  SCR = "p30700_i0"
)
bmrec <- bmrec |>
  rename(all_of(renamer)) |>
  inner_join(wrec, by = join_by(eid)) |>
  mutate(
    sbp_a = rowMeans(cbind(sbp_a1, sbp_a2), na.rm = TRUE),
    sbp_m = rowMeans(cbind(sbp_m1, sbp_m2), na.rm = TRUE),
    dbp_a = rowMeans(cbind(dbp_a1, dbp_a2), na.rm = TRUE),
    dbp_m = rowMeans(cbind(dbp_m1, dbp_m2), na.rm = TRUE),
    SBP = coalesce(sbp_a, sbp_m),
    DBP = coalesce(dbp_a, dbp_m),
    WHR = wstc / hipc,
    bmi = weight / ((height / 100)^2)
  ) |>
  select(
    -c(
      sbp_a1, sbp_a2, sbp_a,
      sbp_m1, sbp_m2, sbp_m,
      dbp_a1, dbp_a2, dbp_a,
      dbp_m1, dbp_m2, dbp_m,
      wstc, hipc,
      weight, height,
      date_entry
    )
  ) |>
  select(eid, sex, age, bmi, everything()) |>
  mutate(
    ALT = if_else(ALT >= 5 & ALT <= 70, ALT, NA),
    AST = if_else(AST >= 15 & AST <= 60, AST, NA),
    CRP = if_else(CRP <= 10, CRP, NA),
    DBP = if_else(DBP >= 50 & DBP <= 120, DBP, NA),
    GGT = if_else(GGT >= 10 & GGT <= 120, GGT, NA),
    HDL = if_else(HDL >= 0.5 & HDL <= 3, HDL, NA),
    LDL = if_else(LDL >= 1 & LDL <= 7, LDL, NA),
    GLU = if_else(GLU >= 3 & GLU <= 15, GLU, NA),
    SBP = if_else(SBP >= 80 & SBP <= 220, SBP, NA),
    SCR = if_else(SCR >= 40 & SCR <= 150, SCR, NA),
    TG = if_else(TG >= 0.5 & TG <= 7, TG, NA),
    WHR = if_else(WHR >= 0.6 & WHR <= 1.4, WHR, NA)
  )
# Saving
write_tsv(
  bmrec,
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv")
)

# Biomarker data at instance 1
bm1i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "BMOut_Instance1_participant.csv"
  ),
  col_types = cols(.default = "n")
)
# Rename columns
renamer <- c(
  dbp_a1 = "p4079_i1_a0",
  dbp_a2 = "p4079_i1_a1",
  dbp_m1 = "p94_i1_a0",
  dbp_m2 = "p94_i1_a1",
  sbp_a1 = "p4080_i1_a0",
  sbp_a2 = "p4080_i1_a1",
  sbp_m1 = "p93_i1_a0",
  sbp_m2 = "p93_i1_a1",
  wstc = "p48_i1",
  hipc = "p49_i1",
  ALT = "p30620_i1",
  AST = "p30650_i1",
  CRP = "p30710_i1",
  GGT = "p30730_i1",
  GLU = "p30740_i1",
  HDL = "p30760_i1",
  LDL = "p30780_i1",
  TG = "p30870_i1",
  SCR = "p30700_i1"
)
bm1i <- bm1i |>
  rename(all_of(renamer)) |>
  inner_join(w1i, by = join_by(eid)) |>
  inner_join(wrec, by = join_by(eid),  suffix = c(".1i", "")) |>
  mutate(
    sbp_a = rowMeans(cbind(sbp_a1, sbp_a2), na.rm = TRUE),
    sbp_m = rowMeans(cbind(sbp_m1, sbp_m2), na.rm = TRUE),
    dbp_a = rowMeans(cbind(dbp_a1, dbp_a2), na.rm = TRUE),
    dbp_m = rowMeans(cbind(dbp_m1, dbp_m2), na.rm = TRUE),
    SBP = coalesce(sbp_a, sbp_m),
    DBP = coalesce(dbp_a, dbp_m),
    WHR = wstc / hipc,
    bmi = weight.1i / ((height / 100)^2),
    age = age.1i
  ) |>
  select(
    -c(
      sbp_a1, sbp_a2, sbp_a,
      sbp_m1, sbp_m2, sbp_m,
      dbp_a1, dbp_a2, dbp_a,
      dbp_m1, dbp_m2, dbp_m,
      wstc, hipc,
      weight, height,
      weight.1i, height.1i,
      date_entry, date_visit,
      age.1i
    )
  ) |>
  select(eid, sex, age, bmi, everything()) |>
  mutate(
    ALT = if_else(ALT >= 5 & ALT <= 70, ALT, NA),
    AST = if_else(AST >= 15 & AST <= 60, AST, NA),
    CRP = if_else(CRP <= 10, CRP, NA),
    DBP = if_else(DBP >= 50 & DBP <= 120, DBP, NA),
    GGT = if_else(GGT >= 10 & GGT <= 120, GGT, NA),
    HDL = if_else(HDL >= 0.5 & HDL <= 3, HDL, NA),
    LDL = if_else(LDL >= 1 & LDL <= 7, LDL, NA),
    GLU = if_else(GLU >= 3 & GLU <= 15, GLU, NA),
    SBP = if_else(SBP >= 80 & SBP <= 220, SBP, NA),
    SCR = if_else(SCR >= 40 & SCR <= 150, SCR, NA),
    TG = if_else(TG >= 0.5 & TG <= 7, TG, NA),
    WHR = if_else(WHR >= 0.6 & WHR <= 1.4, WHR, NA)
  )
# Saving
write_tsv(
  bm1i,
  file.path(projfld, "data", "UKB", "processed_data", "bm1i.tsv")
)

# Imaging data at instance 2
img2i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Imaging_Instance2_participant.csv"
  ),
  col_types = cols(.default = "n")
)
# Rename columns
renamer <- c(
  AndF_perc = "p23247_i2",
  LiverF_perc = "p40061_i2",
  ArmsF_perc = "p23259_i2",
  GynF_perc = "p23264_i2",
  LegsF_perc = "p23276_i2",
  BF_perc = "p23281_i2",
  TrunkF_perc = "p23286_i2"
)
img2i <- img2i |>
  rename(all_of(renamer)) |>
  inner_join(w2i, by = join_by(eid)) |>
  inner_join(wrec, by = join_by(eid),  suffix = c(".2i", "")) |>
  mutate(
    bmi = weight.2i / ((height / 100)^2),
    age = age.2i
  ) |>
  select(
    -c(
      weight, height,
      weight.2i, height.2i,
      date_entry, date_visit,
      age.2i
    )
  ) |>
  select(eid, sex, age, bmi, everything()) |>
  mutate(
    across(-c(eid, sex, age, bmi, LiverF_perc), \(x) x * 100)
  )
# Saving
write_tsv(
  img2i,
  file.path(projfld, "data", "UKB", "processed_data", "img2i.tsv")
)

#--------------------------------------------------------------------------------
# Covariates

# Genetic ancestry data
genancestry <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "GenAncestry.csv"
  ),
  col_types = "ncc"
)
# Rename columns
renamer <- c(ancestry = "p30079")
ancestrytab <- genancestry |>
  rename(all_of(renamer)) |>
  transmute(
    eid,
    ancestry = gsub("^.*\\((.*?)\\)$", "\\1", ancestry),
    ancestry = replace_na(ancestry, "Unknown")
  )
# Saving
write_tsv(
  ancestrytab,
  file.path(
    projfld, "data", "UKB", "processed_data",
    "ancestrytab.tsv"
  )
)

# Covariates to check across instances
tocheck <- c(
  "Diabetes", "Cancer", "Menopause", "CHD",
  "Stroke", "Asthma", "COPD", "Insulin"
)

# Covariate data at recruitment
cvrec <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Covariates_Recruitment_participant.csv"
  )
)
# Rename columns
renamer <- c(
  Diabetes = "p2443_i0",
  Cancer = "p2453_i0",
  VascHeart = "p6150_i0",
  Diseases = "p6152_i0",
  WtCh = "p2306_i0",
  SeriousDx = "p2473_i0",
  CurrSmoke = "p20116_i0",
  TSI = "p22189",
  Menopause = "p2724_i0",
  MedsM = "p6177_i0",
  MedsF = "p6153_i0",
  Edu = "p6138_i0",
  AgeLastBirth = "p2764_i0",
  Pregnant = "p3140_i0"
)
cvrec <- cvrec |>
  rename(all_of(renamer)) |>
  mutate(
    Diabetes = replace_na(1 * (Diabetes == "Yes"), 0),
    Cancer = replace_na(1 * (Cancer == "Yes"), 0),
    CHD = replace_na(
      1 * grepl("Heart attack", VascHeart, fixed = TRUE),
      0
    ),
    Stroke = replace_na(
      1 * grepl("Stroke", VascHeart, fixed = TRUE),
      0
    ),
    Asthma = replace_na(
      1 * grepl("Asthma", Diseases, fixed = TRUE),
      0
    ),
    COPD = replace_na(
      1 * grepl("Emphysema/chronic bronchitis", Diseases, fixed = TRUE),
      0
    ),
    WtG = replace_na(
      1 * (WtCh == "Yes - gained weight"),
      0
    ),
    WtL = replace_na(
      1 * (WtCh == "Yes - lost weight"),
      0
    ),
    SeriousDx = replace_na(
      1 * (
        SeriousDx == paste(
          "Yes", "-",
          "you will be asked about this later by an interviewer"
        )
      ),
      0
    ),
    CurrSmoke = replace_na(
      1 * (CurrSmoke == "Current"),
      0
    ),
    Menopause = replace_na(
      1 * (
        Menopause %in% c(
          "Yes", "Not sure - other reason", "Not sure - had a hysterectomy"
        )
      ),
      0
    ),
    Insulin = replace_na(
      1 * "|"(
        grepl("Insulin", MedsM, fixed = TRUE),
        grepl("Insulin", MedsF, fixed = TRUE)
      ),
      0
    ),
    HTMed = replace_na(
      1 * "|"(
        grepl("Blood pressure medication", MedsM, fixed = TRUE),
        grepl("Blood pressure medication", MedsF, fixed = TRUE)
      ),
      0
    ),
    LipidLower = replace_na(
      1 * "|"(
        grepl("Cholesterol lowering medication", MedsM, fixed = TRUE),
        grepl("Cholesterol lowering medication", MedsF, fixed = TRUE)
      ),
      0
    ),
    CollegeEd = replace_na(
      1 * grepl("College or University degree", Edu)
    ),
    Pregnant = replace_na(1 * (Pregnant == "Yes"), 0),
    AgeLastBirth = as.numeric(AgeLastBirth)
  ) |>
  select(-c(VascHeart, Diseases, WtCh, SeriousDx, MedsM, MedsF, Edu)) |>
  inner_join(wrec[, c("eid", "age")], by = join_by(eid)) |>
  mutate(RecentPreg = replace_na(1 * ((age - AgeLastBirth) < 5), 0)) |>
  select(-c(age, AgeLastBirth)) |>
  drop_na() |>
  filter(Pregnant == 0, RecentPreg == 0) |>
  select(-c(Pregnant, RecentPreg))
# Saving
write_tsv(
  cvrec,
  file.path(projfld, "data", "UKB", "processed_data", "cvrec.tsv")
)

# Covariate data at instance 1
cv1i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Covariates_Instance1_participant.csv"
  )
)
# Rename columns
renamer <- c(
  Diabetes = "p2443_i1",
  Cancer = "p2453_i1",
  VascHeart = "p6150_i1",
  Diseases = "p6152_i1",
  WtCh = "p2306_i1",
  SeriousDx = "p2473_i1",
  CurrSmoke = "p20116_i1",
  TSI = "p22189",
  Menopause = "p2724_i1",
  MedsM = "p6177_i1",
  MedsF = "p6153_i1",
  Edu = "p6138_i1",
  Pregnant = "p3140_i1",
  AgeLastBirth = "p2764_i1"
)
cv1i <- cv1i |>
  rename(all_of(renamer)) |>
  mutate(
    Diabetes = replace_na(1 * (Diabetes == "Yes"), 0),
    Cancer = replace_na(1 * (Cancer == "Yes"), 0),
    CHD = replace_na(
      1 * grepl("Heart attack", VascHeart, fixed = TRUE),
      0
    ),
    Stroke = replace_na(
      1 * grepl("Stroke", VascHeart, fixed = TRUE),
      0
    ),
    Asthma = replace_na(
      1 * grepl("Asthma", Diseases, fixed = TRUE),
      0
    ),
    COPD = replace_na(
      1 * grepl("Emphysema/chronic bronchitis", Diseases, fixed = TRUE),
      0
    ),
    WtG = replace_na(
      1 * (WtCh == "Yes - gained weight"),
      0
    ),
    WtL = replace_na(
      1 * (WtCh == "Yes - lost weight"),
      0
    ),
    SeriousDx = replace_na(
      1 * (
        SeriousDx == paste(
          "Yes", "-",
          "you will be asked about this later by an interviewer"
        )
      ),
      0
    ),
    CurrSmoke = replace_na(
      1 * (CurrSmoke == "Current"),
      0
    ),
    Menopause = replace_na(
      1 * (Menopause == "Yes"),
      0
    ),
    Insulin = replace_na(
      1 * "|"(
        grepl("Insulin", MedsM, fixed = TRUE),
        grepl("Insulin", MedsF, fixed = TRUE)
      ),
      0
    ),
    HTMed = replace_na(
      1 * "|"(
        grepl("Blood pressure medication", MedsM, fixed = TRUE),
        grepl("Blood pressure medication", MedsF, fixed = TRUE)
      ),
      0
    ),
    LipidLower = replace_na(
      1 * "|"(
        grepl("Cholesterol lowering medication", MedsM, fixed = TRUE),
        grepl("Cholesterol lowering medication", MedsF, fixed = TRUE)
      ),
      0
    ),
    CollegeEd = replace_na(
      1 * grepl("College or University degree", Edu)
    ),
    Pregnant = replace_na(1 * (Pregnant == "Yes"), 0),
    AgeLastBirth = as.numeric(AgeLastBirth)
  ) |>
  select(-c(VascHeart, Diseases, WtCh, SeriousDx, MedsM, MedsF, Edu)) |>
  inner_join(w1i[, c("eid", "age")], by = join_by(eid)) |>
  mutate(RecentPreg = replace_na(1 * ((age - AgeLastBirth) < 5), 0)) |>
  select(-c(age, AgeLastBirth)) |>
  drop_na() |>
  filter(Pregnant == 0, RecentPreg == 0) |>
  select(-c(Pregnant, RecentPreg)) |>
  left_join(
    cvrec[, c("eid", tocheck)],
    by = join_by(eid), suffix = c("", ".rec")
  ) |>
  mutate(
    across(
      all_of(tocheck),
      ~ pmax(
        .x,
        get(paste0(cur_column(), ".rec")),
        na.rm = TRUE
      )
    )
  ) |>
  select(-ends_with(".rec"))
# Saving
write_tsv(
  cv1i,
  file.path(projfld, "data", "UKB", "processed_data", "cv1i.tsv")
)

# Covariate data at instance 2
cv2i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Covariates_Instance2_participant.csv"
  )
)
# Rename columns
renamer <- c(
  Diabetes = "p2443_i2",
  Cancer = "p2453_i2",
  VascHeart = "p6150_i2",
  Diseases = "p6152_i2",
  WtCh = "p2306_i2",
  SeriousDx = "p2473_i2",
  CurrSmoke = "p20116_i2",
  TSI = "p22189",
  Menopause = "p2724_i2",
  MedsM = "p6177_i2",
  MedsF = "p6153_i2",
  Edu = "p6138_i2",
  AgeLastBirth = "p2764_i2",
  Pregnant = "p3140_i2"
)
cv2i <- cv2i |>
  rename(all_of(renamer)) |>
  mutate(
    Diabetes = replace_na(1 * (Diabetes == "Yes"), 0),
    Cancer = replace_na(1 * (Cancer == "Yes"), 0),
    CHD = replace_na(
      1 * "|"(
        grepl("Heart attack", VascHeart, fixed = TRUE),
        grepl("Angina", VascHeart, fixed = TRUE)
      ),
      0
    ),
    Stroke = replace_na(
      1 * grepl("Stroke", VascHeart, fixed = TRUE),
      0
    ),
    Asthma = replace_na(
      1 * grepl("Asthma", Diseases, fixed = TRUE),
      0
    ),
    COPD = replace_na(
      1 * grepl("Emphysema/chronic bronchitis", Diseases, fixed = TRUE),
      0
    ),
    WtG = replace_na(
      1 * (WtCh == "Yes - gained weight"),
      0
    ),
    WtL = replace_na(
      1 * (WtCh == "Yes - lost weight"),
      0
    ),
    SeriousDx = replace_na(
      1 * (
        SeriousDx == paste(
          "Yes", "-",
          "you will be asked about this later by an interviewer"
        )
      ),
      0
    ),
    CurrSmoke = replace_na(
      1 * (CurrSmoke == "Current"),
      0
    ),
    Menopause = replace_na(
      1 * (Menopause == "Yes"),
      0
    ),
    Insulin = replace_na(
      1 * "|"(
        grepl("Insulin", MedsM, fixed = TRUE),
        grepl("Insulin", MedsF, fixed = TRUE)
      ),
      0
    ),
    HTMed = replace_na(
      1 * "|"(
        grepl("Blood pressure medication", MedsM, fixed = TRUE),
        grepl("Blood pressure medication", MedsF, fixed = TRUE)
      ),
      0
    ),
    LipidLower = replace_na(
      1 * "|"(
        grepl("Cholesterol lowering medication", MedsM, fixed = TRUE),
        grepl("Cholesterol lowering medication", MedsF, fixed = TRUE)
      ),
      0
    ),
    CollegeEd = replace_na(
      1 * grepl("College or University degree", Edu)
    ),
    Pregnant = replace_na(1 * (Pregnant == "Yes"), 0),
    AgeLastBirth = as.numeric(AgeLastBirth)
  ) |>
  select(-c(VascHeart, Diseases, WtCh, SeriousDx, MedsM, MedsF, Edu)) |>
  inner_join(w2i[, c("eid", "age")], by = join_by(eid)) |>
  mutate(RecentPreg = replace_na(1 * ((age - AgeLastBirth) < 5), 0)) |>
  select(-c(age, AgeLastBirth)) |>
  drop_na() |>
  filter(Pregnant == 0, RecentPreg == 0) |>
  select(-c(Pregnant, RecentPreg)) |>
  left_join(
    cvrec[, c("eid", tocheck)],
    by = join_by(eid), suffix = c("", ".rec")
  ) |>
  left_join(
    cv1i[, c("eid", tocheck)],
    by = join_by(eid), suffix = c("", ".1i")
  ) |>
  mutate(
    across(
      all_of(tocheck),
      ~ pmax(
        .x, get(paste0(cur_column(), ".rec")),
        get(paste0(cur_column(), ".1i")),
        na.rm = TRUE
      )
    )
  ) |>
  select(-ends_with(".rec"), -ends_with(".1i"))
# Saving
write_tsv(
  cv2i,
  file.path(projfld, "data", "UKB", "processed_data", "cv2i.tsv")
)

#-------------------------------------------------------------------------------
# Survival outcomes data

# Administrative censoring date
censordate <- as.Date("2021-12-31")

# Mortality data
mortdat <- read_tsv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "death.txt"
  ),
  col_types = "n---D",
  col_names = c("eid", "date_death"),
  skip = 1,
  locale = locale(date_format = "%d/%m/%Y")
)

# Outcomes data
outdat <- read_tsv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Outcomes_participant.tsv"
  ),
  col_types = "nDDDDD"
)
# Rename columns
renamer <- c(
  date_CAD = "42000-0.0",
  date_Stroke = "42006-0.0",
  date_LiverD1 = "131662-0.0",
  date_LiverD2 = "131664-0.0",
  date_LiverD3 = "131666-0.0"
)
outdat <- outdat |>
  rename(all_of(renamer)) |>
  mutate(
    eid,
    date_LiverD = pmin(
      date_LiverD1,
      date_LiverD2,
      date_LiverD3,
      na.rm = TRUE
    ),
    across(matches("[0-9]$"), ~NULL),
  ) |>
  inner_join(
    wrec[, c("eid", "date_entry")],
    by = "eid"
  ) |>
  left_join(mortdat, by = "eid") |>
  mutate(
    admin_censor = pmin(
      date_death, censordate,
      as.Date(as.numeric(date_entry) + (365.25 * 10)),
      na.rm = TRUE
    )
  ) |>
  select(-date_death) |>
  pivot_longer(
    c(date_CAD, date_Stroke, date_LiverD),
    names_to = c(".value", "outcome"),
    names_sep = "_",
  ) |>
  filter(
    !(outcome == "CAD" & eid %in% cvrec$eid[cvrec$CHD == 1]),
    !(outcome == "Stroke" & eid %in% cvrec$eid[cvrec$Stroke == 1]),
  ) |>
  mutate(
    event = if_else(!is.na(date) & date <= admin_censor, 1, 0),
    date = pmin(date, admin_censor, na.rm = TRUE),
    time = (as.numeric(date) - as.numeric(date_entry)) / 365.25
  ) |>
  filter(time > 0) |>
  select(-c(date_entry, admin_censor, date))
# Saving
write_tsv(
  outdat,
  file.path(projfld, "data", "UKB", "processed_data", "outdat.tsv")
)