library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

list.files(file.path(projfld, "data", "UKB", "raw_data"))

wrec <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Weight_Recruitment_participant.csv"
  )
)

head(wrec)

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

head(wrec)

nrow(wrec)

w1i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Weight_Instance1_participant.csv"
  )
)

head(w1i)

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

w2i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Weight_Instance2_participant.csv"
  )
)

head(w2i)

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

head(w2i)

nrow(w2i)

wgp <- read_tsv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "OnlyWeights_gp.tsv"
  ),
  col_types = "ncDnnn"
)

head(wgp)

wgp_probl <- problems(wgp)

head(wgp_probl)

wgp_probl |>
  count(col)

wgp_probl |>
  count(col, expected, actual)

wgp |>
  group_by(data_provider) |>
  summarise(
    across(
      starts_with("value"),
      list(
        Med = \(x) median(x, na.rm = TRUE),
        q25 = \(x) quantile(x, probs = .25, na.rm = TRUE),
        q75 = \(x) quantile(x, probs = .75, na.rm = TRUE)
      )
    )
  ) |>
  pivot_longer(
    -data_provider,
    names_to = c("value", ".value"),
    names_sep = "_"
  )

wgpc <- wgp |>
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

head(wgpc)

intersect(unique(wgpc$eid), wrec$eid) |> length()

bmrec <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "BMOut_Recruitment_participant.csv"
  )
)

head(bmrec)

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

head(bmrec)

nrow(bmrec)

bm1i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "BMOut_Instance1_participant.csv"
  ),
  col_types = cols(.default = "n")
)

head(bm1i)

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

head(bm1i)

nrow(bm1i)

img2i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Imaging_Instance2_participant.csv"
  ),
  col_types = cols(.default = "n")
)

head(img2i)

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

head(img2i)

nrow(img2i)

genancestry <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "GenAncestry.csv"
  ),
  col_types = "ncc"
)
head(genancestry)

count(genancestry, p30079)

renamer <- c(ancestry = "p30079")
ancestrytab <- genancestry |>
  rename(all_of(renamer)) |>
  transmute(
    eid,
    ancestry = gsub("^.*\\((.*?)\\)$", "\\1", ancestry),
    ancestry = replace_na(ancestry, "Unknown")
  )
head(ancestrytab)

count(ancestrytab, ancestry)

write_tsv(
  ancestrytab,
  file.path(
    projfld, "data", "UKB", "processed_data",
    "ancestrytab.tsv"
  )
)

cvrec <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Covariates_Recruitment_participant.csv"
  )
)

head(cvrec)

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

head(cvrec)

nrow(cvrec)

cv1i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Covariates_Instance1_participant.csv"
  )
)

head(cv1i)

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
  select(-c(Pregnant, RecentPreg))

tocheck <- c(
  "Diabetes", "Cancer", "Menopause", "CHD",
  "Stroke", "Asthma", "COPD", "Insulin"
)

cv1i <- cv1i |>
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

head(cv1i)

nrow(cv1i)

cv2i <- read_csv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Covariates_Instance2_participant.csv"
  )
)

head(cv2i)

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
  select(-c(Pregnant, RecentPreg))

cv2i <- cv2i |>
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

head(cv2i)

nrow(cv2i)

outdat <- read_tsv(
  file.path(
    projfld, "data", "UKB", "raw_data",
    "Outcomes_participant.tsv"
  ),
  col_types = "nDDDDD"
)

head(outdat)

renamer <- c(
  date_CAD = "42000-0.0",
  date_Stroke = "42006-0.0",
  date_LiverD1 = "131662-0.0",
  date_LiverD2 = "131664-0.0",
  date_LiverD3 = "131666-0.0"
)

censordate <- as.Date("2022-05-31")

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
    across(
      starts_with("date_"),
      ~as.integer(!is.na(.x) & .x <= censordate),
      .names = "event_{gsub(\"date_\", \"\", .col)}"
    )
  ) |>
  pivot_longer(
    -eid,
    names_to = c(".value", "outcome"),
    names_sep = "_"
  ) |>
  filter(
    !(outcome == "CAD" & eid %in% cvrec$eid[cvrec$CHD == 1]),
    !(outcome == "Stroke" & eid %in% cvrec$eid[cvrec$Stroke == 1]),
  ) |>
  mutate(
    date = pmin(date, censordate, na.rm = TRUE),
  ) |>
  inner_join(wrec[, c("eid", "date_entry")], by = join_by(eid)) |>
  filter(date >= date_entry) |>
  mutate(
    time = (as.numeric(date) - as.numeric(date_entry)) / 365.25
  ) |>
  select(-c(date_entry, date)) |>
  mutate(
    event = if_else(time > 10, 0, event),
    time = pmin(time, 10)
  )
head(outdat)

bmilongrec <- wrec %>%
  # Joining longitudinal data with data at recruitment
  transmute(
    eid,
    event_dt = date_entry,
    weight = weight
  ) %>%
  drop_na %>%
  bind_rows(wgpc) %>%
  # Only including data up to recruitment
  inner_join(wrec, by = "eid", suffix = c("", ".ref")) %>%
  transmute(
    eid, sex,
    time = as.numeric(event_dt - date_entry) / 365.25,
    bmi0 = weight.ref / ((height / 100)^2),
    bmi = weight / ((height / 100)^2),
    bmi = ifelse(time == 0 & bmi != bmi0, NA, bmi),
    across(c(bmi, bmi0), \(x) round(x, digits = 2))
  ) %>%
  filter(time >= -5, time <= 0) %>%
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
  select(-c(percbase, over6m)) |>
  filter(eid %in% bmrec$eid & eid %in% cvrec$eid)

head(bmilongrec)

length(unique(bmilongrec$eid))

bmilong1i <- bind_rows(
  w1i %>%
    transmute(
      eid,
      event_dt = date_visit,
      weight = weight
    ),
  wrec %>%
    transmute(
      eid,
      event_dt = date_entry,
      weight = weight
    ),
  wgpc
) %>%
  drop_na %>%
  # Only including data up to recruitment
  inner_join(w1i, by = "eid", suffix = c("", ".ref")) %>%
  inner_join(wrec[, c("eid", "sex")], by = join_by(eid)) %>%
  transmute(
    eid, sex,
    time = as.numeric(event_dt - date_visit) / 365.25,
    bmi0 = weight.ref / ((height / 100)^2),
    bmi = weight / ((height / 100)^2),
    bmi = ifelse(time == 0 & bmi != bmi0, NA, bmi),
    across(c(bmi, bmi0), \(x) round(x, digits = 2))
  ) %>%
  filter(time >= -5, time <= 0) %>%
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
  select(-c(percbase, over6m)) |>
  filter(eid %in% bm1i$eid & eid %in% cv1i$eid)

head(bmilong1i)

length(unique(bmilong1i$eid))

bmilong2i <- bind_rows(
  w2i %>%
    transmute(
      eid,
      event_dt = date_visit,
      weight = weight
    ),
  w1i %>%
    transmute(
      eid,
      event_dt = date_visit,
      weight = weight
    ),
  wrec %>%
    transmute(
      eid,
      event_dt = date_entry,
      weight = weight
    ),
  wgpc
) %>%
  drop_na %>%
  # Only including data up to recruitment
  inner_join(w2i, by = "eid", suffix = c("", ".ref")) %>%
  inner_join(wrec[, c("eid", "sex")], by = join_by(eid)) %>%
  transmute(
    eid, sex,
    time = as.numeric(event_dt - date_visit) / 365.25,
    bmi0 = weight.ref / ((height / 100)^2),
    bmi = weight / ((height / 100)^2),
    bmi = ifelse(time == 0 & bmi != bmi0, NA, bmi),
    across(c(bmi, bmi0), \(x) round(x, digits = 2))
  ) %>%
  filter(time >= -5, time <= 0) %>%
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
  select(-c(percbase, over6m)) |>
  filter(eid %in% img2i$eid & eid %in% cv2i$eid)

head(bmilong2i)

length(unique(bmilong2i$eid))

bmrec |>
  pivot_longer(-c(eid, sex, age, bmi)) |>
  inner_join(cvrec, by = join_by(eid)) |>
  group_by(name) |>
  summarise(
    Ntotal = n(),
    Nnonmiss = sum(!is.na(value)),
    Nwithlong = sum(eid %in% bmilongrec$eid),
    Nnonmisswithlong = sum(!is.na(value) & eid %in% bmilongrec$eid)
  ) |>
  arrange(name)

bm1i |>
  pivot_longer(-c(eid, sex, age, bmi)) |>
  inner_join(cv1i, by = join_by(eid)) |>
  group_by(name) |>
  summarise(
    Ntotal = n(),
    Nnonmiss = sum(!is.na(value)),
    Nwithlong = sum(eid %in% bmilong1i$eid),
    Nnonmisswithlong = sum(!is.na(value) & eid %in% bmilong1i$eid)
  ) |>
  arrange(name)

img2i |>
  pivot_longer(-c(eid, sex, age, bmi)) |>
  inner_join(cv2i, by = join_by(eid)) |>
  group_by(name) |>
  summarise(
    Ntotal = n(),
    Nnonmiss = sum(!is.na(value)),
    Nwithlong = sum(eid %in% bmilong2i$eid),
    Nnonmisswithlong = sum(!is.na(value) & eid %in% bmilong2i$eid)
  ) |>
  arrange(name)

bmilongrec |>
  transmute(eid, sex, included = TRUE) |>
  unique() |>
  right_join(bmrec[, c("eid", "sex", "age", "bmi")], by = join_by(eid, sex)) |>
  replace_na(list(included = FALSE)) |>
  inner_join(cvrec[, "eid"], by = join_by(eid)) |>
  group_by(sex, included) |>
  summarise(
    across(
      c(age, bmi),
      list(
        Med = \(x) median(x, na.rm = TRUE),
        q25 = \(x) quantile(x, probs = .25, na.rm = TRUE),
        q75 = \(x) quantile(x, probs = .75, na.rm = TRUE)
      )
    ),
    Ntotal = n(),
    .groups = "drop"
  ) |>
  pivot_longer(
    -c(sex, included, Ntotal),
    names_to = c("variable", ".value"),
    names_sep = "_",
    cols_vary = "slowest"
  )

bmilongrec |>
  inner_join(bmrec[, "eid"], by = join_by(eid)) |>
  inner_join(cvrec[, "eid"], by = join_by(eid)) |>
  group_by(eid) |>
  summarise(
    NMs = n(),
    FUT = min(time) - max(time),
    .groups = "drop"
  ) |>
  summarise(
    across(
      c(NMs, FUT),
      list(
        Med = \(x) median(x, na.rm = TRUE),
        q25 = \(x) quantile(x, probs = .25, na.rm = TRUE),
        q75 = \(x) quantile(x, probs = .75, na.rm = TRUE)
      )
    ),
    Prop3Ms = 100 * mean(NMs >= 3, na.rm = TRUE)
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

bmilongrec |>
  inner_join(bmrec[, "eid"], by = join_by(eid)) |>
  inner_join(cvrec[, "eid"], by = join_by(eid)) |>
  group_by(sex, eid) |>
  summarise(
    NMs = n(),
    FUT = min(time) - max(time),
    .groups = "drop_last"
  ) |>
  summarise(
    across(
      c(NMs, FUT),
      list(
        Med = \(x) median(x, na.rm = TRUE),
        q25 = \(x) quantile(x, probs = 0.25, na.rm = TRUE),
        q75 = \(x) quantile(x, probs = 0.75, na.rm = TRUE)
      )
    )
  ) |>
  mutate(
    across(where(is.numeric), \(x) round(x, digits = 2))
  ) |>
  pivot_longer(
    -sex,
    names_to = c("variable", ".value"),
    names_sep = "_",
    cols_vary = "slowest"
  )

wrec |>
  select(eid, date_entry) |>
  inner_join(w1i, by = join_by(eid)) |>
  mutate(
    time = as.numeric(date_visit - date_entry) / 365.25
  ) |>
  summarise(
    across(
      time,
      list(
        Med = \(x) median(x, na.rm = TRUE),
        q25 = \(x) quantile(x, probs = .25, na.rm = TRUE),
        q75 = \(x) quantile(x, probs = .75, na.rm = TRUE),
        n = \(x) sum(!is.na(x))
      )
    )
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

write_tsv(
  bmrec,
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv")
)
write_tsv(
  cvrec,
  file.path(projfld, "data", "UKB", "processed_data", "cvrec.tsv")
)
write_tsv(
  bmilongrec,
  file.path(projfld, "data", "UKB", "processed_data", "bmilongrec.tsv")
)

bmilongrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmilongrec.tsv"),
  show_col_types = FALSE
)

head(bmilongrec)

write_tsv(
  bm1i,
  file.path(projfld, "data", "UKB", "processed_data", "bm1i.tsv")
)
write_tsv(
  cv1i,
  file.path(projfld, "data", "UKB", "processed_data", "cv1i.tsv")
)
write_tsv(
  bmilong1i,
  file.path(projfld, "data", "UKB", "processed_data", "bmilong1i.tsv")
)

write_tsv(
  img2i,
  file.path(projfld, "data", "UKB", "processed_data", "img2i.tsv")
)
write_tsv(
  cv2i,
  file.path(projfld, "data", "UKB", "processed_data", "cv2i.tsv")
)
write_tsv(
  bmilong2i,
  file.path(projfld, "data", "UKB", "processed_data", "bmilong2i.tsv")
)

write_tsv(
  outdat,
  file.path(projfld, "data", "UKB", "processed_data", "outdat.tsv")
)
