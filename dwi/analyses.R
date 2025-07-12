# Code to perform analyses for the study
# 

# 1) load libraries
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(janitor)

# 2) read in your two CSVs
lh <- read_csv("/Users/a9ws/imaging/datasets/mrgfus/derivatives/fba/analysis/sweetspot/lh_mean-fd.csv") %>% 
  rename(sub_id = `sub-id`)
rh <- read_csv("/Users/a9ws/imaging/datasets/mrgfus/derivatives/fba/analysis/sweetspot/rh_mean-fd.csv") %>% 
  rename(sub_id = `sub-id`)

# 3) read in hemis.txt
#    assuming it's two columns, whitespace‐delimited, no header:
hemis <- read_table(
  "/Users/a9ws/imaging/datasets/mrgfus/derivatives/study_files/fba/hemis.txt",
  col_names = FALSE,       # file has no header row
  col_types = cols(        # both are character columns
    X1 = col_character(),
    X2 = col_character()
  )
) %>%
  # now assign sensible names
  rename(
    sub_id = X1,
    hemi   = X2
  )

# 3) Read clinical scores from Excel
clin_wide <- read_xlsx(
  "/Users/a9ws/imaging/datasets/mrgfus/derivatives/study_files/pheno_full.xlsx",
  sheet = 1
) %>%
  # rename the original Excel names into exactly what you want
  rename(
    sub_id = SUB_ID,
    ses_01 = PRE,
    ses_02 = `6m`,
    ses_03 = `12m`
  ) 

# 4) for each subject, pick the correct row from either lh or rh
result <- hemis %>%
  rowwise() %>%
  do({
    sid <- .$sub_id
    if (.$hemi == "lh") {
      lh %>% filter(sub_id == sid)
    } else {
      rh %>% filter(sub_id == sid)
    }
  }) %>%
  ungroup()

# clean up all the names at once:
result <- result %>%
  clean_names()  
# now you should have: sub_id, hemi, ses_01, ses_02, ses_03

# then pivot:
imaging_long <- result %>%
  pivot_longer(
    cols      = starts_with("ses_"),
    names_to  = "session",
    values_to = "imaging_metric"
  ) %>%
  mutate(session = factor(session, levels = c("ses_01","ses_02","ses_03")))

clin_long <- clin_wide %>%
  # ensure all sessions are numeric
  mutate(across(
    c(ses_01, ses_02, ses_03),
    ~ as.numeric(.)
  )) %>%
  pivot_longer(
    cols     = starts_with("ses_"),
    names_to = "session",
    values_to= "clinical_score"
  ) %>%
  mutate(
    session = factor(session, levels = c("ses_01","ses_02","ses_03"))
  )

# 7) Join imaging and clinical data
df_long <- imaging_long %>%
  inner_join(clin_long, by = c("sub_id","session"))

df_long <- df_long %>%
  group_by(sub_id) %>%
  # age is constant within each subject so group_mean = age itself;
  # better to do grand-mean centering outside the group:
  ungroup() %>%
  mutate(
    age_c = AGE - mean(AGE, na.rm=TRUE)
  )

# treat age and sex as fixed effects
m_cov <- lmer(
  clinical_score ~ imaging_metric   # your main within-subject predictor
                 + age_c          # between-subject
                 + SEX          # between-subject
                 + (1 | sub_id),# random intercept per subject
  data      = df_long,
  na.action = na.exclude
)
summary(m_cov)



