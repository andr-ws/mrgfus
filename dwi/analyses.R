# Code to perform GAMM for the nmap metrics

fd=${fba}/analysis/sweetspot/mean-fd.csv
fc=${fba}/analysis/sweetspot/mean-log_fc.csv
fdc=${fba}/analysis/sweetspot/mean-fdc.csv

library(tidyverse)

# Load data
fd <- read_csv("fd_data.csv") ^^^
clin <- read_csv("pheno_clean.csv")

# Reshape both to long format
fd_long <- fd %>%
  pivot_longer(cols = starts_with("fd_tp"), names_to = "timepoint", values_to = "fibre_density")

clin_long <- clin %>%
  pivot_longer(cols = starts_with("tp"), names_to = "timepoint", values_to = "tremor_score") %>%
  mutate(timepoint = str_replace(timepoint, "tremor_", "fd_"))  # to match timepoint col in fd_long

# Merge
data_long <- left_join(fd_long, clin_long, by = c("subject_id", "timepoint"))

# Clean timepoint label if needed
data_long <- data_long %>%
  mutate(timepoint = case_when(
    timepoint == "fd_tp1" ~ 0,
    timepoint == "fd_tp2" ~ 6,
    timepoint == "fd_tp3" ~ 12,
    TRUE ~ NA_real_
  ))
