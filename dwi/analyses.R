# Code to perform GAMM for the nmap metrics

library(tidyverse)

clin <- read_csv("clinical_data.csv")

clin_long <- clin %>%
  pivot_longer(
    cols = starts_with("tremor_tp"),
    names_to = "timepoint",
    names_pattern = "tremor_(tp\\d)",
    values_to = "tremor_score"
  )

metric_data <- read_csv("combined_metrics_long.csv")

# Merge clinical with metric data
df <- left_join(metric_data, clin_long, by = c("sub-id", "timepoint"))

library(mgcv)

metrics <- c("fdc", "log_fc", "fd")

models <- list()

for (m in metrics) {
  df_sub <- filter(df, metric == m)

  model <- gamm(
    tremor_score ~ s(value) + s(timepoint, k = 3),
    random = list(`sub-id` = ~1),
    data = df_sub
  )

  models[[m]] <- model
}

summary(models[["fdc"]]$gam)
plot(models[["fdc"]]$gam)

summary(models[["log_fc"]]$gam)
plot(models[["log_fc"]]$gam)

summary(models[["fd"]]$gam)
plot(models[["fd"]]$gam)

