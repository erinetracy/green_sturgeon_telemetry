#==============================================================================
# 00_BUILD_MODEL12B_COVARIATES.R
# Builds fish-level covariates for model 12b from clean environmental data
# Requires: gs_daily_environmental_data.csv, gs_multistate_data.RData
# Saves: gs_model12b_covariates.csv, gs_model12b_std_params.csv
#==============================================================================

library(dplyr)
library(tidyr)
library(writexl)

# --- Load clean environmental data ---
daily_env <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_daily_environmental_data.csv"
) %>%
  mutate(date = as.Date(date))

cat("Environmental data loaded:", nrow(daily_env), "days\n")

# --- Load detection data ---
load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_multistate_data.RData")

# safety check
stopifnot(all(ch_mat_nimble %in% c(1, 2, 4, 5, 6, 7)))
cat("Detection history:", nrow(detection_history), "fish\n")

# --- Build arrival dates ---
benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(benicia_date = as.Date(min(first_detection)),
                   .groups = "drop")

occ2_dates <- model1_events %>%
  filter(occasion == 2) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(arr_date_occ2 = as.Date(min(first_detection)),
                   .groups = "drop")

occ3_dates <- model1_events %>%
  filter(occasion == 3) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(arr_date_occ3 = as.Date(min(first_detection)),
                   .groups = "drop")

cat("Benicia dates:", nrow(benicia_dates), "\n")
cat("Occ2 dates:",   nrow(occ2_dates), "\n")
cat("Occ3 dates:",   nrow(occ3_dates), "\n")

# --- Extract fish-level covariates ---
fish_cov_12b <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(occ2_dates,    by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,    by = c("animal_id", "water_year")) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # same-day Rio Vista flow at occ2 (Georgiana routing)
    flow_occ2_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(daily_env$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else daily_env$flow_rv_cfs[idx]
      }
    },
    
    # same-day GES flow at occ3 (SS routing + detection)
    flow_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(daily_env$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else daily_env$flow_ges_cfs[idx]
      }
    },
    
    # 7-day mean temperature before Benicia (incomplete migration)
    temp_7day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- daily_env$temp_c[
          daily_env$date >= (benicia_date - 7) &
            daily_env$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    }
    
  ) %>%
  ungroup()

cat("\nMissing flow_occ2:", sum(is.na(fish_cov_12b$flow_occ2_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_cov_12b$flow_occ3_raw)), "\n")
cat("Missing temp_7day:", sum(is.na(fish_cov_12b$temp_7day_raw)), "\n")

# verify row order
stopifnot(identical(fish_cov_12b$animal_id, detection_history$animal_id))
stopifnot(identical(as.numeric(fish_cov_12b$water_year),
                    as.numeric(detection_history$water_year)))
cat("Row order verified\n")

# --- Standardize ---
flow_occ2_mean_12b <- mean(fish_cov_12b$flow_occ2_raw, na.rm = TRUE)
flow_occ2_sd_12b   <- sd(fish_cov_12b$flow_occ2_raw,   na.rm = TRUE)
flow_occ3_mean_12b <- mean(fish_cov_12b$flow_occ3_raw, na.rm = TRUE)
flow_occ3_sd_12b   <- sd(fish_cov_12b$flow_occ3_raw,   na.rm = TRUE)
temp_7day_mean_12b <- mean(fish_cov_12b$temp_7day_raw, na.rm = TRUE)
temp_7day_sd_12b   <- sd(fish_cov_12b$temp_7day_raw,   na.rm = TRUE)

fish_cov_12b <- fish_cov_12b %>%
  mutate(
    flow_occ2_std = (flow_occ2_raw - flow_occ2_mean_12b) / flow_occ2_sd_12b,
    flow_occ3_std = (flow_occ3_raw - flow_occ3_mean_12b) / flow_occ3_sd_12b,
    temp_7day_std = (temp_7day_raw - temp_7day_mean_12b) / temp_7day_sd_12b
  )

cat("\nflow_occ2 mean:", round(flow_occ2_mean_12b, 0),
    "SD:", round(flow_occ2_sd_12b, 0), "cfs\n")
cat("flow_occ3 mean:", round(flow_occ3_mean_12b, 0),
    "SD:", round(flow_occ3_sd_12b, 0), "cfs\n")
cat("temp_7day mean:", round(temp_7day_mean_12b, 2),
    "SD:", round(temp_7day_sd_12b, 2), "C\n")

# --- Save covariates ---
write.csv(
  fish_cov_12b %>%
    dplyr::select(animal_id, water_year,
                  flow_occ2_raw, flow_occ2_std,
                  flow_occ3_raw, flow_occ3_std,
                  temp_7day_raw, temp_7day_std),
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_model12b_covariates.csv",
  row.names = FALSE
)

# --- Save standardization parameters ---
write.csv(
  data.frame(
    covariate = c("flow_occ2", "flow_occ3", "temp_7day"),
    mean      = c(flow_occ2_mean_12b, flow_occ3_mean_12b, temp_7day_mean_12b),
    sd        = c(flow_occ2_sd_12b,   flow_occ3_sd_12b,   temp_7day_sd_12b)
  ),
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_model12b_std_params.csv",
  row.names = FALSE
)

cat("\nSaved:\n")
cat("  gs_model12b_covariates.csv\n")
cat("  gs_model12b_std_params.csv\n")
cat("Load these directly in model 12b script\n")

