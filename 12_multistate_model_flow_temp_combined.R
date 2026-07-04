#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - 12
# Script: 12_multistate_model_flow_temp_combined.R
# Author: Erin Tracy
# Last updated: June 2026
#
# PURPOSE:
# Combined flow and temperature model testing whether discharge and
# temperature independently predict routing, failure, and detection.
# Flow and temperature are moderately correlated (r = -0.397 for 30-day
# variables) but below the collinearity threshold, justifying joint testing.
#
# MODEL STRUCTURE:
# Routing occ2: psi_geo[i] = ilogit(alpha_geo + beta_geo * flow_occ2[i])
# Routing occ3: psi_ss[i]  = ilogit(alpha_ss  + beta_ss  * flow_occ3[i])
# Failure: phi_fail[i] = ilogit(alpha_fail + beta_fail_f * flow_30day[i]
#                                           + beta_fail_t * temp_7day[i])
# SS detection: p_ss_i[i] = ilogit(alpha_p_ss + beta_p_ss * flow_occ3[i])
#
# KEY QUESTIONS:
# 1. Does temperature add independent information about failure beyond flow?
# 2. Does flow remain significant when temperature is included?
# 3. Do the flow routing effects remain stable with temperature in the model?
#
# COVARIATES:
# flow_occ2:  same-day Rio Vista flow at occ2 arrival (routing)
# flow_occ3:  same-day GES flow at occ3 arrival (routing + SS detection)
# flow_30day: 30-day mean Rio Vista flow before Benicia (failure)
# temp_7day:  7-day mean Rio Vista temperature before Benicia (failure)
#==============================================================================

#==============================================================================
# SECTION 1: LOAD LIBRARIES AND DATA
#==============================================================================

library(nimble)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(dplyr)
library(lubridate)
library(tidyr)
library(zoo)

load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_multistate_data.RData")
load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_09_run1.RData")

# SAFETY CHECK
stopifnot(all(ch_mat_nimble %in% c(1, 2, 4, 5, 6, 7)))
cat("ch_mat_nimble pre-recoding values:",
    sort(unique(as.vector(ch_mat_nimble))), "\n")

#==============================================================================
# SECTION 2: RECODE DETECTION HISTORY TO 5-STATE MODEL
#==============================================================================

nstate <- 5

ch_mat_nimble_5state <- ch_mat_nimble
ch_mat_nimble_5state[ch_mat_nimble == 4] <- 3
ch_mat_nimble_5state[ch_mat_nimble == 5] <- 4
ch_mat_nimble_5state[ch_mat_nimble == 6] <- 5
ch_mat_nimble_5state[ch_mat_nimble == 7] <- 6

ch_mat_nimble <- ch_mat_nimble_5state

cat("Unique values after recoding:\n")
print(table(as.vector(ch_mat_nimble)))
cat("Total fish:", nrow(ch_mat_nimble), "\n")

#==============================================================================
# SECTION 3: BUILD ALL COVARIATES
# Flow: same-day RVB (occ2), same-day GES (occ3), 30-day RVB (failure)
# Temp: 7-day mean RVB before Benicia (failure only)
#==============================================================================

# --- Load flow data ---
rv_flow <- read.csv(
  "C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv"
)

rv_flow_complete <- rv_flow %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_cfs = value) %>%
  filter(!is.na(flow_cfs),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_cfs = zoo::na.approx(flow_cfs, na.rm = FALSE))

geo_flow <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/RioVista_Confluence_Flows.csv"
)

ges_flow_complete <- geo_flow %>%
  mutate(date = as.Date(dateTime, format = "%m/%d/%Y %H:%M")) %>%
  group_by(date) %>%
  dplyr::summarise(flow_GES = mean(GES, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(date),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_GES = zoo::na.approx(flow_GES, na.rm = FALSE))

# --- Load temperature data ---
temp_flow_raw <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/daily_tidalfilter_riovista_temp_flow.csv"
)

temp_flow_clean <- temp_flow_raw %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  dplyr::select(date, temp_c = daily_temp_c) %>%
  filter(!is.na(date),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(as.Date("2006-09-01"),
                             as.Date("2017-09-30"),
                             by = "day")) %>%
  mutate(temp_c = zoo::na.approx(temp_c, na.rm = FALSE))

cat("RV flow days:", nrow(rv_flow_complete),
    "| Missing:", sum(is.na(rv_flow_complete$flow_cfs)), "\n")
cat("GES flow days:", nrow(ges_flow_complete),
    "| Missing:", sum(is.na(ges_flow_complete$flow_GES)), "\n")
cat("Temp days:", nrow(temp_flow_clean),
    "| Missing:", sum(is.na(temp_flow_clean$temp_c)), "\n")

# --- Get arrival dates ---
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

# --- Extract all covariates ---
fish_flow_temp_12 <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(occ2_dates,    by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,    by = c("animal_id", "water_year")) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # Same-day Rio Vista flow at occ2
    flow_occ2_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(rv_flow_complete$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else rv_flow_complete$flow_cfs[idx]
      }
    },
    
    # Same-day GES flow at occ3
    flow_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(ges_flow_complete$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else ges_flow_complete$flow_GES[idx]
      }
    },
    
    # 30-day mean Rio Vista flow before Benicia
    flow_30day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- rv_flow_complete$flow_cfs[
          rv_flow_complete$date >= (benicia_date - 30) &
            rv_flow_complete$date <   benicia_date]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 7-day mean temperature before Benicia
    temp_7day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 7) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    }
    
  ) %>%
  ungroup()

cat("\nMissing flow_occ2:", sum(is.na(fish_flow_temp_12$flow_occ2_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_flow_temp_12$flow_occ3_raw)), "\n")
cat("Missing flow_30day:", sum(is.na(fish_flow_temp_12$flow_30day_raw)), "\n")
cat("Missing temp_7day:", sum(is.na(fish_flow_temp_12$temp_7day_raw)), "\n")

# --- Standardize ---
flow_occ2_mean_12  <- mean(fish_flow_temp_12$flow_occ2_raw,  na.rm = TRUE)
flow_occ2_sd_12    <- sd(fish_flow_temp_12$flow_occ2_raw,    na.rm = TRUE)
flow_occ3_mean_12  <- mean(fish_flow_temp_12$flow_occ3_raw,  na.rm = TRUE)
flow_occ3_sd_12    <- sd(fish_flow_temp_12$flow_occ3_raw,    na.rm = TRUE)
flow_30day_mean_12 <- mean(fish_flow_temp_12$flow_30day_raw, na.rm = TRUE)
flow_30day_sd_12   <- sd(fish_flow_temp_12$flow_30day_raw,   na.rm = TRUE)
temp_7day_mean_12  <- mean(fish_flow_temp_12$temp_7day_raw,  na.rm = TRUE)
temp_7day_sd_12    <- sd(fish_flow_temp_12$temp_7day_raw,    na.rm = TRUE)

fish_flow_temp_12 <- fish_flow_temp_12 %>%
  mutate(
    flow_occ2_std  = (flow_occ2_raw  - flow_occ2_mean_12)  / flow_occ2_sd_12,
    flow_occ3_std  = (flow_occ3_raw  - flow_occ3_mean_12)  / flow_occ3_sd_12,
    flow_30day_std = (flow_30day_raw - flow_30day_mean_12) / flow_30day_sd_12,
    temp_7day_std  = (temp_7day_raw  - temp_7day_mean_12)  / temp_7day_sd_12
  )

cat("\nStandardized covariate summaries:\n")
cat("Flow occ2 mean:", round(flow_occ2_mean_12, 0), "cfs SD:",
    round(flow_occ2_sd_12, 0), "\n")
cat("Flow occ3 mean:", round(flow_occ3_mean_12, 0), "cfs SD:",
    round(flow_occ3_sd_12, 0), "\n")
cat("Flow 30day mean:", round(flow_30day_mean_12, 0), "cfs SD:",
    round(flow_30day_sd_12, 0), "\n")
cat("Temp 7day mean:", round(temp_7day_mean_12, 2), "C SD:",
    round(temp_7day_sd_12, 2), "\n")

# check correlation between 30-day flow and 7-day temp in model 12 dataset
cat("\nCorrelation flow_30day vs temp_7day:\n")
cat("r =", round(cor(fish_flow_temp_12$flow_30day_raw,
                     fish_flow_temp_12$temp_7day_raw,
                     use = "complete.obs"), 3), "\n")

# verify row order
identical(fish_flow_temp_12$animal_id, detection_history$animal_id)
identical(as.numeric(fish_flow_temp_12$water_year),
          as.numeric(detection_history$water_year))

# extract vectors
flow_occ2_12  <- fish_flow_temp_12$flow_occ2_std
flow_occ3_12  <- fish_flow_temp_12$flow_occ3_std
flow_30day_12 <- fish_flow_temp_12$flow_30day_std
temp_7day_12  <- fish_flow_temp_12$temp_7day_std

flow_occ2_12[is.na(flow_occ2_12)]   <- 0
flow_occ3_12[is.na(flow_occ3_12)]   <- 0
flow_30day_12[is.na(flow_30day_12)] <- 0
temp_7day_12[is.na(temp_7day_12)]   <- 0

cat("\nflow_occ2 range:",  round(range(flow_occ2_12),  2), "\n")
cat("flow_occ3 range:",   round(range(flow_occ3_12),  2), "\n")
cat("flow_30day range:",  round(range(flow_30day_12), 2), "\n")
cat("temp_7day range:",   round(range(temp_7day_12),  2), "\n")

#==============================================================================
# SECTION 4: PLACEHOLDER MATRICES AND LIKELIHOOD CHECK
#==============================================================================

S_sac    <- 0.99; S_geo <- 0.95; S_ss <- 0.97
phi_fail <- 0.09; lambda <- 0.99
psi_geo_mean <- 0.124; psi_ss_mean <- 0.374
p_sac2 <- 0.977; p_sac3 <- 0.899; p_sac4 <- 0.910
p_sac5 <- 0.893; p_sac6 <- 0.986
p_geo  <- 0.771; p_ss   <- 0.690

temp_mat <- matrix(0, nrow = nstate + 1, ncol = nstate + 1)
temp_mat[nstate + 1, nstate + 1] <- 1

tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac*(1-psi_geo_mean)*(1-phi_fail),
                S_sac*psi_geo_mean*(1-phi_fail),
                0, (1-S_sac), S_sac*phi_fail, 0)
tr_12[2, 4] <- 1; tr_12[3, 4] <- 1
tr_12[4, 4] <- 1; tr_12[5, 5] <- 1

tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac*(1-psi_ss_mean)*(1-phi_fail),
                0, S_sac*psi_ss_mean*(1-phi_fail),
                (1-S_sac), S_sac*phi_fail, 0)
tr_23[2, 2] <- 1
tr_23[3, 4] <- 1; tr_23[4, 4] <- 1; tr_23[5, 5] <- 1

tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac*(1-phi_fail), 0, 0, (1-S_sac), S_sac*phi_fail, 0)
tr_34[2, ] <- c(S_geo*(1-phi_fail), 0, 0, (1-S_geo), S_geo*phi_fail, 0)
tr_34[3, ] <- c(0, 0, (1-phi_fail), 0, phi_fail, 0)
tr_34[4, 4] <- 1; tr_34[5, 5] <- 1

tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac*(1-phi_fail), 0, 0, (1-S_sac), S_sac*phi_fail, 0)
tr_45[2, 4] <- 1
tr_45[3, ] <- c(S_ss*(1-phi_fail), 0, 0, (1-S_ss), S_ss*phi_fail, 0)
tr_45[4, 4] <- 1; tr_45[5, 5] <- 1

tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac*(1-phi_fail), 0, 0, (1-S_sac), S_sac*phi_fail, 0)
tr_56[2, 4] <- 1; tr_56[3, 4] <- 1
tr_56[4, 4] <- 1; tr_56[5, 5] <- 1

tr_67 <- temp_mat
tr_67[1, ] <- c(lambda*(1-phi_fail), 0, 0, (1-lambda), lambda*phi_fail, 0)
tr_67[2, 4] <- 1; tr_67[3, 4] <- 1
tr_67[4, 4] <- 1; tr_67[5, 5] <- 1

tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

p_mat1 <- temp_mat; p_mat1[1, 1] <- 1
p_mat1[2, nstate+1] <- 1; p_mat1[3, nstate+1] <- 1
p_mat1[4, 4] <- 1; p_mat1[5, 5] <- 1

p_mat2 <- temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, 0, (1-p_sac2))
p_mat2[2, ] <- c(0, p_geo, 0, 0, 0, (1-p_geo))
p_mat2[3, nstate+1] <- 1
p_mat2[4, 4] <- 1; p_mat2[5, 5] <- 1

p_mat3 <- temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, 0, (1-p_sac3))
p_mat3[2, nstate+1] <- 1
p_mat3[3, ] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
p_mat3[4, 4] <- 1; p_mat3[5, 5] <- 1

p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, 0, (1-p_sac4))
p_mat4[2, nstate+1] <- 1
p_mat4[3, ] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
p_mat4[4, 4] <- 1; p_mat4[5, 5] <- 1

p_mat5 <- temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, 0, (1-p_sac5))
p_mat5[2, nstate+1] <- 1; p_mat5[3, nstate+1] <- 1
p_mat5[4, 4] <- 1; p_mat5[5, 5] <- 1

p_mat6 <- temp_mat
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, 0, (1-p_sac6))
p_mat6[2, nstate+1] <- 1; p_mat6[3, nstate+1] <- 1
p_mat6[4, 4] <- 1; p_mat6[5, 5] <- 1

p_mat7 <- temp_mat
p_mat7[1, 1] <- 1
p_mat7[2, nstate+1] <- 1; p_mat7[3, nstate+1] <- 1
p_mat7[4, 4] <- 1; p_mat7[5, 5] <- 1

p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4,
               p_mat5, p_mat6, p_mat7, along = 3)

rel_vec <- c(1, 0, 0, 0, 0, 0)

all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

cat("NaN count:", sum(is.nan(all_ll)), "\n")
cat("Inf count:", sum(is.infinite(all_ll)), "\n")
cat("LL range:", range(all_ll[!is.nan(all_ll) & !is.infinite(all_ll)]), "\n")

#==============================================================================
# SECTION 5: NIMBLE MODEL 12
# Combined flow routing + flow/temperature failure + flow SS detection
#==============================================================================

nimCode_12 <- nimbleCode({
  
  #--- PRIORS ---
  S_sac    ~ dbeta(1, 1)
  S_geo    ~ dbeta(1, 1)
  S_ss     ~ dbeta(1, 1)
  
  p_sac2   ~ dbeta(1, 1)
  p_sac3   ~ dbeta(1, 1)
  p_sac4   ~ dbeta(1, 1)
  p_sac5   ~ dbeta(1, 1)
  p_sac6   ~ dbeta(1, 1)
  p_geo    ~ dbeta(1, 1)
  
  lambda   ~ dbeta(1, 1)
  
  # Routing — same-day gauge-specific flow (from model 09)
  alpha_geo ~ dnorm(0, sd = 1.5)
  alpha_ss  ~ dnorm(0, sd = 1.5)
  beta_geo  ~ dnorm(0, sd = 1.0)
  beta_ss   ~ dnorm(0, sd = 1.0)
  
  # Failure — 30-day flow AND 7-day temperature
  alpha_fail  ~ dnorm(0, sd = 1.5)
  beta_fail_f ~ dnorm(0, sd = 1.0)   # flow slope — expected negative
  beta_fail_t ~ dnorm(0, sd = 1.0)   # temp slope — expected positive
  
  # SS detection — same-day GES flow (from model 09)
  alpha_p_ss ~ dnorm(0, sd = 1.5)
  beta_p_ss  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL PARAMETERS ---
  for(i in 1:nfish){
    
    # Routing — flow only
    psi_geo[i] <- ilogit(alpha_geo + beta_geo * flow_occ2[i])
    psi_ss[i]  <- ilogit(alpha_ss  + beta_ss  * flow_occ3[i])
    
    # Failure — flow + temperature
    phi_fail[i] <- ilogit(alpha_fail +
                            beta_fail_f * flow_30day[i] +
                            beta_fail_t * temp_7day[i])
    
    # SS detection — flow
    p_ss_i[i] <- ilogit(alpha_p_ss + beta_p_ss * flow_occ3[i])
    
  }
  
  #--- FIXED OBSERVATION ARRAY ---
  p_arr[1, 1:6, 1] <- c(1, 0, 0, 0, 0, 0)
  p_arr[2, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 1] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 2] <- c(p_sac2, 0, 0, 0, 0, (1-p_sac2))
  p_arr[2, 1:6, 2] <- c(0, p_geo, 0, 0, 0, (1-p_geo))
  p_arr[3, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 2] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 3] <- c(p_sac3, 0, 0, 0, 0, (1-p_sac3))
  p_arr[2, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 3] <- c(0, 0, 0.5, 0, 0, 0.5)
  p_arr[4, 1:6, 3] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 3] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 4] <- c(p_sac4, 0, 0, 0, 0, (1-p_sac4))
  p_arr[2, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 4] <- c(0, 0, p_sac4, 0, 0, (1-p_sac4))
  p_arr[4, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 4] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 5] <- c(p_sac5, 0, 0, 0, 0, (1-p_sac5))
  p_arr[2, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 5] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 6] <- c(p_sac6, 0, 0, 0, 0, (1-p_sac6))
  p_arr[2, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 6] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 7] <- c(1, 0, 0, 0, 0, 0)
  p_arr[2, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 7] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 7] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  
  #--- LIKELIHOOD ---
  for(i in 1:nfish){
    
    p_arr_i[i, 1, 1:6, 1] <- p_arr[1, 1:6, 1]
    p_arr_i[i, 2, 1:6, 1] <- p_arr[2, 1:6, 1]
    p_arr_i[i, 3, 1:6, 1] <- p_arr[3, 1:6, 1]
    p_arr_i[i, 4, 1:6, 1] <- p_arr[4, 1:6, 1]
    p_arr_i[i, 5, 1:6, 1] <- p_arr[5, 1:6, 1]
    p_arr_i[i, 6, 1:6, 1] <- p_arr[6, 1:6, 1]
    
    p_arr_i[i, 1, 1:6, 2] <- p_arr[1, 1:6, 2]
    p_arr_i[i, 2, 1:6, 2] <- p_arr[2, 1:6, 2]
    p_arr_i[i, 3, 1:6, 2] <- p_arr[3, 1:6, 2]
    p_arr_i[i, 4, 1:6, 2] <- p_arr[4, 1:6, 2]
    p_arr_i[i, 5, 1:6, 2] <- p_arr[5, 1:6, 2]
    p_arr_i[i, 6, 1:6, 2] <- p_arr[6, 1:6, 2]
    
    # Occasion 3 — fish-specific p_ss
    p_arr_i[i, 1, 1:6, 3] <- p_arr[1, 1:6, 3]
    p_arr_i[i, 2, 1:6, 3] <- p_arr[2, 1:6, 3]
    p_arr_i[i, 3, 1:6, 3] <- c(0, 0, p_ss_i[i], 0, 0, (1-p_ss_i[i]))
    p_arr_i[i, 4, 1:6, 3] <- p_arr[4, 1:6, 3]
    p_arr_i[i, 5, 1:6, 3] <- p_arr[5, 1:6, 3]
    p_arr_i[i, 6, 1:6, 3] <- p_arr[6, 1:6, 3]
    
    p_arr_i[i, 1, 1:6, 4] <- p_arr[1, 1:6, 4]
    p_arr_i[i, 2, 1:6, 4] <- p_arr[2, 1:6, 4]
    p_arr_i[i, 3, 1:6, 4] <- p_arr[3, 1:6, 4]
    p_arr_i[i, 4, 1:6, 4] <- p_arr[4, 1:6, 4]
    p_arr_i[i, 5, 1:6, 4] <- p_arr[5, 1:6, 4]
    p_arr_i[i, 6, 1:6, 4] <- p_arr[6, 1:6, 4]
    
    p_arr_i[i, 1, 1:6, 5] <- p_arr[1, 1:6, 5]
    p_arr_i[i, 2, 1:6, 5] <- p_arr[2, 1:6, 5]
    p_arr_i[i, 3, 1:6, 5] <- p_arr[3, 1:6, 5]
    p_arr_i[i, 4, 1:6, 5] <- p_arr[4, 1:6, 5]
    p_arr_i[i, 5, 1:6, 5] <- p_arr[5, 1:6, 5]
    p_arr_i[i, 6, 1:6, 5] <- p_arr[6, 1:6, 5]
    
    p_arr_i[i, 1, 1:6, 6] <- p_arr[1, 1:6, 6]
    p_arr_i[i, 2, 1:6, 6] <- p_arr[2, 1:6, 6]
    p_arr_i[i, 3, 1:6, 6] <- p_arr[3, 1:6, 6]
    p_arr_i[i, 4, 1:6, 6] <- p_arr[4, 1:6, 6]
    p_arr_i[i, 5, 1:6, 6] <- p_arr[5, 1:6, 6]
    p_arr_i[i, 6, 1:6, 6] <- p_arr[6, 1:6, 6]
    
    p_arr_i[i, 1, 1:6, 7] <- p_arr[1, 1:6, 7]
    p_arr_i[i, 2, 1:6, 7] <- p_arr[2, 1:6, 7]
    p_arr_i[i, 3, 1:6, 7] <- p_arr[3, 1:6, 7]
    p_arr_i[i, 4, 1:6, 7] <- p_arr[4, 1:6, 7]
    p_arr_i[i, 5, 1:6, 7] <- p_arr[5, 1:6, 7]
    p_arr_i[i, 6, 1:6, 7] <- p_arr[6, 1:6, 7]
    
    # Transition matrices
    tr_arr_i[i, 1, 1:6, 1] <- c(S_sac*(1-psi_geo[i])*(1-phi_fail[i]),
                                S_sac*psi_geo[i]*(1-phi_fail[i]),
                                0, (1-S_sac), S_sac*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 1] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 2] <- c(S_sac*(1-psi_ss[i])*(1-phi_fail[i]),
                                0,
                                S_sac*psi_ss[i]*(1-phi_fail[i]),
                                (1-S_sac), S_sac*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 2] <- c(0, 1, 0, 0, 0, 0)
    tr_arr_i[i, 3, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 2] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 3] <- c(S_sac*(1-phi_fail[i]), 0, 0,
                                (1-S_sac), S_sac*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 3] <- c(S_geo*(1-phi_fail[i]), 0, 0,
                                (1-S_geo), S_geo*phi_fail[i], 0)
    tr_arr_i[i, 3, 1:6, 3] <- c(0, 0, (1-phi_fail[i]), 0, phi_fail[i], 0)
    tr_arr_i[i, 4, 1:6, 3] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 3] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 4] <- c(S_sac*(1-phi_fail[i]), 0, 0,
                                (1-S_sac), S_sac*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 4] <- c(S_ss*(1-phi_fail[i]), 0, 0,
                                (1-S_ss), S_ss*phi_fail[i], 0)
    tr_arr_i[i, 4, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 4] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 5] <- c(S_sac*(1-phi_fail[i]), 0, 0,
                                (1-S_sac), S_sac*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 5] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 6] <- c(lambda*(1-phi_fail[i]), 0, 0,
                                (1-lambda), lambda*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 6] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
    
    ch_mat[i, 1:7] ~ dDHMMo(
      init         = c(1, 0, 0, 0, 0, 0)[1:6],
      probObs      = p_arr_i[i, 1:6, 1:6, 1:7],
      probTrans    = tr_arr_i[i, 1:6, 1:6, 1:6],
      len          = 7,
      checkRowSums = 0
    )
    
  } #end i loop
})

#==============================================================================
# SECTION 6: BUILD MODEL AND CONFIGURE MCMC
#==============================================================================

inits_12 <- list(
  S_sac       = 0.99, S_geo = 0.95, S_ss = 0.97,
  p_sac2      = 0.977, p_sac3 = 0.899,
  p_sac4      = 0.910, p_sac5 = 0.893, p_sac6 = 0.986,
  p_geo       = 0.771, lambda = 0.986,
  alpha_geo   = -1.95, beta_geo   = 0.0,
  alpha_ss    = -0.51, beta_ss    = 0.0,
  alpha_fail  = -2.30,
  beta_fail_f = 0.0,
  beta_fail_t = 0.0,
  alpha_p_ss  = -0.39, beta_p_ss  = 0.0
)

nimMod_12 <- nimbleModel(
  code      = nimCode_12,
  inits     = inits_12,
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(
    nfish      = nrow(ch_mat_nimble),
    flow_occ2  = flow_occ2_12,
    flow_occ3  = flow_occ3_12,
    flow_30day = flow_30day_12,
    temp_7day  = temp_7day_12
  )
)

nimMod_12$calculate()

inits_fn_12 <- function(){
  list(
    S_sac       = runif(1, 0.95, 1.0),
    S_geo       = runif(1, 0.80, 1.0),
    S_ss        = runif(1, 0.85, 1.0),
    p_sac2      = runif(1, 0.90, 1.0),
    p_sac3      = runif(1, 0.75, 0.95),
    p_sac4      = runif(1, 0.80, 0.95),
    p_sac5      = runif(1, 0.80, 0.95),
    p_sac6      = runif(1, 0.95, 1.0),
    p_geo       = runif(1, 0.50, 0.90),
    lambda      = runif(1, 0.95, 1.0),
    alpha_geo   = rnorm(1, -1.95, 0.3), beta_geo   = rnorm(1, 0, 0.2),
    alpha_ss    = rnorm(1, -0.51, 0.3), beta_ss    = rnorm(1, 0, 0.2),
    alpha_fail  = rnorm(1, -2.30, 0.3),
    beta_fail_f = rnorm(1, 0, 0.2),
    beta_fail_t = rnorm(1, 0, 0.2),
    alpha_p_ss  = rnorm(1, -0.39, 0.3), beta_p_ss = rnorm(1, 0, 0.2)
  )
}

params_12 <- c(
  "S_sac", "S_geo", "S_ss",
  "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
  "p_geo", "lambda",
  "alpha_geo", "beta_geo",
  "alpha_ss",  "beta_ss",
  "alpha_fail", "beta_fail_f", "beta_fail_t",
  "alpha_p_ss", "beta_p_ss"
)

confMCMC_12 <- configureMCMC(nimMod_12, onlySlice = TRUE)
confMCMC_12$addMonitors(params_12)
MCMC_12   <- buildMCMC(confMCMC_12)
CModel_12 <- compileNimble(nimMod_12)
CMCMC_12  <- compileNimble(MCMC_12, project = CModel_12)

#==============================================================================
# SECTION 7: RUN MCMC AND SAVE
#==============================================================================

# Short test run first
# mcmc_test_12 <- runMCMC(CMCMC_12, niter = 1000, nchains = 1,
#                          nburnin = 100, thin = 1,
#                          inits = list(inits_fn_12()),
#                          samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test_12, round = 3)

mcmc_out_12 <- runMCMC(
  CMCMC_12,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn_12(), inits_fn_12(), inits_fn_12()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_12, round = 3)

save(mcmc_out_12,
     flow_occ2_12, flow_occ3_12, flow_30day_12, temp_7day_12,
     fish_flow_temp_12,
     flow_occ2_mean_12, flow_occ2_sd_12,
     flow_occ3_mean_12, flow_occ3_sd_12,
     flow_30day_mean_12, flow_30day_sd_12,
     temp_7day_mean_12, temp_7day_sd_12,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_temp_12_run1.RData")

MCMCtrace(mcmc_out_12,
          params   = c("alpha_geo", "beta_geo",
                       "alpha_ss",  "beta_ss",
                       "alpha_fail", "beta_fail_f", "beta_fail_t",
                       "alpha_p_ss", "beta_p_ss"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_flow_temp_12_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)

cat("Model 12 complete\n")
