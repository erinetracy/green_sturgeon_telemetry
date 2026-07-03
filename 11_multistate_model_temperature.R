#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - 11
# Script: 11_multistate_model_temperature.R
# Author: Erin Tracy
# Last updated: June 2026
#
# PURPOSE:
# Temperature-only sensitivity analysis parallel to flow model (model 09).
# Tests whether pre-migration water temperature predicts routing probability,
# migration failure, and detection probability independently of discharge.
#
# TEMPERATURE COVARIATE STRUCTURE (from temporal scale screening):
# psi_geo:  NO temperature covariate — no meaningful signal in screening
# psi_ss:   30-day mean temperature before Benicia passage
#           (signal strengthens with longer windows, r=0.985 with same-day
#            so single antecedent covariate used for both routing decisions)
# phi_fail: 7-day mean temperature before Benicia passage
#           (strongest predictor, pseudo R2=0.443, AUC=0.901)
# p_ss_i:   same-day GES flow retained from model 09
#           (flow-dependent SS detection kept — not replaced by temperature)
#
# EXPECTED DIRECTIONS (from screening):
# beta_ss_t  < 0: higher temperature reduces SS route use
# beta_fail_t > 0: higher temperature increases migration failure probability
#
# BIOLOGICAL RATIONALE:
# Failed migrants experienced mean 7-day temperature of 19.2C vs 13.7C for
# successful migrants. Warm low-flow conditions (late spring/early summer)
# are associated with migration failure; cool high-flow conditions with
# success. Temperature and flow carry partially independent information
# (r = -0.397) justifying separate model testing before combining.
#
# COMPARISON:
# Model 09: flow-only (routing + failure + SS detection)
# Model 11: temperature routing/failure + flow SS detection (this model)
# Model 12: flow + temperature combined (if model 11 shows significant effects)
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
load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/temp_covariates.RData")

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
print(table(detection_history$status))

#==============================================================================
# SECTION 3: BUILD TEMPERATURE COVARIATES
# 7-day for failure, 30-day for SS routing
# Both anchored to Benicia passage date
# Same-day GES flow retained for SS detection from model 09
#==============================================================================

temp_flow_raw <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/daily_tidalfilter_riovista_temp_flow.csv"
)

temp_flow_clean <- temp_flow_raw %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  dplyr::select(date, temp_c = daily_temp_c, flow_cfs) %>%
  filter(!is.na(date),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(as.Date("2006-09-01"),
                             as.Date("2017-09-30"),
                             by = "day")) %>%
  mutate(
    temp_c   = zoo::na.approx(temp_c,   na.rm = FALSE),
    flow_cfs = zoo::na.approx(flow_cfs, na.rm = FALSE)
  )

cat("Temperature data coverage:\n")
cat("Total days:", nrow(temp_flow_clean), "\n")
cat("Missing temp_c:", sum(is.na(temp_flow_clean$temp_c)), "\n")

# GES flow for same-day SS detection (retained from model 09)
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

# arrival dates
benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(benicia_date = as.Date(min(first_detection)), .groups = "drop")

occ3_dates <- model1_events %>%
  filter(occasion == 3) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(arr_date_occ3 = as.Date(min(first_detection)), .groups = "drop")

# build covariates
fish_temp_11 <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,   by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # 7-day mean temperature before Benicia — for failure
    temp_7day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 7) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 30-day mean temperature before Benicia — for SS routing
    temp_30day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 30) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # same-day GES flow at occ3 — for SS detection (retained from model 09)
    flow_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(ges_flow_complete$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else ges_flow_complete$flow_GES[idx]
      }
    }
    
  ) %>%
  ungroup()

cat("\nCovariate summary:\n")
cat("Missing temp_7day:", sum(is.na(fish_temp_11$temp_7day_raw)), "\n")
cat("Missing temp_30day:", sum(is.na(fish_temp_11$temp_30day_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_temp_11$flow_occ3_raw)), "\n")

# standardize
temp_7day_mean_11  <- mean(fish_temp_11$temp_7day_raw,  na.rm = TRUE)
temp_7day_sd_11    <- sd(fish_temp_11$temp_7day_raw,    na.rm = TRUE)
temp_30day_mean_11 <- mean(fish_temp_11$temp_30day_raw, na.rm = TRUE)
temp_30day_sd_11   <- sd(fish_temp_11$temp_30day_raw,   na.rm = TRUE)
flow_occ3_mean_11  <- mean(fish_temp_11$flow_occ3_raw,  na.rm = TRUE)
flow_occ3_sd_11    <- sd(fish_temp_11$flow_occ3_raw,    na.rm = TRUE)

fish_temp_11 <- fish_temp_11 %>%
  mutate(
    temp_7day_std  = (temp_7day_raw  - temp_7day_mean_11)  / temp_7day_sd_11,
    temp_30day_std = (temp_30day_raw - temp_30day_mean_11) / temp_30day_sd_11,
    flow_occ3_std  = (flow_occ3_raw  - flow_occ3_mean_11)  / flow_occ3_sd_11
  )

cat("\nStandardized covariate summaries:\n")
cat("Mean raw 7-day temp:", round(temp_7day_mean_11, 2), "C\n")
cat("SD raw 7-day temp:",   round(temp_7day_sd_11, 2),   "C\n")
cat("Mean raw 30-day temp:", round(temp_30day_mean_11, 2), "C\n")
cat("SD raw 30-day temp:",   round(temp_30day_sd_11, 2),   "C\n")
cat("Mean raw GES flow:", round(flow_occ3_mean_11, 0), "cfs\n")
cat("SD raw GES flow:",   round(flow_occ3_sd_11, 0),   "cfs\n")

print(summary(fish_temp_11$temp_7day_std))
print(summary(fish_temp_11$temp_30day_std))
print(summary(fish_temp_11$flow_occ3_std))

# verify row order
identical(fish_temp_11$animal_id, detection_history$animal_id)
identical(as.numeric(fish_temp_11$water_year),
          as.numeric(detection_history$water_year))

# extract vectors
temp_7day_11  <- fish_temp_11$temp_7day_std
temp_30day_11 <- fish_temp_11$temp_30day_std
flow_occ3_11  <- fish_temp_11$flow_occ3_std

temp_7day_11[is.na(temp_7day_11)]   <- 0
temp_30day_11[is.na(temp_30day_11)] <- 0
flow_occ3_11[is.na(flow_occ3_11)]   <- 0

cat("\nFinal vectors:\n")
cat("temp_7day range:",  round(range(temp_7day_11),  2), "\n")
cat("temp_30day range:", round(range(temp_30day_11), 2), "\n")
cat("flow_occ3 range:",  round(range(flow_occ3_11),  2), "\n")

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
# SECTION 5: NIMBLE MODEL 11
# Temperature-only routing and failure
# Flow-dependent SS detection retained from model 09
#==============================================================================

nimCode_11 <- nimbleCode({
  
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
  
  # Georgiana routing — no temperature covariate, fixed intercept only
  psi_geo_mu ~ dbeta(1, 1)
  
  # SS routing — 30-day temperature
  alpha_ss_t  ~ dnorm(0, sd = 1.5)
  beta_ss_t   ~ dnorm(0, sd = 1.0)
  
  # Failure — 7-day temperature
  # Expected positive: higher temp increases failure
  alpha_fail_t ~ dnorm(0, sd = 1.5)
  beta_fail_t  ~ dnorm(0, sd = 1.0)
  
  # SS detection — same-day GES flow (retained from model 09)
  alpha_p_ss ~ dnorm(0, sd = 1.5)
  beta_p_ss  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL PARAMETERS ---
  for(i in 1:nfish){
    
    # Georgiana routing — constant (no temp effect)
    psi_geo[i] <- psi_geo_mu
    
    # SS routing — 30-day temperature
    psi_ss[i]   <- ilogit(alpha_ss_t  + beta_ss_t   * temp_30day[i])
    
    # Failure — 7-day temperature
    phi_fail[i] <- ilogit(alpha_fail_t + beta_fail_t * temp_7day[i])
    
    # SS detection — same-day GES flow
    p_ss_i[i]   <- ilogit(alpha_p_ss  + beta_p_ss   * flow_occ3[i])
    
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
  
  # Occasion 3 — SS row overridden per fish
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

inits_11 <- list(
  S_sac        = 0.99, S_geo = 0.95, S_ss = 0.97,
  p_sac2       = 0.977, p_sac3 = 0.899,
  p_sac4       = 0.910, p_sac5 = 0.893, p_sac6 = 0.986,
  p_geo        = 0.771, lambda = 0.986,
  psi_geo_mu   = 0.124,
  alpha_ss_t   = -0.51, beta_ss_t   = 0.0,
  alpha_fail_t = -2.30, beta_fail_t = 0.0,
  alpha_p_ss   = -0.39, beta_p_ss   = 0.0
)

nimMod_11 <- nimbleModel(
  code      = nimCode_11,
  inits     = inits_11,
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(
    nfish      = nrow(ch_mat_nimble),
    temp_7day  = temp_7day_11,
    temp_30day = temp_30day_11,
    flow_occ3  = flow_occ3_11
  )
)

nimMod_11$calculate()

inits_fn_11 <- function(){
  list(
    S_sac        = runif(1, 0.95, 1.0),
    S_geo        = runif(1, 0.80, 1.0),
    S_ss         = runif(1, 0.85, 1.0),
    p_sac2       = runif(1, 0.90, 1.0),
    p_sac3       = runif(1, 0.75, 0.95),
    p_sac4       = runif(1, 0.80, 0.95),
    p_sac5       = runif(1, 0.80, 0.95),
    p_sac6       = runif(1, 0.95, 1.0),
    p_geo        = runif(1, 0.50, 0.90),
    lambda       = runif(1, 0.95, 1.0),
    psi_geo_mu   = runif(1, 0.05, 0.25),
    alpha_ss_t   = rnorm(1, -0.51, 0.3), beta_ss_t   = rnorm(1, 0, 0.2),
    alpha_fail_t = rnorm(1, -2.30, 0.3), beta_fail_t = rnorm(1, 0, 0.2),
    alpha_p_ss   = rnorm(1, -0.39, 0.3), beta_p_ss   = rnorm(1, 0, 0.2)
  )
}

params_11 <- c(
  "S_sac", "S_geo", "S_ss",
  "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
  "p_geo", "lambda",
  "psi_geo_mu",
  "alpha_ss_t",   "beta_ss_t",
  "alpha_fail_t", "beta_fail_t",
  "alpha_p_ss",   "beta_p_ss"
)

confMCMC_11 <- configureMCMC(nimMod_11, onlySlice = TRUE)
confMCMC_11$addMonitors(params_11)
MCMC_11   <- buildMCMC(confMCMC_11)
CModel_11 <- compileNimble(nimMod_11)
CMCMC_11  <- compileNimble(MCMC_11, project = CModel_11)

#==============================================================================
# SECTION 7: RUN MCMC AND SAVE
#==============================================================================

# Short test run first
# mcmc_test_11 <- runMCMC(CMCMC_11, niter = 1000, nchains = 1,
#                          nburnin = 100, thin = 1,
#                          inits = list(inits_fn_11()),
#                          samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test_11, round = 3)

mcmc_out_11 <- runMCMC(
  CMCMC_11,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn_11(), inits_fn_11(), inits_fn_11()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_11, round = 3)

save(mcmc_out_11,
     temp_7day_11, temp_30day_11, flow_occ3_11,
     fish_temp_11,
     temp_7day_mean_11, temp_7day_sd_11,
     temp_30day_mean_11, temp_30day_sd_11,
     flow_occ3_mean_11, flow_occ3_sd_11,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_temp_11_run1.RData")

MCMCtrace(mcmc_out_11,
          params   = c("psi_geo_mu",
                       "alpha_ss_t",   "beta_ss_t",
                       "alpha_fail_t", "beta_fail_t",
                       "alpha_p_ss",   "beta_p_ss"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_temp_11_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)

cat("Model 11 complete\n")

