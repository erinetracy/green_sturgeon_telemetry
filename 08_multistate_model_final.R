#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - 08
# Script: 08_multistate_model_final.R
# Author: Erin Tracy
# Last updated: June 2026
#
# PURPOSE:
# Final model combining gauge-specific same-day junction flow for routing
# and 30-day antecedent Rio Vista flow for failure probability.
#
# BIOLOGICAL JUSTIFICATION FOR MIXED TEMPORAL SCALES:
# Routing decisions reflect instantaneous hydraulic conditions at the
# junction at the moment of the behavioral decision — fish respond to
# local flow signals at Georgiana Slough and Steamboat/Sutter junctions.
# Migration failure reflects sustained hydrological context accumulated
# during weeks of estuarine staging before freshwater entry — fish
# physiological condition and motivation integrate flow conditions over
# longer timescales. Using different temporal scales for routing vs
# failure is therefore biologically justified and supported by model
# comparison: same-day flow best predicted routing (models 06, 07) while
# 30-day antecedent flow produced the strongest and most precise failure
# effect (model 07b: beta_fail = -1.412, CI: -2.203 to -0.753 vs
# model 07: beta_fail = -0.257, CI: -0.590 to 0.027).
#
# MODEL STRUCTURE:
# Routing at occ2: psi_geo[i] = ilogit(alpha_geo + beta_geo * flow_occ2[i])
#   flow_occ2 = same-day Rio Vista discharge at occ2 arrival
# Routing at occ3: psi_ss[i]  = ilogit(alpha_ss  + beta_ss  * flow_occ3[i])
#   flow_occ3 = same-day GES discharge at occ3 arrival
# Failure:    phi_fail[i] = ilogit(alpha_fail + beta_fail * flow_30day[i])
#   flow_30day = 30-day mean Rio Vista discharge before Benicia passage
#
# COMPARISON MODELS:
# Model 05: baseline, no flow covariates
# Model 08: this model — primary result for publication
# Models 06, 07, 07b: supplementary sensitivity analyses
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

cat("Total fish:", nrow(ch_mat_nimble), "\n")
print(table(detection_history$status))

#==============================================================================
# SECTION 3: BUILD FLOW COVARIATES
#
# Three flow covariates:
# 1. flow_occ2: same-day Rio Vista flow at occ2 arrival (Georgiana routing)
# 2. flow_occ3: same-day GES flow at occ3 arrival (SS routing)
# 3. flow_30day: 30-day mean Rio Vista flow before Benicia (failure)
#==============================================================================

# --- Load and clean Rio Vista flow ---
rv_flow <- read.csv("C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv")

rv_flow_complete <- rv_flow %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_cfs = value) %>%
  filter(!is.na(flow_cfs),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_cfs = zoo::na.approx(flow_cfs, na.rm = FALSE))

cat("Rio Vista flow days:", nrow(rv_flow_complete), "\n")
cat("NAs:", sum(is.na(rv_flow_complete$flow_cfs)), "\n")

# --- Load and clean GES flow ---
geo_flow <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/RioVista_Confluence_Flows.csv")

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

cat("GES flow days:", nrow(ges_flow_complete), "\n")
cat("NAs:", sum(is.na(ges_flow_complete$flow_GES)), "\n")

# --- Get arrival dates at occ2 and occ3 ---
occ2_dates <- model1_events %>%
  filter(occasion == 2) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    arr_date_occ2 = as.Date(min(first_detection)),
    .groups = "drop"
  )

occ3_dates <- model1_events %>%
  filter(occasion == 3) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    arr_date_occ3 = as.Date(min(first_detection)),
    .groups = "drop"
  )

# --- Get Benicia passage dates for 30-day window ---
benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    benicia_date = as.Date(min(first_detection)),
    .groups = "drop"
  )

cat("Fish with occ2 dates:", nrow(occ2_dates), "\n")
cat("Fish with occ3 dates:", nrow(occ3_dates), "\n")
cat("Fish with Benicia dates:", nrow(benicia_dates), "\n")

# --- Join and calculate all three flow covariates ---
fish_flow_08 <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(occ2_dates,    by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,    by = c("animal_id", "water_year")) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # 1. Same-day Rio Vista flow at occ2
    flow_occ2_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(rv_flow_complete$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else rv_flow_complete$flow_cfs[idx]
      }
    },
    
    # 2. Same-day GES flow at occ3
    flow_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(ges_flow_complete$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else ges_flow_complete$flow_GES[idx]
      }
    },
    
    # 3. 30-day mean Rio Vista flow before Benicia
    flow_30day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- rv_flow_complete$flow_cfs[
          rv_flow_complete$date >= (benicia_date - 30) &
            rv_flow_complete$date <   benicia_date]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    }
    
  ) %>%
  ungroup()

cat("\nMissing flow_occ2:", sum(is.na(fish_flow_08$flow_occ2_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_flow_08$flow_occ3_raw)), "\n")
cat("Missing flow_30day:", sum(is.na(fish_flow_08$flow_30day_raw)), "\n")

# --- Standardize all three covariates ---
flow_occ2_mean_08 <- mean(fish_flow_08$flow_occ2_raw, na.rm = TRUE)
flow_occ2_sd_08   <- sd(fish_flow_08$flow_occ2_raw,   na.rm = TRUE)
flow_occ3_mean_08 <- mean(fish_flow_08$flow_occ3_raw, na.rm = TRUE)
flow_occ3_sd_08   <- sd(fish_flow_08$flow_occ3_raw,   na.rm = TRUE)
flow_30day_mean_08 <- mean(fish_flow_08$flow_30day_raw, na.rm = TRUE)
flow_30day_sd_08   <- sd(fish_flow_08$flow_30day_raw,   na.rm = TRUE)

fish_flow_08 <- fish_flow_08 %>%
  mutate(
    flow_occ2_std  = (flow_occ2_raw  - flow_occ2_mean_08)  / flow_occ2_sd_08,
    flow_occ3_std  = (flow_occ3_raw  - flow_occ3_mean_08)  / flow_occ3_sd_08,
    flow_30day_std = (flow_30day_raw - flow_30day_mean_08) / flow_30day_sd_08
  )

cat("\nFlow occ2 (Rio Vista same-day):\n")
cat("Mean:", round(flow_occ2_mean_08, 0), "cfs | SD:", round(flow_occ2_sd_08, 0), "cfs\n")
cat("\nFlow occ3 (GES same-day):\n")
cat("Mean:", round(flow_occ3_mean_08, 0), "cfs | SD:", round(flow_occ3_sd_08, 0), "cfs\n")
cat("\nFlow 30-day (Rio Vista antecedent):\n")
cat("Mean:", round(flow_30day_mean_08, 0), "cfs | SD:", round(flow_30day_sd_08, 0), "cfs\n")

# --- Verify row order ---
identical(fish_flow_08$animal_id, detection_history$animal_id)
identical(as.numeric(fish_flow_08$water_year),
          as.numeric(detection_history$water_year))

# --- Extract vectors for NIMBLE ---
flow_occ2  <- fish_flow_08$flow_occ2_std;  flow_occ2[is.na(flow_occ2)]   <- 0
flow_occ3  <- fish_flow_08$flow_occ3_std;  flow_occ3[is.na(flow_occ3)]   <- 0
flow_30day <- fish_flow_08$flow_30day_std; flow_30day[is.na(flow_30day)] <- 0

cat("\nflow_occ2 range:",  round(range(flow_occ2),  2), "\n")
cat("flow_occ3 range:",   round(range(flow_occ3),  2), "\n")
cat("flow_30day range:",  round(range(flow_30day), 2), "\n")

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
# SECTION 5: NIMBLE MODEL
#
# Three flow covariates — different temporal scales for different processes:
#   flow_occ2:  same-day Rio Vista — Georgiana routing decision
#   flow_occ3:  same-day GES      — Steamboat/Sutter routing decision
#   flow_30day: 30-day Rio Vista  — migration failure probability
#==============================================================================

nimCode_08 <- nimbleCode({
  
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
  p_ss     ~ dbeta(1, 1)
  
  lambda   ~ dbeta(1, 1)
  
  # Routing — same-day gauge-specific flow
  alpha_geo ~ dnorm(0, sd = 1.5)
  alpha_ss  ~ dnorm(0, sd = 1.5)
  beta_geo  ~ dnorm(0, sd = 1.0)
  beta_ss   ~ dnorm(0, sd = 1.0)
  
  # Failure — 30-day antecedent flow
  alpha_fail ~ dnorm(0, sd = 1.5)
  beta_fail  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL PARAMETERS ---
  for(i in 1:nfish){
    
    # Routing at occ2 — same-day Rio Vista flow
    psi_geo[i] <- ilogit(alpha_geo + beta_geo * flow_occ2[i])
    
    # Routing at occ3 — same-day GES flow
    psi_ss[i]  <- ilogit(alpha_ss  + beta_ss  * flow_occ3[i])
    
    # Failure — 30-day antecedent Rio Vista flow
    phi_fail[i] <- ilogit(alpha_fail + beta_fail * flow_30day[i])
    
  }
  
  #--- OBSERVATION ARRAY ---
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
  p_arr[3, 1:6, 3] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[4, 1:6, 3] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 3] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 4] <- c(p_sac4, 0, 0, 0, 0, (1-p_sac4))
  p_arr[2, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 4] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[4, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 4] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  
  p_arr[1, 1:6, 5] <- c(p_sac5, 0, 0, 0, 0, (1-p_sac5))
  p_arr[2, 1:6, 5] <- c(0