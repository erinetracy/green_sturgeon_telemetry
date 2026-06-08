#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - GAUGE-SPECIFIC FLOW
# Script: 07_multistate_model_flow_gaugespecific.R
# Author: Erin Tracy
# Last updated: June 2026
#
# PURPOSE:
# Extension of model 06 with two changes:
# 1. Gauge-specific flow covariates:
#    - occ2 routing (Georgiana decision): Rio Vista same-day flow
#    - occ3 routing (SS decision): GES same-day flow
# 2. Flow-dependent migration failure:
#    logit(phi_fail[i]) = alpha_fail + beta_fail * flow_occ2[i]
#    Hypothesis: low flow increases failure probability
#    Expected sign: beta_fail < 0
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
# SECTION 3: BUILD GAUGE-SPECIFIC FLOW COVARIATES
#
# flow_occ2[i]: Rio Vista same-day flow on date fish i arrived at occ2
#               (Georgiana routing decision)
# flow_occ3[i]: GES same-day flow on date fish i arrived at occ3
#               (Steamboat/Sutter routing decision)
# flow_fail[i]: Rio Vista same-day flow on date fish i arrived at occ2
#               (failure probability — same as occ2 flow, migration entry)
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

cat("Fish with occ2 dates:", nrow(occ2_dates), "\n")
cat("Fish with occ3 dates:", nrow(occ3_dates), "\n")

# --- Join to detection history and extract flow ---
fish_flow_07 <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(occ2_dates, by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    # Rio Vista flow on day fish arrived at occ2
    flow_occ2_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(rv_flow_complete$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else rv_flow_complete$flow_cfs[idx]
      }
    },
    # GES flow on day fish arrived at occ3
    flow_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(ges_flow_complete$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else ges_flow_complete$flow_GES[idx]
      }
    }
  ) %>%
  ungroup()

cat("\nMissing flow_occ2:", sum(is.na(fish_flow_07$flow_occ2_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_flow_07$flow_occ3_raw)), "\n")

# --- Standardize both flow covariates ---
flow_occ2_mean <- mean(fish_flow_07$flow_occ2_raw, na.rm = TRUE)
flow_occ2_sd   <- sd(fish_flow_07$flow_occ2_raw,   na.rm = TRUE)
flow_occ3_mean <- mean(fish_flow_07$flow_occ3_raw, na.rm = TRUE)
flow_occ3_sd   <- sd(fish_flow_07$flow_occ3_raw,   na.rm = TRUE)

fish_flow_07 <- fish_flow_07 %>%
  mutate(
    flow_occ2_std = (flow_occ2_raw - flow_occ2_mean) / flow_occ2_sd,
    flow_occ3_std = (flow_occ3_raw - flow_occ3_mean) / flow_occ3_sd
  )

cat("\nFlow occ2 summary (Rio Vista, standardized):\n")
print(summary(fish_flow_07$flow_occ2_std))
cat("\nFlow occ3 summary (GES, standardized):\n")
print(summary(fish_flow_07$flow_occ3_std))

# --- Verify row order matches detection_history ---
identical(fish_flow_07$animal_id,   detection_history$animal_id)
identical(as.numeric(fish_flow_07$water_year),
          as.numeric(detection_history$water_year))

# --- Extract vectors for NIMBLE ---
# Impute 0 (mean) for fish missing dates
flow_occ2 <- fish_flow_07$flow_occ2_std
flow_occ3 <- fish_flow_07$flow_occ3_std
flow_occ2[is.na(flow_occ2)] <- 0
flow_occ3[is.na(flow_occ3)] <- 0

cat("\nfinal flow_occ2 range:", round(range(flow_occ2), 2), "\n")
cat("final flow_occ3 range:", round(range(flow_occ3), 2), "\n")

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
# KEY CHANGES FROM MODEL 06:
# 1. psi_geo uses flow_occ2 (Rio Vista at occ2 arrival)
# 2. psi_ss  uses flow_occ3 (GES at occ3 arrival)
# 3. phi_fail is now individual-level function of flow_occ2:
#    logit(phi_fail[i]) = alpha_fail + beta_fail * flow_occ2[i]
#    Expected: beta_fail < 0 (higher flow = lower failure)
#==============================================================================

nimCode_07 <- nimbleCode({
  
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
  
  # Routing intercepts and slopes
  alpha_geo ~ dnorm(0, sd = 1.5)
  alpha_ss  ~ dnorm(0, sd = 1.5)
  beta_geo  ~ dnorm(0, sd = 1.0)
  beta_ss   ~ dnorm(0, sd = 1.0)
  
  # Failure intercept and slope
  # alpha_fail intercept: log-odds of failure at mean flow
  # beta_fail slope: expected negative (higher flow = lower failure)
  alpha_fail ~ dnorm(0, sd = 1.5)
  beta_fail  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL PARAMETERS ---
  for(i in 1:nfish){
    
    # Routing at occ2 — Rio Vista flow
    psi_geo[i] <- ilogit(alpha_geo + beta_geo * flow_occ2[i])
    
    # Routing at occ3 — GES flow
    psi_ss[i]  <- ilogit(alpha_ss  + beta_ss  * flow_occ3[i])
    
    # Failure probability — Rio Vista flow at occ2
    phi_fail[i] <- ilogit(alpha_fail + beta_fail * flow_occ2[i])
    
  }
  
  #--- OBSERVATION ARRAY ---
  # Identical to model 06
  
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
    
    # Transition 1->2: Georgiana routing — Rio Vista flow
    tr_arr_i[i, 1, 1:6, 1] <- c(S_sac*(1-psi_geo[i])*(1-phi_fail[i]),
                                S_sac*psi_geo[i]*(1-phi_fail[i]),
                                0, (1-S_sac), S_sac*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 1] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
    
    # Transition 2->3: SS routing — GES flow
    tr_arr_i[i, 1, 1:6, 2] <- c(S_sac*(1-psi_ss[i])*(1-phi_fail[i]),
                                0,
                                S_sac*psi_ss[i]*(1-phi_fail[i]),
                                (1-S_sac), S_sac*phi_fail[i], 0)
    tr_arr_i[i, 2, 1:6, 2] <- c(0, 1, 0, 0, 0, 0)
    tr_arr_i[i, 3, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 2] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
    
    # Transitions 3->7: phi_fail[i] now individual-level
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
      probObs      = p_arr[1:6, 1:6, 1:7],
      probTrans    = tr_arr_i[i, 1:6, 1:6, 1:6],
      len          = 7,
      checkRowSums = 0
    )
    
  } #end i loop
})

#==============================================================================
# SECTION 6: BUILD MODEL AND CONFIGURE MCMC
#==============================================================================

inits_07 <- list(
  S_sac      = 0.99, S_geo = 0.95, S_ss = 0.97,
  p_sac2     = 0.977, p_sac3 = 0.899,
  p_sac4     = 0.910, p_sac5 = 0.893, p_sac6 = 0.986,
  p_geo      = 0.771, p_ss   = 0.690,
  lambda     = 0.986,
  alpha_geo  = -1.95,
  alpha_ss   = -0.51,
  beta_geo   = 0.0,
  beta_ss    = 0.0,
  alpha_fail = -2.3,    # ilogit(-2.3) ~ 0.09 = baseline phi_fail
  beta_fail  = 0.0
)

nimMod_07 <- nimbleModel(
  code      = nimCode_07,
  inits     = inits_07,
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(
    nfish      = nrow(ch_mat_nimble),
    flow_occ2  = flow_occ2,
    flow_occ3  = flow_occ3
  )
)
nimMod_07$calculate()

inits_fn_07 <- function(){
  list(
    S_sac      = runif(1, 0.95, 1.0),
    S_geo      = runif(1, 0.80, 1.0),
    S_ss       = runif(1, 0.85, 1.0),
    p_sac2     = runif(1, 0.90, 1.0),
    p_sac3     = runif(1, 0.75, 0.95),
    p_sac4     = runif(1, 0.80, 0.95),
    p_sac5     = runif(1, 0.80, 0.95),
    p_sac6     = runif(1, 0.95, 1.0),
    p_geo      = runif(1, 0.50, 0.90),
    p_ss       = runif(1, 0.50, 0.85),
    lambda     = runif(1, 0.95, 1.0),
    alpha_geo  = rnorm(1, mean = -1.95, sd = 0.3),
    alpha_ss   = rnorm(1, mean = -0.51, sd = 0.3),
    beta_geo   = rnorm(1, mean = 0, sd = 0.2),
    beta_ss    = rnorm(1, mean = 0, sd = 0.2),
    alpha_fail = rnorm(1, mean = -2.3,  sd = 0.3),
    beta_fail  = rnorm(1, mean = 0, sd = 0.2)
  )
}

params_07 <- c(
  "S_sac", "S_geo", "S_ss",
  "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
  "p_geo", "p_ss",
  "lambda",
  "alpha_geo", "alpha_ss", "beta_geo", "beta_ss",
  "alpha_fail", "beta_fail"
)

confMCMC_07 <- configureMCMC(nimMod_07, onlySlice = TRUE)
confMCMC_07$addMonitors(params_07)
MCMC_07   <- buildMCMC(confMCMC_07)
CModel_07 <- compileNimble(nimMod_07)
CMCMC_07  <- compileNimble(MCMC_07, project = CModel_07)

#==============================================================================
# SECTION 7: RUN MCMC AND SAVE
#==============================================================================

# Short test run first — uncomment to check before full run
# mcmc_test_07 <- runMCMC(CMCMC_07, niter = 1000, nchains = 1,
#                          nburnin = 100, thin = 1,
#                          inits = list(inits_fn_07()),
#                          samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test_07, round = 3)

mcmc_out_07 <- runMCMC(
  CMCMC_07,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn_07(), inits_fn_07(), inits_fn_07()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_07, round = 3)

save(mcmc_out_07, flow_occ2, flow_occ3, fish_flow_07,
     flow_occ2_mean, flow_occ2_sd, flow_occ3_mean, flow_occ3_sd,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_07_run1.RData")

MCMCtrace(mcmc_out_07,
          params   = c("alpha_geo", "alpha_ss", "beta_geo", "beta_ss",
                       "alpha_fail", "beta_fail", "lambda"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_flow_07_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)


