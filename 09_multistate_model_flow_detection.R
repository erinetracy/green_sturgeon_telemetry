#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - 09
# Script: 09_multistate_model_flow_detection.R
# Author: Erin Tracy
# Last updated: June 2026
#
# PURPOSE:
# Extension of model 08 testing whether discharge influences detection
# probability at Steamboat/Sutter Slough arrays (p_ss), following
# Perry et al. (2018). beta_p_geo dropped after model 09 test showed
# CI including zero. beta_p_sac3 tested in model 09b and not supported.
#
# FINAL MODEL STRUCTURE:
# Routing occ2: psi_geo[i] = ilogit(alpha_geo + beta_geo * flow_occ2[i])
# Routing occ3: psi_ss[i]  = ilogit(alpha_ss  + beta_ss  * flow_occ3[i])
# Failure:      phi_fail[i] = ilogit(alpha_fail + beta_fail * flow_30day[i])
# SS detection: p_ss_i[i]   = ilogit(alpha_p_ss + beta_p_ss * flow_occ3[i])
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

# SAFETY CHECK
stopifnot(all(ch_mat_nimble %in% c(1, 2, 4, 5, 6, 7)))
cat("ch_mat_nimble pre-recoding values:", 
    sort(unique(as.vector(ch_mat_nimble))), "\n")

#==============================================================================
# SECTION 2: RECODE DETECTION HISTORY TO 5-STATE MODEL
#==============================================================================

nstate <- 5

ch_mat_nimble_5state <- ch_mat_nimble
ch_mat_nimble_5state[ch_mat_nimble == 4] <- 3   # SS: 4->3
ch_mat_nimble_5state[ch_mat_nimble == 5] <- 4   # dead: 5->4
ch_mat_nimble_5state[ch_mat_nimble == 6] <- 5   # failed: 6->5
ch_mat_nimble_5state[ch_mat_nimble == 7] <- 6   # not-detected: 7->6

ch_mat_nimble <- ch_mat_nimble_5state

cat("Unique values after recoding:\n")
print(table(as.vector(ch_mat_nimble)))
cat("Expected: 1 2 3 4 5 6\n")
cat("Total fish:", nrow(ch_mat_nimble), "\n")
print(table(detection_history$status))

#==============================================================================
# SECTION 3: BUILD FLOW COVARIATES
#==============================================================================

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

occ2_dates <- model1_events %>%
  filter(occasion == 2) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(arr_date_occ2 = as.Date(min(first_detection)), .groups = "drop")

occ3_dates <- model1_events %>%
  filter(occasion == 3) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(arr_date_occ3 = as.Date(min(first_detection)), .groups = "drop")

benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(benicia_date = as.Date(min(first_detection)), .groups = "drop")

fish_flow_09 <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(occ2_dates,    by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,    by = c("animal_id", "water_year")) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    flow_occ2_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(rv_flow_complete$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else rv_flow_complete$flow_cfs[idx]
      }
    },
    flow_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(ges_flow_complete$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else ges_flow_complete$flow_GES[idx]
      }
    },
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

cat("Missing flow_occ2:", sum(is.na(fish_flow_09$flow_occ2_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_flow_09$flow_occ3_raw)), "\n")
cat("Missing flow_30day:", sum(is.na(fish_flow_09$flow_30day_raw)), "\n")

flow_occ2_mean_09  <- mean(fish_flow_09$flow_occ2_raw,  na.rm = TRUE)
flow_occ2_sd_09    <- sd(fish_flow_09$flow_occ2_raw,    na.rm = TRUE)
flow_occ3_mean_09  <- mean(fish_flow_09$flow_occ3_raw,  na.rm = TRUE)
flow_occ3_sd_09    <- sd(fish_flow_09$flow_occ3_raw,    na.rm = TRUE)
flow_30day_mean_09 <- mean(fish_flow_09$flow_30day_raw, na.rm = TRUE)
flow_30day_sd_09   <- sd(fish_flow_09$flow_30day_raw,   na.rm = TRUE)

fish_flow_09 <- fish_flow_09 %>%
  mutate(
    flow_occ2_std  = (flow_occ2_raw  - flow_occ2_mean_09)  / flow_occ2_sd_09,
    flow_occ3_std  = (flow_occ3_raw  - flow_occ3_mean_09)  / flow_occ3_sd_09,
    flow_30day_std = (flow_30day_raw - flow_30day_mean_09) / flow_30day_sd_09
  )

identical(fish_flow_09$animal_id, detection_history$animal_id)
identical(as.numeric(fish_flow_09$water_year),
          as.numeric(detection_history$water_year))

flow_occ2_09  <- fish_flow_09$flow_occ2_std
flow_occ3_09  <- fish_flow_09$flow_occ3_std
flow_30day_09 <- fish_flow_09$flow_30day_std
flow_occ2_09[is.na(flow_occ2_09)]   <- 0
flow_occ3_09[is.na(flow_occ3_09)]   <- 0
flow_30day_09[is.na(flow_30day_09)] <- 0

cat("flow_occ2 range:",  round(range(flow_occ2_09),  2), "\n")
cat("flow_occ3 range:",  round(range(flow_occ3_09),  2), "\n")
cat("flow_30day range:", round(range(flow_30day_09), 2), "\n")

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
# SECTION 5: NIMBLE MODEL 09
#==============================================================================

nimCode_09 <- nimbleCode({
  
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
  
  alpha_geo  ~ dnorm(0, sd = 1.5)
  alpha_ss   ~ dnorm(0, sd = 1.5)
  beta_geo   ~ dnorm(0, sd = 1.0)
  beta_ss    ~ dnorm(0, sd = 1.0)
  
  alpha_fail ~ dnorm(0, sd = 1.5)
  beta_fail  ~ dnorm(0, sd = 1.0)
  
  alpha_p_ss ~ dnorm(0, sd = 1.5)
  beta_p_ss  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL PARAMETERS ---
  for(i in 1:nfish){
    psi_geo[i]  <- ilogit(alpha_geo  + beta_geo  * flow_occ2[i])
    psi_ss[i]   <- ilogit(alpha_ss   + beta_ss   * flow_occ3[i])
    phi_fail[i] <- ilogit(alpha_fail + beta_fail * flow_30day[i])
    p_ss_i[i]   <- ilogit(alpha_p_ss + beta_p_ss * flow_occ3[i])
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
  
  # Occasion 3 placeholder — SS row overridden per fish
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

inits_09 <- list(
  S_sac      = 0.99, S_geo = 0.95, S_ss = 0.97,
  p_sac2     = 0.977, p_sac3 = 0.899,
  p_sac4     = 0.910, p_sac5 = 0.893, p_sac6 = 0.986,
  p_geo      = 0.771,
  lambda     = 0.986,
  alpha_geo  = -1.95, beta_geo  = 0.0,
  alpha_ss   = -0.51, beta_ss   = 0.0,
  alpha_fail = -2.3,  beta_fail = 0.0,
  alpha_p_ss = -0.39, beta_p_ss = 0.0
)

nimMod_09 <- nimbleModel(
  code      = nimCode_09,
  inits     = inits_09,
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(
    nfish      = nrow(ch_mat_nimble),
    flow_occ2  = flow_occ2_09,
    flow_occ3  = flow_occ3_09,
    flow_30day = flow_30day_09
  )
)

nimMod_09$calculate()

inits_fn_09 <- function(){
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
    lambda     = runif(1, 0.95, 1.0),
    alpha_geo  = rnorm(1, -1.95, 0.3), beta_geo  = rnorm(1, 0, 0.2),
    alpha_ss   = rnorm(1, -0.51, 0.3), beta_ss   = rnorm(1, 0, 0.2),
    alpha_fail = rnorm(1, -2.3,  0.3), beta_fail = rnorm(1, 0, 0.2),
    alpha_p_ss = rnorm(1, -0.39, 0.3), beta_p_ss = rnorm(1, 0, 0.2)
  )
}

params_09 <- c(
  "S_sac", "S_geo", "S_ss",
  "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
  "p_geo", "lambda",
  "alpha_geo", "beta_geo",
  "alpha_ss",  "beta_ss",
  "alpha_fail", "beta_fail",
  "alpha_p_ss", "beta_p_ss"
)

confMCMC_09 <- configureMCMC(nimMod_09, onlySlice = TRUE)
confMCMC_09$addMonitors(params_09)
MCMC_09   <- buildMCMC(confMCMC_09)
CModel_09 <- compileNimble(nimMod_09)
CMCMC_09  <- compileNimble(MCMC_09, project = CModel_09)

#==============================================================================
# SECTION 7: RUN MCMC AND SAVE
#==============================================================================

# Short test run first
# mcmc_test_09 <- runMCMC(CMCMC_09, niter = 1000, nchains = 1,
#                          nburnin = 100, thin = 1,
#                          inits = list(inits_fn_09()),
#                          samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test_09, round = 3)

mcmc_out_09 <- runMCMC(
  CMCMC_09,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn_09(), inits_fn_09(), inits_fn_09()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_09, round = 3)

save(mcmc_out_09, flow_occ2_09, flow_occ3_09, flow_30day_09,
     fish_flow_09,
     flow_occ2_mean_09, flow_occ2_sd_09,
     flow_occ3_mean_09, flow_occ3_sd_09,
     flow_30day_mean_09, flow_30day_sd_09,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_09_run1.RData")

MCMCtrace(mcmc_out_09,
          params   = c("alpha_geo", "beta_geo",
                       "alpha_ss",  "beta_ss",
                       "alpha_fail", "beta_fail",
                       "alpha_p_ss", "beta_p_ss",
                       "lambda"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_flow_09_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)

cat("Model 09 saved\n")

#==============================================================================
# SECTION 8: POSTERIOR PREDICTIVE CHECKS
# USING NIMBLE'S NATIVE SIMULATE()
# Guarantees simulation exactly matches the fitted model structure
# =============================================================================

library(nimble)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Make sure compiled model objects from model 09 are still in environment ---
# nimMod_09, CModel_09 should exist from the model 09 run
# If not, you'll need to rebuild and recompile before running this

cat("Checking required objects exist...\n")
cat("nimMod_09 exists:", exists("nimMod_09"), "\n")
cat("CModel_09 exists:", exists("CModel_09"), "\n")

# --- Extract posterior samples ---
post_samples <- do.call(rbind, lapply(mcmc_out_09, as.matrix))
n_samples    <- nrow(post_samples)
nfish        <- nrow(ch_mat_nimble)
nocc         <- 7

cat("Posterior samples:", n_samples, "\n")

# Use a subset for computational speed
sample_idx <- seq(1, n_samples, by = 20)
n_ppc      <- length(sample_idx)
cat("Using", n_ppc, "posterior draws for PPC\n")

# --- Storage for simulated detection rates ---
rep_stats <- matrix(NA, nrow = n_ppc, ncol = nocc)

# --- Identify the data node names for ch_mat ---
data_nodes <- CModel_09$getNodeNames(dataOnly = TRUE)
cat("Data nodes found:", length(data_nodes), "\n")
cat("Example node names:", head(data_nodes, 3), "\n")

# --- Run simulation loop using the compiled model ---
cat("\nSimulating", n_ppc, "replicated datasets using NIMBLE simulate()...\n")

for(k in 1:n_ppc){
  s <- sample_idx[k]
  
  # set parameter values in the compiled model
  CModel_09$alpha_geo  <- post_samples[s, "alpha_geo"]
  CModel_09$beta_geo   <- post_samples[s, "beta_geo"]
  CModel_09$alpha_ss   <- post_samples[s, "alpha_ss"]
  CModel_09$beta_ss    <- post_samples[s, "beta_ss"]
  CModel_09$alpha_fail <- post_samples[s, "alpha_fail"]
  CModel_09$beta_fail  <- post_samples[s, "beta_fail"]
  CModel_09$alpha_p_ss <- post_samples[s, "alpha_p_ss"]
  CModel_09$beta_p_ss  <- post_samples[s, "beta_p_ss"]
  CModel_09$S_sac      <- post_samples[s, "S_sac"]
  CModel_09$S_geo      <- post_samples[s, "S_geo"]
  CModel_09$S_ss       <- post_samples[s, "S_ss"]
  CModel_09$p_sac2     <- post_samples[s, "p_sac2"]
  CModel_09$p_sac3     <- post_samples[s, "p_sac3"]
  CModel_09$p_sac4     <- post_samples[s, "p_sac4"]
  CModel_09$p_sac5     <- post_samples[s, "p_sac5"]
  CModel_09$p_sac6     <- post_samples[s, "p_sac6"]
  CModel_09$p_geo      <- post_samples[s, "p_geo"]
  CModel_09$lambda     <- post_samples[s, "lambda"]
  
  # recalculate deterministic nodes (psi_geo, psi_ss, phi_fail, p_ss_i,
  # p_arr_i, tr_arr_i) given the new parameter values
  CModel_09$calculate()
  
  # simulate the data node (ch_mat) from the model's own data-generating
  # process — this uses the actual fitted transition/observation logic
  CModel_09$simulate(data_nodes, includeData = TRUE)
  
  # extract simulated detection histories
  sim_ch <- CModel_09$ch_mat
  
  # compute detection rate by occasion (not state 6 = not detected)
  rep_stats[k, ] <- colMeans(sim_ch != 6)
  
  if(k %% 50 == 0) cat("Completed", k, "of", n_ppc, "simulations\n")
}

cat("Simulation complete\n")

# --- Compute observed test statistics ---
obs_stats <- colMeans(ch_mat_nimble != 6)
cat("\nObserved detection rates by occasion:\n")
print(round(obs_stats, 3))

# --- Bayesian p-values ---
bayesian_pvals <- colMeans(sweep(rep_stats, 2, obs_stats, ">="))

ppc_results <- data.frame(
  occasion     = paste0("occ", 1:nocc),
  observed     = round(obs_stats, 3),
  sim_mean     = round(colMeans(rep_stats), 3),
  sim_lower    = round(apply(rep_stats, 2, quantile, 0.025), 3),
  sim_upper    = round(apply(rep_stats, 2, quantile, 0.975), 3),
  bayesian_p   = round(bayesian_pvals, 3),
  adequate_fit = bayesian_pvals > 0.05 & bayesian_pvals < 0.95
)

cat("\n=== POSTERIOR PREDICTIVE CHECK RESULTS (NIMBLE native simulate) ===\n")
print(ppc_results)

# --- Plot ---
ppc_plot_data <- as.data.frame(rep_stats)
colnames(ppc_plot_data) <- paste0("occ", 1:nocc)
ppc_long <- tidyr::pivot_longer(ppc_plot_data,
                                cols = everything(),
                                names_to = "occasion",
                                values_to = "sim_rate")

obs_df <- data.frame(
  occasion = paste0("occ", 1:nocc),
  obs_rate = obs_stats
)

p_ppc <- ggplot(ppc_long, aes(x = sim_rate)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(data = obs_df,
             aes(xintercept = obs_rate),
             color = "red", linewidth = 1.2) +
  facet_wrap(~ occasion, ncol = 4, scales = "free") +
  labs(
    x        = "Simulated detection rate",
    y        = "Count",
    title    = "Posterior predictive check — model 09",
    subtitle = "Red line = observed detection rate | Simulated using NIMBLE's native data-generating process"
  ) +
  theme_bw(base_size = 11, base_family = "Times New Roman") +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 9, color = "gray40")
  )

print(p_ppc)

ggsave(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/figures/gs_ppc_model09_native.pdf",
  plot = p_ppc, width = 10, height = 6, device = "pdf")

ggsave(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/figures/gs_ppc_model09_native.png",
  plot = p_ppc, width = 10, height = 6, dpi = 300)

write.csv(ppc_results,
          "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/ppc_results_model09_native.csv",
          row.names = FALSE)

cat("\nPosterior predictive check complete and saved\n")
