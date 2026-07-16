#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - 12b
# Script: 12b_multistate_model_flow_temp_final.R
# Author: Erin Tracy
# Last updated: June 2026
#
# PURPOSE:
# Final combined model dropping 30-day flow from failure probability
# following model 12 results showing beta_fail_f collapsed to zero
# when 7-day temperature was included (beta_fail_f = 0.082,
# CI: -0.453 to 0.515). Temperature subsumes flow as the driver
# of migration failure.
#
# FINAL MODEL STRUCTURE:
# psi_geo[i]  = ilogit(alpha_geo  + beta_geo   * flow_occ2[i])
# psi_ss[i]   = ilogit(alpha_ss   + beta_ss    * flow_occ3[i])
# phi_fail[i] = ilogit(alpha_fail + beta_fail_t * temp_7day[i])
# p_ss_i[i]   = ilogit(alpha_p_ss + beta_p_ss  * flow_occ3[i])
#
# BIOLOGICAL INTERPRETATION:
# Flow drives junction routing decisions (instantaneous hydraulic cues)
# Temperature drives migration failure (thermal conditions during staging)
# Flow drives SS detection (acoustic detection efficiency)
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

# SAFETY CHECK, re-load Rdata file before each run
stopifnot(all(ch_mat_nimble %in% c(1, 2, 4, 5, 6, 7)))
cat("ch_mat_nimble pre-recoding values:",
    sort(unique(as.vector(ch_mat_nimble))), "\n")

#==============================================================================
# SECTION 2: RECODE DETECTION HISTORY TO 5-STATE MODEL
#originally contained delta cross but no fish used as final route upstream during 2007-2017
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
# SECTION 3: BUILD COVARIATES
# flow_occ2:  same-day Rio Vista (routing occ2)
# flow_occ3:  same-day GES (routing occ3 + SS detection)
# temp_7day:  7-day mean Rio Vista temp before Benicia (failure)
# NOTE: 30-day flow dropped from failure — not needed when temp included
#==============================================================================
# load pre-built covariates — run gs_build_model12b_covariates.R first
fish_cov_12b <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_model12b_covariates.csv"
)

std_params <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_model12b_std_params.csv"
)

# restore standardization parameters
flow_occ2_mean_12b <- std_params$mean[std_params$covariate == "flow_occ2"]
flow_occ2_sd_12b   <- std_params$sd[std_params$covariate   == "flow_occ2"]
flow_occ3_mean_12b <- std_params$mean[std_params$covariate == "flow_occ3"]
flow_occ3_sd_12b   <- std_params$sd[std_params$covariate   == "flow_occ3"]
temp_7day_mean_12b <- std_params$mean[std_params$covariate == "temp_7day"]
temp_7day_sd_12b   <- std_params$sd[std_params$covariate   == "temp_7day"]

# extract vectors and impute 0 for missing
flow_occ2_12b <- fish_cov_12b$flow_occ2_std
flow_occ3_12b <- fish_cov_12b$flow_occ3_std
temp_7day_12b <- fish_cov_12b$temp_7day_std

flow_occ2_12b[is.na(flow_occ2_12b)] <- 0
flow_occ3_12b[is.na(flow_occ3_12b)] <- 0
temp_7day_12b[is.na(temp_7day_12b)] <- 0

cat("Covariates loaded\n")
cat("n fish:", nrow(fish_cov_12b), "\n")


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
# SECTION 5: NIMBLE MODEL 12b
#==============================================================================

nimCode_12b <- nimbleCode({
  
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
  
  # Routing — same-day gauge-specific flow
  alpha_geo ~ dnorm(0, sd = 1.5)
  alpha_ss  ~ dnorm(0, sd = 1.5)
  beta_geo  ~ dnorm(0, sd = 1.0)
  beta_ss   ~ dnorm(0, sd = 1.0)
  
  # Failure — 7-day temperature only
  alpha_fail  ~ dnorm(0, sd = 1.5)
  beta_fail_t ~ dnorm(0, sd = 1.0)
  
  # SS detection — same-day GES flow
  alpha_p_ss ~ dnorm(0, sd = 1.5)
  beta_p_ss  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL PARAMETERS ---
  for(i in 1:nfish){
    psi_geo[i]  <- ilogit(alpha_geo  + beta_geo   * flow_occ2[i])
    psi_ss[i]   <- ilogit(alpha_ss   + beta_ss    * flow_occ3[i])
    phi_fail[i] <- ilogit(alpha_fail + beta_fail_t * temp_7day[i])
    p_ss_i[i]   <- ilogit(alpha_p_ss + beta_p_ss  * flow_occ3[i])
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

inits_12b <- list(
  S_sac       = 0.99, S_geo = 0.95, S_ss = 0.97,
  p_sac2      = 0.977, p_sac3 = 0.899,
  p_sac4      = 0.910, p_sac5 = 0.893, p_sac6 = 0.986,
  p_geo       = 0.771, lambda = 0.986,
  alpha_geo   = -1.95, beta_geo    = 0.0,
  alpha_ss    = -0.51, beta_ss     = 0.0,
  alpha_fail  = -2.30, beta_fail_t = 0.0,
  alpha_p_ss  = -0.39, beta_p_ss   = 0.0
)

nimMod_12b <- nimbleModel(
  code      = nimCode_12b,
  inits     = inits_12b,
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(
    nfish      = nrow(ch_mat_nimble),
    flow_occ2  = flow_occ2_12b,
    flow_occ3  = flow_occ3_12b,
    temp_7day  = temp_7day_12b
  )
)

nimMod_12b$calculate()

inits_fn_12b <- function(){
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
    alpha_geo   = rnorm(1, -1.95, 0.3), beta_geo    = rnorm(1, 0, 0.2),
    alpha_ss    = rnorm(1, -0.51, 0.3), beta_ss     = rnorm(1, 0, 0.2),
    alpha_fail  = rnorm(1, -2.30, 0.3), beta_fail_t = rnorm(1, 0, 0.2),
    alpha_p_ss  = rnorm(1, -0.39, 0.3), beta_p_ss   = rnorm(1, 0, 0.2)
  )
}

params_12b <- c(
  "S_sac", "S_geo", "S_ss",
  "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
  "p_geo", "lambda",
  "alpha_geo", "beta_geo",
  "alpha_ss",  "beta_ss",
  "alpha_fail", "beta_fail_t",
  "alpha_p_ss", "beta_p_ss"
)

confMCMC_12b <- configureMCMC(nimMod_12b, onlySlice = TRUE)
confMCMC_12b$addMonitors(params_12b)
MCMC_12b   <- buildMCMC(confMCMC_12b)
CModel_12b <- compileNimble(nimMod_12b)
CMCMC_12b  <- compileNimble(MCMC_12b, project = CModel_12b)

#==============================================================================
# SECTION 7: RUN MCMC AND SAVE
#==============================================================================

# Short test run first
# mcmc_test_12b <- runMCMC(CMCMC_12b, niter = 1000, nchains = 1,
#                           nburnin = 100, thin = 1,
#                           inits = list(inits_fn_12b()),
#                           samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test_12b, round = 3)

mcmc_out_12b <- runMCMC(
  CMCMC_12b,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn_12b(), inits_fn_12b(), inits_fn_12b()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_12b, round = 3)

save(mcmc_out_12b,
     flow_occ2_12b, flow_occ3_12b, temp_7day_12b,
     fish_cov_12b,
     flow_occ2_mean_12b, flow_occ2_sd_12b,
     flow_occ3_mean_12b, flow_occ3_sd_12b,
     temp_7day_mean_12b, temp_7day_sd_12b,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_temp_12b_run1.RData")

MCMCtrace(mcmc_out_12b,
          params   = c("alpha_geo", "beta_geo",
                       "alpha_ss",  "beta_ss",
                       "alpha_fail", "beta_fail_t",
                       "alpha_p_ss", "beta_p_ss"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_flow_temp_12b_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)

cat("Model 12b complete\n")


#==============================================================================
# POSTERIOR PREDICTIVE CHECKS — MODEL 12b (FINAL MODEL)
# Using NIMBLE's native simulate() for exact model replication
#==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# confirm compiled model objects exist from model 12b run
cat("CModel_12b exists:", exists("CModel_12b"), "\n")
cat("mcmc_out_12b exists:", exists("mcmc_out_12b"), "\n")

# if starting fresh session reload everything
# load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_multistate_data.RData")
# load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_temp_12b_run1.RData")
# then rerun sections 1-6 of model 12b script to rebuild CModel_12b

# extract posterior samples
post_samples <- do.call(rbind, lapply(mcmc_out_12b, as.matrix))
n_samples    <- nrow(post_samples)
nfish        <- nrow(ch_mat_nimble)
nocc         <- 7

cat("Posterior samples:", n_samples, "\n")

# use every 20th sample for speed
sample_idx <- seq(1, n_samples, by = 20)
n_ppc      <- length(sample_idx)
cat("Using", n_ppc, "posterior draws\n")

# storage
rep_stats <- matrix(NA, nrow = n_ppc, ncol = nocc)

# identify data nodes
data_nodes <- CModel_12b$getNodeNames(dataOnly = TRUE)
cat("Data nodes found:", length(data_nodes), "\n")
cat("Example:", head(data_nodes, 3), "\n")

# run simulation loop
cat("\nSimulating", n_ppc, "replicated datasets...\n")

for(k in 1:n_ppc){
  s <- sample_idx[k]
  
  # set parameter values in compiled model
  CModel_12b$alpha_geo   <- post_samples[s, "alpha_geo"]
  CModel_12b$beta_geo    <- post_samples[s, "beta_geo"]
  CModel_12b$alpha_ss    <- post_samples[s, "alpha_ss"]
  CModel_12b$beta_ss     <- post_samples[s, "beta_ss"]
  CModel_12b$alpha_fail  <- post_samples[s, "alpha_fail"]
  CModel_12b$beta_fail_t <- post_samples[s, "beta_fail_t"]
  CModel_12b$alpha_p_ss  <- post_samples[s, "alpha_p_ss"]
  CModel_12b$beta_p_ss   <- post_samples[s, "beta_p_ss"]
  CModel_12b$S_sac       <- post_samples[s, "S_sac"]
  CModel_12b$S_geo       <- post_samples[s, "S_geo"]
  CModel_12b$S_ss        <- post_samples[s, "S_ss"]
  CModel_12b$p_sac2      <- post_samples[s, "p_sac2"]
  CModel_12b$p_sac3      <- post_samples[s, "p_sac3"]
  CModel_12b$p_sac4      <- post_samples[s, "p_sac4"]
  CModel_12b$p_sac5      <- post_samples[s, "p_sac5"]
  CModel_12b$p_sac6      <- post_samples[s, "p_sac6"]
  CModel_12b$p_geo       <- post_samples[s, "p_geo"]
  CModel_12b$lambda      <- post_samples[s, "lambda"]
  
  # recalculate deterministic nodes
  CModel_12b$calculate()
  
  # simulate detection histories
  CModel_12b$simulate(data_nodes, includeData = TRUE)
  
  # extract simulated ch_mat and compute detection rates
  sim_ch <- CModel_12b$ch_mat
  rep_stats[k, ] <- colMeans(sim_ch != 6)
  
  if(k %% 50 == 0) cat("Completed", k, "of", n_ppc, "\n")
}

cat("Simulation complete\n")

# observed detection rates
obs_stats <- colMeans(ch_mat_nimble != 6)
cat("\nObserved detection rates by occasion:\n")
print(round(obs_stats, 3))

# Bayesian p-values
bayesian_pvals <- colMeans(sweep(rep_stats, 2, obs_stats, ">="))

ppc_results_12b <- data.frame(
  occasion     = paste0("occ", 1:nocc),
  observed     = round(obs_stats, 3),
  sim_mean     = round(colMeans(rep_stats), 3),
  sim_lower    = round(apply(rep_stats, 2, quantile, 0.025), 3),
  sim_upper    = round(apply(rep_stats, 2, quantile, 0.975), 3),
  bayesian_p   = round(bayesian_pvals, 3),
  adequate_fit = bayesian_pvals > 0.05 & bayesian_pvals < 0.95
)

cat("\n=== POSTERIOR PREDICTIVE CHECK RESULTS — MODEL 12b ===\n")
print(ppc_results_12b)

# compare to model 09 PPC results
cat("\n=== COMPARISON TO MODEL 09 PPC ===\n")
cat("Model 09 bayesian p-values: 1.000 0.325 0.043 0.327 0.462 0.558 1.000\n")
cat("Model 12b bayesian p-values:", round(bayesian_pvals, 3), "\n")

# plot
ppc_plot_data <- as.data.frame(rep_stats)
colnames(ppc_plot_data) <- paste0("occ", 1:nocc)

ppc_long <- tidyr::pivot_longer(ppc_plot_data,
                                cols = everything(),
                                names_to  = "occasion",
                                values_to = "sim_rate")

obs_df <- data.frame(
  occasion = paste0("occ", 1:nocc),
  obs_rate = obs_stats
)

p_ppc_12b <- ggplot(ppc_long, aes(x = sim_rate)) +
  geom_histogram(bins = 30, fill = "steelblue",
                 alpha = 0.7, color = "white") +
  geom_vline(data = obs_df,
             aes(xintercept = obs_rate),
             color = "red", linewidth = 1.2) +
  facet_wrap(~ occasion, ncol = 4, scales = "free") +
  labs(
    x        = "Simulated detection rate",
    y        = "Count",
    title    = "Posterior predictive check — final model (12b)",
    subtitle = "Red line = observed | Simulated using NIMBLE native data-generating process"
  ) +
  theme_bw(base_size = 11, base_family = "Times New Roman") +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 9, color = "gray40")
  )

print(p_ppc_12b)

ggsave(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/figures/gs_ppc_model12b_final.pdf",
  plot = p_ppc_12b, width = 10, height = 6, device = "pdf")

ggsave(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/figures/gs_ppc_model12b_final.png",
  plot = p_ppc_12b, width = 10, height = 6, dpi = 300)

write.csv(ppc_results_12b,
          "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/ppc_results_model12b_final.csv",
          row.names = FALSE)

cat("\nPPC complete and saved\n")
