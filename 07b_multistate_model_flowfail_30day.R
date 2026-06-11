#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - 07b
# Script: 07b_multistate_model_flowfail_30day.R
# Author: Erin Tracy
# Last updated: June 2026
#
# PURPOSE:
# Sensitivity analysis to model 07 using 30-day antecedent Rio Vista flow
# before Benicia passage as the single flow covariate for both routing
# decisions and flow-dependent failure. Single flow covariate justified
# by high collinearity (r=0.91) between 30-day Rio Vista and 30-day GES
# windows — fish integrating conditions over 30 days before entering
# freshwater are responding to a single systemic hydrological signal.
#
# DIFFERENCES FROM MODEL 07:
# 1. Flow covariate: 30-day mean before Benicia (vs same-day junction flow)
# 2. Single flow vector for both psi_geo, psi_ss, and phi_fail
# 3. Gauge: Rio Vista only (GES dropped — collinear over 30-day window)
#
# BIOLOGICAL RATIONALE:
# Fish staging in the estuary integrate flow conditions over weeks before
# committing to migrate. The 30-day mean before Benicia represents the
# hydrological context during staging and may predict routing behavior
# if fish are responding to broader flow regime rather than local
# instantaneous hydraulics at the junction.
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
# SECTION 3: BUILD 30-DAY ANTECEDENT FLOW COVARIATE
# Single covariate: 30-day mean Rio Vista flow before Benicia passage
# Used for psi_geo, psi_ss, and phi_fail
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
cat("NAs after interpolation:", sum(is.na(rv_flow_complete$flow_cfs)), "\n")

# --- Get Benicia passage dates ---
benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    benicia_date = as.Date(min(first_detection)),
    .groups = "drop"
  )

cat("Fish with Benicia dates:", nrow(benicia_dates), "\n")

# --- Calculate 30-day mean Rio Vista flow before Benicia ---
fish_flow_07b <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    flow_30day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- rv_flow_complete$flow_cfs[
          rv_flow_complete$date >= (benicia_date - 30) &
            rv_flow_complete$date <   benicia_date]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    n_days_in_window = {
      if(is.na(benicia_date)) 0L else
        sum(rv_flow_complete$date >= (benicia_date - 30) &
              rv_flow_complete$date <   benicia_date &
              !is.na(rv_flow_complete$flow_cfs))
    }
  ) %>%
  ungroup()

cat("Missing flow_30day:", sum(is.na(fish_flow_07b$flow_30day_raw)), "\n")
cat("Window size summary:\n")
print(table(fish_flow_07b$n_days_in_window))

# --- Standardize ---
flow_mean_07b <- mean(fish_flow_07b$flow_30day_raw, na.rm = TRUE)
flow_sd_07b   <- sd(fish_flow_07b$flow_30day_raw,   na.rm = TRUE)

fish_flow_07b <- fish_flow_07b %>%
  mutate(
    flow_30day_std = (flow_30day_raw - flow_mean_07b) / flow_sd_07b
  )

cat("\n30-day antecedent flow summary:\n")
cat("Mean raw flow:", round(flow_mean_07b, 0), "cfs\n")
cat("SD raw flow:",   round(flow_sd_07b, 0),   "cfs\n")
print(summary(fish_flow_07b$flow_30day_std))

# --- Verify row order matches detection_history ---
identical(fish_flow_07b$animal_id, detection_history$animal_id)
identical(as.numeric(fish_flow_07b$water_year),
          as.numeric(detection_history$water_year))

# --- Extract single flow vector for NIMBLE ---
# Same vector used for psi_geo, psi_ss, and phi_fail
flow_30day <- fish_flow_07b$flow_30day_std
flow_30day[is.na(flow_30day)] <- 0

cat("\nFinal flow vector:\n")
cat("n fish:", length(flow_30day), "\n")
cat("NAs imputed to 0:", sum(is.na(fish_flow_07b$flow_30day_std)), "\n")
cat("Range:", round(range(flow_30day), 2), "\n")

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
# Single 30-day Rio Vista flow covariate for routing and failure
#==============================================================================

nimCode_07b <- nimbleCode({
  
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
  
  # Routing intercepts and slopes — 30-day antecedent flow
  alpha_geo  ~ dnorm(0, sd = 1.5)
  alpha_ss   ~ dnorm(0, sd = 1.5)
  beta_geo   ~ dnorm(0, sd = 1.0)
  beta_ss    ~ dnorm(0, sd = 1.0)
  
  # Failure intercept and slope — same 30-day flow
  alpha_fail ~ dnorm(0, sd = 1.5)
  beta_fail  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL PARAMETERS ---
  # Single flow covariate for all three individual-level parameters
  for(i in 1:nfish){
    psi_geo[i]  <- ilogit(alpha_geo  + beta_geo  * flow_30day[i])
    psi_ss[i]   <- ilogit(alpha_ss   + beta_ss   * flow_30day[i])
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

inits_07b <- list(
  S_sac      = 0.99, S_geo = 0.95, S_ss = 0.97,
  p_sac2     = 0.977, p_sac3 = 0.899,
  p_sac4     = 0.910, p_sac5 = 0.893, p_sac6 = 0.986,
  p_geo      = 0.771, p_ss   = 0.690,
  lambda     = 0.986,
  alpha_geo  = -1.95,
  alpha_ss   = -0.51,
  beta_geo   = 0.0,
  beta_ss    = 0.0,
  alpha_fail = -2.3,
  beta_fail  = 0.0
)

nimMod_07b <- nimbleModel(
  code      = nimCode_07b,
  inits     = inits_07b,
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(
    nfish      = nrow(ch_mat_nimble),
    flow_30day = flow_30day
  )
)
nimMod_07b$calculate()

inits_fn_07b <- function(){
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

params_07b <- c(
  "S_sac", "S_geo", "S_ss",
  "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
  "p_geo", "p_ss",
  "lambda",
  "alpha_geo", "alpha_ss", "beta_geo", "beta_ss",
  "alpha_fail", "beta_fail"
)

confMCMC_07b <- configureMCMC(nimMod_07b, onlySlice = TRUE)
confMCMC_07b$addMonitors(params_07b)
MCMC_07b   <- buildMCMC(confMCMC_07b)
CModel_07b <- compileNimble(nimMod_07b)
CMCMC_07b  <- compileNimble(MCMC_07b, project = CModel_07b)

#==============================================================================
# SECTION 7: RUN MCMC AND SAVE
#==============================================================================

# Short test run first — uncomment to verify before full run
# mcmc_test_07b <- runMCMC(CMCMC_07b, niter = 1000, nchains = 1,
#                           nburnin = 100, thin = 1,
#                           inits = list(inits_fn_07b()),
#                           samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test_07b, round = 3)

mcmc_out_07b <- runMCMC(
  CMCMC_07b,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn_07b(), inits_fn_07b(), inits_fn_07b()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_07b, round = 3)

save(mcmc_out_07b, flow_30day, fish_flow_07b,
     flow_mean_07b, flow_sd_07b,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_07b_run1.RData")

MCMCtrace(mcmc_out_07b,
          params   = c("alpha_geo", "alpha_ss", "beta_geo", "beta_ss",
                       "alpha_fail", "beta_fail", "lambda"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_flow_07b_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)
