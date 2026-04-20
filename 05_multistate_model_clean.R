#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL
# Script: 05_multistate_model.R
# Author: Erin Tracy
# Last updated: April 2026
#
# PURPOSE:
# Build a discrete hidden Markov multistate model (dDHMMo) to estimate
# route-specific survival, detection, and routing probabilities for adult
# green sturgeon making upstream spawning migrations in the Sacramento River
# system. Follows Perry et al. (2018) framework using nimbleEcology.
#
# INPUTS:
#   - gs_multistate_data.RData (from 04_multistate_data_prep.R):
#       ch_mat_nimble: 221 x 7 detection history matrix formatted for dDHMMo
#       ch_mat: 221 x 7 detection history matrix (0 = not detected)
#       detection_history: fish-level detection history with metadata
#       nstate: number of states (will be updated to 5 here)
#       fish_info: fish-level metadata
#       exploration_summary: fish-level exploration behavior
#
# OUTPUTS:
#   - mcmc_out_full: MCMC posterior samples
#   - gs_mcmc_full_run3_221fish_5state.RData: saved MCMC results
#   - gs_mcmc_traces_run3.pdf: trace plots for key parameters
#
# MODEL STRUCTURE:
#   States: 1=Sacramento, 2=Georgiana, 3=Steamboat/Sutter,
#           4=Dead (absorbing), 5=Failed migration (absorbing)
#   NOTE: DCC dropped — 0 confirmed DCC upstream migrants after removing
#   post-spawn downstream detections. All previous DCC assignments were
#   artifacts of downstream migration being incorrectly included.
#   Occasions: 7 (Benicia -> Rio Vista -> SR_MOUTH -> SR_BLWSTEAM ->
#              SR_KK345R/SR_FREEPORT -> Upper Sac -> Spawning Ground)
#   Years: 2007-2017 (Model 1, no Yolo Bypass)
#   Fish: 221 (134 up_complete, 83 up_incomplete, 4 incomplete_dead)
#==============================================================================

#==============================================================================
# SECTION 1: LOAD LIBRARIES AND DATA
#==============================================================================

library(nimble)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(dplyr)

load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_multistate_data.RData")

#==============================================================================
# SECTION 2: RECODE DETECTION HISTORY TO 5-STATE MODEL
# Old coding: 1=Sac, 2=Geo, 3=DCC(dropped), 4=SS, 5=Dead, 6=Failed, 7=Not-detected
# New coding: 1=Sac, 2=Geo, 3=SS, 4=Dead, 5=Failed, 6=Not-detected
#==============================================================================

nstate <- 5

# Recode states in ch_mat_nimble
ch_mat_nimble_5state <- ch_mat_nimble
ch_mat_nimble_5state[ch_mat_nimble == 4] <- 3   # SS: 4 -> 3
ch_mat_nimble_5state[ch_mat_nimble == 5] <- 4   # dead: 5 -> 4
ch_mat_nimble_5state[ch_mat_nimble == 6] <- 5   # failed: 6 -> 5
ch_mat_nimble_5state[ch_mat_nimble == 7] <- 6   # not-detected: 7 -> 6

# Use 5-state matrix going forward
ch_mat_nimble <- ch_mat_nimble_5state

#==============================================================================
# SECTION 3: VERIFY DETECTION HISTORY
#==============================================================================

cat("Total fish:", nrow(ch_mat_nimble), "\n")  # should be 221
cat("\nStatus breakdown:\n")
print(table(detection_history$status))
cat("\nDetection rates per occasion:\n")
print(round(colMeans(ch_mat > 0), 3))
cat("\nState counts in 5-state ch_mat_nimble:\n")
print(table(ch_mat_nimble))
# State 3 (old DCC) should be absent - confirmed 0 DCC upstream migrants
cat("\nGeorgiana fish:", apply(ch_mat_nimble, 1, function(x) any(x == 2)) %>% sum(), "\n")
cat("Steamboat/Sutter fish:", apply(ch_mat_nimble, 1, function(x) any(x == 3)) %>% sum(), "\n")

#==============================================================================
# SECTION 4: BUILD TRANSITION AND OBSERVATION MATRICES
# Placeholder parameter values used for likelihood verification only.
# All parameters estimated by NIMBLE MCMC in Section 6.
#==============================================================================

# Placeholder parameter values
S_sac1 <- 0.95; S_sac2 <- 0.95; S_sac3 <- 0.95
S_sac4 <- 0.95; S_sac5 <- 0.95
S_geo  <- 0.95; S_ss   <- 0.95
psi_geo <- 0.15; psi_ss <- 0.30
phi_fail <- 0.10; lambda <- 0.95
p_sac2 <- 0.90; p_sac3 <- 0.90; p_sac4 <- 0.85
p_sac5 <- 0.95; p_sac6 <- 0.95
p_geo  <- 0.80; p_ss   <- 0.80

# Template matrix - all zeros, not-detected state (row 6) stays absorbing
temp_mat <- matrix(0, nrow = nstate + 1, ncol = nstate + 1)
temp_mat[nstate + 1, nstate + 1] <- 1

#------------------------------------------------------------------------------
# TRANSITION MATRICES
# States: 1=Sacramento, 2=Georgiana, 3=Steamboat/Sutter,
#         4=Dead (absorbing), 5=Failed (absorbing), 6=Not-detected
# Rows = from state, Columns = to state, Row sums must = 1
#
# KEY STRUCTURAL DECISIONS:
#
# FAILED MIGRATION (state 5) possible at every transition from state 1:
#   phi_fail estimated from 83 up_incomplete fish.
#
# GEO FISH REJOIN SACRAMENTO AT tr_34 (occ3->occ4):
#   Georgiana exits back to Sacramento upstream of SR_BLWSTEAM.
#
# SS FISH PASS THROUGH AT tr_34 (occ3->occ4), REJOIN AT tr_45 (occ4->occ5):
#   SS can fail between occ3->occ4 (phi_fail applies to SS state at tr_34).
#   SR_KK345R is first receiver above both Steamboat and Sutter rejoin points.
#   SS fish confirmed back on Sacramento at occ5.
#
# MISSED DETECTIONS map to column 6 (not detected), NOT column 4 (dead).
#------------------------------------------------------------------------------

# Transition 1->2: Benicia to Rio Vista junction
# Sac fish can route to Geo, fail, or die
tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac1*(1-psi_geo)*(1-phi_fail),
                S_sac1*psi_geo*(1-phi_fail),
                0,
                (1-S_sac1), S_sac1*phi_fail, 0)
tr_12[2, 4] <- 1  # Geo shouldn't exist yet -> dead
tr_12[3, 4] <- 1  # SS shouldn't exist yet -> dead
tr_12[4, 4] <- 1  # dead stays dead
tr_12[5, 5] <- 1  # failed stays failed

# Transition 2->3: Rio Vista to SR_MOUTH/Steamboat junction
# Sac fish can route to SS, fail, or die
# Geo fish pass through (no receivers at this occasion)
tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac2*(1-psi_ss)*(1-phi_fail),
                0,
                S_sac2*psi_ss*(1-phi_fail),
                (1-S_sac2), S_sac2*phi_fail, 0)
tr_23[2, 2] <- 1  # Geo pass through - still in Georgiana channel
tr_23[3, 4] <- 1  # SS shouldn't exist yet -> dead
tr_23[4, 4] <- 1
tr_23[5, 5] <- 1

# Transition 3->4: SR_MOUTH to SR_BLWSTEAM
# Geo rejoins Sacramento, SS passes through (can fail)
tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac3*(1-phi_fail), 0, 0, (1-S_sac3), S_sac3*phi_fail, 0)
tr_34[2, ] <- c(S_geo*(1-phi_fail),  0, 0, (1-S_geo),  S_geo*phi_fail,  0)
tr_34[3, ] <- c(0, 0, (1-phi_fail), 0, phi_fail, 0)  # SS pass through, can fail
tr_34[4, 4] <- 1
tr_34[5, 5] <- 1

# Transition 4->5: SR_BLWSTEAM to SR_KK345R (SS rejoins Sacramento)
tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac4*(1-phi_fail), 0, 0, (1-S_sac4), S_sac4*phi_fail, 0)
tr_45[2, 4] <- 1  # shouldn't exist - Geo rejoined at occ4
tr_45[3, ] <- c(S_ss*(1-phi_fail), 0, 0, (1-S_ss), S_ss*phi_fail, 0)
tr_45[4, 4] <- 1
tr_45[5, 5] <- 1

# Transition 5->6: SR_FREEPORT to upper Sacramento
tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac5*(1-phi_fail), 0, 0, (1-S_sac5), S_sac5*phi_fail, 0)
tr_56[2, 4] <- 1; tr_56[3, 4] <- 1
tr_56[4, 4] <- 1; tr_56[5, 5] <- 1

# Transition 6->7: upper Sacramento to spawning ground
tr_67 <- temp_mat
tr_67[1, ] <- c(lambda*(1-phi_fail), 0, 0, (1-lambda), lambda*phi_fail, 0)
tr_67[2, 4] <- 1; tr_67[3, 4] <- 1
tr_67[4, 4] <- 1; tr_67[5, 5] <- 1

# Combine into transition array [nstate+1 x nstate+1 x n_transitions]
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

cat("Transition matrix row sums:\n")
for(i in 1:6) cat("tr", i, ":", rowSums(tr_arr[,,i]), "\n")

#------------------------------------------------------------------------------
# OBSERVATION MATRICES
# Rows = true state, Columns = observed state
# Missed detections -> column 6 (not detected)
#------------------------------------------------------------------------------

# Occasion 1: Benicia/Carquinez - certain detection
p_mat1 <- temp_mat
p_mat1[1, 1] <- 1
p_mat1[2, nstate+1] <- 1; p_mat1[3, nstate+1] <- 1
p_mat1[4, 4] <- 1; p_mat1[5, 5] <- 1

# Occasion 2: Rio Vista junction
p_mat2 <- temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, 0, (1-p_sac2))
p_mat2[2, ] <- c(0, p_geo, 0, 0, 0, (1-p_geo))
p_mat2[3, nstate+1] <- 1  # SS shouldn't exist yet
p_mat2[4, 4] <- 1; p_mat2[5, 5] <- 1

# Occasion 3: SR_MOUTH/Steamboat junction
# Geo pass through undetected, SS detected at mouth receivers
p_mat3 <- temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, 0, (1-p_sac3))
p_mat3[2, nstate+1] <- 1  # Geo pass through - undetected
p_mat3[3, ] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
p_mat3[4, 4] <- 1; p_mat3[5, 5] <- 1

# Occasion 4: SR_BLWSTEAM area - Geo rejoined, SS still in channel
p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, 0, (1-p_sac4))
p_mat4[2, nstate+1] <- 1  # shouldn't exist - Geo rejoined
p_mat4[3, ] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
p_mat4[4, 4] <- 1; p_mat4[5, 5] <- 1

# Occasion 5: SR_KK345R + SR_FREEPORT - SS has rejoined Sacramento
p_mat5 <- temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, 0, (1-p_sac5))
p_mat5[2, nstate+1] <- 1; p_mat5[3, nstate+1] <- 1
p_mat5[4, 4] <- 1; p_mat5[5, 5] <- 1

# Occasion 6: upper Sacramento
p_mat6 <- temp_mat
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, 0, (1-p_sac6))
p_mat6[2, nstate+1] <- 1; p_mat6[3, nstate+1] <- 1
p_mat6[4, 4] <- 1; p_mat6[5, 5] <- 1

# Occasion 7: spawning ground - certain detection
p_mat7 <- temp_mat
p_mat7[1, 1] <- 1
p_mat7[2, nstate+1] <- 1; p_mat7[3, nstate+1] <- 1
p_mat7[4, 4] <- 1; p_mat7[5, 5] <- 1

# Combine into observation array [nstate+1 x nstate+1 x n_occasions]
p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)

cat("Observation matrix row sums:\n")
for(i in 1:7) cat("p", i, ":", rowSums(p_arr[,,i]), "\n")

# Initial state vector - all fish start in state 1 (Sacramento mainstem)
rel_vec <- c(1, 0, 0, 0, 0, 0)

#==============================================================================
# SECTION 5: VERIFY MATRICES WITH dDHMMo LIKELIHOOD CHECK
#==============================================================================

all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

cat("NaN count:", sum(is.nan(all_ll)), "\n")
cat("Inf count:", sum(is.infinite(all_ll)), "\n")
cat("LL range:", range(all_ll[!is.nan(all_ll) & !is.infinite(all_ll)]), "\n")

# If NaNs exist investigate which fish
if(sum(is.nan(all_ll)) > 0){
  cat("Fish producing NaN:\n")
  print(detection_history[which(is.nan(all_ll)),
                          c("animal_id", "water_year", "occ_1",
                            "occ_2", "occ_3", "occ_4", "occ_5",
                            "occ_6", "occ_7", "status")])
  print(ch_mat_nimble[which(is.nan(all_ll)), ])
}

#==============================================================================
# SECTION 6: NIMBLE MODEL
# Following Perry et al. (2018) structure
# 5 states: 1=Sacramento, 2=Georgiana, 3=Steamboat/Sutter,
#           4=Dead (absorbing), 5=Failed migration (absorbing)
# Uninformative flat priors dbeta(1,1) on all parameters
# Slice samplers for MCMC following Perry et al.
#==============================================================================

nimCode <- nimbleCode({
  
  #--- PRIORS ---
  S_sac1   ~ dbeta(1, 1)
  S_sac2   ~ dbeta(1, 1)
  S_sac3   ~ dbeta(1, 1)
  S_sac4   ~ dbeta(1, 1)
  S_sac5   ~ dbeta(1, 1)
  S_geo    ~ dbeta(1, 1)
  S_ss     ~ dbeta(1, 1)
  # p_sac1 not estimated - all fish certain to be detected at occ1 by definition
  p_sac2   ~ dbeta(1, 1)
  p_sac3   ~ dbeta(1, 1)
  p_sac4   ~ dbeta(1, 1)
  p_sac5   ~ dbeta(1, 1)
  p_sac6   ~ dbeta(1, 1)
  p_geo    ~ dbeta(1, 1)
  p_ss     ~ dbeta(1, 1)
  lambda   ~ dbeta(1, 1)
  psi_geo  ~ dbeta(1, 1)
  psi_ss   ~ dbeta(1, 1)
  phi_fail ~ dbeta(1, 1)
  
  #--- OBSERVATION ARRAY ---
  # States: 1=Sac, 2=Geo, 3=SS, 4=Dead, 5=Failed, 6=Not-detected
  
  # Occasion 1: Benicia/Carquinez - certain detection
  p_arr[1, 1:6, 1] <- c(1, 0, 0, 0, 0, 0)
  p_arr[2, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 1] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 2: Rio Vista junction
  p_arr[1, 1:6, 2] <- c(p_sac2, 0, 0, 0, 0, (1-p_sac2))
  p_arr[2, 1:6, 2] <- c(0, p_geo, 0, 0, 0, (1-p_geo))
  p_arr[3, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 2] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 3: SR_MOUTH/Steamboat junction
  # Geo pass through undetected, SS detected at mouth receivers
  p_arr[1, 1:6, 3] <- c(p_sac3, 0, 0, 0, 0, (1-p_sac3))
  p_arr[2, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 3] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[4, 1:6, 3] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 3] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 4: SR_BLWSTEAM area - Geo rejoined, SS still in channel
  p_arr[1, 1:6, 4] <- c(p_sac4, 0, 0, 0, 0, (1-p_sac4))
  p_arr[2, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 4] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[4, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 4] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 5: SR_KK345R + SR_FREEPORT - SS has rejoined Sacramento
  p_arr[1, 1:6, 5] <- c(p_sac5, 0, 0, 0, 0, (1-p_sac5))
  p_arr[2, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 5] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 6: upper Sacramento
  p_arr[1, 1:6, 6] <- c(p_sac6, 0, 0, 0, 0, (1-p_sac6))
  p_arr[2, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 6] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 7: spawning ground - certain detection
  p_arr[1, 1:6, 7] <- c(1, 0, 0, 0, 0, 0)
  p_arr[2, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 7] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 7] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  
  #--- TRANSITION ARRAY ---
  
  # Transition 1->2: Benicia to Rio Vista
  tr_arr[1, 1:6, 1] <- c(S_sac1*(1-psi_geo)*(1-phi_fail),
                         S_sac1*psi_geo*(1-phi_fail),
                         0, (1-S_sac1), S_sac1*phi_fail, 0)
  tr_arr[2, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:6, 1] <- c(0, 0, 0, 0, 1, 0)
  tr_arr[6, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  
  # Transition 2->3: Rio Vista to SR_MOUTH/Steamboat
  tr_arr[1, 1:6, 2] <- c(S_sac2*(1-psi_ss)*(1-phi_fail),
                         0, S_sac2*psi_ss*(1-phi_fail),
                         (1-S_sac2), S_sac2*phi_fail, 0)
  tr_arr[2, 1:6, 2] <- c(0, 1, 0, 0, 0, 0)
  tr_arr[3, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:6, 2] <- c(0, 0, 0, 0, 1, 0)
  tr_arr[6, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
  
  # Transition 3->4: SR_MOUTH to SR_BLWSTEAM
  # Geo rejoins Sacramento, SS passes through (can fail)
  tr_arr[1, 1:6, 3] <- c(S_sac3*(1-phi_fail), 0, 0,
                         (1-S_sac3), S_sac3*phi_fail, 0)
  tr_arr[2, 1:6, 3] <- c(S_geo*(1-phi_fail), 0, 0,
                         (1-S_geo), S_geo*phi_fail, 0)
  tr_arr[3, 1:6, 3] <- c(0, 0, (1-phi_fail), 0, phi_fail, 0)
  tr_arr[4, 1:6, 3] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:6, 3] <- c(0, 0, 0, 0, 1, 0)
  tr_arr[6, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  
  # Transition 4->5: SR_BLWSTEAM to SR_KK345R (SS rejoins Sacramento)
  tr_arr[1, 1:6, 4] <- c(S_sac4*(1-phi_fail), 0, 0,
                         (1-S_sac4), S_sac4*phi_fail, 0)
  tr_arr[2, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:6, 4] <- c(S_ss*(1-phi_fail), 0, 0,
                         (1-S_ss), S_ss*phi_fail, 0)
  tr_arr[4, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:6, 4] <- c(0, 0, 0, 0, 1, 0)
  tr_arr[6, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  
  # Transition 5->6: SR_FREEPORT to upper Sacramento
  tr_arr[1, 1:6, 5] <- c(S_sac5*(1-phi_fail), 0, 0,
                         (1-S_sac5), S_sac5*phi_fail, 0)
  tr_arr[2, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:6, 5] <- c(0, 0, 0, 0, 1, 0)
  tr_arr[6, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  
  # Transition 6->7: upper Sacramento to spawning ground
  tr_arr[1, 1:6, 6] <- c(lambda*(1-phi_fail), 0, 0,
                         (1-lambda), lambda*phi_fail, 0)
  tr_arr[2, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:6, 6] <- c(0, 0, 0, 0, 1, 0)
  tr_arr[6, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  
  #--- LIKELIHOOD ---
  for(i in 1:nfish){
    ch_mat[i, 1:7] ~ dDHMMo(
      init         = c(1, 0, 0, 0, 0, 0)[1:6],
      probObs      = p_arr[1:6, 1:6, 1:7],
      probTrans    = tr_arr[1:6, 1:6, 1:6],
      len          = 7,
      checkRowSums = 0
    )
  }
})

#==============================================================================
# SECTION 7: BUILD NIMBLE MODEL AND CONFIGURE MCMC
#==============================================================================

inits <- list(
  S_sac1 = 0.95, S_sac2 = 0.95, S_sac3 = 0.95,
  S_sac4 = 0.95, S_sac5 = 0.95,
  S_geo  = 0.95, S_ss   = 0.95,
  p_sac2 = 0.90, p_sac3 = 0.90,
  p_sac4 = 0.85, p_sac5 = 0.95, p_sac6 = 0.95,
  p_geo  = 0.80, p_ss   = 0.80,
  lambda   = 0.95,
  psi_geo  = 0.10,
  psi_ss   = 0.20,
  phi_fail = 0.10
)

nimMod <- nimbleModel(
  code      = nimCode,
  inits     = inits,
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(nfish = nrow(ch_mat_nimble))
)

nimMod$calculate()

inits_fn <- function(){
  list(
    S_sac1 = runif(1, 0.8, 1), S_sac2 = runif(1, 0.8, 1),
    S_sac3 = runif(1, 0.8, 1), S_sac4 = runif(1, 0.8, 1),
    S_sac5 = runif(1, 0.8, 1),
    S_geo  = runif(1, 0.7, 1), S_ss   = runif(1, 0.7, 1),
    p_sac2 = runif(1, 0.7, 1), p_sac3 = runif(1, 0.7, 1),
    p_sac4 = runif(1, 0.7, 1), p_sac5 = runif(1, 0.8, 1),
    p_sac6 = runif(1, 0.8, 1),
    p_geo  = runif(1, 0.5, 1), p_ss   = runif(1, 0.5, 1),
    lambda   = runif(1, 0.7, 1),
    psi_geo  = runif(1, 0.05, 0.20),
    psi_ss   = runif(1, 0.15, 0.40),
    phi_fail = runif(1, 0.05, 0.20)
  )
}

params <- c("S_sac1", "S_sac2", "S_sac3", "S_sac4", "S_sac5",
            "S_geo", "S_ss",
            "psi_geo", "psi_ss",
            "phi_fail", "lambda",
            "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
            "p_geo", "p_ss")

# Configure MCMC with slice samplers following Perry et al.
confMCMC <- configureMCMC(nimMod, onlySlice = TRUE)
confMCMC$addMonitors(params)
MCMC   <- buildMCMC(confMCMC)
CModel <- compileNimble(nimMod)
CMCMC  <- compileNimble(MCMC, project = CModel)

#==============================================================================
# SECTION 8: RUN MCMC AND SAVE RESULTS
#==============================================================================

# Short test run first to confirm model works before full run
# mcmc_test <- runMCMC(CMCMC, niter = 1000, nchains = 1, nburnin = 100,
#                      thin = 1, inits = list(inits_fn()),
#                      samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test, round = 3)

# Full run - 50000 iterations, 3 chains, 10000 burnin, thin = 10
mcmc_out_full <- runMCMC(
  CMCMC,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn(), inits_fn(), inits_fn()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_full, round = 3)

# Save immediately
save(mcmc_out_full,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_full_run3_221fish_5state.RData")

# Trace plots for key parameters
MCMCtrace(mcmc_out_full,
          params   = c("psi_geo", "psi_ss", "phi_fail", "lambda"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_run3",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)
