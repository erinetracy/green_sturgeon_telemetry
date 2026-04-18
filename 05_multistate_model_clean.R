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
#       ch_mat_nimble: 228 x 7 detection history matrix formatted for dDHMMo
#       ch_mat: 228 x 7 detection history matrix (0 = not detected)
#       detection_history: fish-level detection history with metadata
#       nstate: number of states (6)
#       fish_info: fish-level metadata
#       exploration_summary: fish-level exploration behavior
#
# OUTPUTS:
#   - mcmc_out_full: MCMC posterior samples
#   - gs_mcmc_full_run2.RData: saved MCMC results
#   - gs_mcmc_traces_run2.pdf: trace plots for key parameters
#
# MODEL STRUCTURE:
#   States: 1=Sacramento, 2=Georgiana, 3=DCC, 4=Steamboat/Sutter,
#           5=Dead (absorbing), 6=Failed migration (absorbing)
#   Occasions: 7 (Benicia -> Rio Vista -> SR_MOUTH -> SR_BLWSTEAM ->
#              SR_KK345R/SR_FREEPORT -> Upper Sac -> Spawning Ground)
#   Years: 2007-2017 (Model 1, no Yolo Bypass)
#   Fish: 228 (140 up_complete, 83 up_incomplete, 5 incomplete_dead)
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
# SECTION 2: VERIFY DETECTION HISTORY
#==============================================================================

cat("Total fish:", nrow(ch_mat_nimble), "\n")  # should be 228
cat("\nStatus breakdown:\n")
print(table(detection_history$status))
cat("\nDetection rates per occasion:\n")
print(round(colMeans(ch_mat > 0), 3))
cat("\nState counts in ch_mat_nimble:\n")
print(table(ch_mat_nimble))
cat("\nGeorgiana fish:", apply(ch_mat_nimble, 1, function(x) any(x == 2)) %>% sum(), "\n")
cat("DCC fish:", apply(ch_mat_nimble, 1, function(x) any(x == 3)) %>% sum(), "\n")
cat("Steamboat/Sutter fish:", apply(ch_mat_nimble, 1, function(x) any(x == 4)) %>% sum(), "\n")

#==============================================================================
# SECTION 3: BUILD TRANSITION AND OBSERVATION MATRICES
# Placeholder parameter values used for likelihood verification only.
# All parameters estimated by NIMBLE MCMC in Section 6.
#==============================================================================

# Placeholder parameter values
S_sac1 <- 0.95; S_sac2 <- 0.95; S_sac3 <- 0.95
S_sac4 <- 0.95; S_sac5 <- 0.95
S_geo  <- 0.95; S_dcc  <- 0.95; S_ss   <- 0.95
psi_geo <- 0.15; psi_dcc <- 0.10; psi_ss <- 0.30
phi_fail <- 0.10; lambda <- 0.95
p_sac2 <- 0.90; p_sac3 <- 0.90; p_sac4 <- 0.85
p_sac5 <- 0.95; p_sac6 <- 0.95
p_geo  <- 0.80; p_dcc  <- 0.80; p_ss   <- 0.80

# Template matrix - all zeros, not-detected state stays absorbing
temp_mat <- matrix(0, nrow = nstate + 1, ncol = nstate + 1)
temp_mat[nstate + 1, nstate + 1] <- 1

#------------------------------------------------------------------------------
# TRANSITION MATRICES
# Rows = from state, Columns = to state, Row sums must = 1
#
# KEY STRUCTURAL DECISIONS:
#
# FAILED MIGRATION (state 6) possible at every transition from state 1:
#   Green sturgeon can abort migration at any point after Benicia/Carquinez.
#   phi_fail estimated from up_incomplete fish.
#
# GEO AND DCC FISH REJOIN SACRAMENTO AT tr_34 (occ3->occ4):
#   Georgiana and DCC both exit back to Sacramento upstream of SR_BLWSTEAM.
#   Their last detection was in the Geo/DCC channel (state 2/3) so transition
#   codes them to state 1.
#
# SS FISH PASS THROUGH AT tr_34 (occ3->occ4), REJOIN AT tr_45 (occ4->occ5):
#   SR_BLWSTEAM and SR_BLWSUTTER are BELOW where Steamboat and Sutter Slough
#   rejoin the Sacramento River. SS fish that went all the way through have
#   NOT rejoined yet at occasion 4 — they are still in the channel (state 4
#   pass-through). SR_KK345R (occasion 5) is the first receiver ABOVE both
#   rejoin points, confirming SS fish are back on Sacramento at occ5.
#
# STATES THAT CANNOT EXIST at a given transition are mapped to dead (state 5):
#   Mathematical requirement — all rows must sum to 1. These impossible states
#   are never actually occupied so the choice does not affect the likelihood.
#
# MISSED DETECTIONS map to column 7 (not detected), NOT column 5 (dead):
#   A fish not detected at an occasion is not dead. Mapping to column 5 would
#   incorrectly assign zero probability to fish reappearing at later occasions.
#------------------------------------------------------------------------------

# Transition 1->2: Benicia to Rio Vista junction
tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac1*(1-psi_geo-psi_dcc)*(1-phi_fail),
                S_sac1*psi_geo*(1-phi_fail),
                S_sac1*psi_dcc*(1-phi_fail),
                0, (1-S_sac1), S_sac1*phi_fail, 0)
tr_12[2, 5] <- 1  # Geo shouldn't exist yet
tr_12[3, 5] <- 1  # DCC shouldn't exist yet
tr_12[4, 5] <- 1  # SS shouldn't exist yet
tr_12[5, 5] <- 1
tr_12[6, 6] <- 1

# Transition 2->3: Rio Vista to SR_MOUTH/Steamboat junction
tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac2*(1-psi_ss)*(1-phi_fail), 0, 0,
                S_sac2*psi_ss*(1-phi_fail), (1-S_sac2), S_sac2*phi_fail, 0)
tr_23[2, 2] <- 1  # Geo pass through - still in Georgiana channel
tr_23[3, 3] <- 1  # DCC pass through - still in DCC channel
tr_23[4, 5] <- 1  # SS shouldn't exist yet
tr_23[5, 5] <- 1
tr_23[6, 6] <- 1

# Transition 3->4: SR_MOUTH to SR_BLWSTEAM
# Geo/DCC rejoin Sacramento, SS passes through
tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac3*(1-phi_fail), 0, 0, 0, (1-S_sac3), S_sac3*phi_fail, 0)
tr_34[2, ] <- c(S_geo*(1-phi_fail),  0, 0, 0, (1-S_geo),  S_geo*phi_fail,  0)
tr_34[3, ] <- c(S_dcc*(1-phi_fail),  0, 0, 0, (1-S_dcc),  S_dcc*phi_fail,  0)
tr_34[4, ] <- c(0, 0, 0, (1-phi_fail), 0, phi_fail, 0)  # SS pass through, can fail
tr_34[5, 5] <- 1
tr_34[6, 6] <- 1

# Transition 4->5: SR_BLWSTEAM to SR_KK345R (SS rejoins Sacramento)
tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac4*(1-phi_fail), 0, 0, 0, (1-S_sac4), S_sac4*phi_fail, 0)
tr_45[2, 5] <- 1  # shouldn't exist - Geo rejoined at occ4
tr_45[3, 5] <- 1  # shouldn't exist - DCC rejoined at occ4
tr_45[4, ] <- c(S_ss*(1-phi_fail), 0, 0, 0, (1-S_ss), S_ss*phi_fail, 0)
tr_45[5, 5] <- 1
tr_45[6, 6] <- 1

# Transition 5->6: SR_FREEPORT to upper Sacramento
tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac5*(1-phi_fail), 0, 0, 0, (1-S_sac5), S_sac5*phi_fail, 0)
tr_56[2, 5] <- 1; tr_56[3, 5] <- 1; tr_56[4, 5] <- 1
tr_56[5, 5] <- 1; tr_56[6, 6] <- 1

# Transition 6->7: upper Sacramento to spawning ground
tr_67 <- temp_mat
tr_67[1, ] <- c(lambda*(1-phi_fail), 0, 0, 0, (1-lambda), lambda*phi_fail, 0)
tr_67[2, 5] <- 1; tr_67[3, 5] <- 1; tr_67[4, 5] <- 1
tr_67[5, 5] <- 1; tr_67[6, 6] <- 1

# Combine into transition array [nstate+1 x nstate+1 x n_transitions]
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

cat("Transition matrix row sums:\n")
for(i in 1:6) cat("tr", i, ":", rowSums(tr_arr[,,i]), "\n")

#------------------------------------------------------------------------------
# OBSERVATION MATRICES
# Rows = true state, Columns = observed state
# Missed detections -> column 7 (not detected)
#------------------------------------------------------------------------------

# Occasion 1: Benicia/Carquinez - certain detection
p_mat1 <- temp_mat
p_mat1[1, 1] <- 1
p_mat1[2, nstate+1] <- 1; p_mat1[3, nstate+1] <- 1; p_mat1[4, nstate+1] <- 1
p_mat1[5, 5] <- 1; p_mat1[6, 6] <- 1

# Occasion 2: Rio Vista junction
p_mat2 <- temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, 0, 0, (1-p_sac2))
p_mat2[2, ] <- c(0, p_geo, 0, 0, 0, 0, (1-p_geo))
p_mat2[3, ] <- c(0, 0, p_dcc, 0, 0, 0, (1-p_dcc))
p_mat2[4, nstate+1] <- 1  # SS shouldn't exist yet
p_mat2[5, 5] <- 1; p_mat2[6, 6] <- 1

# Occasion 3: SR_MOUTH/Steamboat junction
# Geo and DCC pass through undetected, SS detected at mouth receivers
p_mat3 <- temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, 0, 0, (1-p_sac3))
p_mat3[2, nstate+1] <- 1  # Geo pass through - undetected
p_mat3[3, nstate+1] <- 1  # DCC pass through - undetected
p_mat3[4, ] <- c(0, 0, 0, p_ss, 0, 0, (1-p_ss))
p_mat3[5, 5] <- 1; p_mat3[6, 6] <- 1

# Occasion 4: SR_BLWSTEAM area
# Geo/DCC rejoined Sacramento, SS still in channel
p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, 0, 0, (1-p_sac4))
p_mat4[2, nstate+1] <- 1  # shouldn't exist - Geo rejoined
p_mat4[3, nstate+1] <- 1  # shouldn't exist - DCC rejoined
p_mat4[4, ] <- c(0, 0, 0, p_ss, 0, 0, (1-p_ss))
p_mat4[5, 5] <- 1; p_mat4[6, 6] <- 1

# Occasion 5: SR_KK345R + SR_FREEPORT - SS has rejoined Sacramento
p_mat5 <- temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, 0, 0, (1-p_sac5))
p_mat5[2, nstate+1] <- 1; p_mat5[3, nstate+1] <- 1; p_mat5[4, nstate+1] <- 1
p_mat5[5, 5] <- 1; p_mat5[6, 6] <- 1

# Occasion 6: upper Sacramento
p_mat6 <- temp_mat
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, 0, 0, (1-p_sac6))
p_mat6[2, nstate+1] <- 1; p_mat6[3, nstate+1] <- 1; p_mat6[4, nstate+1] <- 1
p_mat6[5, 5] <- 1; p_mat6[6, 6] <- 1

# Occasion 7: spawning ground - certain detection
p_mat7 <- temp_mat
p_mat7[1, 1] <- 1
p_mat7[2, nstate+1] <- 1; p_mat7[3, nstate+1] <- 1; p_mat7[4, nstate+1] <- 1
p_mat7[5, 5] <- 1; p_mat7[6, 6] <- 1

# Combine into observation array [nstate+1 x nstate+1 x n_occasions]
p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)

cat("Observation matrix row sums:\n")
for(i in 1:7) cat("p", i, ":", rowSums(p_arr[,,i]), "\n")

# Initial state vector - all fish start in state 1 (Sacramento mainstem)
rel_vec <- c(1, 0, 0, 0, 0, 0, 0)

#==============================================================================
# SECTION 4: VERIFY MATRICES WITH dDHMMo LIKELIHOOD CHECK
#==============================================================================

all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

cat("NaN count:", sum(is.nan(all_ll)), "\n")
cat("Inf count:", sum(is.infinite(all_ll)), "\n")
cat("LL range:", range(all_ll[!is.nan(all_ll) & !is.infinite(all_ll)]), "\n")

#==============================================================================
# SECTION 5: NIMBLE MODEL
# Following Perry et al. (2018) structure
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
  S_dcc    ~ dbeta(1, 1)
  S_ss     ~ dbeta(1, 1)
  # p_sac1 not estimated - all fish certain to be detected at occ1 by definition
  p_sac2   ~ dbeta(1, 1)
  p_sac3   ~ dbeta(1, 1)
  p_sac4   ~ dbeta(1, 1)
  p_sac5   ~ dbeta(1, 1)
  p_sac6   ~ dbeta(1, 1)
  p_geo    ~ dbeta(1, 1)
  p_dcc    ~ dbeta(1, 1)
  p_ss     ~ dbeta(1, 1)
  lambda   ~ dbeta(1, 1)
  psi_geo  ~ dbeta(1, 1)
  psi_dcc  ~ dbeta(1, 1)
  psi_ss   ~ dbeta(1, 1)
  phi_fail ~ dbeta(1, 1)
  
  #--- OBSERVATION ARRAY ---
  # Occasion 1: Benicia/Carquinez - certain detection
  p_arr[1, 1:7, 1] <- c(1, 0, 0, 0, 0, 0, 0)
  p_arr[2, 1:7, 1] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[3, 1:7, 1] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[4, 1:7, 1] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[5, 1:7, 1] <- c(0, 0, 0, 0, 1, 0, 0)
  p_arr[6, 1:7, 1] <- c(0, 0, 0, 0, 0, 1, 0)
  p_arr[7, 1:7, 1] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Occasion 2: Rio Vista junction
  p_arr[1, 1:7, 2] <- c(p_sac2, 0, 0, 0, 0, 0, (1-p_sac2))
  p_arr[2, 1:7, 2] <- c(0, p_geo, 0, 0, 0, 0, (1-p_geo))
  p_arr[3, 1:7, 2] <- c(0, 0, p_dcc, 0, 0, 0, (1-p_dcc))
  p_arr[4, 1:7, 2] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[5, 1:7, 2] <- c(0, 0, 0, 0, 1, 0, 0)
  p_arr[6, 1:7, 2] <- c(0, 0, 0, 0, 0, 1, 0)
  p_arr[7, 1:7, 2] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Occasion 3: SR_MOUTH/Steamboat junction
  p_arr[1, 1:7, 3] <- c(p_sac3, 0, 0, 0, 0, 0, (1-p_sac3))
  p_arr[2, 1:7, 3] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[3, 1:7, 3] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[4, 1:7, 3] <- c(0, 0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[5, 1:7, 3] <- c(0, 0, 0, 0, 1, 0, 0)
  p_arr[6, 1:7, 3] <- c(0, 0, 0, 0, 0, 1, 0)
  p_arr[7, 1:7, 3] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Occasion 4: SR_BLWSTEAM area - Geo/DCC rejoined, SS still in channel
  p_arr[1, 1:7, 4] <- c(p_sac4, 0, 0, 0, 0, 0, (1-p_sac4))
  p_arr[2, 1:7, 4] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[3, 1:7, 4] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[4, 1:7, 4] <- c(0, 0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[5, 1:7, 4] <- c(0, 0, 0, 0, 1, 0, 0)
  p_arr[6, 1:7, 4] <- c(0, 0, 0, 0, 0, 1, 0)
  p_arr[7, 1:7, 4] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Occasion 5: SR_KK345R + SR_FREEPORT - SS has rejoined Sacramento
  p_arr[1, 1:7, 5] <- c(p_sac5, 0, 0, 0, 0, 0, (1-p_sac5))
  p_arr[2, 1:7, 5] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[3, 1:7, 5] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[4, 1:7, 5] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[5, 1:7, 5] <- c(0, 0, 0, 0, 1, 0, 0)
  p_arr[6, 1:7, 5] <- c(0, 0, 0, 0, 0, 1, 0)
  p_arr[7, 1:7, 5] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Occasion 6: upper Sacramento
  p_arr[1, 1:7, 6] <- c(p_sac6, 0, 0, 0, 0, 0, (1-p_sac6))
  p_arr[2, 1:7, 6] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[3, 1:7, 6] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[4, 1:7, 6] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[5, 1:7, 6] <- c(0, 0, 0, 0, 1, 0, 0)
  p_arr[6, 1:7, 6] <- c(0, 0, 0, 0, 0, 1, 0)
  p_arr[7, 1:7, 6] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Occasion 7: spawning ground - certain detection
  p_arr[1, 1:7, 7] <- c(1, 0, 0, 0, 0, 0, 0)
  p_arr[2, 1:7, 7] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[3, 1:7, 7] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[4, 1:7, 7] <- c(0, 0, 0, 0, 0, 0, 1)
  p_arr[5, 1:7, 7] <- c(0, 0, 0, 0, 1, 0, 0)
  p_arr[6, 1:7, 7] <- c(0, 0, 0, 0, 0, 1, 0)
  p_arr[7, 1:7, 7] <- c(0, 0, 0, 0, 0, 0, 1)
  
  #--- TRANSITION ARRAY ---
  # Transition 1->2: Benicia to Rio Vista
  tr_arr[1, 1:7, 1] <- c(S_sac1*(1-psi_geo-psi_dcc)*(1-phi_fail),
                         S_sac1*psi_geo*(1-phi_fail),
                         S_sac1*psi_dcc*(1-phi_fail),
                         0, (1-S_sac1), S_sac1*phi_fail, 0)
  tr_arr[2, 1:7, 1] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:7, 1] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:7, 1] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:7, 1] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[6, 1:7, 1] <- c(0, 0, 0, 0, 0, 1, 0)
  tr_arr[7, 1:7, 1] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Transition 2->3: Rio Vista to SR_MOUTH/Steamboat
  tr_arr[1, 1:7, 2] <- c(S_sac2*(1-psi_ss)*(1-phi_fail),
                         0, 0, S_sac2*psi_ss*(1-phi_fail),
                         (1-S_sac2), S_sac2*phi_fail, 0)
  tr_arr[2, 1:7, 2] <- c(0, 1, 0, 0, 0, 0, 0)
  tr_arr[3, 1:7, 2] <- c(0, 0, 1, 0, 0, 0, 0)
  tr_arr[4, 1:7, 2] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:7, 2] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[6, 1:7, 2] <- c(0, 0, 0, 0, 0, 1, 0)
  tr_arr[7, 1:7, 2] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Transition 3->4: SR_MOUTH to SR_BLWSTEAM
  # Geo/DCC rejoin Sacramento, SS passes through (can fail)
  tr_arr[1, 1:7, 3] <- c(S_sac3*(1-phi_fail), 0, 0, 0,
                         (1-S_sac3), S_sac3*phi_fail, 0)
  tr_arr[2, 1:7, 3] <- c(S_geo*(1-phi_fail), 0, 0, 0,
                         (1-S_geo), S_geo*phi_fail, 0)
  tr_arr[3, 1:7, 3] <- c(S_dcc*(1-phi_fail), 0, 0, 0,
                         (1-S_dcc), S_dcc*phi_fail, 0)
  tr_arr[4, 1:7, 3] <- c(0, 0, 0, (1-phi_fail), 0, phi_fail, 0)
  tr_arr[5, 1:7, 3] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[6, 1:7, 3] <- c(0, 0, 0, 0, 0, 1, 0)
  tr_arr[7, 1:7, 3] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Transition 4->5: SR_BLWSTEAM to SR_KK345R (SS rejoins Sacramento)
  tr_arr[1, 1:7, 4] <- c(S_sac4*(1-phi_fail), 0, 0, 0,
                         (1-S_sac4), S_sac4*phi_fail, 0)
  tr_arr[2, 1:7, 4] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:7, 4] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:7, 4] <- c(S_ss*(1-phi_fail), 0, 0, 0,
                         (1-S_ss), S_ss*phi_fail, 0)
  tr_arr[5, 1:7, 4] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[6, 1:7, 4] <- c(0, 0, 0, 0, 0, 1, 0)
  tr_arr[7, 1:7, 4] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Transition 5->6: SR_FREEPORT to upper Sacramento
  tr_arr[1, 1:7, 5] <- c(S_sac5*(1-phi_fail), 0, 0, 0,
                         (1-S_sac5), S_sac5*phi_fail, 0)
  tr_arr[2, 1:7, 5] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:7, 5] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:7, 5] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:7, 5] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[6, 1:7, 5] <- c(0, 0, 0, 0, 0, 1, 0)
  tr_arr[7, 1:7, 5] <- c(0, 0, 0, 0, 0, 0, 1)
  
  # Transition 6->7: upper Sacramento to spawning ground
  tr_arr[1, 1:7, 6] <- c(lambda*(1-phi_fail), 0, 0, 0,
                         (1-lambda), lambda*phi_fail, 0)
  tr_arr[2, 1:7, 6] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[3, 1:7, 6] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[4, 1:7, 6] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[5, 1:7, 6] <- c(0, 0, 0, 0, 1, 0, 0)
  tr_arr[6, 1:7, 6] <- c(0, 0, 0, 0, 0, 1, 0)
  tr_arr[7, 1:7, 6] <- c(0, 0, 0, 0, 0, 0, 1)
  
  #--- LIKELIHOOD ---
  for(i in 1:nfish){
    ch_mat[i, 1:7] ~ dDHMMo(
      init         = c(1, 0, 0, 0, 0, 0, 0)[1:7],
      probObs      = p_arr[1:7, 1:7, 1:7],
      probTrans    = tr_arr[1:7, 1:7, 1:6],
      len          = 7,
      checkRowSums = 0
    )
  }
})

#==============================================================================
# SECTION 6: BUILD NIMBLE MODEL AND CONFIGURE MCMC
#==============================================================================

inits <- list(
  S_sac1 = 0.95, S_sac2 = 0.95, S_sac3 = 0.95,
  S_sac4 = 0.95, S_sac5 = 0.95,
  S_geo  = 0.95, S_dcc  = 0.95, S_ss   = 0.95,
  p_sac2 = 0.90, p_sac3 = 0.90,
  p_sac4 = 0.85, p_sac5 = 0.95, p_sac6 = 0.95,
  p_geo  = 0.80, p_dcc  = 0.80, p_ss   = 0.80,
  lambda   = 0.95,
  psi_geo  = 0.10,
  psi_dcc  = 0.05,
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
    S_geo  = runif(1, 0.7, 1), S_dcc  = runif(1, 0.7, 1),
    S_ss   = runif(1, 0.7, 1),
    p_sac2 = runif(1, 0.7, 1), p_sac3 = runif(1, 0.7, 1),
    p_sac4 = runif(1, 0.7, 1), p_sac5 = runif(1, 0.8, 1),
    p_sac6 = runif(1, 0.8, 1),
    p_geo  = runif(1, 0.5, 1), p_dcc  = runif(1, 0.5, 1),
    p_ss   = runif(1, 0.5, 1),
    lambda   = runif(1, 0.7, 1),
    psi_geo  = runif(1, 0.05, 0.15),
    psi_dcc  = runif(1, 0.02, 0.08),
    psi_ss   = runif(1, 0.15, 0.35),
    phi_fail = runif(1, 0.05, 0.20)
  )
}

params <- c("S_sac1", "S_sac2", "S_sac3", "S_sac4", "S_sac5",
            "S_geo", "S_dcc", "S_ss",
            "psi_geo", "psi_dcc", "psi_ss",
            "phi_fail", "lambda",
            "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
            "p_geo", "p_dcc", "p_ss")

# Configure MCMC with slice samplers following Perry et al.
# Slice samplers handle constrained probability parameters better
# than default Metropolis-Hastings
confMCMC <- configureMCMC(nimMod, onlySlice = TRUE)
confMCMC$addMonitors(params)
MCMC   <- buildMCMC(confMCMC)
CModel <- compileNimble(nimMod)
CMCMC  <- compileNimble(MCMC, project = CModel)

#==============================================================================
# SECTION 7: RUN MCMC AND SAVE RESULTS
#==============================================================================

# Run short test first to confirm model works
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
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_full_run2_228fish.RData")

# Trace plots for key parameters
MCMCtrace(mcmc_out_full,
          params   = c("psi_geo", "psi_dcc", "psi_ss", "phi_fail", "lambda"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_run2",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)

