#mulit state gs MODEL
library(nimble)
library(dplyr)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(lubridate)
library(abind)

# Build Multi-state model
#==============================================================================
# BUILD TRANSITION AND OBSERVATION MATRICES
# States: 1=Sac, 2=Geo, 3=DCC, 4=SS, 5=Dead, 6=Failed
# 7 occasions, so 6 transition matrices (between occasions)
# and 7 observation matrices (one per occasion)
#==============================================================================

nstate <- 6

# Placeholder parameter values for testing matrix structure
# These will be estimated by NIMBLE later
S_sac1  <- 0.95  # survival Sac occ1->2 (Benicia to Rio Vista)
S_sac2  <- 0.95  # survival Sac occ2->3 (Rio Vista to SR_MOUTH)
S_sac3  <- 0.95  # survival Sac occ3->4 (SR_MOUTH to SR_BLWSTEAM)
S_sac4  <- 0.95  # survival Sac occ4->5 (SR_BLWSTEAM to SR_FREEPORT)
S_sac5  <- 0.95  # survival Sac occ5->6 (SR_FREEPORT to upper Sac)
S_geo   <- 0.95  # survival through Georgiana
S_dcc   <- 0.95  # survival through DCC
S_ss    <- 0.95  # survival through Steamboat/Sutter
psi_geo <- 0.15  # routing prob to Georgiana at occ2
psi_dcc <- 0.10  # routing prob to DCC at occ2
psi_ss  <- 0.30  # routing prob to Steamboat/Sutter at occ3
lambda  <- 0.95  # survival occ6->7 (upper Sac to spawning ground)

# Detection probabilities
p_sac1  <- 0.99  # detection at Benicia/Carquinez (occ1) - set high, near certain
p_sac2  <- 0.90  # detection at Rio Vista (occ2)
p_sac3  <- 0.90  # detection at SR_MOUTH (occ3)
p_sac4  <- 0.85  # detection at SR_BLWSTEAM (occ4)
p_sac5  <- 0.95  # detection at SR_FREEPORT (occ5)
p_sac6  <- 0.99  # detection at upper Sac (occ6) - observed 100%
p_geo   <- 0.80  # detection at Georgiana receivers
p_dcc   <- 0.80  # detection at DCC receivers
p_ss    <- 0.80  # detection at Steamboat/Sutter receivers

#------------------------------------------------------------------------------
# TRANSITION MATRICES
# Rows = from state, Columns = to state
# Row sums must = 1
# States 5 (dead) and 6 (failed) are absorbing
#------------------------------------------------------------------------------

# Template - all zeros, absorbing state stays absorbing
temp_mat <- matrix(0, nrow = nstate + 1, ncol = nstate + 1)
temp_mat[nstate + 1, nstate + 1] <- 1

# Transition 1->2: Benicia to Rio Vista junction
# Sac fish: survive and stay Sac, or enter Geo, or enter DCC, or die
# No other states exist yet
tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac1*(1-psi_geo-psi_dcc), S_sac1*psi_geo, S_sac1*psi_dcc, 0, (1-S_sac1), 0, 0)
# Check row sums
rowSums(tr_12)

# Transition 2->3: Rio Vista to SR_MOUTH/Steamboat junction
# Sac fish: survive and stay Sac or enter SS, or die
# Geo fish: pass through (stay Geo with prob 1, no survival penalty yet)
# DCC fish: pass through (stay DCC with prob 1, no survival penalty yet)
# Failed fish: absorbing
tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac2*(1-psi_ss), 0, 0, S_sac2*psi_ss, (1-S_sac2), 0, 0)
tr_23[2, 2] <- 1  # Geo fish pass through
tr_23[3, 3] <- 1  # DCC fish pass through
tr_23[6, 6] <- 1  # Failed fish absorbing
rowSums(tr_23)

# Transition 3->4: SR_MOUTH to SR_BLWSTEAM (Geo/DCC rejoin)
# Sac fish: survive
# Geo fish: rejoin Sac with survival S_geo
# DCC fish: rejoin Sac with survival S_dcc
# SS fish: pass through (still in Steamboat/Sutter)
# Failed fish: absorbing
tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac3, 0, 0, 0, (1-S_sac3), 0, 0)
tr_34[2, ] <- c(S_geo,  0, 0, 0, (1-S_geo),  0, 0)  # Geo rejoins Sac
tr_34[3, ] <- c(S_dcc,  0, 0, 0, (1-S_dcc),  0, 0)  # DCC rejoins Sac
tr_34[4, 4] <- 1  # SS fish pass through
tr_34[6, 6] <- 1  # Failed fish absorbing
rowSums(tr_34)

# Transition 4->5: SR_BLWSTEAM to SR_FREEPORT (SS rejoins)
# Sac fish: survive
# SS fish: rejoin Sac with survival S_ss
# Failed fish: absorbing
tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac4, 0, 0, 0, (1-S_sac4), 0, 0)
tr_45[4, ] <- c(S_ss,   0, 0, 0, (1-S_ss),   0, 0)  # SS rejoins Sac
tr_45[6, 6] <- 1  # Failed fish absorbing
rowSums(tr_45)

# Transition 5->6: SR_FREEPORT to upper Sac
# All fish now in Sac state 1
tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac5, 0, 0, 0, (1-S_sac5), 0, 0)
tr_56[6, 6] <- 1  # Failed fish absorbing
rowSums(tr_56)

# Transition 6->7: upper Sac to spawning ground
# lambda = survival/probability of reaching spawning ground
tr_67 <- temp_mat
tr_67[1, ] <- c(lambda, 0, 0, 0, (1-lambda), 0, 0)
tr_67[6, 6] <- 1  # Failed fish absorbing
rowSums(tr_67)

# Combine into transition array
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)
dim(tr_arr)  # should be 7 x 7 x 6

#------------------------------------------------------------------------------
# OBSERVATION MATRICES
# Rows = true state, Columns = observed state
# p_obs[true, observed]
#------------------------------------------------------------------------------

# Occasion 1: Benicia/Carquinez - all fish in Sac state, detection = 1
p_mat1 <- temp_mat
p_mat1[1, 1] <- 1  # certain detection at start

# Occasion 2: Rio Vista junction
# Sac fish detected at Rio Vista with prob p_sac2
# Geo fish detected at Georgiana with prob p_geo
# DCC fish detected at DCC receivers with prob p_dcc
# SS, Dead, Failed not present at this occasion
p_mat2 <- temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, (1-p_sac2), 0, 0)
p_mat2[2, ] <- c(0, p_geo, 0, 0, (1-p_geo),   0, 0)
p_mat2[3, ] <- c(0, 0, p_dcc, 0, (1-p_dcc),   0, 0)

# Occasion 3: SR_MOUTH/Steamboat junction
# Sac fish detected at SR_MOUTH with prob p_sac3
# Geo fish: pass through, no receivers here -> not detected (prob=1 to unobserved)
# DCC fish: pass through, no receivers here -> not detected
# SS fish detected at Steamboat/Sutter receivers with prob p_ss
p_mat3 <- temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, (1-p_sac3), 0, 0)
p_mat3[2, nstate+1] <- 1  # Geo pass through - undetected
p_mat3[3, nstate+1] <- 1  # DCC pass through - undetected
p_mat3[4, ] <- c(0, 0, 0, p_ss, (1-p_ss), 0, 0)

# Occasion 4: SR_BLWSTEAM (Geo/DCC rejoined)
# All fish now in Sac state 1 (Geo/DCC rejoined in transition 3->4)
# SS still in Steamboat/Sutter - pass through, no detection
p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, (1-p_sac4), 0, 0)
p_mat4[4, nstate+1] <- 1  # SS pass through - undetected

# Occasion 5: SR_FREEPORT (SS rejoined)
# All fish now in Sac state 1
p_mat5 <- temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, (1-p_sac5), 0, 0)

# Occasion 6: upper Sac
p_mat6 <- temp_mat
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, (1-p_sac6), 0, 0)

# Occasion 7: spawning ground - detection = 1 if alive
p_mat7 <- temp_mat
p_mat7[1, 1] <- 1  # certain detection at spawning ground

# Combine into observation array
p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)
dim(p_arr)  # should be 7 x 7 x 7

# Initial state vector - all fish start in Sac state 1
rel_vec <- c(1, 0, 0, 0, 0, 0, 0)

# Test on a single fish using Perry's dDHMMo function
dDHMMo(ch_mat[1, ],
       init = rel_vec,
       probObs = p_arr,
       probTrans = tr_arr,
       len = 7,
       checkRowSums = FALSE,
       log = TRUE)