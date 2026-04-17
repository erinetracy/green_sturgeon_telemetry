#mulit state gs MODEL
library(nimble)
library(dplyr)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(lubridate)

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
p_sac6  <- 0.95  # detection at upper Sac (occ6) - observed 100%
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


# Recode 0 (not detected) to nstate+1 (7) for dDHMMo
ch_mat_nimble <- ch_mat
ch_mat_nimble[ch_mat_nimble == 0] <- nstate + 1

# Check
table(ch_mat_nimble)

# Test again
dDHMMo(ch_mat_nimble[1, ],
       init = rel_vec,
       probObs = p_arr,
       probTrans = tr_arr,
       len = 7,
       checkRowSums = FALSE,
       log = TRUE)

# Test on up_incomplete fish (should have state 6 in detection history)
inc_rows <- which(detection_history$status == "up_incomplete")
dDHMMo(ch_mat_nimble[inc_rows[1], ],
       init = rel_vec,
       probObs = p_arr,
       probTrans = tr_arr,
       len = 7,
       checkRowSums = FALSE,
       log = TRUE)

# Test on incomplete_dead fish
dead_rows <- which(detection_history$status == "incomplete_dead")
dDHMMo(ch_mat_nimble[dead_rows[1], ],
       init = rel_vec,
       probObs = p_arr,
       probTrans = tr_arr,
       len = 7,
       checkRowSums = FALSE,
       log = TRUE)

# Test on all fish - check none return NaN or -Inf
all_ll <- apply(ch_mat_nimble, 1, function(x) 
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

sum(is.nan(all_ll))   # should be 0
sum(is.infinite(all_ll))  # should be 0
range(all_ll)


#there are alot of problems when trying to run the model 
# Find which fish are returning NaN or -Inf
problem_rows <- which(is.nan(all_ll) | is.infinite(all_ll))
length(problem_rows)

# Look at their detection histories
ch_mat_nimble[problem_rows[1:10], ]

# Check their status
detection_history$status[problem_rows[1:10]]

# What states are in the problem rows
apply(ch_mat_nimble[problem_rows, ], 1, table)

# Fix 1: tr_23 needs state 4 (SS) row - SS doesn't exist yet at occ2->3
# but we need a row so it doesn't sum to 0
# Fix 2: tr_12 needs state 6 (failed) row
tr_12[6, 6] <- 1  # failed fish absorbing from start

# Fix 3: p_mat4 needs state 4 (SS) observation row
# SS fish are passing through at occ4 - undetected
p_mat4[4, nstate+1] <- 1  # already there - check this

# Fix 4: Check tr_23 has SS row
tr_23[4, 4] <- 1  # SS doesn't exist yet but need row to sum to 1

# Check all transition matrix row sums
for(i in 1:6){
  cat("Transition matrix", i, "row sums:\n")
  print(rowSums(tr_arr[,,i]))
}


#try to fix 0s
# Fix all transition matrices to ensure every row sums to 1
# Rule: if a state cannot exist at a given transition, 
# it transitions to dead (state 5) with prob 1
# This is a mathematical necessity - biologically these states 
# simply shouldn't exist at those occasions

# Transition 1->2
tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac1*(1-psi_geo-psi_dcc), S_sac1*psi_geo, S_sac1*psi_dcc, 0, (1-S_sac1), 0, 0)
tr_12[2, 5] <- 1  # Geo shouldn't exist yet -> dead
tr_12[3, 5] <- 1  # DCC shouldn't exist yet -> dead
tr_12[4, 5] <- 1  # SS shouldn't exist yet -> dead
tr_12[5, 5] <- 1  # dead stays dead
tr_12[6, 6] <- 1  # failed stays failed

# Transition 2->3
tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac2*(1-psi_ss), 0, 0, S_sac2*psi_ss, (1-S_sac2), 0, 0)
tr_23[2, 2] <- 1  # Geo pass through
tr_23[3, 3] <- 1  # DCC pass through
tr_23[4, 5] <- 1  # SS shouldn't exist yet -> dead
tr_23[5, 5] <- 1  # dead stays dead
tr_23[6, 6] <- 1  # failed stays failed

# Transition 3->4
tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac3, 0, 0, 0, (1-S_sac3), 0, 0)
tr_34[2, ] <- c(S_geo, 0, 0, 0, (1-S_geo), 0, 0)
tr_34[3, ] <- c(S_dcc, 0, 0, 0, (1-S_dcc), 0, 0)
tr_34[4, 4] <- 1  # SS pass through
tr_34[5, 5] <- 1  # dead stays dead
tr_34[6, 6] <- 1  # failed stays failed

# Transition 4->5
tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac4, 0, 0, 0, (1-S_sac4), 0, 0)
tr_45[2, 5] <- 1  # Geo shouldn't exist -> dead
tr_45[3, 5] <- 1  # DCC shouldn't exist -> dead
tr_45[4, ] <- c(S_ss, 0, 0, 0, (1-S_ss), 0, 0)
tr_45[5, 5] <- 1  # dead stays dead
tr_45[6, 6] <- 1  # failed stays failed

# Transition 5->6
tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac5, 0, 0, 0, (1-S_sac5), 0, 0)
tr_56[2, 5] <- 1  # Geo shouldn't exist -> dead
tr_56[3, 5] <- 1  # DCC shouldn't exist -> dead
tr_56[4, 5] <- 1  # SS shouldn't exist -> dead
tr_56[5, 5] <- 1  # dead stays dead
tr_56[6, 6] <- 1  # failed stays failed

# Transition 6->7
tr_67 <- temp_mat
tr_67[1, ] <- c(lambda, 0, 0, 0, (1-lambda), 0, 0)
tr_67[2, 5] <- 1  # Geo shouldn't exist -> dead
tr_67[3, 5] <- 1  # DCC shouldn't exist -> dead
tr_67[4, 5] <- 1  # SS shouldn't exist -> dead
tr_67[5, 5] <- 1  # dead stays dead
tr_67[6, 6] <- 1  # failed stays failed

# Rebuild transition array
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

# Verify all row sums = 1
for(i in 1:6){
  cat("Transition matrix", i, "row sums:\n")
  print(rowSums(tr_arr[,,i]))
}

# Fix observation matrices - every row must sum to 1
# States that cannot be observed at an occasion 
# get probability 1 of being "not detected" (state 7)

p_mat1 <- temp_mat
p_mat1[1, 1] <- 1
p_mat1[2, nstate+1] <- 1
p_mat1[3, nstate+1] <- 1
p_mat1[4, nstate+1] <- 1
p_mat1[5, 5] <- 1
p_mat1[6, 6] <- 1

p_mat2 <- temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, (1-p_sac2), 0, 0)
p_mat2[2, ] <- c(0, p_geo, 0, 0, (1-p_geo), 0, 0)
p_mat2[3, ] <- c(0, 0, p_dcc, 0, (1-p_dcc), 0, 0)
p_mat2[4, nstate+1] <- 1
p_mat2[5, 5] <- 1
p_mat2[6, 6] <- 1

p_mat3 <- temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, (1-p_sac3), 0, 0)
p_mat3[2, nstate+1] <- 1
p_mat3[3, nstate+1] <- 1
p_mat3[4, ] <- c(0, 0, 0, p_ss, (1-p_ss), 0, 0)
p_mat3[5, 5] <- 1
p_mat3[6, 6] <- 1

p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, (1-p_sac4), 0, 0)
p_mat4[2, nstate+1] <- 1
p_mat4[3, nstate+1] <- 1
p_mat4[4, nstate+1] <- 1
p_mat4[5, 5] <- 1
p_mat4[6, 6] <- 1

p_mat5 <- temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, (1-p_sac5), 0, 0)
p_mat5[2, nstate+1] <- 1
p_mat5[3, nstate+1] <- 1
p_mat5[4, nstate+1] <- 1
p_mat5[5, 5] <- 1
p_mat5[6, 6] <- 1

p_mat6 <- temp_mat
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, (1-p_sac6), 0, 0)
p_mat6[2, nstate+1] <- 1
p_mat6[3, nstate+1] <- 1
p_mat6[4, nstate+1] <- 1
p_mat6[5, 5] <- 1
p_mat6[6, 6] <- 1

p_mat7 <- temp_mat
p_mat7[1, 1] <- 1
p_mat7[2, nstate+1] <- 1
p_mat7[3, nstate+1] <- 1
p_mat7[4, nstate+1] <- 1
p_mat7[5, 5] <- 1
p_mat7[6, 6] <- 1

# Rebuild observation array
p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)

# Test all fish
all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

sum(is.nan(all_ll))
sum(is.infinite(all_ll))
range(all_ll)

# Look at first problem fish detection history
problem_rows <- which(is.nan(all_ll))
ch_mat_nimble[problem_rows[1], ]

# Test it individually with log=FALSE to see the actual probability
dDHMMo(ch_mat_nimble[problem_rows[1], ],
       init = rel_vec,
       probObs = p_arr,
       probTrans = tr_arr,
       len = 7,
       checkRowSums = FALSE,
       log = FALSE)

# Also check what states these fish have
table(ch_mat_nimble[problem_rows, ])

# Get animal_ids of the 50 problem fish
problem_fish <- detection_history[ch_mat_nimble[,3] == 4 & ch_mat_nimble[,4] == 1,
                                  c("animal_id", "water_year")]

# Look at their actual detections at occasions 3 and 4
model1_events %>%
  inner_join(problem_fish, by = c("animal_id", "water_year")) %>%
  filter(occasion %in% c(3, 4)) %>%
  distinct(animal_id, water_year, occasion, location, receiver_group, state, 
           first_detection, last_detection) %>%
  arrange(animal_id, water_year, first_detection)


# Rebuild detection history using LAST detection state at each occasion, still resulted in issues because of back and forth movement
# Rebuild using FIRST detection state at each occasion, made the problem worse
#actually just recode routes but still Use max state detection history to address back and forth movement

# Fix tr_34 - SS fish rejoin Sac here (not pass-through)
tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac3, 0, 0, 0, (1-S_sac3), 0, 0)
tr_34[2, ] <- c(S_geo,  0, 0, 0, (1-S_geo),  0, 0)  # Geo rejoins
tr_34[3, ] <- c(S_dcc,  0, 0, 0, (1-S_dcc),  0, 0)  # DCC rejoins
tr_34[4, ] <- c(S_ss,   0, 0, 0, (1-S_ss),   0, 0)  # SS rejoins here
tr_34[5, 5] <- 1
tr_34[6, 6] <- 1

# Fix tr_45 - SS no longer exists at this transition
# all fish now in state 1
tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac4, 0, 0, 0, (1-S_sac4), 0, 0)
tr_45[2, 5] <- 1  # shouldnt exist
tr_45[3, 5] <- 1  # shouldnt exist
tr_45[4, 5] <- 1  # shouldnt exist - SS already rejoined
tr_45[5, 5] <- 1
tr_45[6, 6] <- 1

# Also fix p_mat4 - SS fish now observed as Sac state 1 at occasion 4
# Remove the SS pass-through undetected row
p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, (1-p_sac4), 0, 0)
p_mat4[2, nstate+1] <- 1  # shouldnt exist
p_mat4[3, nstate+1] <- 1  # shouldnt exist
p_mat4[4, nstate+1] <- 1  # shouldnt exist - SS already rejoined
p_mat4[5, 5] <- 1
p_mat4[6, 6] <- 1

# Rebuild arrays
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)
p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)

# Use max state detection history (back to this approach)
# since last detection was better than first detection
detection_history <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
  arrange(animal_id, water_year, occasion, first_detection) %>%
  group_by(animal_id, water_year, occasion) %>%
  summarise(state = max(state), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = occasion,
    values_from = state,
    names_prefix = "occ_",
    values_fill = 0
  ) %>%
  select(animal_id, water_year, occ_1, occ_2, occ_3, occ_4, occ_5, occ_6, occ_7) %>%
  left_join(
    migration_status %>% select(animal_id, water_year, status),
    by = c("animal_id", "water_year")
  ) %>%
  rowwise() %>%
  mutate(
    last_occ = max(which(c(occ_1, occ_2, occ_3, occ_4, occ_5, occ_6, occ_7) > 0)),
    absorbing_state = case_when(
      status == "incomplete_dead" ~ 5,
      status == "up_incomplete"   ~ 6,
      TRUE ~ NA_real_
    ),
    occ_1 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 1 > last_occ, absorbing_state, occ_1),
    occ_2 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 2 > last_occ, absorbing_state, occ_2),
    occ_3 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 3 > last_occ, absorbing_state, occ_3),
    occ_4 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 4 > last_occ, absorbing_state, occ_4),
    occ_5 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 5 > last_occ, absorbing_state, occ_5),
    occ_6 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 6 > last_occ, absorbing_state, occ_6),
    occ_7 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 7 > last_occ, absorbing_state, occ_7)
  ) %>%
  ungroup() %>%
  select(-last_occ, -absorbing_state)

ch_mat <- detection_history %>%
  select(occ_1:occ_7) %>%
  as.matrix()

ch_mat_nimble <- ch_mat
ch_mat_nimble[ch_mat_nimble == 0] <- nstate + 1

# Test all fish
all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

sum(is.nan(all_ll))
sum(is.infinite(all_ll))
range(all_ll)

p_arr[,,3]
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, (1-p_sac3), 0, 0)
# Fix ALL observation matrices - missed detections go to column 7 not column 5
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, 0, 0, (1-p_sac2))
p_mat2[2, ] <- c(0, p_geo, 0, 0, 0, 0, (1-p_geo))
p_mat2[3, ] <- c(0, 0, p_dcc, 0, 0, 0, (1-p_dcc))

p_mat3[1, ] <- c(p_sac3, 0, 0, 0, 0, 0, (1-p_sac3))
p_mat3[4, ] <- c(0, 0, 0, p_ss, 0, 0, (1-p_ss))

p_mat4[1, ] <- c(p_sac4, 0, 0, 0, 0, 0, (1-p_sac4))
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, 0, 0, (1-p_sac5))
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, 0, 0, (1-p_sac6))

# Rebuild p_arr and test
p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)

all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

sum(is.nan(all_ll))
sum(is.infinite(all_ll))
range(all_ll)


# Add failure probability to all transitions from state 1
phi_fail <- 0.10  # placeholder - will be estimated by NIMBLE

# Transition 1->2
tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac1*(1-psi_geo-psi_dcc)*(1-phi_fail), 
                S_sac1*psi_geo*(1-phi_fail), 
                S_sac1*psi_dcc*(1-phi_fail), 
                0, 
                (1-S_sac1), 
                S_sac1*phi_fail, 
                0)
tr_12[2, 5] <- 1
tr_12[3, 5] <- 1
tr_12[4, 5] <- 1
tr_12[5, 5] <- 1
tr_12[6, 6] <- 1
rowSums(tr_12)

# Transition 2->3
tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac2*(1-psi_ss)*(1-phi_fail), 
                0, 
                0, 
                S_sac2*psi_ss*(1-phi_fail), 
                (1-S_sac2), 
                S_sac2*phi_fail, 
                0)
tr_23[2, 2] <- 1
tr_23[3, 3] <- 1
tr_23[4, 5] <- 1
tr_23[5, 5] <- 1
tr_23[6, 6] <- 1
rowSums(tr_23)

# Transition 3->4
tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac3*(1-phi_fail), 0, 0, 0, (1-S_sac3), S_sac3*phi_fail, 0)
tr_34[2, ] <- c(S_geo*(1-phi_fail),  0, 0, 0, (1-S_geo),  S_geo*phi_fail,  0)
tr_34[3, ] <- c(S_dcc*(1-phi_fail),  0, 0, 0, (1-S_dcc),  S_dcc*phi_fail,  0)
tr_34[4, ] <- c(S_ss*(1-phi_fail),   0, 0, 0, (1-S_ss),   S_ss*phi_fail,   0)
tr_34[5, 5] <- 1
tr_34[6, 6] <- 1
rowSums(tr_34)

# Transition 4->5
tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac4*(1-phi_fail), 0, 0, 0, (1-S_sac4), S_sac4*phi_fail, 0)
tr_45[2, 5] <- 1
tr_45[3, 5] <- 1
tr_45[4, 5] <- 1
tr_45[5, 5] <- 1
tr_45[6, 6] <- 1
rowSums(tr_45)

# Transition 5->6
tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac5*(1-phi_fail), 0, 0, 0, (1-S_sac5), S_sac5*phi_fail, 0)
tr_56[2, 5] <- 1
tr_56[3, 5] <- 1
tr_56[4, 5] <- 1
tr_56[5, 5] <- 1
tr_56[6, 6] <- 1
rowSums(tr_56)

# Transition 6->7
tr_67 <- temp_mat
tr_67[1, ] <- c(lambda*(1-phi_fail), 0, 0, 0, (1-lambda), lambda*phi_fail, 0)
tr_67[2, 5] <- 1
tr_67[3, 5] <- 1
tr_67[4, 5] <- 1
tr_67[5, 5] <- 1
tr_67[6, 6] <- 1
rowSums(tr_67)

# Rebuild transition array
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

# Test all fish
all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

sum(is.nan(all_ll))
sum(is.infinite(all_ll))
range(all_ll)


# Export raw detections for problem fish to examine
problem_fish_data <- model1_events %>%
  inner_join(
    detection_history[problem_rows, c("animal_id", "water_year")],
    by = c("animal_id", "water_year")
  ) %>%
  select(animal_id, water_year, location, receiver_group, 
         state, occasion, first_detection, last_detection) %>%
  arrange(animal_id, water_year, first_detection)

# How many unique fish are we talking about?
problem_fish_data %>%
  distinct(animal_id, water_year) %>%
  nrow()

# Write to csv so you can share it
write.csv(problem_fish_data, 
          "problem_fish_detections.csv", 
          row.names = FALSE)

# Create exploration variable
# IMPORTANT BIOLOGICAL NOTE:
# Fish detected at steamboat_sutter receivers (state 4) that subsequently
# appear at SR_BLWSTEAM or SR_BLWSUTTER on the Sacramento mainstem
# definitively entered Steamboat/Sutter Slough then TURNED AROUND and
# exited back to Sacramento mainstem. SR_BLWSTEAM is BELOW the point
# where Steamboat Slough rejoins the Sacramento River - so detection
# there after SS detection means the fish reversed course back to Sac.
# These fish are correctly coded as explored_ss = TRUE even though
# their final committed route (last detection) is Sacramento (state 1).
# This exploration behavior is preserved separately for IBM analysis
# and to test whether exploratory behavior is flow-dependent.

exploration_summary <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
  arrange(animal_id, water_year, first_detection) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    explored_geo = any(state == 2),
    explored_dcc = any(state == 3),
    explored_ss  = any(state == 4),
    explored_any = any(state %in% c(2, 3, 4)),
    n_alt_detections = sum(state %in% c(2, 3, 4)),
    .groups = "drop"
  )

# Overall exploration summary
exploration_summary %>%
  dplyr::summarise(
    total_fish = n(),
    explored_any = sum(explored_any),
    explored_geo = sum(explored_geo),
    explored_dcc = sum(explored_dcc),
    explored_ss  = sum(explored_ss),
    pct_explored_any = round(mean(explored_any) * 100, 1)
  )

# Fish that explored SS but ultimately took Sac (turned around)
exploration_summary %>%
  left_join(
    detection_history %>% 
      dplyr::select(animal_id, water_year, occ_3) %>%
      mutate(water_year = as.numeric(water_year)),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(committed_ss = occ_3 == 4) %>%
  dplyr::summarise(
    explored_ss_total = sum(explored_ss),
    committed_ss_total = sum(committed_ss, na.rm = TRUE),
    explored_but_turned_around = sum(explored_ss & !committed_ss, na.rm = TRUE)
  )

exploration_summary %>%
  dplyr::summarise(
    total_fish = n(),
    explored_any = sum(explored_any),
    explored_geo = sum(explored_geo),
    explored_dcc = sum(explored_dcc),
    explored_ss  = sum(explored_ss),
    pct_explored_any = round(sum(explored_any) / n() * 100, 1)
  )

# Rebuild detection history using LAST detection at each occasion
# This captures the route the fish ultimately committed to
detection_history <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
  arrange(animal_id, water_year, occasion, first_detection) %>%
  group_by(animal_id, water_year, occasion) %>%
  summarise(state = last(state), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = occasion,
    values_from = state,
    names_prefix = "occ_",
    values_fill = 0
  ) %>%
  select(animal_id, water_year, occ_1, occ_2, occ_3, occ_4, occ_5, occ_6, occ_7) %>%
  left_join(
    migration_status %>% select(animal_id, water_year, status),
    by = c("animal_id", "water_year")
  ) %>%
  # Join exploration variable
  left_join(exploration_summary, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    last_occ = max(which(c(occ_1, occ_2, occ_3, occ_4, occ_5, occ_6, occ_7) > 0)),
    absorbing_state = case_when(
      status == "incomplete_dead" ~ 5,
      status == "up_incomplete"   ~ 6,
      TRUE ~ NA_real_
    ),
    occ_1 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 1 > last_occ, absorbing_state, occ_1),
    occ_2 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 2 > last_occ, absorbing_state, occ_2),
    occ_3 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 3 > last_occ, absorbing_state, occ_3),
    occ_4 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 4 > last_occ, absorbing_state, occ_4),
    occ_5 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 5 > last_occ, absorbing_state, occ_5),
    occ_6 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 6 > last_occ, absorbing_state, occ_6),
    occ_7 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 7 > last_occ, absorbing_state, occ_7)
  ) %>%
  ungroup() %>%
  select(-last_occ, -absorbing_state)

# Rebuild ch_mat
ch_mat <- detection_history %>%
  select(occ_1:occ_7) %>%
  as.matrix()

ch_mat_nimble <- ch_mat
ch_mat_nimble[ch_mat_nimble == 0] <- nstate + 1

# Check NaN with last detection approach
all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

sum(is.nan(all_ll))
sum(is.infinite(all_ll))

# Check how many fish explored alternative routes
table(exploration_summary$explored_any)
table(exploration_summary$explored_geo)
table(exploration_summary$explored_ss)
table(exploration_summary$explored_dcc)

