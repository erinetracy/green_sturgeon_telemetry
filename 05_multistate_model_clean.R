#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL
# Script: 05_multistate_model.R
# Author: Erin Tracy
# Last updated: March 2026
#
# PURPOSE:
# Build a discrete hidden Markov multistate model (dDHMMo) to estimate
# route-specific survival, detection, and routing probabilities for adult
# green sturgeon making upstream spawning migrations in the Sacramento River
# system. Follows Perry et al. framework using nimbleEcology.
#
# INPUTS (from 04_multistate_data_prep.R):
#   - filtered_events: detection events for up_complete, up_incomplete,
#                      incomplete_dead fish
#   - migration_status: fish-level migration classification
#   - model1_events: filtered to 2007-2017 with occasions and states assigned
#
# OUTPUTS:
#   - detection_history: 222 x 7 detection history matrix
#   - exploration_summary: fish-level exploration behavior variable
#   - ch_mat_nimble: matrix formatted for dDHMMo
#   - tr_arr, p_arr: transition and observation arrays
#
# MODEL STRUCTURE:
#   States: 1=Sacramento, 2=Georgiana, 3=DCC, 4=Steamboat/Sutter,
#           5=Dead (absorbing), 6=Failed migration (absorbing)
#   Occasions: 7 (Benicia -> Rio Vista -> SR_MOUTH -> SR_BLWSTEAM ->
#              SR_KK345R/SR_FREEPORT -> Upper Sac -> Spawning Ground)
#   Years: 2007-2017 (Model 1, no Yolo Bypass)
#   Fish: 222 (138 up_complete, 77 up_incomplete, 7 incomplete_dead)
#==============================================================================

library(nimble)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(dplyr)
library(lubridate)

#==============================================================================
# SECTION 1: FILTER TO MODEL 1 FISH
#==============================================================================

# Filter to model years and relevant migration statuses
# DECISION: Include up_complete, up_incomplete, and incomplete_dead
# RATIONALE: up_incomplete fish provide information about where migrations
#   fail (coded as state 6 absorbing). incomplete_dead fish provide survival
#   information (coded as state 5 absorbing). Excluding them would bias
#   survival estimates upward.
# NOTE: Model 1 is 2007-2017 without Yolo Bypass because Georgiana,
#   DCC, and Steamboat/Sutter receivers were removed after 2017, making
#   the full route model impossible beyond that year.

model1_events <- filtered_events %>%
  filter(status %in% c("up_complete", "incomplete_dead", "up_incomplete"),
         water_year >= 2007, water_year <= 2017)

#==============================================================================
# SECTION 2: RECEIVER GROUP RELABELING FOR DCC
#==============================================================================

# DECISION: Rename mok_deltacross to DCC and restrict to SR_DCC and SR_DCC2 only
# RATIONALE: The mok_deltacross receiver group contained receivers at multiple
#   locations including Mok River confluences where fish could branch onto the
#   Mok River WITHOUT entering the Delta Cross Channel. Only SR_DCC and SR_DCC2
#   are physically located inside the DCC channel itself. Detection there
#   unambiguously confirms DCC entry.
# TESTED AND REJECTED: Using all mok_deltacross receivers - would have included
#   fish that may have been on Mok River, not DCC specifically.

model1_events <- model1_events %>%
  mutate(receiver_group = case_when(
    receiver_group == "mok_deltacross" ~ "DCC",
    TRUE ~ receiver_group
  )) %>%
  mutate(receiver_group = case_when(
    receiver_group == "DCC" &
      !location %in% c("SR_DCC", "SR_DCC2") ~ NA_character_,
    TRUE ~ receiver_group
  ))

#==============================================================================
# SECTION 3: ASSIGN OCCASIONS
#==============================================================================

# Occasions are defined by receiver location representing key decision points
# along the upstream migration route.
#
# OCCASION STRUCTURE:
# Occ 1: Benicia + Carquinez - migration start, all fish enter here
# Occ 2: Rio Vista (Sac) OR Georgiana OR SR_DCC/SR_DCC2 - first junction
# Occ 3: SR_MOUTH (Sac) OR Steamboat/Sutter MOUTH receivers - second junction
# Occ 4: SR_BLWSTEAM area (Sac) - Geo/DCC rejoined, SS still passing through channel
# Occ 5: SR_KK345R + SR_FREEPORT area (Sac) - SS has now rejoined Sacramento
# Occ 6: SR_BLWCHIBEND + SR_BUTTEBR (Sac) - upper Sacramento
# Occ 7: Spawning ground - terminal state
#
# KEY DECISIONS ON RECEIVER ASSIGNMENT:
#
# DECKER ISLAND EXCLUDED from occasion 2:
#   Fish can split into interior delta side channels AT Decker, making it
#   ambiguous as a route indicator. Rio Vista is north of this split and
#   confirms fish stayed on Sacramento mainstem past the interior delta entry.
#
# RIO VISTA used as first decision point (not Decker):
#   SR_RV receivers at ~38.15-38.17N are past the Decker/interior delta
#   split and before the Steamboat/Sutter junction (~38.18N).
#
# STEAMBOAT/SUTTER MOUTH receivers at occasion 3 (not occasion 4):
#   STEAMBOATSL_MOUTH1/2 at ~38.18N are the entry point into Steamboat Slough
#   from the Sacramento River. Fish detected here have entered the slough.
#
# SR_BLWSTEAM, SR_BLWGEORGIA, SR_BLWSUTTER at occasion 4:
#   BLW = BELOW. These Sacramento mainstem receivers are DOWNSTREAM of where
#   each respective slough rejoins the Sacramento River.
#   - SR_BLWGEORGIA: below Georgiana rejoin — Geo fish detected here only if
#     they turned around. With last(state) coding, those fish are already
#     coded as Sac (state 1) since their last detection was back on Sac.
#   - SR_BLWSTEAM: below Steamboat rejoin — SS fish detected here only if
#     they turned around. Same logic — last(state) codes them as Sac.
#   - SR_BLWSUTTER: below Sutter rejoin — same logic as above.
#   Fish that went ALL THE WAY THROUGH Steamboat/Sutter WITHOUT turning
#   around are still in the channel at occasion 4 (state 4, pass-through).
#   They have NOT yet rejoined the Sacramento at SR_BLWSTEAM.
#
# SR_KK345R at occasion 5:
#   This receiver is ABOVE where both Steamboat AND Sutter Slough rejoin
#   the Sacramento River. Fish detected here are confirmed back on the
#   Sacramento mainstem regardless of which route they took. This is where
#   SS fish that went all the way through are first confirmed to have rejoined.
#
# THEREFORE:
#   tr_34: Geo/DCC rejoin Sacramento, SS fish PASS THROUGH (still in channel)
#   tr_45: SS fish rejoin Sacramento (confirmed at SR_KK345R and above)

model1_events <- model1_events %>%
  mutate(occasion = case_when(
    
    # Occasion 1: migration start
    receiver_group %in% c("benicia", "carquinez") ~ 1,
    
    # Occasion 2: first junction
    # Georgiana and DCC = entered delta route
    # Rio Vista Sacramento receivers = stayed on mainstem past delta entry
    receiver_group %in% c("georgiana", "DCC") ~ 2,
    receiver_group == "sacramento" &
      location %in% c("SR_RV10_7L", "SR_RV125L", "RIOVISTABR01", "SR_RV127L",
                      "RIOVISTABR02", "RIOVISTABR03", "SR_RV169L", "SR_RV169R") ~ 2,
    
    # Occasion 3: second junction
    # Steamboat/Sutter MOUTH receivers only = entered SS from Sacramento
    # SR_MOUTH = stayed on Sacramento past SS entry point
    receiver_group == "steamboat_sutter" &
      location %in% c("STEAMBOATSL_MOUTH1", "STEAMBOATSL_MOUTH2") ~ 3,
    receiver_group == "sacramento" &
      location %in% c("SR_MOUTH_2", "SR_MOUTH", "SR_RV150R") ~ 3,
    
    # Occasion 4: SR_BLWSTEAM area
    # Geo/DCC fish have rejoined Sacramento — their last detection was in
    #   Geo/DCC channel so tr_34 transitions them to state 1 here.
    # SS fish still in channel — pass-through in tr_34, not yet rejoined.
    #   SR_BLWSTEAM and SR_BLWSUTTER are BELOW the SS rejoin points so
    #   only SS fish that TURNED AROUND would ping these Sac receivers.
    #   With last(state) coding those fish are already coded state 1 (Sac).
    # Interior SS receivers = fish still traveling through Steamboat/Sutter
    receiver_group == "steamboat_sutter" &
      location %in% c("SUTSTMSLS1", "SUTSTMSLS2", "SUTSTMSLS",
                      "SUTSL_BLWMINERSL", "STEAMBOATSL",
                      "STMSL_MARI2", "STMSL_MARI1B", "STMSL_MARI",
                      "SUTSLOUGH", "SUTSLOUGH2", "SUTSLOUGH1") ~ 4,
    receiver_group == "sacramento" &
      location %in% c("SR_KK240L", "SR_RYDE", "SR_BLWGEORGIA2", "SR_BLWGEORGIA",
                      "SR_KK250L", "SR_DCCSOUTH", "SR_DCCSOUTH2", "SR_KK269L",
                      "SR_DCCNORTH", "SR_BLWSTEAM2", "SR_BLWSTEAM",
                      "SR_BLWSUTTER", "SR_BLWSUTTER2") ~ 4,
    
    # Occasion 5: SR_KK345R + SR_FREEPORT area
    # SR_KK345R is the first receiver ABOVE both Steamboat AND Sutter rejoin points.
    # SS fish confirmed back on Sacramento mainstem from SR_KK345R onwards.
    # tr_45 transitions SS fish (state 4) to state 1 (Sac) at this occasion.
    receiver_group == "sacramento" &
      location %in% c("SR_KK345R",
                      "SR_GB437R", "SR_GB447R", "SR_FREEPORT", "SR_FREEPORT_1",
                      "SR_GB470L", "SR_FREEPORTDIV_W", "SR_FREEPORTDIV_S",
                      "SR_FREEPORTDIV_N", "SR_ABVFREEPORTW", "SR_GB479R",
                      "SR_ABVFREEPORTE", "SR_ABVFREEPORTE-2", "SR_GB502L",
                      "SR_GB503R", "SR_CCFBRW", "SR_CCFBRW_RT", "SR_CCFBRE",
                      "SR_CCFBR", "SR_TOWERBRE", "SR_TOWERBRE_RT", "SR_TOWERBRW",
                      "SR_TOWERBRW_RT", "SR_ISTBRIDGE_E", "SR_ISTBRIDGE_W",
                      "SR_ISTBRIDGE_E-1", "SR_SRWTPD_VPS_SW", "SR_SRWTPD_VPS_SE",
                      "SR_SRWTPD_VPS_WSW", "SR_SRWTPD_VPS_ESE", "SR_SRWTPD_S",
                      "SR_SRWTPD_VPS_W", "SR_SRWTPD_E", "SR_SRWTPD_N-1",
                      "SR_SRWTPD_VPS_WNW", "SR_SRWTPD_VPS_ENE", "SR_SRWTPD_OLD",
                      "SR_SRWTPD_VPS_NW", "SR_BLWAMERICANW", "SR_BLWAMERICANE",
                      "SR_ABVAMERICANW", "SR_ABVAMERICANE", "SR_RIVERVIEW_MARIW",
                      "SR_RIVERVIEW_MARIE") ~ 5,
    
    # Occasion 6: upper Sacramento (Yolo rejoins here in Model 2)
    receiver_group == "sacramento" &
      location %in% c("SR_RM69-1", "SR_RM69-2", "SR_RM69-5", "SR_RM69-7",
                      "SR_RM69-8", "SR_RM69-10", "SR_RM69-12", "SR_EH699R",
                      "SR_RM69-13", "SR_RM69-14", "SR_RM69-15", "SR_RM72-1",
                      "SR_RM72-2", "SR_RM72-3", "SR_RM72-4", "SR_RM72-5",
                      "SR_RM72-6", "SR_RM85-1", "SR_RM85-2", "SR_ABVFEATHER",
                      "SR_KL856R", "SR_RM85-3", "SR_RM85-4", "SR_RM85-5",
                      "SR_ABVFEATHER1", "SR_ABVFEATHER1_RT", "SR_ABVFEATHER2",
                      "SR_BLWFEATHERR2", "SR_BLWFEATHERR", "SR_FEATHER2_RT",
                      "SR_KNIGHTSBR_W", "SR_KNIGHTSBRS_RT", "SR_KNIGHTSBR",
                      "SR_KNIGHTSBRN_RT", "SR_KNIGHTSLAND", "SR_RM92-1",
                      "SR_RM92-2", "SR_RM92-3", "SR_KL919R", "SR_BLWCHIBEND_W2",
                      "SR_BLWCHIBEND_E", "SR_BLWCHIBEND_W", "SR_ABVTISDALE_E2",
                      "SR_ABVTISDALE_E1", "SR_ABVTISDALE_E1_RT", "SR_MERIDIANBR2",
                      "SR_MERIDIANBR", "SR_BLWWARDSE", "SR_BLWWARDS",
                      "SR_ABVCOLUSABR1", "SR_ABVCOLUSABR1_RT", "SR_ABVCOLUSABR2",
                      "SR_ABVCOLUSABR2_RT", "SR_BLWBUTTE1", "SR_BLWBUTTE2",
                      "SR_BUTTEBR_E", "SR_BUTTEBR", "SR_BLWORD2", "SR_BLWORD1",
                      "SR_ORDBEND", "SR_ABVORDBR2", "SR_ABVCHICOCK") ~ 6,
    
    # Occasion 7: spawning ground - terminal state
    receiver_group == "spawning_ground" ~ 7,
    
    TRUE ~ NA_real_
  ))

#==============================================================================
# SECTION 4: ASSIGN STATES
#==============================================================================

model1_events <- model1_events %>%
  mutate(state = case_when(
    receiver_group %in% c("benicia", "carquinez", "sacramento", "spawning_ground") ~ 1,
    receiver_group == "georgiana"        ~ 2,
    receiver_group == "DCC"              ~ 3,
    receiver_group == "steamboat_sutter" ~ 4,
    TRUE ~ NA_real_
  ))

#==============================================================================
# SECTION 5: EXPLORATION BEHAVIOR VARIABLE
#==============================================================================

# BIOLOGICAL FINDING: Many green sturgeon exhibit extended staging behavior
# at junctions, making repeated exploratory forays into alternative routes
# over periods of days to weeks before committing to a migration route.
# 47.3% of fish (105/222) explored at least one alternative route.
#
# IMPORTANT: Fish detected at steamboat_sutter receivers that subsequently
# appear at SR_BLWSTEAM on the Sacramento mainstem definitively entered
# Steamboat/Sutter Slough then TURNED AROUND. SR_BLWSTEAM is BELOW where
# Steamboat Slough rejoins the Sacramento - detection there after SS detection
# means the fish reversed course. These fish are coded explored_ss = TRUE
# even though their final committed route is Sacramento (state 1).
#
# This variable is preserved for IBM analysis and to test whether
# exploratory behavior itself is flow-dependent.

exploration_summary <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
  arrange(animal_id, water_year, first_detection) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    explored_geo     = any(state == 2),   # ever detected at Georgiana receivers
    explored_dcc     = any(state == 3),   # ever detected at SR_DCC/SR_DCC2
    explored_ss      = any(state == 4),   # ever detected at SS receivers
    explored_any     = any(state %in% c(2, 3, 4)),
    n_alt_detections = sum(state %in% c(2, 3, 4)),
    .groups = "drop"
  )

# Summary: 105/222 fish (47.3%) explored alternative routes
# explored_geo: 27 fish, explored_dcc: 13 fish, explored_ss: 82 fish

#==============================================================================
# SECTION 6: BUILD DETECTION HISTORY MATRIX
#==============================================================================

# DECISION: Use last(state) at each occasion for detection history
# RATIONALE: Captures the route the fish ULTIMATELY COMMITTED TO rather
#   than exploratory forays. For the IBM we want P(fish takes route | flow)
#   which is best represented by the final route choice, not exploration.
#   Exploration behavior is preserved separately in exploration_summary.
#
# TESTED AND REJECTED:
#   - max(state): only 1 fish had within-occasion Sac + alternative conflict.
#     max(state) would have coded these as alternative route even if fish
#     ultimately returned to Sac. Not appropriate for last-committed-route goal.
#   - first(state): made the NaN problem worse (46 NaN vs 20 with last).
#     First detection is before staging exploration completes.
#
# NOTE ON ABSORBING STATES:
#   up_incomplete fish: state 6 (failed) filled forward after last detection
#   incomplete_dead fish: state 5 (dead) filled forward after last detection
#   Specific water years assigned to dead fish based on individual
#   detection history review (not all years for multi-year tagged fish).

detection_history <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
  arrange(animal_id, water_year, occasion, first_detection) %>%
  group_by(animal_id, water_year, occasion) %>%
  dplyr::summarise(state = last(state), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = occasion,
    values_from = state,
    names_prefix = "occ_",
    values_fill = 0   # 0 = not detected, recoded to nstate+1 for dDHMMo
  ) %>%
  dplyr::select(animal_id, water_year,
                occ_1, occ_2, occ_3, occ_4, occ_5, occ_6, occ_7) %>%
  left_join(
    migration_status %>% dplyr::select(animal_id, water_year, status),
    by = c("animal_id", "water_year")
  ) %>%
  left_join(exploration_summary, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    # Find last occasion where fish was detected (state > 0)
    last_occ = max(which(c(occ_1, occ_2, occ_3, occ_4, occ_5, occ_6, occ_7) > 0)),
    # Assign absorbing state for non-complete fish
    absorbing_state = case_when(
      status == "incomplete_dead" ~ 5,
      status == "up_incomplete"   ~ 6,
      TRUE ~ NA_real_
    ),
    # Fill absorbing state forward after last detection
    occ_1 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 1 > last_occ, absorbing_state, occ_1),
    occ_2 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 2 > last_occ, absorbing_state, occ_2),
    occ_3 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 3 > last_occ, absorbing_state, occ_3),
    occ_4 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 4 > last_occ, absorbing_state, occ_4),
    occ_5 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 5 > last_occ, absorbing_state, occ_5),
    occ_6 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 6 > last_occ, absorbing_state, occ_6),
    occ_7 = ifelse(status %in% c("incomplete_dead","up_incomplete") & 7 > last_occ, absorbing_state, occ_7)
  ) %>%
  ungroup() %>%
  dplyr::select(-last_occ, -absorbing_state)

# Store fish metadata separately
fish_info <- detection_history %>%
  dplyr::select(animal_id, water_year, status, explored_geo, explored_dcc,
                explored_ss, explored_any, n_alt_detections)

# Convert to numeric matrix for dDHMMo
# Recode 0 (not detected) to nstate+1 (7) as required by dDHMMo
ch_mat <- detection_history %>%
  dplyr::select(occ_1:occ_7) %>%
  as.matrix()

nstate <- 6
ch_mat_nimble <- ch_mat
ch_mat_nimble[ch_mat_nimble == 0] <- nstate + 1

# Verify matrix
cat("Matrix dimensions:", dim(ch_mat_nimble), "\n")
cat("Status breakdown:\n")
print(table(detection_history$status))
cat("Detection rates per occasion:\n")
print(colMeans(ch_mat > 0))
cat("State counts across matrix:\n")
print(table(ch_mat_nimble))

#==============================================================================
# SECTION 7: TRANSITION AND OBSERVATION MATRICES
#==============================================================================

# Placeholder parameter values for testing matrix structure
# All parameters will be estimated by NIMBLE MCMC
S_sac1  <- 0.95  # survival Sac occ1->2 (Benicia to Rio Vista)
S_sac2  <- 0.95  # survival Sac occ2->3 (Rio Vista to SR_MOUTH)
S_sac3  <- 0.95  # survival Sac occ3->4 (SR_MOUTH to SR_BLWSTEAM)
S_sac4  <- 0.95  # survival Sac occ4->5 (SR_BLWSTEAM to SR_FREEPORT)
S_sac5  <- 0.95  # survival Sac occ5->6 (SR_FREEPORT to upper Sac)
S_geo   <- 0.95  # survival through Georgiana Slough
S_dcc   <- 0.95  # survival through Delta Cross Channel
S_ss    <- 0.95  # survival through Steamboat/Sutter
psi_geo <- 0.15  # routing probability to Georgiana at occ2
psi_dcc <- 0.10  # routing probability to DCC at occ2
psi_ss  <- 0.30  # routing probability to Steamboat/Sutter at occ3
phi_fail <- 0.10 # probability of failed migration at each occasion
lambda  <- 0.95  # survival occ6->7 (upper Sac to spawning ground)

# Detection probabilities (placeholder values)
p_sac1  <- 0.99  # Benicia/Carquinez (occ1) - near certain
p_sac2  <- 0.90  # Rio Vista (occ2)
p_sac3  <- 0.90  # SR_MOUTH (occ3)
p_sac4  <- 0.85  # SR_BLWSTEAM area (occ4)
p_sac5  <- 0.95  # SR_FREEPORT area (occ5)
p_sac6  <- 0.95  # upper Sac (occ6)
p_geo   <- 0.80  # Georgiana receivers
p_dcc   <- 0.80  # DCC receivers (SR_DCC, SR_DCC2)
p_ss    <- 0.80  # Steamboat/Sutter receivers

# Template matrix - all zeros, absorbing state (row 7) stays absorbing
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
#   phi_fail estimated from 77 up_incomplete fish.
#
# GEO AND DCC FISH REJOIN SACRAMENTO AT tr_34 (occ3->occ4):
#   Geographic review confirms Georgiana and DCC both exit back to Sacramento
#   upstream of the SR_BLWSTEAM receiver cluster. Their last detection was
#   in the Geo/DCC channel (state 2/3) so transition codes them to state 1.
#
# SS FISH PASS THROUGH AT tr_34 (occ3->occ4), REJOIN AT tr_45 (occ4->occ5):
#   SR_BLWSTEAM and SR_BLWSUTTER are BELOW where Steamboat and Sutter Slough
#   rejoin the Sacramento River. SS fish that went all the way through have
#   NOT rejoined yet at occasion 4 — they are still in the channel (state 4
#   pass-through). SR_KK345R (occasion 5) is the first receiver ABOVE both
#   rejoin points, so SS fish are confirmed back on Sacramento at occ5.
#   NOTE: SS fish that TURNED AROUND are already coded as Sac (state 1) via
#   last(state) detection history — their turnaround is captured in
#   explored_ss = TRUE in exploration_summary.
#
# STATES THAT CANNOT EXIST at a given transition are mapped to dead (state 5):
#   Mathematical requirement - all rows must sum to 1. These impossible states
#   are never actually occupied so the choice does not affect the likelihood.
#
# MISSED DETECTIONS map to column 7 (not detected), NOT column 5 (dead):
#   A fish not detected at an occasion is not dead. Column 7 accounts for
#   imperfect detection probability p. Mapping to column 5 would incorrectly
#   assign zero probability to fish reappearing at later occasions.
#------------------------------------------------------------------------------

# Transition 1->2: Benicia to Rio Vista junction
# Sac fish can route to Geo, DCC, fail, or die. No other states exist yet.
tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac1*(1-psi_geo-psi_dcc)*(1-phi_fail),
                S_sac1*psi_geo*(1-phi_fail),
                S_sac1*psi_dcc*(1-phi_fail),
                0,
                (1-S_sac1),
                S_sac1*phi_fail,
                0)
tr_12[2, 5] <- 1  # Geo shouldn't exist yet -> dead (mathematical necessity)
tr_12[3, 5] <- 1  # DCC shouldn't exist yet -> dead
tr_12[4, 5] <- 1  # SS shouldn't exist yet -> dead
tr_12[5, 5] <- 1  # dead stays dead
tr_12[6, 6] <- 1  # failed stays failed

# Transition 2->3: Rio Vista to SR_MOUTH/Steamboat junction
# Sac fish can route to SS, fail, or die.
# Geo and DCC fish pass through (no receivers at this occasion).
tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac2*(1-psi_ss)*(1-phi_fail),
                0,
                0,
                S_sac2*psi_ss*(1-phi_fail),
                (1-S_sac2),
                S_sac2*phi_fail,
                0)
tr_23[2, 2] <- 1  # Geo pass through - still in Georgiana channel
tr_23[3, 3] <- 1  # DCC pass through - still in DCC channel
tr_23[4, 5] <- 1  # SS shouldn't exist yet -> dead
tr_23[5, 5] <- 1
tr_23[6, 6] <- 1

# Transition 3->4: SR_MOUTH to SR_BLWSTEAM
# Geo and DCC fish rejoin Sacramento here (exited their channels upstream of SR_BLWSTEAM)
# SS fish PASS THROUGH - still in Steamboat/Sutter channel, not yet rejoined
# Sac fish can fail or die
tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac3*(1-phi_fail), 0, 0, 0, (1-S_sac3), S_sac3*phi_fail, 0)
tr_34[2, ] <- c(S_geo*(1-phi_fail),  0, 0, 0, (1-S_geo),  S_geo*phi_fail,  0)  # Geo rejoins Sac
tr_34[3, ] <- c(S_dcc*(1-phi_fail),  0, 0, 0, (1-S_dcc),  S_dcc*phi_fail,  0)  # DCC rejoins Sac
tr_34[4, 4] <- 1  # SS pass through - still in channel, rejoins at occ5
tr_34[5, 5] <- 1
tr_34[6, 6] <- 1

# Transition 4->5: SR_BLWSTEAM to SR_KK345R (SS rejoins)
# SR_KK345R is first receiver above both Steamboat and Sutter rejoin points
# SS fish finally back on Sacramento mainstem here
tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac4*(1-phi_fail), 0, 0, 0, (1-S_sac4), S_sac4*phi_fail, 0)
tr_45[2, 5] <- 1  # shouldn't exist - Geo rejoined at occ4
tr_45[3, 5] <- 1  # shouldn't exist - DCC rejoined at occ4
tr_45[4, ] <- c(S_ss*(1-phi_fail), 0, 0, 0, (1-S_ss), S_ss*phi_fail, 0)  # SS rejoins Sac
tr_45[5, 5] <- 1
tr_45[6, 6] <- 1

# Transition 5->6: SR_FREEPORT to upper Sac
tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac5*(1-phi_fail), 0, 0, 0, (1-S_sac5), S_sac5*phi_fail, 0)
tr_56[2, 5] <- 1
tr_56[3, 5] <- 1
tr_56[4, 5] <- 1
tr_56[5, 5] <- 1
tr_56[6, 6] <- 1

# Transition 6->7: upper Sac to spawning ground
# lambda = probability of reaching spawning ground given survival
tr_67 <- temp_mat
tr_67[1, ] <- c(lambda*(1-phi_fail), 0, 0, 0, (1-lambda), lambda*phi_fail, 0)
tr_67[2, 5] <- 1
tr_67[3, 5] <- 1
tr_67[4, 5] <- 1
tr_67[5, 5] <- 1
tr_67[6, 6] <- 1

# Combine into transition array [nstate+1 x nstate+1 x n_transitions]
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

# Verify all row sums = 1
for(i in 1:6){
  cat("Transition matrix", i, "row sums:", rowSums(tr_arr[,,i]), "\n")
}

#------------------------------------------------------------------------------
# OBSERVATION MATRICES
# Rows = true state, Columns = observed state
# Missed detections -> column 7 (not detected), NOT column 5 (dead)
#------------------------------------------------------------------------------

# Occasion 1: Benicia/Carquinez - all fish start in state 1, detection = 1
p_mat1 <- temp_mat
p_mat1[1, 1] <- 1
p_mat1[2, nstate+1] <- 1  # shouldn't exist at occ1
p_mat1[3, nstate+1] <- 1
p_mat1[4, nstate+1] <- 1
p_mat1[5, 5] <- 1
p_mat1[6, 6] <- 1

# Occasion 2: Rio Vista junction
# Sac fish detected at Rio Vista with prob p_sac2
# Geo fish detected at Georgiana with prob p_geo
# DCC fish detected at SR_DCC/SR_DCC2 with prob p_dcc
p_mat2 <- temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, 0, 0, (1-p_sac2))
p_mat2[2, ] <- c(0, p_geo, 0, 0, 0, 0, (1-p_geo))
p_mat2[3, ] <- c(0, 0, p_dcc, 0, 0, 0, (1-p_dcc))
p_mat2[4, nstate+1] <- 1  # SS shouldn't exist yet
p_mat2[5, 5] <- 1
p_mat2[6, 6] <- 1

# Occasion 3: SR_MOUTH/Steamboat junction
# Geo and DCC fish pass through - no receivers, undetected
# SS fish detected at SS mouth receivers with prob p_ss
p_mat3 <- temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, 0, 0, (1-p_sac3))
p_mat3[2, nstate+1] <- 1  # Geo pass through - undetected
p_mat3[3, nstate+1] <- 1  # DCC pass through - undetected
p_mat3[4, ] <- c(0, 0, 0, p_ss, 0, 0, (1-p_ss))
p_mat3[5, 5] <- 1
p_mat3[6, 6] <- 1

# Occasion 4: SR_BLWSTEAM area
# Geo/DCC fish have rejoined Sacramento — detected as state 1
# SS fish still in channel — pass-through, detected as state 4 at interior
#   SS receivers or undetected at Sac receivers (they are not on Sac yet)
p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, 0, 0, (1-p_sac4))
p_mat4[2, nstate+1] <- 1  # shouldn't exist - Geo rejoined
p_mat4[3, nstate+1] <- 1  # shouldn't exist - DCC rejoined
p_mat4[4, ] <- c(0, 0, 0, p_ss, 0, 0, (1-p_ss))  # SS in channel, detected at interior SS receivers
p_mat4[5, 5] <- 1
p_mat4[6, 6] <- 1

# Occasion 5: SR_KK345R + SR_FREEPORT area - SS has rejoined Sacramento
# All fish now state 1 after tr_45 transitions SS to Sac
p_mat5 <- temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, 0, 0, (1-p_sac5))
p_mat5[2, nstate+1] <- 1  # shouldn't exist
p_mat5[3, nstate+1] <- 1  # shouldn't exist
p_mat5[4, nstate+1] <- 1  # shouldn't exist - SS rejoined at occ5 transition
p_mat5[5, 5] <- 1
p_mat5[6, 6] <- 1

# Occasion 6: upper Sacramento
p_mat6 <- temp_mat
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, 0, 0, (1-p_sac6))
p_mat6[2, nstate+1] <- 1
p_mat6[3, nstate+1] <- 1
p_mat6[4, nstate+1] <- 1
p_mat6[5, 5] <- 1
p_mat6[6, 6] <- 1

# Occasion 7: spawning ground - detection = 1 if fish reached here
p_mat7 <- temp_mat
p_mat7[1, 1] <- 1
p_mat7[2, nstate+1] <- 1
p_mat7[3, nstate+1] <- 1
p_mat7[4, nstate+1] <- 1
p_mat7[5, 5] <- 1
p_mat7[6, 6] <- 1

# Combine into observation array [nstate+1 x nstate+1 x n_occasions]
p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)

# Verify all row sums = 1
for(i in 1:7){
  cat("Observation matrix", i, "row sums:", rowSums(p_arr[,,i]), "\n")
}

# Initial state vector - all fish start in state 1 (Sacramento)
rel_vec <- c(1, 0, 0, 0, 0, 0, 0)

#==============================================================================
# SECTION 8: VERIFY MATRICES WITH dDHMMo
#==============================================================================

# Test likelihood on all fish
all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

cat("NaN count:", sum(is.nan(all_ll)), "\n")
cat("Inf count:", sum(is.infinite(all_ll)), "\n")
cat("LL range:", range(all_ll[!is.nan(all_ll) & !is.infinite(all_ll)]), "\n")

# NOTE: Some fish may still produce NaN due to temporal reversal staging
# behavior (fish visiting occ3 receivers before occ2 receivers in real time).
# These represent genuine biological complexity - fish staging at junctions
# and making back-and-forth movements over days to weeks.

# Fix problem 1: SS at occ3 then Sac at occ4 - recode occ3 to Sac
# These fish turned around, exploration preserved in explored_ss
ch_mat_nimble[ch_mat_nimble[,3] == 4 & ch_mat_nimble[,4] == 1, 3] <- 1

# Fix problem 2: Geo/DCC at occ2 then Sac at occ3 - recode occ2 to Sac
# These fish turned around, exploration preserved in explored_geo/dcc
ch_mat_nimble[(ch_mat_nimble[,2] %in% c(2,3)) & ch_mat_nimble[,3] == 1, 2] <- 1

# Fix problem 3: Sac at occ3 then SS at occ4 - recode occ4 to not-detected
# Fish entered SS after SR_MOUTH, cannot be accommodated by model structure
ch_mat_nimble[ch_mat_nimble[,3] == 1 & ch_mat_nimble[,4] == 4, 4] <- nstate + 1

# Fix tr_34 - SS fish can fail between occ3 and occ4
tr_34[4, ] <- c(0, 0, 0, (1-phi_fail), 0, phi_fail, 0)

# Rebuild transition array
tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

# Retest
all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

cat("NaN count:", sum(is.nan(all_ll)), "\n")
cat("Inf count:", sum(is.infinite(all_ll)), "\n")
cat("LL range:", range(all_ll[!is.nan(all_ll) & !is.infinite(all_ll)]), "\n")

#==============================================================================
# SECTION 9: SAVE OBJECTS FOR NIMBLE MODEL
#==============================================================================

save(detection_history,
     ch_mat,
     ch_mat_nimble,
     fish_info,
     exploration_summary,
     tr_arr,
     p_arr,
     rel_vec,
     nstate,
     model1_events,
     migration_status,
     file = "gs_multistate_model1_objects.RData")

cat("All objects saved to gs_multistate_model1_objects.RData\n")

#==============================================================================
# SECTION 10: NIMBLE MODEL
#==============================================================================

# Bundle data for NIMBLE
nimble_data <- list(
  y = ch_mat_nimble  # detection history matrix [nfish x nocc]
)

nimble_constants <- list(
  nfish    = nrow(ch_mat_nimble),
  nocc     = 7,
  nstate   = nstate + 1,  # includes unobserved state (7)
  init     = rel_vec
)

#------------------------------------------------------------------------------
# NIMBLE MODEL CODE
#------------------------------------------------------------------------------

gs_model <- nimbleCode({
  
  #--- PRIORS ---
  
  # Survival probabilities
  S_sac1 ~ dbeta(1, 1)
  S_sac2 ~ dbeta(1, 1)
  S_sac3 ~ dbeta(1, 1)
  S_sac4 ~ dbeta(1, 1)
  S_sac5 ~ dbeta(1, 1)
  S_geo  ~ dbeta(1, 1)
  S_dcc  ~ dbeta(1, 1)
  S_ss   ~ dbeta(1, 1)
  
  # Routing probabilities - using Dirichlet-like structure
  # At occ2: fish on Sac can go to Geo, DCC, or stay Sac
  psi_geo ~ dbeta(1, 1)
  psi_dcc ~ dbeta(1, 1)
  
  # At occ3: fish on Sac can go to SS or stay Sac
  psi_ss  ~ dbeta(1, 1)
  
  # Failure probability (same at each occasion)
  phi_fail ~ dbeta(1, 1)
  
  # Survival to spawning ground
  lambda ~ dbeta(1, 1)
  
  # Detection probabilities
  p_sac1 ~ dbeta(1, 1)
  p_sac2 ~ dbeta(1, 1)
  p_sac3 ~ dbeta(1, 1)
  p_sac4 ~ dbeta(1, 1)
  p_sac5 ~ dbeta(1, 1)
  p_sac6 ~ dbeta(1, 1)
  p_geo  ~ dbeta(1, 1)
  p_dcc  ~ dbeta(1, 1)
  p_ss   ~ dbeta(1, 1)
  
  #--- TRANSITION MATRICES ---
  
  # Transition 1->2
  T[1,1,1] <- S_sac1*(1-psi_geo-psi_dcc)*(1-phi_fail)
  T[1,2,1] <- S_sac1*psi_geo*(1-phi_fail)
  T[1,3,1] <- S_sac1*psi_dcc*(1-phi_fail)
  T[1,4,1] <- 0
  T[1,5,1] <- (1-S_sac1)
  T[1,6,1] <- S_sac1*phi_fail
  T[1,7,1] <- 0
  T[2,5,1] <- 1; T[2,1,1] <- 0; T[2,2,1] <- 0; T[2,3,1] <- 0; T[2,4,1] <- 0; T[2,6,1] <- 0; T[2,7,1] <- 0
  T[3,5,1] <- 1; T[3,1,1] <- 0; T[3,2,1] <- 0; T[3,3,1] <- 0; T[3,4,1] <- 0; T[3,6,1] <- 0; T[3,7,1] <- 0
  T[4,5,1] <- 1; T[4,1,1] <- 0; T[4,2,1] <- 0; T[4,3,1] <- 0; T[4,4,1] <- 0; T[4,6,1] <- 0; T[4,7,1] <- 0
  T[5,5,1] <- 1; T[5,1,1] <- 0; T[5,2,1] <- 0; T[5,3,1] <- 0; T[5,4,1] <- 0; T[5,6,1] <- 0; T[5,7,1] <- 0
  T[6,6,1] <- 1; T[6,1,1] <- 0; T[6,2,1] <- 0; T[6,3,1] <- 0; T[6,4,1] <- 0; T[6,5,1] <- 0; T[6,7,1] <- 0
  T[7,7,1] <- 1; T[7,1,1] <- 0; T[7,2,1] <- 0; T[7,3,1] <- 0; T[7,4,1] <- 0; T[7,5,1] <- 0; T[7,6,1] <- 0
  
  # Transition 2->3
  T[1,1,2] <- S_sac2*(1-psi_ss)*(1-phi_fail)
  T[1,2,2] <- 0
  T[1,3,2] <- 0
  T[1,4,2] <- S_sac2*psi_ss*(1-phi_fail)
  T[1,5,2] <- (1-S_sac2)
  T[1,6,2] <- S_sac2*phi_fail
  T[1,7,2] <- 0
  T[2,2,2] <- 1; T[2,1,2] <- 0; T[2,3,2] <- 0; T[2,4,2] <- 0; T[2,5,2] <- 0; T[2,6,2] <- 0; T[2,7,2] <- 0
  T[3,3,2] <- 1; T[3,1,2] <- 0; T[3,2,2] <- 0; T[3,4,2] <- 0; T[3,5,2] <- 0; T[3,6,2] <- 0; T[3,7,2] <- 0
  T[4,5,2] <- 1; T[4,1,2] <- 0; T[4,2,2] <- 0; T[4,3,2] <- 0; T[4,4,2] <- 0; T[4,6,2] <- 0; T[4,7,2] <- 0
  T[5,5,2] <- 1; T[5,1,2] <- 0; T[5,2,2] <- 0; T[5,3,2] <- 0; T[5,4,2] <- 0; T[5,6,2] <- 0; T[5,7,2] <- 0
  T[6,6,2] <- 1; T[6,1,2] <- 0; T[6,2,2] <- 0; T[6,3,2] <- 0; T[6,4,2] <- 0; T[6,5,2] <- 0; T[6,7,2] <- 0
  T[7,7,2] <- 1; T[7,1,2] <- 0; T[7,2,2] <- 0; T[7,3,2] <- 0; T[7,4,2] <- 0; T[7,5,2] <- 0; T[7,6,2] <- 0
  
  # Transition 3->4
  T[1,1,3] <- S_sac3*(1-phi_fail)
  T[1,2,3] <- 0; T[1,3,3] <- 0; T[1,4,3] <- 0
  T[1,5,3] <- (1-S_sac3)
  T[1,6,3] <- S_sac3*phi_fail
  T[1,7,3] <- 0
  T[2,1,3] <- S_geo*(1-phi_fail); T[2,2,3] <- 0; T[2,3,3] <- 0; T[2,4,3] <- 0; T[2,5,3] <- (1-S_geo); T[2,6,3] <- S_geo*phi_fail; T[2,7,3] <- 0
  T[3,1,3] <- S_dcc*(1-phi_fail); T[3,2,3] <- 0; T[3,3,3] <- 0; T[3,4,3] <- 0; T[3,5,3] <- (1-S_dcc); T[3,6,3] <- S_dcc*phi_fail; T[3,7,3] <- 0
  T[4,4,3] <- (1-phi_fail); T[4,1,3] <- 0; T[4,2,3] <- 0; T[4,3,3] <- 0; T[4,5,3] <- 0; T[4,6,3] <- phi_fail; T[4,7,3] <- 0
  T[5,5,3] <- 1; T[5,1,3] <- 0; T[5,2,3] <- 0; T[5,3,3] <- 0; T[5,4,3] <- 0; T[5,6,3] <- 0; T[5,7,3] <- 0
  T[6,6,3] <- 1; T[6,1,3] <- 0; T[6,2,3] <- 0; T[6,3,3] <- 0; T[6,4,3] <- 0; T[6,5,3] <- 0; T[6,7,3] <- 0
  T[7,7,3] <- 1; T[7,1,3] <- 0; T[7,2,3] <- 0; T[7,3,3] <- 0; T[7,4,3] <- 0; T[7,5,3] <- 0; T[7,6,3] <- 0
  
  # Transition 4->5
  T[1,1,4] <- S_sac4*(1-phi_fail)
  T[1,2,4] <- 0; T[1,3,4] <- 0; T[1,4,4] <- 0
  T[1,5,4] <- (1-S_sac4)
  T[1,6,4] <- S_sac4*phi_fail
  T[1,7,4] <- 0
  T[2,5,4] <- 1; T[2,1,4] <- 0; T[2,2,4] <- 0; T[2,3,4] <- 0; T[2,4,4] <- 0; T[2,6,4] <- 0; T[2,7,4] <- 0
  T[3,5,4] <- 1; T[3,1,4] <- 0; T[3,2,4] <- 0; T[3,3,4] <- 0; T[3,4,4] <- 0; T[3,6,4] <- 0; T[3,7,4] <- 0
  T[4,1,4] <- S_ss*(1-phi_fail); T[4,2,4] <- 0; T[4,3,4] <- 0; T[4,4,4] <- 0; T[4,5,4] <- (1-S_ss); T[4,6,4] <- S_ss*phi_fail; T[4,7,4] <- 0
  T[5,5,4] <- 1; T[5,1,4] <- 0; T[5,2,4] <- 0; T[5,3,4] <- 0; T[5,4,4] <- 0; T[5,6,4] <- 0; T[5,7,4] <- 0
  T[6,6,4] <- 1; T[6,1,4] <- 0; T[6,2,4] <- 0; T[6,3,4] <- 0; T[6,4,4] <- 0; T[6,5,4] <- 0; T[6,7,4] <- 0
  T[7,7,4] <- 1; T[7,1,4] <- 0; T[7,2,4] <- 0; T[7,3,4] <- 0; T[7,4,4] <- 0; T[7,5,4] <- 0; T[7,6,4] <- 0
  
  # Transition 5->6
  T[1,1,5] <- S_sac5*(1-phi_fail)
  T[1,2,5] <- 0; T[1,3,5] <- 0; T[1,4,5] <- 0
  T[1,5,5] <- (1-S_sac5)
  T[1,6,5] <- S_sac5*phi_fail
  T[1,7,5] <- 0
  T[2,5,5] <- 1; T[2,1,5] <- 0; T[2,2,5] <- 0; T[2,3,5] <- 0; T[2,4,5] <- 0; T[2,6,5] <- 0; T[2,7,5] <- 0
  T[3,5,5] <- 1; T[3,1,5] <- 0; T[3,2,5] <- 0; T[3,3,5] <- 0; T[3,4,5] <- 0; T[3,6,5] <- 0; T[3,7,5] <- 0
  T[4,5,5] <- 1; T[4,1,5] <- 0; T[4,2,5] <- 0; T[4,3,5] <- 0; T[4,4,5] <- 0; T[4,6,5] <- 0; T[4,7,5] <- 0
  T[5,5,5] <- 1; T[5,1,5] <- 0; T[5,2,5] <- 0; T[5,3,5] <- 0; T[5,4,5] <- 0; T[5,6,5] <- 0; T[5,7,5] <- 0
  T[6,6,5] <- 1; T[6,1,5] <- 0; T[6,2,5] <- 0; T[6,3,5] <- 0; T[6,4,5] <- 0; T[6,5,5] <- 0; T[6,7,5] <- 0
  T[7,7,5] <- 1; T[7,1,5] <- 0; T[7,2,5] <- 0; T[7,3,5] <- 0; T[7,4,5] <- 0; T[7,5,5] <- 0; T[7,6,5] <- 0
  
  # Transition 6->7
  T[1,1,6] <- lambda*(1-phi_fail)
  T[1,2,6] <- 0; T[1,3,6] <- 0; T[1,4,6] <- 0
  T[1,5,6] <- (1-lambda)
  T[1,6,6] <- lambda*phi_fail
  T[1,7,6] <- 0
  T[2,5,6] <- 1; T[2,1,6] <- 0; T[2,2,6] <- 0; T[2,3,6] <- 0; T[2,4,6] <- 0; T[2,6,6] <- 0; T[2,7,6] <- 0
  T[3,5,6] <- 1; T[3,1,6] <- 0; T[3,2,6] <- 0; T[3,3,6] <- 0; T[3,4,6] <- 0; T[3,6,6] <- 0; T[3,7,6] <- 0
  T[4,5,6] <- 1; T[4,1,6] <- 0; T[4,2,6] <- 0; T[4,3,6] <- 0; T[4,4,6] <- 0; T[4,6,6] <- 0; T[4,7,6] <- 0
  T[5,5,6] <- 1; T[5,1,6] <- 0; T[5,2,6] <- 0; T[5,3,6] <- 0; T[5,4,6] <- 0; T[5,6,6] <- 0; T[5,7,6] <- 0
  T[6,6,6] <- 1; T[6,1,6] <- 0; T[6,2,6] <- 0; T[6,3,6] <- 0; T[6,4,6] <- 0; T[6,5,6] <- 0; T[6,7,6] <- 0
  T[7,7,6] <- 1; T[7,1,6] <- 0; T[7,2,6] <- 0; T[7,3,6] <- 0; T[7,4,6] <- 0; T[7,5,6] <- 0; T[7,6,6] <- 0
  
  #--- OBSERVATION MATRICES ---
  
  # Occasion 1
  O[1,1,1] <- 1;     O[1,2,1] <- 0; O[1,3,1] <- 0; O[1,4,1] <- 0; O[1,5,1] <- 0; O[1,6,1] <- 0; O[1,7,1] <- 0
  O[2,7,1] <- 1;     O[2,1,1] <- 0; O[2,2,1] <- 0; O[2,3,1] <- 0; O[2,4,1] <- 0; O[2,5,1] <- 0; O[2,6,1] <- 0
  O[3,7,1] <- 1;     O[3,1,1] <- 0; O[3,2,1] <- 0; O[3,3,1] <- 0; O[3,4,1] <- 0; O[3,5,1] <- 0; O[3,6,1] <- 0
  O[4,7,1] <- 1;     O[4,1,1] <- 0; O[4,2,1] <- 0; O[4,3,1] <- 0; O[4,4,1] <- 0; O[4,5,1] <- 0; O[4,6,1] <- 0
  O[5,5,1] <- 1;     O[5,1,1] <- 0; O[5,2,1] <- 0; O[5,3,1] <- 0; O[5,4,1] <- 0; O[5,6,1] <- 0; O[5,7,1] <- 0
  O[6,6,1] <- 1;     O[6,1,1] <- 0; O[6,2,1] <- 0; O[6,3,1] <- 0; O[6,4,1] <- 0; O[6,5,1] <- 0; O[6,7,1] <- 0
  O[7,7,1] <- 1;     O[7,1,1] <- 0; O[7,2,1] <- 0; O[7,3,1] <- 0; O[7,4,1] <- 0; O[7,5,1] <- 0; O[7,6,1] <- 0
  
  # Occasion 2
  O[1,1,2] <- p_sac2; O[1,2,2] <- 0; O[1,3,2] <- 0; O[1,4,2] <- 0; O[1,5,2] <- 0; O[1,6,2] <- 0; O[1,7,2] <- (1-p_sac2)
  O[2,2,2] <- p_geo;  O[2,1,2] <- 0; O[2,3,2] <- 0; O[2,4,2] <- 0; O[2,5,2] <- 0; O[2,6,2] <- 0; O[2,7,2] <- (1-p_geo)
  O[3,3,2] <- p_dcc;  O[3,1,2] <- 0; O[3,2,2] <- 0; O[3,4,2] <- 0; O[3,5,2] <- 0; O[3,6,2] <- 0; O[3,7,2] <- (1-p_dcc)
  O[4,7,2] <- 1;      O[4,1,2] <- 0; O[4,2,2] <- 0; O[4,3,2] <- 0; O[4,4,2] <- 0; O[4,5,2] <- 0; O[4,6,2] <- 0
  O[5,5,2] <- 1;      O[5,1,2] <- 0; O[5,2,2] <- 0; O[5,3,2] <- 0; O[5,4,2] <- 0; O[5,6,2] <- 0; O[5,7,2] <- 0
  O[6,6,2] <- 1;      O[6,1,2] <- 0; O[6,2,2] <- 0; O[6,3,2] <- 0; O[6,4,2] <- 0; O[6,5,2] <- 0; O[6,7,2] <- 0
  O[7,7,2] <- 1;      O[7,1,2] <- 0; O[7,2,2] <- 0; O[7,3,2] <- 0; O[7,4,2] <- 0; O[7,5,2] <- 0; O[7,6,2] <- 0
  
  # Occasion 3
  O[1,1,3] <- p_sac3; O[1,2,3] <- 0; O[1,3,3] <- 0; O[1,4,3] <- 0; O[1,5,3] <- 0; O[1,6,3] <- 0; O[1,7,3] <- (1-p_sac3)
  O[2,7,3] <- 1;      O[2,1,3] <- 0; O[2,2,3] <- 0; O[2,3,3] <- 0; O[2,4,3] <- 0; O[2,5,3] <- 0; O[2,6,3] <- 0
  O[3,7,3] <- 1;      O[3,1,3] <- 0; O[3,2,3] <- 0; O[3,3,3] <- 0; O[3,4,3] <- 0; O[3,5,3] <- 0; O[3,6,3] <- 0
  O[4,4,3] <- p_ss;   O[4,1,3] <- 0; O[4,2,3] <- 0; O[4,3,3] <- 0; O[4,5,3] <- 0; O[4,6,3] <- 0; O[4,7,3] <- (1-p_ss)
  O[5,5,3] <- 1;      O[5,1,3] <- 0; O[5,2,3] <- 0; O[5,3,3] <- 0; O[5,4,3] <- 0; O[5,6,3] <- 0; O[5,7,3] <- 0
  O[6,6,3] <- 1;      O[6,1,3] <- 0; O[6,2,3] <- 0; O[6,3,3] <- 0; O[6,4,3] <- 0; O[6,5,3] <- 0; O[6,7,3] <- 0
  O[7,7,3] <- 1;      O[7,1,3] <- 0; O[7,2,3] <- 0; O[7,3,3] <- 0; O[7,4,3] <- 0; O[7,5,3] <- 0; O[7,6,3] <- 0
  
  # Occasion 4
  O[1,1,4] <- p_sac4; O[1,2,4] <- 0; O[1,3,4] <- 0; O[1,4,4] <- 0; O[1,5,4] <- 0; O[1,6,4] <- 0; O[1,7,4] <- (1-p_sac4)
  O[2,7,4] <- 1;      O[2,1,4] <- 0; O[2,2,4] <- 0; O[2,3,4] <- 0; O[2,4,4] <- 0; O[2,5,4] <- 0; O[2,6,4] <- 0
  O[3,7,4] <- 1;      O[3,1,4] <- 0; O[3,2,4] <- 0; O[3,3,4] <- 0; O[3,4,4] <- 0; O[3,5,4] <- 0; O[3,6,4] <- 0
  O[4,4,4] <- p_ss;   O[4,1,4] <- 0; O[4,2,4] <- 0; O[4,3,4] <- 0; O[4,5,4] <- 0; O[4,6,4] <- 0; O[4,7,4] <- (1-p_ss)
  O[5,5,4] <- 1;      O[5,1,4] <- 0; O[5,2,4] <- 0; O[5,3,4] <- 0; O[5,4,4] <- 0; O[5,6,4] <- 0; O[5,7,4] <- 0
  O[6,6,4] <- 1;      O[6,1,4] <- 0; O[6,2,4] <- 0; O[6,3,4] <- 0; O[6,4,4] <- 0; O[6,5,4] <- 0; O[6,7,4] <- 0
  O[7,7,4] <- 1;      O[7,1,4] <- 0; O[7,2,4] <- 0; O[7,3,4] <- 0; O[7,4,4] <- 0; O[7,5,4] <- 0; O[7,6,4] <- 0
  
  # Occasion 5
  O[1,1,5] <- p_sac5; O[1,2,5] <- 0; O[1,3,5] <- 0; O[1,4,5] <- 0; O[1,5,5] <- 0; O[1,6,5] <- 0; O[1,7,5] <- (1-p_sac5)
  O[2,7,5] <- 1;      O[2,1,5] <- 0; O[2,2,5] <- 0; O[2,3,5] <- 0; O[2,4,5] <- 0; O[2,5,5] <- 0; O[2,6,5] <- 0
  O[3,7,5] <- 1;      O[3,1,5] <- 0; O[3,2,5] <- 0; O[3,3,5] <- 0; O[3,4,5] <- 0; O[3,5,5] <- 0; O[3,6,5] <- 0
  O[4,7,5] <- 1;      O[4,1,5] <- 0; O[4,2,5] <- 0; O[4,3,5] <- 0; O[4,4,5] <- 0; O[4,5,5] <- 0; O[4,6,5] <- 0
  O[5,5,5] <- 1;      O[5,1,5] <- 0; O[5,2,5] <- 0; O[5,3,5] <- 0; O[5,4,5] <- 0; O[5,6,5] <- 0; O[5,7,5] <- 0
  O[6,6,5] <- 1;      O[6,1,5] <- 0; O[6,2,5] <- 0; O[6,3,5] <- 0; O[6,4,5] <- 0; O[6,5,5] <- 0; O[6,7,5] <- 0
  O[7,7,5] <- 1;      O[7,1,5] <- 0; O[7,2,5] <- 0; O[7,3,5] <- 0; O[7,4,5] <- 0; O[7,5,5] <- 0; O[7,6,5] <- 0
  
  # Occasion 6
  O[1,1,6] <- p_sac6; O[1,2,6] <- 0; O[1,3,6] <- 0; O[1,4,6] <- 0; O[1,5,6] <- 0; O[1,6,6] <- 0; O[1,7,6] <- (1-p_sac6)
  O[2,7,6] <- 1;      O[2,1,6] <- 0; O[2,2,6] <- 0; O[2,3,6] <- 0; O[2,4,6] <- 0; O[2,5,6] <- 0; O[2,6,6] <- 0
  O[3,7,6] <- 1;      O[3,1,6] <- 0; O[3,2,6] <- 0; O[3,3,6] <- 0; O[3,4,6] <- 0; O[3,5,6] <- 0; O[3,6,6] <- 0
  O[4,7,6] <- 1;      O[4,1,6] <- 0; O[4,2,6] <- 0; O[4,3,6] <- 0; O[4,4,6] <- 0; O[4,5,6] <- 0; O[4,6,6] <- 0
  O[5,5,6] <- 1;      O[5,1,6] <- 0; O[5,2,6] <- 0; O[5,3,6] <- 0; O[5,4,6] <- 0; O[5,6,6] <- 0; O[5,7,6] <- 0
  O[6,6,6] <- 1;      O[6,1,6] <- 0; O[6,2,6] <- 0; O[6,3,6] <- 0; O[6,4,6] <- 0; O[6,5,6] <- 0; O[6,7,6] <- 0
  O[7,7,6] <- 1;      O[7,1,6] <- 0; O[7,2,6] <- 0; O[7,3,6] <- 0; O[7,4,6] <- 0; O[7,5,6] <- 0; O[7,6,6] <- 0
  
  # Occasion 7
  O[1,1,7] <- 1;      O[1,2,7] <- 0; O[1,3,7] <- 0; O[1,4,7] <- 0; O[1,5,7] <- 0; O[1,6,7] <- 0; O[1,7,7] <- 0
  O[2,7,7] <- 1;      O[2,1,7] <- 0; O[2,2,7] <- 0; O[2,3,7] <- 0; O[2,4,7] <- 0; O[2,5,7] <- 0; O[2,6,7] <- 0
  O[3,7,7] <- 1;      O[3,1,7] <- 0; O[3,2,7] <- 0; O[3,3,7] <- 0; O[3,4,7] <- 0; O[3,5,7] <- 0; O[3,6,7] <- 0
  O[4,7,7] <- 1;      O[4,1,7] <- 0; O[4,2,7] <- 0; O[4,3,7] <- 0; O[4,4,7] <- 0; O[4,5,7] <- 0; O[4,6,7] <- 0
  O[5,5,7] <- 1;      O[5,1,7] <- 0; O[5,2,7] <- 0; O[5,3,7] <- 0; O[5,4,7] <- 0; O[5,6,7] <- 0; O[5,7,7] <- 0
  O[6,6,7] <- 1;      O[6,1,7] <- 0; O[6,2,7] <- 0; O[6,3,7] <- 0; O[6,4,7] <- 0; O[6,5,7] <- 0; O[6,7,7] <- 0
  O[7,7,7] <- 1;      O[7,1,7] <- 0; O[7,2,7] <- 0; O[7,3,7] <- 0; O[7,4,7] <- 0; O[7,5,7] <- 0; O[7,6,7] <- 0
  
  #--- LIKELIHOOD ---
  for(i in 1:nfish){
    y[i, 1:nocc] ~ dDHMMo(
      init     = init[1:nstate],
      probObs  = O[1:nstate, 1:nstate, 1:nocc],
      probTrans = T[1:nstate, 1:nstate, 1:(nocc-1)],
      len      = nocc,
      checkRowSums = 0
    )
  }
})

#------------------------------------------------------------------------------
# INITIAL VALUES
#------------------------------------------------------------------------------

inits <- function(){
  list(
    S_sac1 = runif(1, 0.8, 1),
    S_sac2 = runif(1, 0.8, 1),
    S_sac3 = runif(1, 0.8, 1),
    S_sac4 = runif(1, 0.8, 1),
    S_sac5 = runif(1, 0.8, 1),
    S_geo  = runif(1, 0.8, 1),
    S_dcc  = runif(1, 0.8, 1),
    S_ss   = runif(1, 0.8, 1),
    psi_geo  = runif(1, 0.05, 0.3),
    psi_dcc  = runif(1, 0.02, 0.2),
    psi_ss   = runif(1, 0.1, 0.5),
    phi_fail = runif(1, 0.05, 0.2),
    lambda   = runif(1, 0.7, 1),
    p_sac1 = runif(1, 0.8, 1),
    p_sac2 = runif(1, 0.7, 1),
    p_sac3 = runif(1, 0.7, 1),
    p_sac4 = runif(1, 0.7, 1),
    p_sac5 = runif(1, 0.7, 1),
    p_sac6 = runif(1, 0.8, 1),
    p_geo  = runif(1, 0.5, 1),
    p_dcc  = runif(1, 0.5, 1),
    p_ss   = runif(1, 0.5, 1)
  )
}

#------------------------------------------------------------------------------
# PARAMETERS TO MONITOR
#------------------------------------------------------------------------------

params <- c("S_sac1", "S_sac2", "S_sac3", "S_sac4", "S_sac5",
            "S_geo", "S_dcc", "S_ss",
            "psi_geo", "psi_dcc", "psi_ss",
            "phi_fail", "lambda",
            "p_sac1", "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
            "p_geo", "p_dcc", "p_ss")

#------------------------------------------------------------------------------
# BUILD AND RUN MODEL
#------------------------------------------------------------------------------

gs_nimble <- nimbleModel(
  code      = gs_model,
  constants = nimble_constants,
  data      = nimble_data,
  inits     = inits()
)

# Check model is properly configured
gs_nimble$calculate()

# Configure and build MCMC
gs_mcmc_conf <- configureMCMC(gs_nimble, monitors = params)
gs_mcmc      <- buildMCMC(gs_mcmc_conf)

# Compile
gs_nimble_c <- compileNimble(gs_nimble)
gs_mcmc_c   <- compileNimble(gs_mcmc, project = gs_nimble)

# Run MCMC - start with short run to check
# Increase niter to 50000+ for final run
samples <- runMCMC(
  gs_mcmc_c,
  nchains  = 3,
  niter    = 10000,
  nburnin  = 2000,
  thin     = 5,
  inits    = list(inits(), inits(), inits()),
  samplesAsCodaMCMC = TRUE
)

# Quick diagnostics
MCMCsummary(samples, round = 3)
MCMCtrace(samples, params = c("psi_geo", "psi_dcc", "psi_ss"), pdf = FALSE)