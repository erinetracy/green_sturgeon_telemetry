#==============================================================================
# GREEN STURGEON MULTISTATE MODEL - DATA PREPARATION
# Script: 04_multistate_data_prep.R
# Author: Erin Tracy
# Last updated: April 2026
#
# PURPOSE:
# Prepare detection history matrix for green sturgeon upstream migration
# multistate model. Filters events to model period, assigns occasions and
# states, builds detection history matrix for use in NIMBLE model.
#
# INPUTS:
#   - events_with_receivergroups_032026.csv: cleaned detection events with
#     receiver group labels (from 01_receiver_cleaning.R and 02_detection_cleaning.R)
#   - arc_receivers_update.csv: receiver metadata with group labels
#   - migration_status: object from 03_migration_status.R
#
# OUTPUTS:
#   - filtered_events: all upstream migration events (up_complete,
#     up_incomplete, incomplete_dead)
#   - model1_events: filtered to 2007-2017 with occasions and states assigned
#   - detection_history: 222 x 7 detection history matrix with fish metadata
#   - ch_mat: numeric matrix for NIMBLE (0 = not detected)
#   - ch_mat_nimble: numeric matrix for dDHMMo (0 recoded to nstate+1)
#   - fish_info: fish-level metadata (animal_id, water_year, status)
#   - gs_multistate_data.RData: all objects saved for modeling script
#
# MODEL STRUCTURE:
#   States: 1=Sacramento, 2=Georgiana, 3=DCC, 4=Steamboat/Sutter,
#           5=Death (absorbing), 6=Failed migration (absorbing)
#   Occasions: 7 (Benicia -> Rio Vista -> SR_MOUTH -> SR_BLWSTEAM ->
#              SR_KK345R/SR_FREEPORT -> Upper Sac -> Spawning Ground)
#   Years: 2007-2017 (Model 1, no Yolo Bypass)

#mulit state gs DATA PREPARATION
library(dplyr)
library(lubridate)
library(tidyr)
library(nimbleEcology)
library(abind)

#==============================================================================
# SECTION 1: LOAD DATA
#==============================================================================

events <- read.csv ("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_032026.csv")
migration_status <- read.csv ("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/migratory_status_03162026.csv")

#write.csv(filtered_events, "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_upmigration_032126.csv")
#==============================================================================
#==============================================================================
# SECTION 2: FILTER TO MODEL 1 YEARS AND ASSIGN RECEIVER GROUP LABELS
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
# Step 1: Filter events to model years and status types


# Build explicit inclusion list from migration_status
# This ensures only the correct fish x water_year x status combinations
# are included, avoiding join artifacts from the events dataset
include_list <- migration_status %>%
  filter(
    status %in% c("up_complete", "up_incomplete", "incomplete_dead"),
    water_year >= 2007,
    water_year <= 2017
  ) %>%
  dplyr::select(animal_id, water_year, status)

model1_events <- events %>%
  dplyr::select(-any_of("status")) %>%
  inner_join(include_list, by = c("animal_id", "water_year"))

# CRITICAL FIX: Remove post-spawning detections
# Use exact POSIXct timestamp not just date to catch same-day
# spawning ground + downstream detections
first_sg_date <- model1_events %>%
  filter(receiver_group == "spawning_ground") %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    first_sg_date = min(as.POSIXct(first_detection)),  # timestamp not date
    .groups = "drop"
  )

model1_events <- model1_events %>%
  left_join(first_sg_date, by = c("animal_id", "water_year")) %>%
  filter(
    is.na(first_sg_date) |
      as.POSIXct(first_detection) <= first_sg_date
  ) %>%
  dplyr::select(-first_sg_date)

# Verify
cat("Detections after fix:", nrow(model1_events), "\n")

# DECISION: Rename mok_deltacross to DCC and restrict to SR_DCC and SR_DCC2 only
# RATIONALE: The mok_deltacross receiver group contained receivers at multiple
#   locations including Mok River confluences where fish could branch onto the
#   Mok River WITHOUT entering the Delta Cross Channel. Only SR_DCC and SR_DCC2
#   are physically located inside the DCC channel itself. Detection there
#   unambiguously confirms DCC entry.
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
# Occ 4: SR_BLWSTEAM area (Sac) - Geo/DCC rejoined, SS still passing through
# Occ 5: SR_KK345R + SR_FREEPORT area (Sac) - SS has now rejoined Sacramento
# Occ 6: SR_BLWCHIBEND + SR_BUTTEBR (Sac) - upper Sacramento
# Occ 7: Spawning ground - terminal state
#
# KEY DECISIONS:
# - Decker Island EXCLUDED from occasion 2: fish can split into interior delta
#   side channels at Decker making it an ambiguous route indicator
# - Rio Vista used as first decision point (SR_RV receivers at ~38.15-38.17N)
# - STEAMBOATSL_MOUTH1/2 at occasion 3: entry point into Steamboat Slough
# - SR_BLWSTEAM/SR_BLWGEORGIA/SR_BLWSUTTER at occasion 4:
#   BLW = BELOW. These are downstream of where each slough rejoins Sacramento.
#   Geo/DCC fish rejoined upstream of here. SS fish still in channel (pass-through).
#   SS fish that TURNED AROUND would be detected here but with last(state)
#   coding they are already coded as Sac (state 1).
# - SR_KK345R moved to occasion 5: first receiver ABOVE both Steamboat AND
#   Sutter Slough rejoin points. SS fish confirmed back on Sacramento here.

model1_events <- model1_events %>%
  mutate(occasion = case_when(
    
    # Occasion 1: migration start
    receiver_group %in% c("benicia", "carquinez") ~ 1,
    
    # occ2: Rio Vista junction receivers
    # SR_RV169L and SR_RV169R removed — physically upstream of SR_MOUTH
    # despite appearing lower in lat/long due to Sacramento River bend
    receiver_group == "sacramento" &
      location %in% c("RIOVISTABR01", "RIOVISTABR02", "RIOVISTABR03",
                      "SR_RV10_7L", "SR_RV125L", "SR_RV127L") ~ 2,
    
    # occ3: SR_MOUTH junction receivers
    # SR_RV169L and SR_RV169R added here — confirmed upstream of SR_MOUTH
    # by physical inspection, river bend causes misleading lat/long position
    receiver_group == "sacramento" &
      location %in% c("SR_MOUTH", "SR_MOUTH_2", "SR_RV150R",
                      "SR_RV169L", "SR_RV169R") ~ 3,
    
    # Occasion 4: SR_BLWSTEAM area
    # Geo/DCC fish have rejoined Sacramento (exited upstream of SR_BLWSTEAM)
    # SS fish still in channel - pass-through (SR_BLWSTEAM is BELOW SS rejoin)
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
    # SR_KK345R is ABOVE both Steamboat and Sutter Slough rejoin points
    # SS fish confirmed back on Sacramento mainstem from SR_KK345R onwards
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

# BIOLOGICAL FINDING: 47.3% of fish (105/222) exhibit extended staging behavior
# at junctions, making repeated exploratory forays into alternative routes
# over periods of days to weeks before committing to a migration route.
#
# IMPORTANT: Fish detected at steamboat_sutter receivers that subsequently
# appear at SR_BLWSTEAM on the Sacramento mainstem definitively entered
# Steamboat/Sutter Slough then TURNED AROUND. SR_BLWSTEAM is BELOW where
# Steamboat Slough rejoins the Sacramento - detection there after SS detection
# means the fish reversed course. These fish are correctly coded as
# explored_ss = TRUE even though their final committed route is Sacramento.
# This exploration behavior is preserved separately for IBM analysis.

exploration_summary <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
  arrange(animal_id, water_year, first_detection) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    explored_geo     = any(state == 2),
    explored_dcc     = any(state == 3),
    explored_ss      = any(state == 4),
    explored_any     = any(state %in% c(2, 3, 4)),
    n_alt_detections = sum(state %in% c(2, 3, 4)),
    .groups = "drop"
  )

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
#   - max(state): only 1 fish had within-occasion Sac + alternative conflict
#   - first(state): made the NaN problem worse (46 NaN vs 20 with last)
#
# ABSORBING STATE ASSIGNMENT:
#   up_incomplete fish: state 6 (failed) filled forward after last detection
#   incomplete_dead fish: state 5 (dead) filled forward after last detection

detection_history <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
  arrange(animal_id, water_year, occasion, first_detection) %>%
  group_by(animal_id, water_year, occasion) %>%
  dplyr::summarise(state = last(state), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = occasion,
    values_from = state,
    names_prefix = "occ_",
    values_fill = 0
  ) %>%
  dplyr::select(animal_id, water_year,
                occ_1, occ_2, occ_3, occ_4, occ_5, occ_6, occ_7) %>%
  # Use status from model1_events directly - already filtered to correct statuses
  left_join(
    model1_events %>% 
      distinct(animal_id, water_year, status),
    by = c("animal_id", "water_year")
  ) %>%
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
  dplyr::select(-last_occ, -absorbing_state)

cat("Total fish:", nrow(detection_history), "\n")
print(table(detection_history$status))

#==============================================================================
# SECTION 7: CONVERT TO MATRIX AND FIX IMPOSSIBLE DETECTION SEQUENCES
#==============================================================================

# Store fish metadata separately
fish_info <- detection_history %>%
  dplyr::select(animal_id, water_year, status,
                explored_geo, explored_dcc, explored_ss,
                explored_any, n_alt_detections)

# Convert to numeric matrix
# 0 = not detected at that occasion
ch_mat <- detection_history %>%
  dplyr::select(occ_1:occ_7) %>%
  as.matrix()

# Recode 0 (not detected) to nstate+1 (7) as required by dDHMMo
nstate <- 6
ch_mat_nimble <- ch_mat
ch_mat_nimble[ch_mat_nimble == 0] <- nstate + 1

# Fix impossible detection sequences caused by staging behavior
# (fish visiting receivers from multiple occasions in real time)
#
# FIX 1: SS at occ3 then Sac at occ4 - fish entered SS then turned around
# Recode occ3 to Sac - exploration preserved in explored_ss = TRUE
ch_mat_nimble[ch_mat_nimble[,3] == 4 & ch_mat_nimble[,4] == 1, 3] <- 1

# FIX 2: SS at occ3 then failed (state 6) at occ4 - same turnaround logic
# Handled by allowing SS->failed transition in tr_34 (phi_fail on SS row)
# No recoding needed - handled in transition matrix

# FIX 3: Geo/DCC at occ2 then Sac at occ3 - fish entered Geo/DCC then turned around
# Recode occ2 to Sac - exploration preserved in explored_geo/dcc = TRUE
ch_mat_nimble[(ch_mat_nimble[,2] %in% c(2,3)) & ch_mat_nimble[,3] == 1, 2] <- 1

# FIX 4: Sac at occ3 then SS at occ4 - late SS entry after SR_MOUTH
# Cannot be accommodated by model structure - recode occ4 to not-detected
ch_mat_nimble[ch_mat_nimble[,3] == 1 & ch_mat_nimble[,4] == 4, 4] <- nstate + 1

#==============================================================================
# SECTION 8: VERIFY
#==============================================================================

cat("Total fish:", nrow(ch_mat_nimble), "\n")  # should be 222
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
# SECTION 9: SAVE ALL OBJECTS FOR MODELING SCRIPT
#==============================================================================
save(detection_history,
     ch_mat,
     ch_mat_nimble,
     nstate,
     fish_info,
     exploration_summary,
     model1_events,
     migration_status,
     events,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_multistate_data.RData")
cat("All objects saved\n")



















#Fixing problem that last occasion label is detecting downstream migration
# Get first spawning ground detection date for each fish x water_year
first_sg_date <- model1_events %>%
  filter(receiver_group == "spawning_ground") %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    first_sg_date = as.POSIXct(min(first_detection)),
    .groups = "drop"
  )

cat("Fish with spawning ground detections:", nrow(first_sg_date), "\n")

# Now check scope of problem
affected_fish <- model1_events %>%
  left_join(first_sg_date, by = c("animal_id", "water_year")) %>%
  filter(
    !is.na(first_sg_date),
    as.POSIXct(first_detection) > first_sg_date
  ) %>%
  distinct(animal_id, water_year, receiver_group) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    post_spawn_groups = paste(unique(receiver_group), collapse = ", "),
    .groups = "drop"
  )

cat("Fish with post-spawning detections:", nrow(affected_fish), "\n")
print(affected_fish, n = 50)

#checking all fish detection history 

# Build route string for each fish showing movement through time
# Only using detections before spawning ground arrival (already fixed)

route_summary <- model1_events %>%
  filter(!is.na(receiver_group)) %>%
  arrange(animal_id, water_year, first_detection) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    # Collapsed route - removes consecutive repeats at same receiver group
    route_string = paste(rle(receiver_group)$values, collapse = " -> "),
    # Key dates
    first_benicia = as.Date(min(first_detection[receiver_group %in% 
                                                  c("benicia", "carquinez")],
                                na.rm = TRUE)),
    first_sg = as.Date(min(first_detection[receiver_group == "spawning_ground"],
                           na.rm = TRUE)),
    last_det = as.Date(max(last_detection)),
    # Migration duration
    migration_days = as.integer(first_sg - first_benicia),
    # Was DCC or Geo detected before spawning ground?
    has_geo_upstream = any(receiver_group == "georgiana" & 
                             as.Date(first_detection) < first_sg),
    has_dcc_upstream = any(receiver_group %in% c("DCC", "mok_deltacross") & 
                             as.Date(first_detection) < first_sg),
    has_ss_upstream  = any(receiver_group == "steamboat_sutter" & 
                             as.Date(first_detection) < first_sg),
    .groups = "drop"
  ) %>%
  left_join(
    detection_history %>% dplyr::select(animal_id, water_year, 
                                        occ_1:occ_7, status),
    by = c("animal_id", "water_year")
  )

# Quick summary
cat("Route string examples:\n")
route_summary %>%
  dplyr::select(animal_id, water_year, status, route_string, 
                migration_days) %>%
  print(n = 20)

# Check for any fish where route goes backwards in time
# (downstream detections after upstream detections)
cat("\nFish with suspiciously long migrations (>60 days):\n")
route_summary %>%
  filter(migration_days > 60, status == "up_complete") %>%
  dplyr::select(animal_id, water_year, migration_days, route_string) %>%
  arrange(desc(migration_days))

# Check Georgiana fish routes
cat("\nGeorgiana fish routes:\n")
route_summary %>%
  filter(has_geo_upstream) %>%
  dplyr::select(animal_id, water_year, route_string, 
                migration_days, occ_2) %>%
  print(n = 20)

# Check SS fish routes
cat("\nSteamboat/Sutter fish routes:\n")
route_summary %>%
  filter(has_ss_upstream) %>%
  dplyr::select(animal_id, water_year, route_string,
                migration_days, occ_3) %>%
  print(n = 72)

# Check detection history state assignments match route strings
cat("\nState 2 (Geo) fish in detection history:\n")
detection_history %>%
  filter(occ_2 == 2) %>%
  left_join(route_summary %>% 
              dplyr::select(animal_id, water_year, route_string),
            by = c("animal_id", "water_year")) %>%
  dplyr::select(animal_id, water_year, occ_1:occ_4, 
                route_string) %>%
  print(n = 20)

# Print full route summary with dates
route_summary_full <- model1_events %>%
  filter(!is.na(receiver_group)) %>%
  arrange(animal_id, water_year, first_detection) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    # Route with dates - shows receiver group and first detection date
    route_with_dates = paste(
      paste0(rle(receiver_group)$values, 
             " (", format(as.Date(
               sapply(rle(receiver_group)$values, function(rg) {
                 min(first_detection[receiver_group == rg])
               })), "%m/%d/%y"), ")"),
      collapse = " -> "),
    # Clean route without dates
    route_string = paste(rle(receiver_group)$values, collapse = " -> "),
    # Key dates
    first_benicia  = as.Date(min(first_detection[receiver_group %in%
                                                   c("benicia","carquinez")], na.rm = TRUE)),
    first_rv       = as.Date(min(first_detection[receiver_group == "sacramento" &
                                                   location %in% c("SR_RV10_7L","SR_RV125L",
                                                                   "RIOVISTABR01","SR_RV127L","RIOVISTABR02",
                                                                   "RIOVISTABR03","SR_RV169L","SR_RV169R")],
                                 na.rm = TRUE)),
    first_geo      = as.Date(min(first_detection[receiver_group == "georgiana"],
                                 na.rm = TRUE)),
    first_ss       = as.Date(min(first_detection[receiver_group == "steamboat_sutter"],
                                 na.rm = TRUE)),
    first_sg       = as.Date(min(first_detection[receiver_group == "spawning_ground"],
                                 na.rm = TRUE)),
    migration_days = as.integer(first_sg - first_benicia),
    .groups = "drop"
  ) %>%
  left_join(
    detection_history %>% dplyr::select(animal_id, water_year, 
                                        occ_1:occ_7, status),
    by = c("animal_id", "water_year")
  )

# Also print full table in R - wide format so set width high
options(width = 300)
route_summary_full %>%
  dplyr::select(animal_id, water_year, status, 
                first_benicia, first_rv, first_geo, 
                first_ss, first_sg, migration_days,
                occ_1, occ_2, occ_3, route_string) %>%
  arrange(water_year, first_benicia) %>%
  print(n = Inf)
