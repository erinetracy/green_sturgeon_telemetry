#mulit state gs DATA PREPARATION
library(nimble)
library(dplyr)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(lubridate)
library(abind)


#starting with figuring out receiver group locations to match up with states 
events %>%
  filter(status == "up_complete") %>%
  group_by(receiver_group) %>%
  summarise(
    mean_lat = mean(mean_latitude, na.rm = TRUE),
    min_lat = min(mean_latitude, na.rm = TRUE),
    max_lat = max(mean_latitude, na.rm = TRUE),
    n_locations = n_distinct(location)
  ) %>%
  arrange(mean_lat)

filtered_events <- events %>%
  inner_join(
    migration_status %>% filter(status %in% c("up_complete","up_incomplete")) %>%
      select(animal_id, water_year, status),
    by = c("animal_id","water_year"))

filtered_events %>%
  filter(status == "up_complete", receiver_group == "sacramento") %>%
  distinct(location, mean_latitude, mean_longitude) %>%
  arrange(mean_latitude) %>%
  as.data.frame()


#try to id receiver points before confluences for route transitions
library(leaflet)

filtered_events %>%
  filter(status == "up_complete", receiver_group == "sacramento") %>%
  distinct(location, mean_latitude, mean_longitude) %>%
  leaflet() %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~mean_longitude,
    lat = ~mean_latitude,
    popup = ~paste(location, "<br>", 
                   "Lat:", mean_latitude, "<br>", 
                   "Lon:", mean_longitude),
    radius = 5,
    color = "blue"
  )

#separating sac receiver group into multiple groups 
filtered_events <- filtered_events %>%
  mutate(occasion = case_when(
    receiver_group %in% c("benicia", "carquinez") ~ 1,
    receiver_group == "sacramento" & 
      location %in% c("DECKER_ISSE", "DECKER_ISCL3", "DECKER_ISCR2", "DECKER_ISSW") ~ 2,
    receiver_group %in% c("georgiana", "mok_deltacross") ~ 2,
    receiver_group == "sacramento" & 
      location %in% c("SR_MOUTH", "SR_MOUTH_2", "SR_RV150R") ~ 3,
    receiver_group %in% c("steamboat_sutter", "yolo_bypass") ~ 3,
    receiver_group == "sacramento" & 
      location %in% c("SR_BLWSTEAM", "SR_BLWSTEAM2") ~ 4,
    receiver_group == "sacramento" & 
      location %in% c("SR_FREEPORT", "SR_FREEPORT_1", "SR_GB470L", 
                      "SR_FREEPORTDIV_W", "SR_FREEPORTDIV_S", "SR_FREEPORTDIV_N") ~ 5,
    receiver_group == "sacramento" & 
      location %in% c("SR_BLWCHIBEND_W2", "SR_BLWCHIBEND_E", "SR_BLWCHIBEND_W",) ~ 6,
    receiver_group == "sacramento" & 
      location %in% c("SR_MERIDIANBR", "SR_MERIDIANBR2", "SR_BUTTEBR",
                      "BUTTEBR2-V", "BUTTEBR3-V", "SR_BUTTEBR") ~ 7,
    receiver_group == "spawning_ground" ~ 7,
    TRUE ~ NA_real_
  ))
##############################################################################################
#looking at receiver coverage and fish detections

#this one is old
#receiver_metadata <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/full_receiver_pull_OTNmoorings_Cleaned_120825.csv")
#this one is current
receiver_metadata <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/arc_receivers_update.csv")
#i think events has outdated receiver_group labels so i updated from receiver_metadata

#how many fish do we have complete and incomplete 
filtered_events <- events %>%
  inner_join(
    migration_status %>% filter(status %in% c("up_complete","up_incomplete", "incomplete_dead")) %>%
      select(animal_id, water_year, status),
    by = c("animal_id","water_year"))

filtered_events <- filtered_events %>%
  select(-status.y) %>%
  rename(status = status.x)

filtered_events %>%
  filter (status =="up_complete") %>%
  group_by(receiver_group) %>%
  summarise(n_fish = n_distinct(animal_id)) %>%
  arrange(desc(n_fish))

write.csv(filtered_events, "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_upmigration_032126.csv")

#this shows the NA and blank detections are outside study area and dont need consideration
events %>%
  filter(receiver_group == "" | is.na(receiver_group)) %>%
  distinct(location, receiver_group) 

#lookin at detection per water_year
events %>%
  filter(receiver_group %in% c("benicia", "carquinez", "sacramento", 
                               "georgiana", "mok_deltacross", "steamboat_sutter",
                               "yolo_bypass", "spawning_ground")) %>%
  group_by(water_year, receiver_group) %>%
  summarise(n_fish = n_distinct(animal_id), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = receiver_group, values_from = n_fish, values_fill = 0) %>%
  arrange(water_year)

#receiver coverage by water year
receiver_metadata %>%
  mutate(year = lubridate::year(lubridate::mdy_hm(startdatetime))) %>%
  filter(receiver_group %in% c("benicia", "carquinez", "sacramento",
                               "georgiana", "mok_deltacross", "steamboat_sutter",
                               "yolo_bypass", "spawning_ground")) %>%
  mutate(receiver_group = case_when(
    receiver_group == "benicia"          ~ "B",
    receiver_group == "carquinez"        ~ "C",
    receiver_group == "sacramento"       ~ "Sac",
    receiver_group == "georgiana"        ~ "G",
    receiver_group == "mok_deltacross"   ~ "DCC",
    receiver_group == "steamboat_sutter" ~ "SS",
    receiver_group == "yolo_bypass"      ~ "Y",
    receiver_group == "spawning_ground"  ~ "SG"
  )) %>%
  group_by(receiver_group, year) %>%
  summarise(
    n = n(),
    good = sum(gap_days < 7, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = receiver_group,
    values_from = c(n, good),
    values_fill = 0
  ) %>%
  arrange(year) %>%
  print(n=50)

#==============================================================================
#==============================================================================
# GREEN STURGEON MULTISTATE MODEL - DATA PREPARATION
# Model 1: 2007-2017
# States: 1=Sacramento, 2=Georgiana, 3=DCC, 4=Steamboat/Sutter
#         5=Death, 6=Failed migration
# Occasions: 7 (Benicia -> Rio Vista -> SR_MOUTH -> SR_BLWSTEAM -> 
#            SR_FREEPORT -> upper Sac -> Spawning ground)
#==============================================================================

# Step 1: Filter events to model years and status types
model1_events <- filtered_events %>%
  filter(status %in% c("up_complete", "incomplete_dead", "up_incomplete"),
         water_year >= 2007, water_year <= 2017)

# Step 2: Rename mok_deltacross to DCC, restrict to SR_DCC and SR_DCC2 only
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

# Step 3: Assign occasions
model1_events <- model1_events %>%
  mutate(occasion = case_when(
    receiver_group %in% c("benicia", "carquinez") ~ 1,
    receiver_group %in% c("georgiana", "DCC") ~ 2,
    receiver_group == "sacramento" & 
      location %in% c("SR_RV10_7L", "SR_RV125L", "RIOVISTABR01", "SR_RV127L",
                      "RIOVISTABR02", "RIOVISTABR03", "SR_RV169L", "SR_RV169R") ~ 2,
    receiver_group == "steamboat_sutter" ~ 3,
    receiver_group == "sacramento" & 
      location %in% c("SR_MOUTH_2", "SR_MOUTH", "SR_RV150R") ~ 3,
    receiver_group == "sacramento" & 
      location %in% c("SR_KK240L", "SR_RYDE", "SR_BLWGEORGIA2", "SR_BLWGEORGIA",
                      "SR_KK250L", "SR_DCCSOUTH", "SR_DCCSOUTH2", "SR_KK269L",
                      "SR_DCCNORTH", "SR_BLWSTEAM2", "SR_BLWSTEAM",
                      "SR_BLWSUTTER", "SR_BLWSUTTER2", "SR_KK345R") ~ 4,
    receiver_group == "sacramento" & 
      location %in% c("SR_GB437R", "SR_GB447R", "SR_FREEPORT", "SR_FREEPORT_1",
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
    receiver_group == "spawning_ground" ~ 7,
    TRUE ~ NA_real_
  ))

# Step 4: Assign states
model1_events <- model1_events %>%
  mutate(state = case_when(
    receiver_group %in% c("benicia", "carquinez", "sacramento", "spawning_ground") ~ 1,
    receiver_group == "georgiana"        ~ 2,
    receiver_group == "DCC"              ~ 3,
    receiver_group == "steamboat_sutter" ~ 4,
    TRUE ~ NA_real_
  ))

# Step 5: Build detection history matrix
# max(state) prioritizes alternative route detections over Sacramento
detection_history <- model1_events %>%
  filter(!is.na(occasion), !is.na(state)) %>%
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

# Step 6: Convert to numeric matrix for NIMBLE
fish_info <- detection_history %>%
  select(animal_id, water_year, status)

ch_mat <- detection_history %>%
  select(occ_1:occ_7) %>%
  as.matrix()

# Step 7: Verify
nrow(ch_mat)                    # should be 222
table(detection_history$status) # check status counts
colMeans(ch_mat > 0)            # detection rates per occasion
apply(ch_mat, 1, function(x) any(x == 2)) %>% sum()  # Georgiana fish
apply(ch_mat, 1, function(x) any(x == 3)) %>% sum()  # DCC fish
apply(ch_mat, 1, function(x) any(x == 4)) %>% sum()  # Steamboat/Sutter fish

########################################################

# Actually build Multi-state model
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