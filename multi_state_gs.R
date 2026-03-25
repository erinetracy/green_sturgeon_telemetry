#mulit state gs
library(nimble)
library(dplyr)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(lubridate)


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
# BUILD DETECTION HISTORY MATRIX FOR GREEN STURGEON MULTISTATE MODEL
# Model 1: 2007-2017, up_complete fish only
# States: 1=Sacramento, 2=Georgiana, 3=Delta Cross Channel, 
#         4=Steamboat/Sutter, 5=Death/not detected
# Occasions: 7 total (see below)
#==============================================================================


#i think im going to do two models 2007-2014 (without yolo) and 2015-2017 (all receiver coverage)
# Step 1: Filter to up_complete fish for model 1 years (2007-2017)
model1_events <- filtered_events %>%
  filter(status %in% c( "up_complete", "incomplete_dead", "up_incomplete"),
         water_year >= 2007, water_year <= 2017)

# Step 2: Check sacramento receiver locations and latitudes
# Need to split sacramento into occasions by location name
# Occasion 2: Decker Island receivers (below Georgiana junction)
# Occasion 4: SR_BLWSTEAM receivers (confirms Geo/DC rejoined)
# Occasion 5: SR_FREEPORT receivers (confirms Steamboat/Sutter rejoined)
# Occasion 6: SR_BLWCHIBEND + SR_BUTTEBR (confirms Yolo rejoined)
model1_events %>%
  filter(receiver_group == "sacramento") %>%
  distinct(location, mean_latitude) %>%
  arrange(mean_latitude)

# this should all receivers to maximize detection 
# Assign occasions to model1_events
# Occasion 1: benicia + carquinez (migration start)
# Occasion 2: Rio vista (state 1=Sac) OR georgiana/DCC (state 2/3)
# Occasion 3: SR_MOUTH area (state 1=Sac) OR steamboat_sutter (state 4)
# Occasion 4: SR_BLWSTEAM (state 1) - confirms Geo/DCC rejoined sac
# Occasion 5: SR_FREEPORT area (state 1) - confirms Steamboat/Sutter rejoined
# Occasion 6: SR_BLWCHIBEND + SR_BUTTEBR (state 1) - upper sac
# Occasion 7: spawning_ground (state 1) - terminal state

# Rename mok_deltacross to DCC in model1_events
# and restrict DCC detections to SR_DCC and SR_DCC2 only
model1_events <- model1_events %>%
  mutate(receiver_group = case_when(
    receiver_group == "mok_deltacross" ~ "DCC",
    TRUE ~ receiver_group
  )) %>%
  # Set DCC receiver group to NA for mok receivers that arent SR_DCC or SR_DCC2
  mutate(receiver_group = case_when(
    receiver_group == "DCC" & 
      !location %in% c("SR_DCC", "SR_DCC2") ~ NA_character_,
    TRUE ~ receiver_group
  ))


model1_events <- model1_events %>%
  mutate(occasion = case_when(
    # Occasion 1: migration start - benicia and carquinez only
    receiver_group %in% c("benicia", "carquinez") ~ 1,
    
    # Occasion 2: first decision point
    # Rio Vista = stayed sac, georgiana/mok = took alternative route
    # Decker/Horseshoe excluded - ambiguous due to side channels
    receiver_group %in% c("georgiana", "DCC") ~ 2,
    receiver_group == "sacramento" & 
      location %in% c("SR_RV10_7L", "SR_RV125L", "RIOVISTABR01", "SR_RV127L",
                      "RIOVISTABR02", "RIOVISTABR03", "SR_RV169L", "SR_RV169R") ~ 2,
    
    # Occasion 3: second decision point
    # SR_MOUTH = stayed sac, steamboat_sutter = took that route
    receiver_group == "steamboat_sutter" ~ 3,
    receiver_group == "sacramento" & 
      location %in% c("SR_MOUTH_2", "SR_MOUTH", "SR_RV150R") ~ 3,
    
    # Occasion 4: SR_BLWSTEAM area - confirms Geo/DC rejoined
    receiver_group == "sacramento" & 
      location %in% c("SR_KK240L", "SR_RYDE", "SR_BLWGEORGIA2", "SR_BLWGEORGIA",
                      "SR_KK250L", "SR_DCCSOUTH", "SR_DCCSOUTH2", "SR_KK269L",
                      "SR_DCCNORTH", "SR_BLWSTEAM2", "SR_BLWSTEAM",
                      "SR_BLWSUTTER", "SR_BLWSUTTER2", "SR_KK345R") ~ 4,
    
    # Occasion 5: SR_FREEPORT area - confirms Steamboat/Sutter rejoined
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
    
    # Occasion 6: upper sac - Yolo rejoined
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

# Check detections per occasion
model1_events %>%
  filter(!is.na(occasion)) %>%
  count(occasion, receiver_group) %>%
  arrange(occasion)


# create state column
model1_events <- model1_events %>%
  mutate(state = case_when(
    receiver_group %in% c("benicia", "carquinez", "sacramento", "spawning_ground") ~ 1,
    receiver_group == "georgiana"        ~ 2,
    receiver_group == "DCC"   ~ 3,
    receiver_group == "steamboat_sutter" ~ 4,
    TRUE ~ NA_real_
  ))


#==============================================================================
# BUILD DETECTION HISTORY MATRIX
# One row per fish per water year
# Columns = occasions 1-7
# Values = state observed (1-4) or 0 if not detected at that occasion
# State 1 = Sacramento, 2 = Georgiana, 3 = Delta Cross, 4 = Steamboat/Sutter
# State 5 = Death (incomplete_dead fish)
# State 6 = Failed migration (up_incomplete fish)
#==============================================================================

# Rebuild detection history with corrected occasions
# Using max(state) to prioritize alternative route detections over Sacramento
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

# Check
nrow(detection_history)
table(detection_history$status)

# Check state counts
ch_mat <- detection_history %>%
  select(occ_1:occ_7) %>%
  as.matrix()

apply(ch_mat, 1, function(x) any(x == 2)) %>% sum()  # Georgiana
apply(ch_mat, 1, function(x) any(x == 3)) %>% sum()  # DCC
apply(ch_mat, 1, function(x) any(x == 4)) %>% sum()  # Steamboat/Sutter
table(ch_mat)

# 0s are misdetections (model account for this), making sure there aren't too many
# What proportion of fish were detected at each occasion
colMeans(ch_mat > 0)
#    occ_1     occ_2     occ_3     occ_4     occ_5     occ_6     occ_7 
#1.0000000 0.9819820 0.8738739 0.8513514 0.9414414 1.0000000 1.0000000 


#investigating fish detections per year
# How many unique fish hit each state per water year
state_summary <- detection_history %>%
  group_by(water_year) %>%
  summarise(
    n_fish = n(),
    n_sac    = sum(apply(select(cur_data(), occ_1:occ_7), 1, function(x) any(x == 1))),
    n_geo    = sum(apply(select(cur_data(), occ_1:occ_7), 1, function(x) any(x == 2))),
    n_dcc    = sum(apply(select(cur_data(), occ_1:occ_7), 1, function(x) any(x == 3))),
    n_ss     = sum(apply(select(cur_data(), occ_1:occ_7), 1, function(x) any(x == 4))),
    n_dead   = sum(apply(select(cur_data(), occ_1:occ_7), 1, function(x) any(x == 5))),
    n_failed = sum(apply(select(cur_data(), occ_1:occ_7), 1, function(x) any(x == 6)))
  ) %>%
  arrange(water_year)

print(state_summary)

model1_events %>%
  filter(receiver_group %in% c("interior_delta", "DCC", "georgiana")) %>%
  distinct(animal_id, water_year) %>%
  count(water_year) %>%
  arrange(water_year)


#map of receivers
library(leaflet)

# Get representative coordinates for each receiver group
route_points <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez", "sacramento", 
                               "interior_delta", "georgiana", "DCC",
                               "steamboat_sutter", "spawning_ground")) %>%
  mutate(route = case_when(
    receiver_group %in% c("georgiana", "interior_delta", "DCC") ~ "Interior Delta",
    receiver_group == "steamboat_sutter" ~ "Steamboat/Sutter",
    receiver_group == "spawning_ground" ~ "Spawning Ground",
    receiver_group %in% c("benicia", "carquinez") ~ "Start",
    receiver_group == "sacramento" ~ "Sacramento"
  )) %>%
  group_by(location, route, receiver_group) %>%
  summarise(lat = mean(mean_latitude), lng = mean(mean_longitude), .groups = "drop")

# Color by route
pal <- colorFactor(
  palette = c("blue", "orange", "green", "red", "purple"),
  domain = c("Sacramento", "Interior Delta", "Steamboat/Sutter", "Start", "Spawning Ground")
)

leaflet(route_points) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~lng, lat = ~lat,
    color = ~pal(route),
    radius = 5,
    popup = ~paste(location, "<br>", route),
    label = ~route
  ) %>%
  addLegend(pal = pal, values = ~route, title = "Route")


