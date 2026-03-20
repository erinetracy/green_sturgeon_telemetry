#mulit state gs
library(nimble)
library(dplyr)
library(nimbleEcology)
library(abind)
library(MCMCvis)


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
  filter(status =="up_complete")

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
#looking at receiver coverage 
receiver_metadata <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/full_receiver_pull_OTNmoorings_Cleaned_120825.csv")

# Get coverage per location per year
good_coverage_yr <- receiver_metadata %>%
  mutate(year = lubridate::year(startdatetime)) %>%
  group_by(relatedcatalogitem, year) %>%
  summarise(
    max_gap = max(gap_days, na.rm = TRUE),
    good_coverage = max_gap < 7,
    .groups = "drop"
  ) %>%
  rename(location = relatedcatalogitem)

# Check
table(good_coverage_yr$good_coverage)

# For each occasion's receiver locations, what years had good coverage?
occasion_coverage <- good_coverage_yr %>%
  left_join(
    filtered_events %>% 
      distinct(location, receiver_group),
    by = "location"
  ) %>%
  filter(!is.na(receiver_group)) %>%
  group_by(receiver_group, year) %>%
  summarise(
    n_locations = n(),
    n_good = sum(good_coverage),
    pct_good = round(n_good/n_locations * 100, 1),
    all_good = all(good_coverage),
    .groups = "drop"
  ) %>%
  arrange(receiver_group, year)

# View it
occasion_coverage %>%
  print(n = Inf)

#so the way this code is working is showing how many receivers are present and a percent
# of how many had gaps of more than 7 days. Curious how to integrate into detection probability

