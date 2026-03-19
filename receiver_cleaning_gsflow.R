#receiver cleaning green sturgeon flow 

library(dplyr)
library(stringr)
library(geosphere)
library(lubridate)
library(ggmap)
library(ggplot2)

# newest data labeled with routes in arc
labeled_receiver <- read.csv("C:/Users/eetracy/Desktop/ST_telemetry/arc_receivers_update.csv")

#original dataset of receivers I pulled, can also add in cleaned later
org_receiver <- read.csv("C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/obis_moorings_OTN_receivers_filtered.csv")


#new pull of receivers 12/06/2025
receiver <- read.csv("C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/full_receiver_pull_OTNmoorings_120725.csv")
unique(receiver$latitude) #1300
unique(receiver$longitude) #1425
unique(receiver$fieldnumber) #1483
unique(receiver$relatedcatalogitem)#1180
unique(org_receiver$location) #1049

setdiff(unique(receiver$relatedcatalogitem), unique(org_receiver$location))


#start trying to clean up these receivers

#previous cleaning script for large scale green sturgeon flow in gs_analysis


#######################################################################################
#ID 35 receivers that had made big movements (>81 meters)

# --- Candidate 2: median + k * MAD (robust) on non-zero moves ---
# Use MAD because distribution is heavy-tailed; use log-space if needed.
movement_threshold_m <- 81   # meters to flag "moved a lot"

receiver_summary <- receiver %>%
  # Create combined ID: actual receiver (fieldnumber) × project (relatedcatalogitem)
  mutate(
    receiver_unit = paste(fieldnumber, relatedcatalogitem, sep = "_"),
    startdatetime_parsed = mdy_hms(startdatetime, tz = "UTC")
  ) %>%
  
  group_by(receiver_unit) %>%
  arrange(startdatetime_parsed) %>%
  
  # Compute movement between deployments
  mutate(
    dist_m = distHaversine(
      cbind(longitude, latitude),
      cbind(lag(longitude), lag(latitude))
    )
  ) %>%
  
  summarise(
    fieldnumber        = first(fieldnumber),
    relatedcatalogitem = first(relatedcatalogitem),
    
    mean_latitude      = mean(latitude, na.rm = TRUE),
    mean_longitude     = mean(longitude, na.rm = TRUE),
    
    max_move_m         = max(dist_m, na.rm = TRUE),
    moved_flag         = max_move_m > movement_threshold_m,
    
    first_start = format(min(startdatetime_parsed, na.rm = TRUE), "%m/%d/%Y %H:%M:%S"),
    last_start  = format(max(startdatetime_parsed, na.rm = TRUE), "%m/%d/%Y %H:%M:%S"),
    
    .groups = "drop"
  )
receiver_summary$max_move_m[is.infinite(receiver_summary$max_move_m)] <- 0


#############################################################################################
#joining this summary dataset with original to add info
receiver <- receiver %>%
  mutate(receiver_unit = paste(fieldnumber, relatedcatalogitem, sep = "_"))

receiver_joined <- receiver %>%
  left_join(
    receiver_summary %>%
      dplyr::select(receiver_unit, mean_latitude, mean_longitude, moved_flag, max_move_m),
    by = "receiver_unit"
  )

##### Filter for CA lat and longs
receiver_CA <- receiver_joined %>%
  filter(latitude >= 32 & latitude <= 42,
         longitude <= -114 & longitude >= -124) %>%
  mutate(
    startdatetime = mdy_hms(startdatetime, tz = "UTC"),
    enddatetime   = mdy_hms(enddatetime, tz = "UTC")
  )

# Then run the ggplot script using receiver_CA

# ---- Step 2: Summarize deployments per receiver unit ----
receiver_CA_summary <- receiver_CA %>%
  arrange(receiver_unit, startdatetime) %>%
  group_by(receiver_unit) %>%
  summarise(
    first_start = min(startdatetime, na.rm = TRUE),
    last_end = max(enddatetime, na.rm = TRUE),
    total_days = as.numeric(difftime(max(enddatetime, na.rm = TRUE),
                                     min(startdatetime, na.rm = TRUE), units = "days")),
    n_deployments = n(),
    moved_any = any(moved_flag),
    .groups = "drop"
  )

# ---- Step 3: Calculate gaps between consecutive deployments ----
#prev_end is the previous enddatetime to figure out gaps between the startdatetime and end of each sequence for a unique receiver
receiver_CA_gaps <- receiver_CA %>%
  arrange(receiver_unit, startdatetime) %>%
  group_by(receiver_unit) %>%
  mutate(
    prev_end = lag(enddatetime),
    gap_days = as.numeric(difftime(startdatetime, prev_end, units = "days")),
    gap_days = ifelse(is.na(gap_days), 0, gap_days)   # Replace NA with 0
  ) %>%
  ungroup()

receiver_CA_gaps <- receiver_CA_gaps %>%
  mutate(
    startdatetime = as.POSIXct(startdatetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    enddatetime   = as.POSIXct(enddatetime,   format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    operational_days = as.numeric(difftime(enddatetime, startdatetime, units = "days"))
  ) %>%
  group_by(receiver_unit) %>%
  mutate(
    total_active_days = sum(operational_days, na.rm = TRUE)
  ) %>%
  ungroup()

write.csv (receiver_CA_gaps, "C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/full_receiver_pull_OTNmoorings_Cleaned_120825.csv")
receiver_CA_gaps <- read.csv("C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/full_receiver_pull_OTNmoorings_Cleaned_120825.csv")
#################################

################################################
#trying to get points of receivers on map showing opertion time with gaps

library(dplyr)
library(ggmap)
library(ggplot2)
library(lubridate)

# ---- Step 0: Make sure start/end datetime are POSIXct ----
receiver_operational <- receiver_operational %>%
  mutate(
    startdatetime = mdy_hms(startdatetime, tz = "UTC"),
    enddatetime   = mdy_hms(enddatetime, tz = "UTC")
  )

# ---- Step 1: Summarize receiver across all years ----
receiver_summary_all_years <- receiver_operational %>%
  group_by(receiver_unit) %>%
  summarise(
    mean_latitude = mean(mean_latitude, na.rm = TRUE),
    mean_longitude = mean(mean_longitude, na.rm = TRUE),
    max_gap_overall = max(max_gap, na.rm = TRUE),
    moved_any = any(moved_any),
    first_start = min(startdatetime, na.rm = TRUE),
    last_end = max(enddatetime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_active_days = as.numeric(difftime(last_end, first_start, units = "days")),
    gap_class = case_when(
      max_gap_overall == 0 ~ "No gap",
      max_gap_overall <= 7 ~ "Small gap (≤7d)",
      max_gap_overall > 7  ~ "Large gap (>7d)"
    )
  )

# ---- Step 2: Define bounding box for California Delta / Golden Gate to Hamilton City ----
bbox <- c(
  left = -122.6,  
  bottom = 37.75,  
  right = -121.2,  
  top = 39.5      
)

# ---- Step 3: Get base map ----
base_map <- get_stadiamap(bbox = bbox, zoom = 10, maptype = "stamen_terrain")

# ---- Step 4: Plot ----
ggmap(base_map, darken = c(0.4, "white")) +
  geom_point(
    data = receiver_summary_all_years,
    aes(
      x = mean_longitude,
      y = mean_latitude,
      color = gap_class,
      shape = moved_any,
      size = total_active_days
    ),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "No gap" = "steelblue",
      "Small gap (≤7d)" = "gold",
      "Large gap (>7d)" = "red"
    ),
    name = "Worst Gap Across Years"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("Stable", "Moved"),
    name = "Receiver Moved"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "Total Active Days"
  ) +
  labs(
    title = "Telemetry Receiver Operational Continuity (All Years)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  coord_quickmap()



###########################
#gap over all years active

# ---- Step 1: Summarize by receiver_unit ----
receiver_summary_all_years <- receiver_CA_gaps %>%
  group_by(receiver_unit) %>%
  summarise(
    mean_latitude = mean(latitude, na.rm = TRUE),
    mean_longitude = mean(longitude, na.rm = TRUE),
    max_gap = max(gap_days, na.rm = TRUE),          # biggest gap across all years
    total_active_days = sum(operational_days, na.rm = TRUE),  # total deployment days
    moved_any = any(moved_flag),                   # flag if moved at any point
    .groups = "drop"
  ) %>%
  mutate(
    gap_class = case_when(
      max_gap == 0 ~ "No gap",
      max_gap <= 7 ~ "Small gap (≤7d)",
      max_gap > 7  ~ "Large gap (>7d)"
    )
  )

# ---- Step 2: Define bounding box for CA Delta / upstream region ----
bbox <- c(
  left = -122.5,
  bottom = 37.9,
  right = -121.25,
  top = 39.0
)

# ---- Step 3: Get base map ----
base_map <- get_stadiamap(bbox = bbox, zoom = 10, maptype = "stamen_terrain")

# ---- Step 4: Generate map ----
map_all_years <- ggmap(base_map, darken = c(0.4, "white")) +
  geom_point(data = receiver_summary_all_years,
             aes(
               x = mean_longitude,
               y = mean_latitude,
               color = gap_class,
               shape = moved_any,
               size = total_active_days
             ),
             alpha = 0.7) +
  scale_color_manual(
    values = c(
      "No gap" = "steelblue",
      "Small gap (≤7d)" = "gold",
      "Large gap (>7d)" = "red"
    ),
    name = "Max Gap"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("Stable", "Moved"),
    name = "Receiver Moved"
  ) +
  scale_size_continuous(range = c(2, 8), name = "Total Active Days") +
  labs(
    title = "Telemetry Receiver Deployment Summary",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  coord_quickmap()

map_all_years

graphics.off()


############################
#receiver gaps
# Top 3 biggest gaps for each
labeled_receiver %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(receiver_group) %>%
  arrange(desc(gap_days)) %>%
  slice(1:10) %>%
  dplyr::select(receiver_group, gap_days, prev_end)

#benicia has gaps in 2014 and 2016 carquinez in 2013, 2009, 2007 so maybe we set that the fish has to pass one of these receiver groups

# Total operational time summary for all receiver groups
labeled_receiver <- labeled_receiver %>%
  mutate(
    startdatetime = format(mdy_hm(startdatetime, tz = "UTC"), "%Y-%m-%d %H:%M:%S"),
    enddatetime = format(mdy_hm(enddatetime, tz = "UTC"), "%Y-%m-%d %H:%M:%S"),
    prev_end = format(mdy_hm(prev_end, tz = "UTC"), "%Y-%m-%d %H:%M:%S")
  )
# Step 2: For each receiver_group, find periods with NO coverage (gaps > 7 days)
coverage_gaps <- labeled_receiver %>%
  filter(receiver_group %in% c("golden_gate", "bay", "carquinez", "benicia", "chipps",
                               "sacramento", "interior_delta", "mok_deltacross", "georgiana", 
                               "steamboat_sutter", "yolo_bypass", "miner", "SJ_exterior_delta", 
                               "feather", "ocean")) %>%
  arrange(receiver_group, startdatetime) %>%
  group_by(receiver_group) %>%
  mutate(
    next_start = lead(startdatetime),
    gap_to_next = as.numeric(difftime(next_start, enddatetime, units = "days"))
  ) %>%
  # Find overlaps or small gaps (coverage is maintained)
  summarise(
    earliest_coverage = min(startdatetime, na.rm = TRUE),
    latest_coverage = max(enddatetime, na.rm = TRUE),
    total_span_days = as.numeric(difftime(latest_coverage, earliest_coverage, units = "days")),
    gaps_over_7days = sum(gap_to_next > 7, na.rm = TRUE),
    max_gap_days = max(gap_to_next, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(gaps_over_7days))

coverage_gaps


library(ggplot2)
library(dplyr)

# Convert to datetime
labeled_receiver <- labeled_receiver %>%
  mutate(
    startdatetime = ymd_hms(startdatetime),
    enddatetime = ymd_hms(enddatetime)
  )
# Create timeline plot
ggplot(labeled_receiver %>%
         filter(receiver_group %in% c("golden_gate", "bay", "carquinez", "benicia", "chipps",
                                      "sacramento", "interior_delta", "mok_deltacross", "georgiana", 
                                      "steamboat_sutter", "yolo_bypass", "miner", "SJ_exterior_delta", 
                                      "feather", "ocean")),
       aes(y = receiver_group, color = receiver_group)) +
  geom_segment(aes(x = startdatetime, xend = enddatetime, 
                   yend = receiver_group), 
               size = 3, alpha = 0.7) +
  labs(title = "Receiver Coverage Timeline by Location",
       x = "Date",
       y = "Receiver Group") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10)) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")
####################################################################################

# Get unique receivers in interior_delta with their coverage times
interior_delta <- labeled_receiver %>%
  filter(receiver_group == "interior_delta") %>%
  group_by(relatedcatalogitem) %>%
  summarise(
    earliest_start = min(startdatetime, na.rm = TRUE),
    latest_end = max(enddatetime, na.rm = TRUE),
    total_operational_days = sum(as.numeric(difftime(enddatetime, startdatetime, units = "days")), na.rm = TRUE),
    number_deployments = n(),
    .groups = "drop"
  ) %>%
  arrange(earliest_start)


#relabel to split out steamboat and sutter slough as different receiver groups
labeled_receiver %>%
  filter(receiver_group == "steamboat_sutter") %>%
  distinct(relatedcatalogitem)

steamboat_items <- c("STEAMBOATSL_MOUTH1", "STMSL_MARI2", "STEAMBOATSL_MOUTH2", 
                     "STEAMBOATSL", "STMSL_MARI", "STMSL_MARI1B", 
                     "STEAMBOATSLABVCACHE_1_1_69", "STEAMBOATSLABVCACHE_1_2_69")

events <- events %>%
  mutate(receiver_group = case_when(
    receiver_group == "steamboat_sutter" & location %in% steamboat_items ~ "steamboat",
    receiver_group == "steamboat_sutter" ~ "sutter",
    TRUE ~ receiver_group
  ))

events <- events %>%
  mutate(receiver_group = case_when(
    receiver_group == "sacramento" & group == 24 ~ "spawning_ground",
    TRUE ~ receiver_group
  ))


write.csv(events, "C:/Users/eetracy/Desktop/ST_telemetry/events_with_receivergroups_021826.csv", row.names = FALSE)
