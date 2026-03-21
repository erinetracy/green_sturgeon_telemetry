#detection cleaning gsflow 01242026
#GS flow data cleaning and analysis 7/2/25
# Core data manipulation
library(data.table)
library(lubridate)
library(tidyverse)      # includes dplyr, ggplot2, tidyr, readr etc.

# Database
library(RPostgres)
library(DBI)
library(dbplyr)

# Spatial
library(mapview)
library(ggmap)
library(dismo)

# Statistics/clustering
library(doParallel)
library(factoextra)
library(cluster)

# Other
library(glatos)
library(reshape2)
library(stringdist)
library(corrplot)
library(utils)          # base R, don't need to load explicitly
#Starting with OTN code to download detections from the PATH database
#Must connect to campus VPN if downloading remotely from database and have UCD credentials
#full script in gs_data_cleaning.R


readRenviron("./.Renviron")
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  host = Sys.getenv("dbhost"),
  user = Sys.getenv("dbuser"),
  password = Sys.getenv("dbpasswd"),
  dbname = "pathnode"
)

#THIS WORKS PERFECT FOR RUNNING DETECTIONS THROUGH GLATOS
query <- "
select 
det.collectioncode AS detectedby,
anm.collectioncode AS collectioncode,
anm.catalognumber AS catalognumber,
anm.scientificname AS scientificname,
anm.commonname AS commonname,
det.datelastmodified AS datelastmodified,
det.collectioncode :: text AS detectedby,
det.collectioncode :: text AS receiver_group,
CASE WHEN det.collectornumber ::text ~~ '%(%' :: text THEN (det.station || '(' :: text) || split_part(det.collectornumber :: text, '(' :: text, 2) ELSE det.station END AS station,
det.collectornumber :: text AS receiver,
det.bottom_depth,
det.receiver_depth,
det.fieldnumber :: text AS tagname,
CASE WHEN det.fieldnumber :: text ~~ 'A%' :: text THEN (
  split_part(det.fieldnumber :: text, '-' :: text, 1) || '-' :: text
) || split_part(det.fieldnumber :: text, '-' :: text, 2) ELSE NULL :: text END AS codespace,
det.sensorname :: text AS sensorname,
det.sensorraw :: text AS sensorraw,
CASE WHEN det.sensorraw IS NULL THEN 'pinger' :: text ELSE det.sensortype :: text END AS sensortype,
det.sensorvalue :: numeric AS sensorvalue,
det.sensorunit :: text AS sensorunit,
det.datecollected,
'UTC' :: text AS timezone,
round(det.longitude :: numeric, 5) AS longitude,
round(det.latitude :: numeric, 5) AS latitude,
det.the_geom,
det.yearcollected :: text AS yearcollected,
det.monthcollected :: text AS monthcollected,
det.daycollected :: text AS daycollected,
det.julianday,
det.timeofday,
local_area,
det.notes,
citation,
(det.collectioncode :: text || '-' :: text) || det.catalognumber :: text AS unqdetecid 
from  (SELECT * FROM obis.otn_detections_early
       UNION ALL
       SELECT * FROM obis.otn_detections_2007
       UNION ALL
       SELECT * FROM obis.otn_detections_2008
       UNION ALL
       SELECT * FROM obis.otn_detections_2009
       UNION ALL
       SELECT * FROM obis.otn_detections_2010
       UNION ALL
       SELECT * FROM obis.otn_detections_2011
       UNION ALL
       SELECT * FROM obis.otn_detections_2012
       UNION ALL
       SELECT * FROM obis.otn_detections_2013
       UNION ALL
       SELECT * FROM obis.otn_detections_2014
       UNION ALL
       SELECT * FROM obis.otn_detections_2015
       UNION ALL
       SELECT * FROM obis.otn_detections_2016
       UNION ALL
       SELECT * FROM obis.otn_detections_2017
       UNION ALL
       SELECT * FROM obis.otn_detections_2018
       UNION ALL
       SELECT * FROM obis.otn_detections_2019
       UNION ALL
       SELECT * FROM obis.otn_detections_2020
       UNION ALL
       SELECT * FROM obis.otn_detections_2021
       UNION ALL
       SELECT * FROM obis.otn_detections_2022
       UNION ALL
       SELECT * FROM obis.otn_detections_2023
       UNION ALL
       SELECT * FROM obis.otn_detections_2024
       UNION ALL
       SELECT * FROM obis.otn_detections_2025
) det
left JOIN obis.otn_animals anm ON det.relatedcatalogitem = anm.catalognumber 
left join obis.otn_resources res on det.collectioncode = res.collectioncode
where det.relationshiptype = 'ANIMAL'
and LOWER(anm.commonname) = 'green sturgeon'
order by det.datecollected"

# Execute the query and fetch the results
gs_detections_original <- dbGetQuery(con, query)
write.csv(gs_detections_original, "green_sturgeon_detection_OTN_012926.csv", row.names = FALSE)
# This results in full detections
gs_detections_original <- read.csv("C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/green_sturgeon_detection_OTN_012926.csv")


#obis otn_animals metadata to join with detections/receiver data 
query <- "
SELECT *
FROM obis.otn_animals
WHERE LOWER(commonname) = 'green sturgeon'
"
otn_animals <- dbGetQuery(con, query)
write.csv(otn_animals, "green_sturgeon_animalmetadata_OTN_020726.csv", row.names = FALSE)

animal_lookup <- otn_animals %>%
  dplyr::select(catalognumber, age, lifestage, length, lengthtype, sex) %>%
  distinct(catalognumber, .keep_all = TRUE)

events <- events %>%
  left_join(
    animal_lookup,
    by = c("animal_id" = "catalognumber")
  )


#############################################################################################
#GLATOS cleaning

#Need to rename columns to run through GLATOS
colnames(gs_detections_original)[c(3, 8, 10, 13, 14, 18, 19, 20, 22, 23)] <- c("animal_id", "glatos_array", "receiver_sn", "transmitter_id", "transmitter_codespace", "sensor_value", "sensor_unit", "detection_timestamp_utc", "deploy_long", "deploy_lat")
#making sure timezone write class
class(gs_detections_original$detection_timestamp_utc)
gs_detections_original$detection_timestamp_utc <-
  as.POSIXct(gs_detections_original$detection_timestamp_utc,format = "%Y-%m-%d %H:%M:%S",
             tz = "UTC")

#False detection filter
#Write the filtered data to a new det_filtered object
#Doesn't delete rows, adds new column if detection was filtered out
detections_filtered <- false_detections(gs_detections_original, tf=3600, show_plot=TRUE)
#chapter 3 results: The filter identified 53505 (0.78%) of 6858330 detections as potentially false
#re-run 020126 The filter identified 55540 (0.8%) of 6965102 detections as potentially false

# Filter based on the column if you're happy with it.
detections_filtered <- detections_filtered[detections_filtered$passed_filter == 1,]
nrow(detections_filtered) # Check that its smaller than before

# Reduce Detections to Detection Events ####
events <- detection_events(detections_filtered,
                           location_col = 'station', # combines events across different receivers in a single array
                           time_sep=Inf)

## trying to join events with labeled_receiver so I can get the group names I assigned
#joining by lat and long didnt really work because many receivers have very similar lat and longs
receiver_lookup <- labeled_receiver %>%
  distinct(relatedcatalogitem, receiver_group)

receiver_lookup %>%
  count(relatedcatalogitem) %>%
  filter(n > 1)

events <- events %>%
  left_join(
    receiver_lookup,
    by = c("location" = "relatedcatalogitem")
  )

###summarize 
animals_per_receiver <- events %>%
  filter(!is.na(receiver_group)) %>%
  group_by(receiver_group) %>%
  summarise(
    n_unique_animals = n_distinct(animal_id),
    .groups = "drop"
  )
# there are 698 unique fish ids in this dataset
#547 are detected passing golden gate 448 in Sacramento river, 411 at benicia

write.csv(events, "C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/events_with_receivergroups_020226.csv", row.names = FALSE)
events <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/events_with_receivergroups_020226.csv")

#Make Abacus plots to look for dropped tags, deaths, and duplicates 
# Ensure the date column is converted to Date format
events_plot <- events$last_detection <- as.Date(events$last_detection)
# Define the output PDF file
output_file <- "C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/abacus_plots_location_020226.pdf"

# Open a PDF device for multi-page output
pdf(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Get unique animal IDs
unique_animals <- unique(events_plot$animal_id)

# Loop through each animal_id
for (animal in unique_animals) {
  # Filter data for the current animal
  animal_data <- events_plot %>%
    filter(animal_id == animal)
  
  # Create the abacus plot for the current animal
  p <- ggplot(animal_data, aes(x = last_detection, y = mean_latitude)) +
    geom_point(size = 3, color = "blue") +
    labs(
      title = paste("Abacus Plot for GS:", animal),
      x = "Last Detection",
      y = "Mean Latitude"
    ) +
    scale_x_date(
      date_breaks = "1 week",          # Display labels at weekly intervals
      date_labels = "%b %d '%y"       # Format as "Month Day 'Year" (e.g., Jan 01 '25)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels for readability
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  # Print the plot to the current page in the PDF
  print(p)
}
# Close the PDF device
dev.off()
# Message to indicate the PDF is ready
cat("Abacus plots have been saved to", output_file, "\n")

#Filter
# Filtering animals based on abacus plots
events$animal_id <- as.character(events$animal_id)
# List of animal_ids to exclude

# Remove spaces
events$animal_id <- trimws(events$animal_id)
any(grepl("^\\s|\\s$", events$animal_id))

#animals excluded based on dropped tag potential (didn't analyze this time)
animal_ids_to_exclude <- c()
events <- events[!(events$animal_id %in% animal_ids_to_exclude), ]

#now exclude 1 of the tags for double tagged fish
animal_ids_to_exclude <- c("UCDHIST-GS0576-2011-05-06", "UCDHIST-GS0775-2012-04-23", "UCDHIST-GS0787-2012-05-07", "UCDHIST-GS0937-2011-06-08", "UCDHIST-GS0939-2011-06-08")
events <- events[!(events$animal_id %in% animal_ids_to_exclude), ]

# Then I did a bunch of exploring the receivers and abacus plots for time to make the groups
# Receivers were grouped to ensure no time gaps and confluences were monitored
# Can be found in gs_analysis.R code the file is receiver_group_final.csv

###############################################################################################################
#matrix presence/absence cluster

#adding in receiver group number column
#group <- read.csv("C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/receiver_group_final.csv")

#joined_data <- events %>%
#  inner_join(group %>% dplyr::select(location, group, area), by = "location")

events <- events %>%
  mutate(year = year(first_detection), 
         month =month(first_detection))
events <- events %>%
  mutate(water_year = ifelse(month > 9, year + 1, year))

########################################################################################################
receiver_metadata <- read.csv("")

##################
#relabel receiver groups from arcgis metadata
# Fix conflicts in correct_groups before relabeling events
# Create correct location to receiver_group lookup from receiver_metadata
# Step 1: Build correct_groups from receiver_metadata
correct_groups <- receiver_metadata %>%
  distinct(relatedcatalogitem, receiver_group) %>%
  filter(!is.na(relatedcatalogitem), relatedcatalogitem != "",
         !is.na(receiver_group), receiver_group != "") %>%
  rename(location = relatedcatalogitem,
         receiver_group_correct = receiver_group) %>%
  # Fix conflicts
  mutate(receiver_group_correct = case_when(
    location == "SR_FEATHER2_RT" ~ "sacramento",
    location == "RICHBR_22_2015" ~ "bay",
    TRUE ~ receiver_group_correct
  )) %>%
  # Remove duplicates
  distinct(location, .keep_all = TRUE)

events <- events %>%
  select(-receiver_group) %>%
  left_join(correct_groups %>% rename(receiver_group = receiver_group_correct), 
            by = "location") %>%
  mutate(receiver_group = ifelse(group == 24, "spawning_ground", receiver_group))

# Verify
unique(events$receiver_group)

###################################################
#sort into up and down migrants 
#Battalie 2024 used Approximate Coordinates (RM 200 / 322 km) which is approximately hamilton city 
# slightly south at lat -122 and long 39.7
#do I still use this as a downstream migration point?
# our lowest point for group 24 is -121.9745, 39.73131

events <- events %>%
  mutate(receiver_group = case_when(
    group == 24 ~ "spawning_ground",
    TRUE ~ receiver_group
  ))

unique(events$receiver_group)

#######################################################################################################
#summary stats
# Step 1: ID fish that passed benicia OR carquinez, excluding juveniles/sub-adults
migration_status <- events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  distinct(animal_id, water_year) %>%
  left_join(
    events %>%
      mutate(first_detection = as.POSIXct(first_detection)) %>%
      arrange(animal_id, water_year, first_detection) %>%
      group_by(animal_id, water_year) %>%
      summarise(
        lifestage = first(lifestage),
        # Collapsed route removing consecutive repeats at same receiver group
        route_groups = paste(rle(receiver_group)$values, collapse = " -> "),
        # First and last receiver group visited
        first_group = first(receiver_group),
        last_group = last(receiver_group),
        # Earliest detection timestamp at each key receiver group
        gg  = suppressWarnings(min(first_detection[receiver_group == "golden_gate"], na.rm = TRUE)),
        bay = suppressWarnings(min(first_detection[receiver_group == "bay"], na.rm = TRUE)),
        bc  = suppressWarnings(min(first_detection[receiver_group %in% c("benicia", "carquinez")], na.rm = TRUE)),
        sac = suppressWarnings(min(first_detection[receiver_group == "sacramento"], na.rm = TRUE)),
        sg  = suppressWarnings(min(first_detection[receiver_group == "spawning_ground"], na.rm = TRUE)),
        .groups = "drop"
      ),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    # Boolean flags for whether fish was detected at each key receiver group
    has_gg  = !is.na(gg) & !is.infinite(gg),
    has_bay = !is.na(bay) & !is.infinite(bay),
    has_bc  = !is.na(bc) & !is.infinite(bc),
    has_sac = !is.na(sac) & !is.infinite(sac),
    has_sg  = !is.na(sg) & !is.infinite(sg),
    # Ocean side receiver (gg or bay) was detected before bc - indicates fish came from ocean
    ocean_before_bc = (has_gg & gg < bc) | (has_bay & bay < bc),
    # Spawning ground detected before bc - indicates fish is moving downstream
    sg_before_bc = has_sg & sg < bc,
    status = case_when(
      # DOWNSTREAM COMPLETE: fish detected at sg then bc then gg - full downstream migration
      sg_before_bc & has_bc & has_gg & bc < gg          ~ "down_complete",
      # DOWNSTREAM COMPLETE: fish came from sac direction and made it to gg
      # counts as down_complete since fish may have spawned downstream of spawning_ground receivers
      has_sac & sac < bc & has_gg & bc < gg             ~ "down_complete",
      # DOWNSTREAM INCOMPLETE: fish detected at sg then bc but never reached gg
      # fish started downstream migration but did not make it to ocean
      sg_before_bc & has_bc & !has_gg                   ~ "down_incomplete",
      # UPSTREAM COMPLETE: fish came from ocean side (gg or bay) then bc then sg
      # full upstream spawning migration confirmed
      ocean_before_bc & has_sg & bc < sg                ~ "up_complete",
      # UPSTREAM INCOMPLETE: fish came from ocean side (gg or bay) then bc
      # but never reached spawning ground - failed or incomplete migration
      ocean_before_bc & !has_sg                         ~ "up_incomplete",
      # BAD: everything else - ambiguous movement, delta resident, or data issues
      TRUE                                              ~ "bad"
    )
  )

migration_status <- migration_status %>%
  mutate(status = case_when(
    animal_id == "CDFWA15-1219838-2016-09-08" & water_year == 2017 ~ "bad",
    animal_id == "UCDHIST-GS0823-2012-07-03" & water_year == 2021 ~ "bad",
    TRUE ~ status
  ))
dead_fish <- c("UCDHIST-GS0276-2005-08-20", "UCDHIST-GS0488-2011-08-10",
               "UCDHIST-GS0512-2012-04-20", "UCDHIST-GS0516-2012-04-26",
               "UCDHIST-GS0634-2011-07-08", "UCDHIST-GS0806-2012-07-01",
               "UCDHIST-GS0814-2012-07-01", "UCDHIST-GS0821-2012-07-02",
               "UCDHIST-GS0823-2012-07-03", "CDFWA15-1306970-2018-12-18")
migration_status %>% count(status)

# Remove status column if it already exists before joining
events <- events %>%
  select(-any_of("status")) %>%
  left_join(migration_status %>% 
              dplyr::select(animal_id, water_year, status),
            by = c("animal_id", "water_year"))

# Remove the duplicate column and rename the one you want to keep
events <- events %>%
  dplyr::select(-status.y) %>%
  rename(status = status.x)

events <- events %>% dplyr::select(-event)

# updated receiver_groups that were not labled 
write.csv (events, "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_032026.csv", row.names = FALSE)

events <- read.csv ("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_032026.csv")
head(events)

