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

