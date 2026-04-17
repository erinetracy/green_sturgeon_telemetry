#flow covariate creation for gs multistate model

library(nimble)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(dplyr)
library(lubridate)
library(tidyr)

#tidally filtered rio vista gauge data 
rv_flow <- read.csv("C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv")
rv_temp <- read.csv("C:/Users/eetracy/Desktop/flow_data/Rio_Vista_daily_temp.csv")

# Clean flow data - keep only what we need
rv_flow_clean <- rv_flow %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_cfs = value) %>%
  filter(!is.na(flow_cfs),
         date >= as.Date("2006-10-01"),  # start of WY2007
         date <= as.Date("2017-09-30")) %>%  # end of WY2017
  arrange(date)

#not sure if we want to do this yet

# Standardize flow (mean = 0, sd = 1) following Perry
# Standardizing improves MCMC mixing and makes beta coefficients interpretable
flow_mean <- mean(rv_flow_clean$flow_cfs)
flow_sd   <- sd(rv_flow_clean$flow_cfs)

rv_flow_clean <- rv_flow_clean %>%
  mutate(flow_std = (flow_cfs - flow_mean) / flow_sd)

# Create study date index - integer day number from start of study
study_start <- as.Date("2006-10-01")
rv_flow_clean <- rv_flow_clean %>%
  mutate(study_day = as.integer(date - study_start) + 1)

cat("Flow data days:", nrow(rv_flow_clean), "\n")
cat("Any missing:", sum(is.na(rv_flow_clean$flow_std)), "\n")

# Get first detection date at occ2 and occ3 for each fish
# These are the dates the fish arrived at each decision point
junction_dates <- model1_events %>%
  filter(occasion %in% c(2, 3),
         !is.na(state)) %>%
  arrange(animal_id, water_year, occasion, first_detection) %>%
  group_by(animal_id, water_year, occasion) %>%
  dplyr::summarise(
    arrival_date = as.Date(min(first_detection)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = occasion,
    values_from = arrival_date,
    names_prefix = "arrival_occ"
  )

head(junction_dates)

# Join flow to each fish at their junction arrival dates
fish_flow <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(junction_dates, by = c("animal_id", "water_year")) %>%
  left_join(
    rv_flow_clean %>% dplyr::select(date, flow_occ2 = flow_std),
    by = c("arrival_occ2" = "date")
  ) %>%
  left_join(
    rv_flow_clean %>% dplyr::select(date, flow_occ3 = flow_std),
    by = c("arrival_occ3" = "date")
  )

# Check
cat("Missing flow_occ2:", sum(is.na(fish_flow$flow_occ2)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_flow$flow_occ3)), "\n")

# For fish with no occ2 detection (some failed fish) impute mean flow (0 after standardizing)
fish_flow <- fish_flow %>%
  mutate(
    flow_occ2 = ifelse(is.na(flow_occ2), 0, flow_occ2),
    flow_occ3 = ifelse(is.na(flow_occ3), 0, flow_occ3)
  )

# Quick plot to see flow distribution
par(mfrow = c(1,2))
hist(fish_flow$flow_occ2, main = "Flow at occ2 (Geo/DCC junction)", 
     xlab = "Standardized flow")
hist(fish_flow$flow_occ3, main = "Flow at occ3 (SS junction)", 
     xlab = "Standardized flow")

#alot of NAs
# Check what's actually driving the NAs first
fish_flow %>%
  mutate(has_occ2 = !is.na(arrival_occ2),
         has_occ3 = !is.na(arrival_occ3)) %>%
  left_join(detection_history %>% 
              dplyr::select(animal_id, water_year, status),
            by = c("animal_id", "water_year")) %>%
  count(status, has_occ2, has_occ3) %>%
  arrange(status)

# Don't impute - keep NAs and pass them to NIMBLE
# NIMBLE only evaluates routing probabilities for fish
# that transition through state 1 at each occasion
# Fish with missing flow will have their routing probability
# evaluated at 0 (mean flow) which is acceptable since
# the likelihood contribution from those fish doesn't
# depend on the routing covariate

# Instead of imputing 0, keep NAs and let NIMBLE handle it
# by passing a separate indicator for which fish have flow data
fish_flow <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(junction_dates, by = c("animal_id", "water_year")) %>%
  left_join(
    rv_flow_clean %>% dplyr::select(date, flow_occ2 = flow_std),
    by = c("arrival_occ2" = "date")
  ) %>%
  left_join(
    rv_flow_clean %>% dplyr::select(date, flow_occ3 = flow_std),
    by = c("arrival_occ3" = "date")
  ) %>%
  mutate(
    # Replace NA with 0 only as a computational placeholder
    # These fish don't influence routing parameter estimates
    # because they are not in state 1 at the relevant transition
    flow_occ2 = ifelse(is.na(flow_occ2), 0, flow_occ2),
    flow_occ3 = ifelse(is.na(flow_occ3), 0, flow_occ3)
  )

# Verify order matches detection_history
identical(fish_flow$animal_id, detection_history$animal_id)
identical(fish_flow$water_year, detection_history$water_year)

# Extract as vectors for NIMBLE
flow_occ2 <- fish_flow$flow_occ2
flow_occ3 <- fish_flow$flow_occ3

# Check distribution of non-imputed values only
par(mfrow = c(1,2))
hist(flow_occ2[flow_occ2 != 0 | !is.na(fish_flow$arrival_occ2)],
     main = "Flow at occ2 junction (non-imputed)",
     xlab = "Standardized flow (cfs)")
hist(flow_occ3[flow_occ3 != 0 | !is.na(fish_flow$arrival_occ3)],
     main = "Flow at occ3 junction (non-imputed)", 
     xlab = "Standardized flow (cfs)")

# Plot only fish that actually have real flow values
par(mfrow = c(1,2))

real_flow_occ2 <- fish_flow$flow_occ2[!is.na(fish_flow$arrival_occ2)]
real_flow_occ3 <- fish_flow$flow_occ3[!is.na(fish_flow$arrival_occ3)]

hist(real_flow_occ2,
     main = "Flow at occ2 (n=156 fish)",
     xlab = "Standardized flow (cfs)")

hist(real_flow_occ3,
     main = "Flow at occ3 (n=106 fish)",
     xlab = "Standardized flow (cfs)")

# Summary stats
cat("occ2 flow - mean:", round(mean(real_flow_occ2), 2),
    "sd:", round(sd(real_flow_occ2), 2),
    "range:", round(range(real_flow_occ2), 2), "\n")
cat("occ3 flow - mean:", round(mean(real_flow_occ3), 2),
    "sd:", round(sd(real_flow_occ3), 2),
    "range:", round(range(real_flow_occ3), 2), "\n")




















###########################################################################################################
#exploring flow data
#this may be for secondary analysis of migration start time influenced by flow

library(ggplot2)
library(dplyr)
library(lubridate)

# Get date each fish passed Benicia/Carquinez (migration start)
entry_dates <- model1_events %>%
  filter(occasion == 1) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    entry_date = as.Date(min(first_detection)),
    .groups = "drop"
  )

# Get date each fish passed Rio Vista (occ2 - first decision point)
rv_dates <- model1_events %>%
  filter(occasion == 2, receiver_group == "sacramento") %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    rv_date = as.Date(min(first_detection)),
    .groups = "drop"
  )

# Join entry and Rio Vista dates
fish_dates <- entry_dates %>%
  left_join(rv_dates, by = c("animal_id", "water_year")) %>%
  left_join(
    detection_history %>% dplyr::select(animal_id, water_year, status,
                                        occ_2, occ_3),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    # Final route taken
    route = case_when(
      occ_2 == 2 ~ "Georgiana",
      occ_2 == 3 ~ "DCC",
      occ_3 == 4 ~ "Steamboat/Sutter",
      TRUE ~ "Sacramento"
    )
  )

head(fish_dates)

# For each fish create a window of flow from 30 days before
# to 30 days after Rio Vista passage
# This centers everything on the routing decision

library(tidyr)

flow_windows <- fish_dates %>%
  filter(!is.na(rv_date)) %>%
  rowwise() %>%
  mutate(
    flow_window = list(
      rv_flow_clean %>%
        filter(date >= rv_date - 30,
               date <= rv_date + 30) %>%
        mutate(
          days_relative = as.integer(date - rv_date)
        )
    )
  ) %>%
  tidyr::unnest(flow_window) %>%
  ungroup()


# Plot 1: Hydrograph centered on Rio Vista passage date
# Each line is one fish, colored by route taken
ggplot(flow_windows, aes(x = days_relative, y = flow_cfs,
                         group = paste(animal_id, water_year),
                         color = route)) +
  geom_line(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -7, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = -3, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = c("Sacramento" = "steelblue",
                                "Georgiana" = "orange",
                                "Steamboat/Sutter" = "forestgreen",
                                "DCC" = "red")) +
  facet_wrap(~water_year, scales = "free_y") +
  labs(title = "Flow at Rio Vista relative to fish passage (day 0)",
       x = "Days relative to Rio Vista passage",
       y = "Flow (cfs)",
       color = "Route taken") +
  theme_minimal()


# Calculate flow metrics for each fish
fish_flow_metrics <- fish_dates %>%
  filter(!is.na(rv_date)) %>%
  rowwise() %>%
  mutate(
    # Flow on day of Rio Vista passage
    flow_day0 = {
      idx <- which(rv_flow_clean$date == rv_date)
      if(length(idx) == 0) NA_real_ else rv_flow_clean$flow_cfs[idx]
    },
    # Mean flow 3 days prior
    flow_3day = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= rv_date - 3 &
        rv_flow_clean$date < rv_date], na.rm = TRUE),
    # Mean flow 7 days prior
    flow_7day = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= rv_date - 7 &
        rv_flow_clean$date < rv_date], na.rm = TRUE),
    # Mean flow 30 days prior
    flow_30day = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= rv_date - 30 &
        rv_flow_clean$date < rv_date], na.rm = TRUE)
  ) %>%
  ungroup()

# Check how many NAs in flow_day0
cat("Missing flow_day0:", sum(is.na(fish_flow_metrics$flow_day0)), "\n")

# Compare flow metrics by route
fish_flow_metrics %>%
  tidyr::pivot_longer(
    cols = c(flow_day0, flow_3day, flow_7day, flow_30day),
    names_to = "metric",
    values_to = "flow"
  ) %>%
  ggplot(aes(x = route, y = flow, fill = route)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y",
             labeller = labeller(metric = c(
               flow_day0 = "Day of passage",
               flow_3day = "3-day mean prior",
               flow_7day = "7-day mean prior",
               flow_30day = "30-day mean prior"
             ))) +
  scale_fill_manual(values = c("Sacramento" = "steelblue",
                               "Georgiana" = "orange",
                               "Steamboat/Sutter" = "forestgreen",
                               "DCC" = "red")) +
  labs(title = "Flow at different temporal scales by route taken",
       x = "Route", y = "Flow (cfs)") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Show the full migration season hydrograph per water year
# with vertical lines showing when fish passed Rio Vista
# colored by route

ggplot() +
  # Full season hydrograph
  geom_line(data = rv_flow_clean %>%
              mutate(water_year = ifelse(month(date) >= 10,
                                         year(date) + 1, year(date))) %>%
              filter(water_year %in% 2007:2017,
                     month(date) %in% c(1,2,3,4,5,6,7)),
            aes(x = date, y = flow_cfs), color = "gray60", linewidth = 0.5) +
  # Fish passage events
  geom_rug(data = fish_dates %>% filter(!is.na(rv_date)),
           aes(x = rv_date, color = route), sides = "b", alpha = 0.7) +
  scale_color_manual(values = c("Sacramento" = "steelblue",
                                "Georgiana" = "orange",
                                "Steamboat/Sutter" = "forestgreen",
                                "DCC" = "red")) +
  facet_wrap(~water_year, scales = "free") +
  labs(title = "Migration season hydrograph with fish passage timing",
       x = "Date", y = "Flow (cfs)", color = "Route taken") +
  theme_minimal()


# Bin fish by flow quantile and look at proportion 
# taking each route
fish_flow_metrics %>%
  filter(!is.na(flow_day0)) %>%
  mutate(
    flow_quantile = ntile(flow_day0, 4),
    flow_label = case_when(
      flow_quantile == 1 ~ "Q1 Low flow",
      flow_quantile == 2 ~ "Q2 Med-low",
      flow_quantile == 3 ~ "Q3 Med-high",
      flow_quantile == 4 ~ "Q4 High flow"
    )
  ) %>%
  count(flow_label, route) %>%
  group_by(flow_label) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = flow_label, y = prop, fill = route)) +
  geom_col() +
  scale_fill_manual(values = c("Sacramento" = "steelblue",
                               "Georgiana" = "orange",
                               "Steamboat/Sutter" = "forestgreen",
                               "DCC" = "red")) +
  labs(title = "Proportion of fish taking each route by flow quartile",
       x = "Flow quartile at day of Rio Vista passage",
       y = "Proportion of fish",
       fill = "Route") +
  theme_minimal()


#################################################################
#looking at start of migration and flow/temp
# Get Benicia passage date for each fish (migration start)
benicia_dates <- model1_events %>%
  filter(occasion == 1) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    benicia_date = as.Date(min(first_detection)),
    .groups = "drop"
  ) %>%
  left_join(
    detection_history %>% dplyr::select(animal_id, water_year, status),
    by = c("animal_id", "water_year")
  )

# Clean temp data
rv_temp_clean <- rv_temp %>%
  mutate(date = as.Date(date)) %>%
  dplyr::select(date, daily_temp) %>%
  filter(!is.na(daily_temp),
         date >= as.Date("2006-10-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date)

# Calculate flow and temp metrics at migration start
migration_metrics <- benicia_dates %>%
  rowwise() %>%
  mutate(
    # Flow metrics
    flow_day0 = {
      idx <- which(rv_flow_clean$date == benicia_date)
      if(length(idx) == 0) NA_real_ else rv_flow_clean$flow_cfs[idx]
    },
    flow_3day = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= benicia_date - 3 &
        rv_flow_clean$date < benicia_date], na.rm = TRUE),
    flow_7day = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= benicia_date - 7 &
        rv_flow_clean$date < benicia_date], na.rm = TRUE),
    flow_30day = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= benicia_date - 30 &
        rv_flow_clean$date < benicia_date], na.rm = TRUE),
    # Temp metrics
    temp_day0 = {
      idx <- which(rv_temp_clean$date == benicia_date)
      if(length(idx) == 0) NA_real_ else rv_temp_clean$daily_temp[idx]
    },
    temp_3day = mean(rv_temp_clean$daily_temp[
      rv_temp_clean$date >= benicia_date - 3 &
        rv_temp_clean$date < benicia_date], na.rm = TRUE),
    temp_7day = mean(rv_temp_clean$daily_temp[
      rv_temp_clean$date >= benicia_date - 7 &
        rv_temp_clean$date < benicia_date], na.rm = TRUE),
    temp_30day = mean(rv_temp_clean$daily_temp[
      rv_temp_clean$date >= benicia_date - 30 &
        rv_temp_clean$date < benicia_date], na.rm = TRUE)
  ) %>%
  ungroup()

# Check
cat("Missing flow day0:", sum(is.na(migration_metrics$flow_day0)), "\n")
cat("Missing temp day0:", sum(is.na(migration_metrics$temp_day0)), "\n")
head(migration_metrics)

# Plot 1: Flow at different temporal scales at migration start
migration_metrics %>%
  tidyr::pivot_longer(
    cols = c(flow_day0, flow_3day, flow_7day, flow_30day),
    names_to = "metric",
    values_to = "flow"
  ) %>%
  mutate(metric = factor(metric, 
                         levels = c("flow_day0", "flow_3day", 
                                    "flow_7day", "flow_30day"),
                         labels = c("Day of", "3-day prior", 
                                    "7-day prior", "30-day prior"))) %>%
  ggplot(aes(x = flow)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  facet_wrap(~metric, scales = "free") +
  labs(title = "Flow at Benicia passage (migration start)",
       x = "Flow (cfs)", y = "Number of fish") +
  theme_minimal()

# Plot 2: Temp at different temporal scales at migration start
migration_metrics %>%
  tidyr::pivot_longer(
    cols = c(temp_day0, temp_3day, temp_7day, temp_30day),
    names_to = "metric",
    values_to = "temp"
  ) %>%
  mutate(metric = factor(metric,
                         levels = c("temp_day0", "temp_3day",
                                    "temp_7day", "temp_30day"),
                         labels = c("Day of", "3-day prior",
                                    "7-day prior", "30-day prior"))) %>%
  ggplot(aes(x = temp)) +
  geom_histogram(bins = 30, fill = "tomato", color = "white", alpha = 0.8) +
  facet_wrap(~metric, scales = "free") +
  labs(title = "Temperature at Benicia passage (migration start)",
       x = "Temperature (°F)", y = "Number of fish") +
  theme_minimal()

# Plot 3: Flow vs temp scatter colored by water year
ggplot(migration_metrics, aes(x = flow_day0, y = temp_day0,
                              color = factor(water_year))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Flow vs temperature at migration start",
       x = "Flow day of (cfs)",
       y = "Temperature day of (°F)",
       color = "Water year") +
  theme_minimal()

# Plot 4: Migration start date by water year with flow context
ggplot(migration_metrics, aes(x = benicia_date, y = flow_day0,
                              color = factor(water_year))) +
  geom_point(size = 3, alpha = 0.7) +
  facet_wrap(~water_year, scales = "free_x") +
  labs(title = "Flow at migration start by water year",
       x = "Migration start date",
       y = "Flow (cfs)",
       color = "Water year") +
  theme_minimal()

########################################
#rate of change
# Calculate rate of change of flow (flow pulse detection)
rv_flow_clean <- rv_flow_clean %>%
  arrange(date) %>%
  mutate(
    flow_change_1day = flow_cfs - lag(flow_cfs, 1),
    flow_change_3day = flow_cfs - lag(flow_cfs, 3),
    flow_change_7day = flow_cfs - lag(flow_cfs, 7),
    # Days above a threshold flow (e.g. 20000 cfs)
    above_threshold = flow_cfs > 20000,
    # Cumulative days above threshold in prior 30 days
    days_above_30 = zoo::rollsum(above_threshold, 30, 
                                 fill = NA, align = "right")
  )

# Join to migration start dates
migration_metrics <- migration_metrics %>%
  left_join(
    rv_flow_clean %>% 
      dplyr::select(date, flow_change_1day, flow_change_3day,
                    flow_change_7day, days_above_30),
    by = c("benicia_date" = "date")
  ) %>%
  left_join(
    rv_temp_clean %>% rename(benicia_date = date),
    by = "benicia_date"
  )

# Plot rate of change vs migration timing
ggplot(migration_metrics, 
       aes(x = flow_change_7day, y = benicia_date,
           color = factor(water_year))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Flow rate of change (7-day) at migration start",
       x = "Change in flow 7 days before passage (cfs)",
       y = "Migration start date",
       color = "Water year") +
  theme_minimal()

# Correlation matrix of all flow/temp metrics
library(GGally)
migration_metrics %>%
  dplyr::select(flow_day0, flow_7day, flow_30day,
                flow_change_7day, days_above_30,
                temp_day0, temp_7day) %>%
  filter(complete.cases(.)) %>%
  ggpairs(title = "Correlation between hydrological metrics at migration start")


##########################################################
#looking at more years to include wet years
# Check what years are available in your full dataset
events %>%
  mutate(water_year = ifelse(month(as.Date(first_detection)) >= 10,
                             year(as.Date(first_detection)) + 1,
                             year(as.Date(first_detection)))) %>%
  filter(status %in% c("up_complete", "up_incomplete", "incomplete_dead")) %>%
  count(water_year) %>%
  arrange(water_year)

# Unique fish per water year in extended dataset
filtered_events %>%
  filter(status %in% c("up_complete", "up_incomplete", "incomplete_dead"),
         water_year >= 2007, water_year <= 2021) %>%
  distinct(animal_id, water_year) %>%
  count(water_year) %>%
  arrange(water_year)

# Check Benicia/Carquinez receiver coverage by year
events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  mutate(water_year = ifelse(month(as.Date(first_detection)) >= 10,
                             year(as.Date(first_detection)) + 1,
                             year(as.Date(first_detection)))) %>%
  distinct(animal_id, water_year) %>%
  count(water_year) %>%
  arrange(water_year)





# Calculate flow metrics at multiple temporal scales at Benicia passage
# Using Rio Vista gauge, looking BEFORE Benicia passage

benicia_flow <- benicia_dates %>%
  rowwise() %>%
  mutate(
    flow_3day   = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= benicia_date - 3 &
        rv_flow_clean$date < benicia_date], na.rm = TRUE),
    flow_7day   = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= benicia_date - 7 &
        rv_flow_clean$date < benicia_date], na.rm = TRUE),
    flow_15day  = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= benicia_date - 15 &
        rv_flow_clean$date < benicia_date], na.rm = TRUE),
    flow_30day  = mean(rv_flow_clean$flow_cfs[
      rv_flow_clean$date >= benicia_date - 30 &
        rv_flow_clean$date < benicia_date], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(
    detection_history %>% dplyr::select(animal_id, water_year, occ_2, occ_3),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    route = case_when(
      occ_2 == 2 ~ "Georgiana",
      occ_2 == 3 ~ "DCC",
      occ_3 == 4 ~ "Steamboat/Sutter",
      TRUE ~ "Sacramento"
    )
  )

# ANOVA for each temporal scale
cat("=== 3-day mean flow ===\n")
summary(aov(flow_3day ~ route, data = benicia_flow))

cat("\n=== 7-day mean flow ===\n")
summary(aov(flow_7day ~ route, data = benicia_flow))

cat("\n=== 15-day mean flow ===\n")
summary(aov(flow_15day ~ route, data = benicia_flow))

cat("\n=== 30-day mean flow ===\n")
summary(aov(flow_30day ~ route, data = benicia_flow))

# Compare F-statistics and p-values across temporal scales
anova_results <- data.frame(
  scale = c("3-day", "7-day", "15-day", "30-day"),
  F_stat = c(
    summary(aov(flow_3day  ~ route, data = benicia_flow))[[1]]$`F value`[1],
    summary(aov(flow_7day  ~ route, data = benicia_flow))[[1]]$`F value`[1],
    summary(aov(flow_15day ~ route, data = benicia_flow))[[1]]$`F value`[1],
    summary(aov(flow_30day ~ route, data = benicia_flow))[[1]]$`F value`[1]
  )
)
# Check p-values
cat("30-day ANOVA p-value:\n")
summary(aov(flow_30day ~ route, data = benicia_flow))

###################################################################


