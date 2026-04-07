#flow covariate creation for gs multistate model
#tidally filtered rio vista gauge data 
rv_flow <- read.csv("C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv")

# Clean flow data - keep only what we need
rv_flow_clean <- rv_flow %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_cfs = value) %>%
  filter(!is.na(flow_cfs),
         date >= as.Date("2006-10-01"),  # start of WY2007
         date <= as.Date("2017-09-30")) %>%  # end of WY2017
  arrange(date)

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

# Check if log transform helps
hist(log(rv_flow_clean$flow_cfs), 
     main = "Log flow distribution",
     xlab = "log(flow cfs)")