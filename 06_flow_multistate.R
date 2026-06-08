#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - WITH FLOW COVARIATE
# Script: 06_multistate_model_flow.R
# Author: Erin Tracy
# Last updated: April 2026
#
# PURPOSE:
# Extend the baseline 5-state multistate model (05_multistate_model.R) to
# include river discharge (flow) as a covariate on routing probabilities.
# Also simplifies survival structure from reach-specific to route-specific
# following the data support available with 221 fish.
#
# FLOW COVARIATE APPROACH (following Perry et al. 2018):
# Perry assigned each fish the flow on the date it first arrived at each
# junction (First Passage Time flow, FPT.flow). We follow this approach
# exactly, using the USGS Rio Vista tidally-filtered daily discharge gauge
# (station 11455420) as our flow variable.
#
# Two junction-specific flow covariates:
#   flow_occ2[i]: flow at Rio Vista on date fish i first arrived at occ2
#                 (Rio Vista junction) -> predictor for psi_geo[i]
#   flow_occ3[i]: flow at Rio Vista on date fish i first arrived at occ3
#                 (SR_MOUTH/SS junction) -> predictor for psi_ss[i]
#
# NOTE ON ALTERNATIVE FLOW METRIC (30-day mean, separate analysis):
# We also explored 30-day mean flow before Benicia passage as a predictor.
# ANOVA comparing temporal scales showed 30-day mean best separates route
# choices (F=2.13, p=0.097) vs 3-day (F=1.51), 7-day (F=1.54), 15-day
# (F=1.64). This metric captures antecedent hydrological conditions and
# may be more appropriate for modeling MIGRATION INITIATION (whether/when
# fish start migrating from the ocean) rather than route selection. The
# 30-day mean approach will be used in a separate migration initiation
# analysis (possibly a survival/time-to-event model) using the extended
# 2007-2021 Benicia dataset.
#
# KEY DIFFERENCES FROM PERRY 2018:
# 1. Perry used Freeport gauge; we use Rio Vista (closest gauge to junctions)
# 2. Perry modeled travel time explicitly to get exact arrival dates;
#    we use observed first detection dates at each junction occasion
# 3. Perry had fish size as a covariate on survival; we do not (tagging
#    to migration lag too long for size at tagging to be meaningful)
# 4. Perry had DCC gate open/closed covariate; we do not have gate records
# 5. Perry had release group random effects; our fish are wild (no releases)
#
# SIMPLIFIED SURVIVAL STRUCTURE:
# Baseline model estimated reach-specific Sacramento survival (S_sac1-S_sac5)
# but all estimates were essentially identical (0.986-0.993) suggesting the
# data cannot distinguish between reach-specific survival rates. This is
# expected given only 4 confirmed mortality events in 221 fish over 11 years.
# We therefore use a single S_sac parameter for all Sacramento reaches,
# keeping S_geo and S_ss separate as they represent biologically different
# habitats and movement corridors.
#
# LOGISTIC REGRESSION ON ROUTING:
# Following Perry et al. we model routing probabilities as individual-level
# functions of flow using a logistic regression link:
#   logit(psi_geo[i]) = alpha_geo + beta_geo * flow_occ2[i]
#   logit(psi_ss[i])  = alpha_ss  + beta_ss  * flow_occ3[i]
# where flow_occ2[i] is flow at Rio Vista on the day fish i arrived at the
# Rio Vista junction (occ2), and flow_occ3[i] is flow on the day fish i
# arrived at the SR_MOUTH/SS junction (occ3). Both are standardized to
# mean=0, sd=1.
#
# PRIORS ON FLOW EFFECTS (following Perry 2018):
# Perry used Student-t priors with scale sigma=2.5 (tau=0.16) on flow slopes
# following Gelman's recommendation for logistic regression coefficients.
# We use Normal(0, sd=1.5) priors which are similar in behavior and more
# natural in NIMBLE. On the logit scale, sd=1.5 allows routing probabilities
# to vary substantially with flow while remaining weakly informative.
#
# INPUTS:
#   - gs_multistate_data.RData (from 04_multistate_data_prep.R)
#   - daily_tidalfilter_riovista.csv: USGS tidally-filtered daily discharge
#
# OUTPUTS:
#   - gs_mcmc_flow_run1.RData: MCMC posterior samples
#   - gs_mcmc_traces_flow_run1.pdf: trace plots
#
# MODEL STRUCTURE:
#   States: 1=Sacramento, 2=Georgiana, 3=Steamboat/Sutter,
#           4=Dead (absorbing), 5=Failed migration (absorbing)
#   Fish: 221 (134 up_complete, 83 up_incomplete, 4 incomplete_dead)
#==============================================================================

#==============================================================================
# SECTION 1: LOAD LIBRARIES AND DATA
#==============================================================================

library(nimble)
library(nimbleEcology)
library(abind)
library(MCMCvis)
library(dplyr)
library(lubridate)
library(tidyr)

load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_multistate_data.RData")

#==============================================================================
# SECTION 2: RECODE DETECTION HISTORY TO 5-STATE MODEL
# Same recoding as baseline model — DCC dropped (0 confirmed upstream migrants)
# Old coding: 1=Sac, 2=Geo, 3=DCC(dropped), 4=SS, 5=Dead, 6=Failed, 7=Not-detected
# New coding: 1=Sac, 2=Geo, 3=SS, 4=Dead, 5=Failed, 6=Not-detected
#==============================================================================

nstate <- 5

ch_mat_nimble_5state <- ch_mat_nimble
ch_mat_nimble_5state[ch_mat_nimble == 4] <- 3   # SS: 4 -> 3
ch_mat_nimble_5state[ch_mat_nimble == 5] <- 4   # dead: 5 -> 4
ch_mat_nimble_5state[ch_mat_nimble == 6] <- 5   # failed: 6 -> 5
ch_mat_nimble_5state[ch_mat_nimble == 7] <- 6   # not-detected: 7 -> 6

ch_mat_nimble <- ch_mat_nimble_5state

cat("Total fish:", nrow(ch_mat_nimble), "\n")
print(table(detection_history$status))

#==============================================================================
# SECTION 3: BUILD FLOW COVARIATES
#
# APPROACH (following Perry et al. 2018 FPT.flow method):
# Assign each fish the flow at Rio Vista on the date it first arrived at
# each junction occasion. This directly links the routing decision to the
# hydrological conditions the fish experienced at the moment of decision.
#
# Perry's exact implementation used cumulative travel time to assign each
# fish its flow at each occasion: FPT.flow[FPT.index[i,t]] where FPT.index
# maps fish i's arrival time at occasion t to the daily flow vector.
# We replicate this using observed first detection dates at each occasion.
#
# TWO JUNCTION-SPECIFIC FLOW COVARIATES:
#
# flow_occ2[i]: flow at Rio Vista on date fish i first arrived at occ2
#   - occ2 is the Rio Vista junction where fish choose Georgiana vs Sacramento
#   - This is the flow at the time of the Georgiana routing decision
#   - Used as predictor for psi_geo[i]
#   - Fish not detected at occ2 (failed migrants that never reached Rio Vista)
#     get flow_occ2 = 0 (mean flow) as a placeholder — these fish do not
#     contribute to the psi_geo likelihood anyway
#
# flow_occ3[i]: flow at Rio Vista on date fish i first arrived at occ3
#   - occ3 is the SR_MOUTH/SS junction where fish choose SS vs Sacramento
#   - This is the flow at the time of the SS routing decision
#   - Used as predictor for psi_ss[i]
#   - Same imputation logic for fish that never reached occ3
#
# NOTE: We use the first detection at the COMMITTED ROUTE receivers, not
# the first detection at any occ2/occ3 receiver. This is consistent with
# last(state) detection history construction and avoids assigning flow
# from exploratory staging behavior to the routing decision.
#==============================================================================

# Load Rio Vista flow data
rv_flow <- read.csv("C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv")

# Clean flow data - keep only study period
rv_flow_complete <- rv_flow %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_cfs = value) %>%
  filter(!is.na(flow_cfs),
         date >= as.Date("2006-10-01"),   # start of WY2007
         date <= as.Date("2017-09-30")) %>%  # end of WY2017
  arrange(date)

cat("Flow data days:", nrow(rv_flow_complete), "\n")
cat("Missing flow values:", sum(is.na(rv_flow_complete$flow_cfs)), "\n")

# CHECK: identify any gaps in the daily flow record
# These are dates where the gauge has no data — fish arriving on gap dates
# will get NA flow values and need to be handled by interpolation below
rv_flow_complete %>%
  mutate(gap_after = as.numeric(lead(date) - date)) %>%
  filter(gap_after > 1) %>%
  dplyr::select(date, flow_cfs, gap_after) %>%
  mutate(gap_start = date + 1,
         gap_end   = date + gap_after - 1) %>%
  print()

# Fill gaps with linear interpolation
# Rationale: gaps are short (confirmed above) so linear interpolation
# between adjacent known values is a reasonable approximation
# Only 1 fish (GS0118 WY2007, arrival 2007-05-21) falls in a gap
library(zoo)
rv_flow_complete <- rv_flow_clean %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_cfs = zoo::na.approx(flow_cfs, na.rm = FALSE))

cat("Days before interpolation:", nrow(rv_flow_clean), "\n")
cat("Days after interpolation:", nrow(rv_flow_complete), "\n")
cat("Remaining NAs:", sum(is.na(rv_flow_complete$flow_cfs)), "\n")

# Then use rv_flow_complete (not rv_flow_clean) for all subsequent joins

geo_flow <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/RioVista_Confluence_Flows.csv")
library(lubridate)

geo_flow_daily <- geo_flow %>%
  mutate(date = as.Date(dateTime, format = "%m/%d/%Y %H:%M")) %>%
  group_by(date) %>%
  dplyr::summarise(
    flow_GES        = mean(GES,        na.rm = TRUE),
    flow_RIO        = mean(RIO,        na.rm = TRUE),
    flow_RYI_MB_GES = mean(RYI.MB.GES, na.rm = TRUE),
    flow_RYI_MB_SOI = mean(RYI.MB.SOI, na.rm = TRUE),
    flow_SXS        = mean(SXS,        na.rm = TRUE),
    flow_SOI        = mean(SOI,        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(date),
         date >= as.Date("2006-10-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date)

cat("Daily flow rows:", nrow(geo_flow_daily), "\n")
cat("Date range:", as.character(range(geo_flow_daily$date)), "\n")
cat("NAs in RYI_MB_GES:", sum(is.na(geo_flow_daily$flow_RYI_MB_GES)), "\n")
cat("NAs in GES:", sum(is.na(geo_flow_daily$flow_GES)), "\n")

# Interpolate gaps in GES
geo_flow_daily <- geo_flow_daily %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_GES = zoo::na.approx(flow_GES, na.rm = FALSE))

cat("NAs after interpolation:", sum(is.na(geo_flow_daily$flow_GES)), "\n")

# Check flow ranges look reasonable
summary(geo_flow_daily$flow_GES)


# Get first detection date at each junction for each fish
# Using the committed route state to be consistent with detection history
# i.e. for a Geo fish at occ2, use first Geo detection date
#      for a Sac fish at occ2, use first Rio Vista Sacramento detection date
junction_dates <- model1_events %>%
  filter(occasion %in% c(2, 3), !is.na(state)) %>%
  arrange(animal_id, water_year, occasion, first_detection) %>%
  group_by(animal_id, water_year, occasion) %>%
  dplyr::summarise(
    # Get arrival date at the committed route state
    # last(state) matches the detection history coding
    committed_state = last(state),
    arrival_date = as.Date(min(first_detection[state == last(state)])),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = occasion,
    values_from = c(arrival_date, committed_state),
    names_sep   = "_occ"
  )

# Join flow on date of arrival at each junction
# Following Perry: flow[i] = daily flow on date fish i arrived at occasion t
fish_flow <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(junction_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # flow_occ2: flow on date fish arrived at Rio Vista junction
    # This is the flow conditions when Georgiana routing decision was made
    flow_occ2_raw = {
      idx <- which(rv_flow_complete$date == arrival_date_occ2)
      if(length(idx) == 0) NA_real_ else rv_flow_complete$flow_cfs[idx]
    },
    
    # flow_occ3: flow on date fish arrived at SR_MOUTH/SS junction
    # This is the flow conditions when SS routing decision was made
    flow_occ3_raw = {
      idx <- which(rv_flow_complete$date == arrival_date_occ3)
      if(length(idx) == 0) NA_real_ else rv_flow_complete$flow_cfs[idx]
    }
    
  ) %>%
  ungroup()

cat("Missing flow_occ2:", sum(is.na(fish_flow$flow_occ2_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_flow$flow_occ3_raw)), "\n")

# Standardize flow following Perry et al.
# Using same mean and SD for both covariates so coefficients are comparable
# Standardizing to mean=0, sd=1 so:
#   flow_std = 0  -> average flow conditions
#   flow_std > 0  -> above average (wet year/high discharge)
#   flow_std < 0  -> below average (dry year/low discharge)
# This also improves MCMC mixing and makes beta coefficients directly
# interpretable as change in log-odds per SD of flow

# Compute scaling using all non-missing flow values from both occasions
all_flow_vals <- c(fish_flow$flow_occ2_raw, fish_flow$flow_occ3_raw)
flow_mean <- mean(all_flow_vals, na.rm = TRUE)
flow_sd   <- sd(all_flow_vals, na.rm = TRUE)

fish_flow <- fish_flow %>%
  mutate(
    flow_occ2_std = (flow_occ2_raw - flow_mean) / flow_sd,
    flow_occ3_std = (flow_occ3_raw - flow_mean) / flow_sd
  )

cat("\nFlow covariate summary:\n")
cat("Scaling mean:", round(flow_mean, 0), "cfs\n")
cat("Scaling SD:", round(flow_sd, 0), "cfs\n")
cat("flow_occ2 range:", round(range(fish_flow$flow_occ2_std, na.rm = TRUE), 2), "\n")
cat("flow_occ3 range:", round(range(fish_flow$flow_occ3_std, na.rm = TRUE), 2), "\n")

# Verify row order matches detection_history before extracting vectors
identical(fish_flow$animal_id, detection_history$animal_id)
identical(as.numeric(fish_flow$water_year), as.numeric(detection_history$water_year))

# Extract flow vectors for NIMBLE
# Impute 0 (mean flow) for fish missing junction dates
# Rationale: fish that never reached a junction do not contribute to the
# routing likelihood for that junction, so the flow value does not affect
# parameter estimates. Using 0 (mean) is a neutral placeholder.
flow_occ2 <- fish_flow$flow_occ2_std
flow_occ3 <- fish_flow$flow_occ3_std
flow_occ2[is.na(flow_occ2)] <- 0
flow_occ3[is.na(flow_occ3)] <- 0

cat("\nFinal flow vectors:\n")
cat("flow_occ2 - n fish:", length(flow_occ2),
    "| NAs imputed:", sum(is.na(fish_flow$flow_occ2_std)), "\n")
cat("flow_occ3 - n fish:", length(flow_occ3),
    "| NAs imputed:", sum(is.na(fish_flow$flow_occ3_std)), "\n")

# Distribution of real (non-imputed) flow values
par(mfrow = c(1,2))
hist(flow_occ2[!is.na(fish_flow$flow_occ2_std)],
     main = "Flow at occ2 junction\n(Georgiana decision, non-imputed)",
     xlab = "Standardized flow (cfs)",
     col = "steelblue")
abline(v = 0, lty = 2, col = "red")
hist(flow_occ3[!is.na(fish_flow$flow_occ3_std)],
     main = "Flow at occ3 junction\n(SS decision, non-imputed)",
     xlab = "Standardized flow (cfs)",
     col = "forestgreen")
abline(v = 0, lty = 2, col = "red")

#==============================================================================
# SECTION 4: VERIFY FLOW COVARIATE MAKES BIOLOGICAL SENSE
# Check that flow at each junction varies meaningfully across water years
# and that observed patterns match exploratory analysis:
#   - SS use DECREASES with flow (confirmed by route quartile plot)
#   - Geo is relatively flow-independent
#==============================================================================

# Mean flow by water year at each junction
cat("Mean flow at occ2 (Geo decision) by water year:\n")
fish_flow %>%
  filter(!is.na(flow_occ2_std)) %>%
  left_join(detection_history %>%
              dplyr::select(animal_id, water_year),
            by = c("animal_id", "water_year")) %>%
  group_by(water_year) %>%
  dplyr::summarise(
    n_fish = n(),
    mean_flow_std = round(mean(flow_occ2_std), 2),
    mean_flow_cfs = round(mean(flow_occ2_raw), 0),
    .groups = "drop"
  ) %>%
  arrange(water_year) %>%
  print()

# Flow by route taken - check expected direction
fish_flow %>%
  left_join(
    detection_history %>% dplyr::select(animal_id, water_year, occ_2, occ_3),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    route = case_when(
      occ_2 == 2 ~ "Georgiana",
      occ_3 == 3 ~ "Steamboat/Sutter",
      TRUE       ~ "Sacramento"
    )
  ) %>%
  group_by(route) %>%
  dplyr::summarise(
    n              = n(),
    mean_flow_occ2 = round(mean(flow_occ2_std, na.rm = TRUE), 2),
    mean_flow_occ3 = round(mean(flow_occ3_std, na.rm = TRUE), 2),
    .groups        = "drop"
  )

#==============================================================================
# SECTION 5: BUILD PLACEHOLDER TRANSITION AND OBSERVATION MATRICES
# Used for likelihood check only — in NIMBLE model these are fish-specific
#==============================================================================

# Placeholder parameter values
S_sac    <- 0.99; S_geo <- 0.95; S_ss <- 0.97
phi_fail <- 0.09; lambda <- 0.99

# Population-mean routing at mean flow for placeholder matrices
# In NIMBLE model these will be individual-level functions of flow
psi_geo_mean <- 0.124   # baseline estimate
psi_ss_mean  <- 0.374   # baseline estimate

# Detection probabilities from baseline model
p_sac2 <- 0.977; p_sac3 <- 0.899; p_sac4 <- 0.910
p_sac5 <- 0.893; p_sac6 <- 0.986
p_geo  <- 0.771; p_ss   <- 0.690

temp_mat <- matrix(0, nrow = nstate + 1, ncol = nstate + 1)
temp_mat[nstate + 1, nstate + 1] <- 1

tr_12 <- temp_mat
tr_12[1, ] <- c(S_sac*(1-psi_geo_mean)*(1-phi_fail),
                S_sac*psi_geo_mean*(1-phi_fail),
                0, (1-S_sac), S_sac*phi_fail, 0)
tr_12[2, 4] <- 1; tr_12[3, 4] <- 1
tr_12[4, 4] <- 1; tr_12[5, 5] <- 1

tr_23 <- temp_mat
tr_23[1, ] <- c(S_sac*(1-psi_ss_mean)*(1-phi_fail),
                0, S_sac*psi_ss_mean*(1-phi_fail),
                (1-S_sac), S_sac*phi_fail, 0)
tr_23[2, 2] <- 1
tr_23[3, 4] <- 1; tr_23[4, 4] <- 1; tr_23[5, 5] <- 1

tr_34 <- temp_mat
tr_34[1, ] <- c(S_sac*(1-phi_fail), 0, 0, (1-S_sac), S_sac*phi_fail, 0)
tr_34[2, ] <- c(S_geo*(1-phi_fail), 0, 0, (1-S_geo), S_geo*phi_fail, 0)
tr_34[3, ] <- c(0, 0, (1-phi_fail), 0, phi_fail, 0)
tr_34[4, 4] <- 1; tr_34[5, 5] <- 1

tr_45 <- temp_mat
tr_45[1, ] <- c(S_sac*(1-phi_fail), 0, 0, (1-S_sac), S_sac*phi_fail, 0)
tr_45[2, 4] <- 1
tr_45[3, ] <- c(S_ss*(1-phi_fail), 0, 0, (1-S_ss), S_ss*phi_fail, 0)
tr_45[4, 4] <- 1; tr_45[5, 5] <- 1

tr_56 <- temp_mat
tr_56[1, ] <- c(S_sac*(1-phi_fail), 0, 0, (1-S_sac), S_sac*phi_fail, 0)
tr_56[2, 4] <- 1; tr_56[3, 4] <- 1
tr_56[4, 4] <- 1; tr_56[5, 5] <- 1

tr_67 <- temp_mat
tr_67[1, ] <- c(lambda*(1-phi_fail), 0, 0, (1-lambda), lambda*phi_fail, 0)
tr_67[2, 4] <- 1; tr_67[3, 4] <- 1
tr_67[4, 4] <- 1; tr_67[5, 5] <- 1

tr_arr <- abind(tr_12, tr_23, tr_34, tr_45, tr_56, tr_67, along = 3)

cat("Transition matrix row sums:\n")
for(i in 1:6) cat("tr", i, ":", rowSums(tr_arr[,,i]), "\n")

# Observation matrices - identical to baseline model
p_mat1 <- temp_mat
p_mat1[1, 1] <- 1
p_mat1[2, nstate+1] <- 1; p_mat1[3, nstate+1] <- 1
p_mat1[4, 4] <- 1; p_mat1[5, 5] <- 1

p_mat2 <- temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, 0, (1-p_sac2))
p_mat2[2, ] <- c(0, p_geo, 0, 0, 0, (1-p_geo))
p_mat2[3, nstate+1] <- 1
p_mat2[4, 4] <- 1; p_mat2[5, 5] <- 1

p_mat3 <- temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, 0, (1-p_sac3))
p_mat3[2, nstate+1] <- 1
p_mat3[3, ] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
p_mat3[4, 4] <- 1; p_mat3[5, 5] <- 1

p_mat4 <- temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, 0, (1-p_sac4))
p_mat4[2, nstate+1] <- 1
p_mat4[3, ] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
p_mat4[4, 4] <- 1; p_mat4[5, 5] <- 1

p_mat5 <- temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, 0, (1-p_sac5))
p_mat5[2, nstate+1] <- 1; p_mat5[3, nstate+1] <- 1
p_mat5[4, 4] <- 1; p_mat5[5, 5] <- 1

p_mat6 <- temp_mat
p_mat6[1, ] <- c(p_sac6, 0, 0, 0, 0, (1-p_sac6))
p_mat6[2, nstate+1] <- 1; p_mat6[3, nstate+1] <- 1
p_mat6[4, 4] <- 1; p_mat6[5, 5] <- 1

p_mat7 <- temp_mat
p_mat7[1, 1] <- 1
p_mat7[2, nstate+1] <- 1; p_mat7[3, nstate+1] <- 1
p_mat7[4, 4] <- 1; p_mat7[5, 5] <- 1

p_arr <- abind(p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, p_mat7, along = 3)

cat("Observation matrix row sums:\n")
for(i in 1:7) cat("p", i, ":", rowSums(p_arr[,,i]), "\n")

rel_vec <- c(1, 0, 0, 0, 0, 0)

#==============================================================================
# SECTION 6: LIKELIHOOD CHECK WITH PLACEHOLDER VALUES
#==============================================================================

all_ll <- apply(ch_mat_nimble, 1, function(x)
  dDHMMo(x, init = rel_vec, probObs = p_arr, probTrans = tr_arr,
         len = 7, checkRowSums = FALSE, log = TRUE))

cat("NaN count:", sum(is.nan(all_ll)), "\n")
cat("Inf count:", sum(is.infinite(all_ll)), "\n")
cat("LL range:", range(all_ll[!is.nan(all_ll) & !is.infinite(all_ll)]), "\n")

if(sum(is.nan(all_ll)) > 0){
  cat("Fish producing NaN:\n")
  print(detection_history[which(is.nan(all_ll)),
                          c("animal_id", "water_year", "occ_1",
                            "occ_2", "occ_3", "occ_4", "occ_5",
                            "occ_6", "occ_7", "status")])
  print(ch_mat_nimble[which(is.nan(all_ll)), ])
}

#==============================================================================
# SECTION 7: NIMBLE MODEL WITH FLOW COVARIATE
#
# STRUCTURE OF FLOW MODEL:
#
# For each fish i the routing probabilities are individual-level functions
# of the flow each fish experienced at the time of its routing decision:
#   logit(psi_geo[i]) = alpha_geo + beta_geo * flow_occ2[i]
#   logit(psi_ss[i])  = alpha_ss  + beta_ss  * flow_occ3[i]
#
# alpha parameters are intercepts = log-odds of taking each route at
# mean flow conditions (flow_std = 0). Back-transforming with ilogit()
# gives the routing probability at mean flow, comparable to the fixed
# psi estimates from the baseline model.
#
# beta parameters are flow slopes:
#   beta > 0: higher flow -> higher probability of taking that route
#   beta < 0: higher flow -> lower probability of taking that route
#
# From exploratory analysis (route selection by flow quartile plot) we expect:
#   beta_ss  < 0 (SS use decreases with increasing flow)
#   beta_geo ~ 0 (Geo relatively flow-independent across quartiles)
#
# TRANSITION MATRICES ARE FISH-SPECIFIC:
# Because psi_geo[i] and psi_ss[i] vary by fish, the transition matrices
# at occ1->2 (Geo decision) and occ2->3 (SS decision) must be computed
# inside the individual fish loop in NIMBLE. This is computationally
# more expensive than the baseline model but is necessary to capture
# individual-level flow effects following Perry's approach.
#
# PRIORS:
# alpha (intercepts): Normal(0, sd=1.5) on logit scale
#   -> routing probability between ~5% and 95% a priori at mean flow
#   -> weakly informative, allows data to dominate
#   -> initialized at baseline model log-odds estimates for efficiency
# beta (slopes): Normal(0, sd=1.0) on logit scale
#   -> one SD change in flow changes log-odds by ~1 unit a priori
#   -> allows moderate to strong flow effects while regularizing
#   -> comparable to Perry's Student-t(0, tau=0.16, df=7) but more
#      natural in NIMBLE
# All other parameters use flat dbeta(1,1) priors as in baseline model
#==============================================================================


nimCode_flow <- nimbleCode({
  
  #--- PRIORS ---
  
  # Simplified survival - single parameter per route
  # Data cannot distinguish reach-specific survival (all 0.986-0.993 in baseline)
  S_sac    ~ dbeta(1, 1)   # all Sacramento reaches combined
  S_geo    ~ dbeta(1, 1)   # Georgiana Slough
  S_ss     ~ dbeta(1, 1)   # Steamboat/Sutter Sloughs
  
  # Detection probabilities - same structure as baseline model
  # p_sac1 not estimated (certain detection at occ1 by definition)
  p_sac2   ~ dbeta(1, 1)
  p_sac3   ~ dbeta(1, 1)
  p_sac4   ~ dbeta(1, 1)
  p_sac5   ~ dbeta(1, 1)
  p_sac6   ~ dbeta(1, 1)
  p_geo    ~ dbeta(1, 1)
  p_ss     ~ dbeta(1, 1)
  
  # Migration failure and spawning ground arrival
  phi_fail ~ dbeta(1, 1)
  lambda   ~ dbeta(1, 1)
  
  # Flow covariate intercepts (alpha) on logit scale
  # Represents routing probability at mean flow (flow_std = 0)
  # Normal(0, sd=1.5): routing prob between ~5-95% at mean flow a priori
  alpha_geo ~ dnorm(0, sd = 1.5)
  alpha_ss  ~ dnorm(0, sd = 1.5)
  
  # Flow covariate slopes (beta) on logit scale
  # Represents change in log-odds of taking route per 1 SD increase in flow
  # Normal(0, sd=1.0): allows moderate-strong effects, regularizing
  # Expected signs from exploratory analysis:
  #   beta_ss  < 0 (SS use decreases with flow)
  #   beta_geo ~ 0 (Geo relatively flow-independent)
  beta_geo ~ dnorm(0, sd = 1.0)
  beta_ss  ~ dnorm(0, sd = 1.0)
  
  #--- INDIVIDUAL-LEVEL ROUTING PROBABILITIES ---
  # Each fish gets its own routing probability based on flow at its junction
  # This is the key addition relative to the baseline model
  for(i in 1:nfish){
    
    # Georgiana routing probability for fish i
    # Uses flow at Rio Vista on day fish i arrived at occ2 (Rio Vista junction)
    # Following Perry: FPT.flow[FPT.index[i, occ2]]
    psi_geo[i] <- ilogit(alpha_geo + beta_geo * flow_occ2[i])
    
    # Steamboat/Sutter routing probability for fish i
    # Uses flow at Rio Vista on day fish i arrived at occ3 (SR_MOUTH junction)
    # Following Perry: FPT.flow[FPT.index[i, occ3]]
    psi_ss[i]  <- ilogit(alpha_ss  + beta_ss  * flow_occ3[i])
    
  } #end i loop for individual routing probabilities
  
  #--- OBSERVATION ARRAY ---
  # Identical structure to baseline model - detection is not flow-dependent
  # States: 1=Sac, 2=Geo, 3=SS, 4=Dead, 5=Failed, 6=Not-detected
  
  # Occasion 1: Benicia/Carquinez - certain detection
  p_arr[1, 1:6, 1] <- c(1, 0, 0, 0, 0, 0)
  p_arr[2, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 1] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 2: Rio Vista junction
  p_arr[1, 1:6, 2] <- c(p_sac2, 0, 0, 0, 0, (1-p_sac2))
  p_arr[2, 1:6, 2] <- c(0, p_geo, 0, 0, 0, (1-p_geo))
  p_arr[3, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 2] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 3: SR_MOUTH/Steamboat junction
  # Geo fish pass through undetected at occ3 (no receivers in Geo channel)
  p_arr[1, 1:6, 3] <- c(p_sac3, 0, 0, 0, 0, (1-p_sac3))
  p_arr[2, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 3] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[4, 1:6, 3] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 3] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 4: SR_BLWSTEAM area - Geo rejoined Sac, SS still in channel
  p_arr[1, 1:6, 4] <- c(p_sac4, 0, 0, 0, 0, (1-p_sac4))
  p_arr[2, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 4] <- c(0, 0, p_ss, 0, 0, (1-p_ss))
  p_arr[4, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 4] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 5: SR_KK345R + SR_FREEPORT - SS has rejoined Sacramento
  p_arr[1, 1:6, 5] <- c(p_sac5, 0, 0, 0, 0, (1-p_sac5))
  p_arr[2, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 5] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 6: upper Sacramento
  p_arr[1, 1:6, 6] <- c(p_sac6, 0, 0, 0, 0, (1-p_sac6))
  p_arr[2, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 6] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
  
  # Occasion 7: spawning ground - certain detection
  p_arr[1, 1:6, 7] <- c(1, 0, 0, 0, 0, 0)
  p_arr[2, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  p_arr[3, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  p_arr[4, 1:6, 7] <- c(0, 0, 0, 1, 0, 0)
  p_arr[5, 1:6, 7] <- c(0, 0, 0, 0, 1, 0)
  p_arr[6, 1:6, 7] <- c(0, 0, 0, 0, 0, 1)
  
  #--- LIKELIHOOD ---
  # tr_arr_i is now a 4D array [fish, from_state, to_state, transition]
  # This correctly gives each fish its own transition matrix
  # and eliminates the multiple definitions warning
  
  for(i in 1:nfish){
    
    # Transition 1->2: Benicia to Rio Vista (Georgiana routing decision)
    # Row 1 uses fish-specific psi_geo[i] — this is the key individual-level effect
    tr_arr_i[i, 1, 1:6, 1] <- c(S_sac*(1-psi_geo[i])*(1-phi_fail),
                                S_sac*psi_geo[i]*(1-phi_fail),
                                0, (1-S_sac), S_sac*phi_fail, 0)
    tr_arr_i[i, 2, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 1] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 1] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 1] <- c(0, 0, 0, 0, 0, 1)
    
    # Transition 2->3: Rio Vista to SR_MOUTH (SS routing decision)
    # Row 1 uses fish-specific psi_ss[i]
    tr_arr_i[i, 1, 1:6, 2] <- c(S_sac*(1-psi_ss[i])*(1-phi_fail),
                                0,
                                S_sac*psi_ss[i]*(1-phi_fail),
                                (1-S_sac), S_sac*phi_fail, 0)
    tr_arr_i[i, 2, 1:6, 2] <- c(0, 1, 0, 0, 0, 0)
    tr_arr_i[i, 3, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 2] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 2] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 2] <- c(0, 0, 0, 0, 0, 1)
    
    # Transitions 3->7: identical across fish but still indexed by i
    # to maintain correct DAG structure
    tr_arr_i[i, 1, 1:6, 3] <- c(S_sac*(1-phi_fail), 0, 0,
                                (1-S_sac), S_sac*phi_fail, 0)
    tr_arr_i[i, 2, 1:6, 3] <- c(S_geo*(1-phi_fail), 0, 0,
                                (1-S_geo), S_geo*phi_fail, 0)
    tr_arr_i[i, 3, 1:6, 3] <- c(0, 0, (1-phi_fail), 0, phi_fail, 0)
    tr_arr_i[i, 4, 1:6, 3] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 3] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 3] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 4] <- c(S_sac*(1-phi_fail), 0, 0,
                                (1-S_sac), S_sac*phi_fail, 0)
    tr_arr_i[i, 2, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 4] <- c(S_ss*(1-phi_fail), 0, 0,
                                (1-S_ss), S_ss*phi_fail, 0)
    tr_arr_i[i, 4, 1:6, 4] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 4] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 4] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 5] <- c(S_sac*(1-phi_fail), 0, 0,
                                (1-S_sac), S_sac*phi_fail, 0)
    tr_arr_i[i, 2, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 5] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 5] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 5] <- c(0, 0, 0, 0, 0, 1)
    
    tr_arr_i[i, 1, 1:6, 6] <- c(lambda*(1-phi_fail), 0, 0,
                                (1-lambda), lambda*phi_fail, 0)
    tr_arr_i[i, 2, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 3, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 4, 1:6, 6] <- c(0, 0, 0, 1, 0, 0)
    tr_arr_i[i, 5, 1:6, 6] <- c(0, 0, 0, 0, 1, 0)
    tr_arr_i[i, 6, 1:6, 6] <- c(0, 0, 0, 0, 0, 1)
    
    # Likelihood for fish i — now uses fish i's own slice of the 4D array
    ch_mat[i, 1:7] ~ dDHMMo(
      init         = c(1, 0, 0, 0, 0, 0)[1:6],
      probObs      = p_arr[1:6, 1:6, 1:7],
      probTrans    = tr_arr_i[i, 1:6, 1:6, 1:6],
      len          = 7,
      checkRowSums = 0
    )
    
  } #end i loop
})

#==============================================================================
# SECTION 8: BUILD NIMBLE MODEL AND CONFIGURE MCMC
#==============================================================================

inits_flow <- list(
  S_sac    = 0.99, S_geo = 0.95, S_ss = 0.97,
  p_sac2   = 0.977, p_sac3 = 0.899,
  p_sac4   = 0.910, p_sac5 = 0.893, p_sac6 = 0.986,
  p_geo    = 0.771, p_ss   = 0.690,
  phi_fail = 0.091, lambda = 0.986,
  # Initialize intercepts at baseline model log-odds estimates
  # ilogit(-1.95) = 0.124 ~ baseline psi_geo
  # ilogit(-0.51) = 0.375 ~ baseline psi_ss
  alpha_geo = -1.95,
  alpha_ss  = -0.51,
  # Initialize slopes at 0 — no flow effect assumption
  # MCMC will update toward data-supported values
  beta_geo  = 0.0,
  beta_ss   = 0.0
)

nimMod_flow <- nimbleModel(
  code      = nimCode_flow,
  inits     = inits_flow,
  data      = list(ch_mat = ch_mat_nimble,
                   flow_occ2 = flow_occ2,   # flow at Rio Vista on date of occ2 arrival
                   flow_occ3 = flow_occ3), # flow at Rio Vista on date of occ3 arrival
  constants = list(
    nfish     = nrow(ch_mat_nimble))
  )


nimMod_flow$calculate()

inits_fn_flow <- function(){
  list(
    S_sac    = runif(1, 0.95, 1.0),
    S_geo    = runif(1, 0.80, 1.0),
    S_ss     = runif(1, 0.85, 1.0),
    p_sac2   = runif(1, 0.90, 1.0),
    p_sac3   = runif(1, 0.75, 0.95),
    p_sac4   = runif(1, 0.80, 0.95),
    p_sac5   = runif(1, 0.80, 0.95),
    p_sac6   = runif(1, 0.95, 1.0),
    p_geo    = runif(1, 0.50, 0.90),
    p_ss     = runif(1, 0.50, 0.85),
    phi_fail = runif(1, 0.05, 0.15),
    lambda   = runif(1, 0.95, 1.0),
    # Intercepts near baseline estimates with small noise
    alpha_geo = rnorm(1, mean = -1.95, sd = 0.3),
    alpha_ss  = rnorm(1, mean = -0.51, sd = 0.3),
    # Slopes near zero
    beta_geo  = rnorm(1, mean = 0, sd = 0.2),
    beta_ss   = rnorm(1, mean = 0, sd = 0.2)
  )
}

params_flow <- c(
  # Survival
  "S_sac", "S_geo", "S_ss",
  # Detection
  "p_sac2", "p_sac3", "p_sac4", "p_sac5", "p_sac6",
  "p_geo", "p_ss",
  # Failure and spawning arrival
  "phi_fail", "lambda",
  # Flow model parameters - KEY RESULTS
  "alpha_geo", "alpha_ss",   # intercepts: routing probability at mean flow
  "beta_geo", "beta_ss"      # slopes: change in log-odds per SD of flow
)

# Slice samplers following Perry et al.
confMCMC_flow <- configureMCMC(nimMod_flow, onlySlice = TRUE)
confMCMC_flow$addMonitors(params_flow)
MCMC_flow   <- buildMCMC(confMCMC_flow)
CModel_flow <- compileNimble(nimMod_flow)
CMCMC_flow  <- compileNimble(MCMC_flow, project = CModel_flow)

#==============================================================================
# SECTION 9: RUN MCMC AND SAVE RESULTS
#==============================================================================
mcmc_out_flow <- runMCMC(
  CMCMC_flow,
  niter   = 50000,
  nchains = 3,
  nburnin = 10000,
  thin    = 10,
  inits   = list(inits_fn_flow(), inits_fn_flow(), inits_fn_flow()),
  samplesAsCodaMCMC = TRUE
)

MCMCsummary(mcmc_out_flow, round = 3)


save(mcmc_out_flow, flow_occ2, flow_occ3, fish_flow,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_run1.RData")

MCMCtrace(mcmc_out_flow,
          params   = c("alpha_geo", "alpha_ss", "beta_geo", "beta_ss",
                       "phi_fail", "lambda"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_flow_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)

#==============================================================================
# SECTION 10: INTERPRET FLOW RESULTS
# Run after MCMC completes
#==============================================================================

flow_summary <- MCMCsummary(mcmc_out_flow, round = 3,
                            params = c("alpha_geo", "alpha_ss",
                                       "beta_geo", "beta_ss"))
print(flow_summary)

# Back-transform intercepts to routing probabilities at mean flow
# alpha -> ilogit(alpha) = psi at mean flow (flow_std = 0)
# Compare to baseline model fixed estimates
cat("\nRouting probabilities at mean flow (flow_std = 0):\n")
alpha_geo_post <- MCMCpstr(mcmc_out_flow, params = "alpha_geo",
                           func = mean)[[1]]
alpha_ss_post  <- MCMCpstr(mcmc_out_flow, params = "alpha_ss",
                           func = mean)[[1]]
beta_geo_post  <- MCMCpstr(mcmc_out_flow, params = "beta_geo",
                           func = mean)[[1]]
beta_ss_post   <- MCMCpstr(mcmc_out_flow, params = "beta_ss",
                           func = mean)[[1]]

cat("psi_geo at mean flow:", round(plogis(alpha_geo_post), 3),
    "(baseline: 0.124)\n")
cat("psi_ss at mean flow:", round(plogis(alpha_ss_post), 3),
    "(baseline: 0.374)\n")
cat("beta_geo:", round(beta_geo_post, 3),
    "-> flow effect on Geo routing\n")
cat("beta_ss:", round(beta_ss_post, 3),
    "-> flow effect on SS routing (expected negative)\n")

# Predict routing probabilities across range of observed flows
# Back-transform to probability scale for plotting
flow_range <- seq(min(flow_occ2[flow_occ2 != 0]),
                  max(flow_occ2),
                  length.out = 100)

psi_geo_pred <- plogis(alpha_geo_post + beta_geo_post * flow_range)
psi_ss_pred  <- plogis(alpha_ss_post  + beta_ss_post  * flow_range)

# Convert standardized flow back to cfs for x-axis
flow_range_cfs <- flow_range * flow_sd + flow_mean

par(mfrow = c(1, 2))

plot(flow_range_cfs, psi_ss_pred,
     type = "l", col = "forestgreen", lwd = 2,
     ylim = c(0, 0.8),
     xlab = "30-day mean flow at junction (cfs)",
     ylab = "Routing probability",
     main = "Steamboat/Sutter routing vs flow")
abline(h = 0.374, lty = 2, col = "gray50")
legend("topright", c("Flow model (posterior mean)", "Baseline estimate"),
       lty = c(1, 2), col = c("forestgreen", "gray50"), cex = 0.8)

plot(flow_range_cfs, psi_geo_pred,
     type = "l", col = "orange", lwd = 2,
     ylim = c(0, 0.4),
     xlab = "30-day mean flow at junction (cfs)",
     ylab = "Routing probability",
     main = "Georgiana routing vs flow")
abline(h = 0.124, lty = 2, col = "gray50")
legend("topright", c("Flow model (posterior mean)", "Baseline estimate"),
       lty = c(1, 2), col = c("orange", "gray50"), cex = 0.8)


# Predictions at observed flow percentiles
flow_quantiles_std <- quantile(flow_occ2[flow_occ2 != 0], 
                               probs = c(0.25, 0.50, 0.75))
flow_quantiles_cfs <- flow_quantiles_std * flow_sd + flow_mean

cat("Predictions at observed flow percentiles:\n")
cat("Flow (cfs):", round(flow_quantiles_cfs), "\n")
cat("psi_geo:   ", round(plogis(alpha_geo_post + beta_geo_post * flow_quantiles_std), 3), "\n")
cat("psi_ss:    ", round(plogis(alpha_ss_post  + beta_ss_post  * flow_quantiles_std), 3), "\n")
