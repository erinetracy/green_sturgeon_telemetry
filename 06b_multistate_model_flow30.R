#==============================================================================
# GREEN STURGEON UPSTREAM MIGRATION MULTISTATE MODEL - 30-DAY ANTECEDENT FLOW
# Script: 06b_multistate_model_flow_30day.R
# Author: Erin Tracy
# Last updated: April 2026
#
# PURPOSE:
# Parallel version of 06_multistate_model_flow.R using 30-day mean antecedent
# flow before Benicia passage as the flow covariate instead of same-day
# junction flow. Intended as a sensitivity analysis and model comparison.
#
# COMPARISON WITH 06_multistate_model_flow.R:
# The same-day junction flow model (06) asks:
#   "Does the flow at the junction at the moment of the routing decision
#    predict which route the fish takes?"
#   -> Follows Perry et al. 2018 FPT.flow approach exactly
#   -> Captures local hydraulic conditions at the decision point
#
# This 30-day antecedent flow model (06b) asks:
#   "Does the broader hydrological context in the 30 days before the fish
#    entered the freshwater system predict which route it takes?"
#   -> Captures systemic flow conditions driving migration behavior
#   -> Single flow value per fish (not junction-specific)
#   -> More closely related to migration initiation conditions
#
# BIOLOGICAL RATIONALE FOR 30-DAY MEAN:
# ANOVA comparing temporal scales of flow as predictors of route choice showed
# 30-day mean flow before Benicia passage produced the strongest separation
# among route groups (F=2.13, p=0.097) compared to shorter windows:
#   3-day:  F=1.51, p=0.21
#   7-day:  F=1.54, p=0.21
#   15-day: F=1.64, p=0.18
#   30-day: F=2.13, p=0.097 <- strongest signal
# Note: p=0.097 is not conventionally significant and the ANOVA did not
# account for detection probability. The 30-day window is chosen as the
# best available temporal scale based on this exploratory analysis.
#
# The biological interpretation is that fish staging in the estuary integrate
# flow conditions over weeks before committing to a route, rather than
# responding to instantaneous flow pulses at the junction. Fish that entered
# during sustained high-flow periods may have different route selection
# behavior than fish entering during low-flow periods regardless of same-day
# junction conditions.
#
# IMPORTANT CONCEPTUAL NOTE:
# 30-day antecedent flow before Benicia is primarily a MIGRATION INITIATION
# variable — it describes conditions when fish committed to beginning their
# upstream migration. Using it as a routing predictor conflates two separate
# behavioral decisions: (1) when to start migrating and (2) which route to
# take. The cleaner use of this metric is in a separate migration initiation
# analysis (time-to-event model) using the extended 2007-2021 Benicia dataset.
# This model is run here as a sensitivity analysis only.
#
# KEY DIFFERENCES FROM PERRY 2018:
# 1. Perry used Freeport same-day FPT.flow; we use 30-day mean before Benicia
# 2. This is a single flow value per fish not junction-specific
# 3. All other differences same as 06_multistate_model_flow.R
#
# SIMPLIFIED SURVIVAL STRUCTURE:
# Same as baseline and 06 — single S_sac for all Sacramento reaches,
# separate S_geo and S_ss. Justified by all reach-specific estimates
# being 0.986-0.993 in baseline model.
#
# LOGISTIC REGRESSION ON ROUTING:
# Both routing probabilities use the same 30-day antecedent flow covariate:
#   logit(psi_geo[i]) = alpha_geo + beta_geo * flow_30day[i]
#   logit(psi_ss[i])  = alpha_ss  + beta_ss  * flow_30day[i]
# where flow_30day[i] is the mean Rio Vista discharge in the 30 days before
# fish i first passed Benicia/Carquinez, standardized to mean=0, sd=1.
#
# PRIORS: Same as 06_multistate_model_flow.R
#   alpha: Normal(0, sd=1.5) on logit scale
#   beta:  Normal(0, sd=1.0) on logit scale
#
# INPUTS:
#   - gs_multistate_data.RData (from 04_multistate_data_prep.R)
#   - daily_tidalfilter_riovista.csv: USGS tidally-filtered daily discharge
#   - gs_mcmc_flow_run1.RData: same-day flow model results for comparison
#
# OUTPUTS:
#   - gs_mcmc_flow_30day_run1.RData: MCMC posterior samples
#   - gs_mcmc_traces_flow_30day_run1.pdf: trace plots
#
# MODEL STRUCTURE:
#   States: 1=Sacramento, 2=Georgiana, 3=Steamboat/Sutter,
#           4=Dead (absorbing), 5=Failed migration (absorbing)
#   Fish: 221 (134 up_complete, 83 up_incomplete, 4 incomplete_dead)
#
# NOTE ON flow_occ2/flow_occ3 AS DATA:
# These are fixed observed covariates with no likelihood distribution.
# They are passed as data (not constants) so they can be updated via
# nimMod$setData() without recompiling — useful for sensitivity analyses
# comparing different flow metrics. NIMBLE treats undeclared data nodes
# as fixed observed values when they only appear in deterministic (<-) 
# expressions, which is the case here.
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
# SECTION 3: BUILD 30-DAY ANTECEDENT FLOW COVARIATE
#
# APPROACH:
# Assign each fish the mean Rio Vista discharge in the 30 days BEFORE it
# first passed Benicia/Carquinez (migration entry point). This is a single
# flow value per fish representing the hydrological context during the period
# when the fish was staging in the estuary and committing to migrate.
#
# This differs from 06_multistate_model_flow.R which used same-day flow at
# each junction (two values per fish, junction-specific). Here we use one
# value per fish applied to both routing decisions (psi_geo and psi_ss).
#
# RATIONALE FOR SAME COVARIATE FOR BOTH JUNCTIONS:
# Unlike same-day junction flow (which should logically differ between the
# Geo junction at occ2 and SS junction at occ3), the 30-day antecedent
# window before Benicia is the same for both routing decisions for a given
# fish. This is a simplification but is appropriate given:
# 1. The 30-day window predates both junctions for all fish
# 2. It represents the migration initiation context not junction-specific cues
# 3. Using the same covariate for both keeps the model parsimonious
#
# WINDOW DEFINITION:
# 30 days = days [-30, -1] before date of first Benicia/Carquinez detection
# Day 0 (Benicia passage date) NOT included — fish already committed at this
# point. We want the conditions that drove the decision to start migrating.
#
# IMPUTATION:
# Fish with no Benicia date get flow_30day = 0 (standardized mean)
# NIMBLE requires fully numeric constants — NA not permitted
#==============================================================================

# Load Rio Vista flow data
rv_flow <- read.csv("C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv")

# Clean flow data — need extra buffer before WY2007 for 30-day windows
rv_flow_clean <- rv_flow %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_cfs = value) %>%
  filter(!is.na(flow_cfs),
         date >= as.Date("2006-09-01"),   # buffer for early WY2007 migrants
         date <= as.Date("2017-09-30")) %>%
  arrange(date)

cat("Flow data days:", nrow(rv_flow_clean), "\n")
cat("Missing flow values:", sum(is.na(rv_flow_clean$flow_cfs)), "\n")

# Check for gaps in flow record
rv_flow_clean %>%
  mutate(gap_after = as.numeric(lead(date) - date)) %>%
  filter(gap_after > 1) %>%
  dplyr::select(date, flow_cfs, gap_after) %>%
  mutate(gap_start = date + 1, gap_end = date + gap_after - 1) %>%
  print()

# Fill gaps with linear interpolation
library(zoo)
rv_flow_complete <- rv_flow_clean %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_cfs = zoo::na.approx(flow_cfs, na.rm = FALSE))

cat("Days after interpolation:", nrow(rv_flow_complete), "\n")
cat("Remaining NAs:", sum(is.na(rv_flow_complete$flow_cfs)), "\n")

# Get Benicia passage date for each fish — anchor for 30-day window
benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    benicia_date = as.Date(min(first_detection)),
    .groups = "drop"
  )

cat("Fish with Benicia dates:", nrow(benicia_dates), "\n")

# Calculate 30-day mean flow before Benicia passage
fish_flow <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    flow_30day_raw = {
      if(is.na(benicia_date)) {
        NA_real_
      } else {
        window_flow <- rv_flow_complete$flow_cfs[
          rv_flow_complete$date >= (benicia_date - 30) &
            rv_flow_complete$date <   benicia_date
        ]
        if(length(window_flow) == 0) NA_real_ else mean(window_flow, na.rm = TRUE)
      }
    },
    n_days_in_window = {
      if(is.na(benicia_date)) 0L else
        sum(rv_flow_complete$date >= (benicia_date - 30) &
              rv_flow_complete$date <   benicia_date &
              !is.na(rv_flow_complete$flow_cfs))
    }
  ) %>%
  ungroup()

cat("Missing flow_30day:", sum(is.na(fish_flow$flow_30day_raw)), "\n")
cat("Window size summary:\n")
print(table(fish_flow$n_days_in_window))

# Standardize to mean=0, sd=1
flow_mean <- mean(fish_flow$flow_30day_raw, na.rm = TRUE)
flow_sd   <- sd(fish_flow$flow_30day_raw,   na.rm = TRUE)

fish_flow <- fish_flow %>%
  mutate(flow_30day_std = (flow_30day_raw - flow_mean) / flow_sd)

cat("\n30-day antecedent flow summary:\n")
cat("Mean raw flow:", round(flow_mean, 0), "cfs\n")
cat("SD raw flow:", round(flow_sd, 0), "cfs\n")
cat("Standardized range:",
    round(range(fish_flow$flow_30day_std, na.rm = TRUE), 2), "\n")

# Verify row order matches detection_history
identical(fish_flow$animal_id, detection_history$animal_id)
identical(as.numeric(fish_flow$water_year),
          as.numeric(detection_history$water_year))

# Extract single flow vector for NIMBLE
# Used for both psi_geo and psi_ss — same antecedent conditions for both
flow_30day <- fish_flow$flow_30day_std
flow_30day[is.na(flow_30day)] <- 0   # impute mean for fish missing Benicia date

cat("\nFinal flow vector:\n")
cat("n fish:", length(flow_30day), "\n")
cat("NAs imputed to 0:", sum(is.na(fish_flow$flow_30day_std)), "\n")
cat("Range:", round(range(flow_30day), 2), "\n")

hist(flow_30day[!is.na(fish_flow$flow_30day_std)],
     main = "30-day mean flow before Benicia (standardized, non-imputed)",
     xlab = "Standardized flow",
     col  = "steelblue")
abline(v = 0, lty = 2, col = "red")

#==============================================================================
# SECTION 4: VERIFY FLOW COVARIATE MAKES BIOLOGICAL SENSE
# Check that 30-day antecedent flow varies meaningfully across water years
# and that observed patterns match exploratory analysis:
#   - SS use DECREASES with flow (confirmed by route quartile plot)
#   - Geo is relatively flow-independent
# Also compare to same-day junction flow from 06 to understand differences
#==============================================================================

# Mean 30-day flow by water year - 2016-2017 should be highest (wet years)
cat("Mean 30-day antecedent flow by water year:\n")
fish_flow %>%
  filter(!is.na(flow_30day_std)) %>%
  group_by(water_year) %>%
  dplyr::summarise(
    n_fish        = n(),
    mean_flow_std = round(mean(flow_30day_std), 2),
    mean_flow_cfs = round(mean(flow_30day_raw), 0),
    .groups = "drop"
  ) %>%
  arrange(water_year) %>%
  print()

# Flow by route taken - check expected direction
# SS fish should show LOWER antecedent flow than average
fish_flow %>%
  left_join(
    detection_history %>% dplyr::select(animal_id, water_year, occ_2, occ_3),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    route = case_when(
      occ_2 == 2 ~ "Georgiana",
      occ_3 == 4 ~ "Steamboat/Sutter",
      TRUE       ~ "Sacramento"
    )
  ) %>%
  group_by(route) %>%
  dplyr::summarise(
    n                = n(),
    mean_flow_30day  = round(mean(flow_30day_std, na.rm = TRUE), 2),
    sd_flow_30day    = round(sd(flow_30day_std,   na.rm = TRUE), 2),
    .groups          = "drop"
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
  # Each fish gets its own routing probability based on 30-day antecedent flow
  # Same flow value used for both psi_geo and psi_ss — single covariate per fish
  # This differs from 06 where junction-specific flow was used
  for(i in 1:nfish){
    
    # Georgiana routing probability for fish i
    # Uses 30-day mean flow before fish i passed Benicia
    psi_geo[i] <- ilogit(alpha_geo + beta_geo * flow_30day[i])
    
    # Steamboat/Sutter routing probability for fish i
    # Same antecedent flow covariate — represents migration initiation context
    psi_ss[i]  <- ilogit(alpha_ss  + beta_ss  * flow_30day[i])
    
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
  data      = list(ch_mat = ch_mat_nimble),
  constants = list(
    nfish      = nrow(ch_mat_nimble),
    flow_30day = flow_30day
  )
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

# Short test run first
# mcmc_test_flow <- runMCMC(CMCMC_flow, niter = 1000, nchains = 1,
#                           nburnin = 100, thin = 1,
#                           inits = list(inits_fn_flow()),
#                           samplesAsCodaMCMC = TRUE)
# MCMCsummary(mcmc_test_flow, round = 3)

# Full run
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

save(mcmc_out_flow, flow_30day, fish_flow,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_30day_run1.RData")

MCMCtrace(mcmc_out_flow,
          params   = c("alpha_geo", "alpha_ss", "beta_geo", "beta_ss",
                       "phi_fail", "lambda"),
          pdf      = TRUE,
          filename = "gs_mcmc_traces_flow_30day_run1",
          ind      = TRUE,
          Rhat     = TRUE,
          n.eff    = TRUE)

#==============================================================================
# SECTION 10: INTERPRET AND COMPARE RESULTS
# Compare 30-day antecedent flow model (this script) with
# same-day junction flow model (06_multistate_model_flow.R)
# Run after both MCMC runs are complete
#==============================================================================

# Load same-day flow model results for comparison
load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_run1.RData")
mcmc_sameday <- mcmc_out_flow

# Load 30-day model results
load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_30day_run1.RData")
mcmc_30day <- mcmc_out_flow
