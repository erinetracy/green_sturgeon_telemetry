#==============================================================================
# TEMPERATURE AND FLOW COVARIATE EXTRACTION
# For multistate model temperature sensitivity analysis
# Data: daily_tidalfilter_riovista_temp_flow.csv
# Starts: 2007-05-30 (incomplete WY2007)
# Temperature in Celsius used throughout
#==============================================================================

library(dplyr)
library(tidyr)
library(zoo)
library(lubridate)

#==============================================================================
# SECTION 1: LOAD AND CLEAN TEMPERATURE DATA
#==============================================================================

temp_flow_raw <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/daily_tidalfilter_riovista_temp_flow.csv"
)

# check structure
cat("Columns:", paste(names(temp_flow_raw), collapse = ", "), "\n")
cat("Rows:", nrow(temp_flow_raw), "\n")
cat("Date range:", range(temp_flow_raw$date), "\n")

#clean
temp_flow_clean <- temp_flow_raw %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  dplyr::select(date, temp_c = daily_temp_c, flow_cfs) %>%
  filter(!is.na(date),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(as.Date("2006-09-01"),
                             as.Date("2017-09-30"),
                             by = "day")) %>%
  mutate(
    temp_c   = zoo::na.approx(temp_c,   na.rm = FALSE),
    flow_cfs = zoo::na.approx(flow_cfs, na.rm = FALSE)
  )

cat("Total days:", nrow(temp_flow_clean), "\n")
cat("Missing temp_c:", sum(is.na(temp_flow_clean$temp_c)), "\n")
cat("Missing flow_cfs:", sum(is.na(temp_flow_clean$flow_cfs)), "\n")
cat("Date range:", as.character(range(temp_flow_clean$date)), "\n")


# check WY2007 coverage specifically
wy2007_temp <- temp_flow_clean %>%
  filter(date >= as.Date("2006-10-01"),
         date <= as.Date("2007-09-30"))
cat("\nWY2007 temperature coverage:\n")
cat("Total days:", nrow(wy2007_temp), "\n")
cat("Days with temp data:", sum(!is.na(wy2007_temp$temp_c)), "\n")
cat("Days missing:", sum(is.na(wy2007_temp$temp_c)), "\n")
cat("First temp date in WY2007:", 
    as.character(min(wy2007_temp$date[!is.na(wy2007_temp$temp_c)])), "\n")

#load detection events
load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_multistate_data.RData")

# confirm key objects loaded
cat("Objects loaded:\n")
cat("detection_history:", nrow(detection_history), "rows\n")
cat("model1_events:", nrow(model1_events), "rows\n")

# now rebuild benicia dates
benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    benicia_date = as.Date(min(first_detection)),
    .groups = "drop"
  )

# check WY2007 fish coverage
benicia_dates %>%
  filter(water_year == 2007) %>%
  arrange(benicia_date) %>%
  mutate(
    has_same_day_temp  = benicia_date >= as.Date("2007-05-30"),
    has_30day_window   = benicia_date >= as.Date("2007-06-29")
  ) %>%
  print()

#rebuild
occ2_dates <- model1_events %>%
  filter(occasion == 2) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    arr_date_occ2 = as.Date(min(first_detection)),
    .groups = "drop"
  )

occ3_dates <- model1_events %>%
  filter(occasion == 3) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(
    arr_date_occ3 = as.Date(min(first_detection)),
    .groups = "drop"
  )

cat("occ2 dates:", nrow(occ2_dates), "\n")
cat("occ3 dates:", nrow(occ3_dates), "\n")

#==============================================================================
# SECTION 2: GET ARRIVAL DATES
# Reuse from model 09 if already in environment,
# otherwise rebuild here
#==============================================================================

# benicia dates for 30-day window
if(!exists("benicia_dates")){
  benicia_dates <- model1_events %>%
    filter(receiver_group %in% c("benicia", "carquinez")) %>%
    group_by(animal_id, water_year) %>%
    dplyr::summarise(
      benicia_date = as.Date(min(first_detection)),
      .groups = "drop"
    )
}

# occ2 arrival dates for same-day covariate
if(!exists("occ2_dates")){
  occ2_dates <- model1_events %>%
    filter(occasion == 2) %>%
    group_by(animal_id, water_year) %>%
    dplyr::summarise(
      arr_date_occ2 = as.Date(min(first_detection)),
      .groups = "drop"
    )
}

# occ3 arrival dates for same-day covariate
if(!exists("occ3_dates")){
  occ3_dates <- model1_events %>%
    filter(occasion == 3) %>%
    group_by(animal_id, water_year) %>%
    dplyr::summarise(
      arr_date_occ3 = as.Date(min(first_detection)),
      .groups = "drop"
    )
}

cat("Benicia dates:", nrow(benicia_dates), "\n")
cat("Occ2 dates:", nrow(occ2_dates), "\n")
cat("Occ3 dates:", nrow(occ3_dates), "\n")

#==============================================================================
# SECTION 3: EXTRACT TEMPERATURE COVARIATES
# Same structure as flow:
# temp_occ2_raw:    same-day mean temp at Rio Vista on occ2 arrival date
# temp_occ3_raw:    same-day mean temp at Rio Vista on occ3 arrival date
# temp_30day_raw:   30-day mean temp before Benicia passage
#==============================================================================

fish_temp <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(occ2_dates,    by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,    by = c("animal_id", "water_year")) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # Same-day temperature at occ2 arrival (Rio Vista)
    temp_occ2_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(temp_flow_clean$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else temp_flow_clean$temp_c[idx]
      }
    },
    
    # Same-day temperature at occ3 arrival (Rio Vista)
    # Note: using Rio Vista temp for both since no gauge-specific
    # temp data at GES — same rationale as 30-day flow
    temp_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(temp_flow_clean$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else temp_flow_clean$temp_c[idx]
      }
    },
    
    # 30-day mean temperature before Benicia passage
    temp_30day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 30) &
            temp_flow_clean$date <   benicia_date]
        # exclude NAs from window (WY2007 early fish)
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # count days with temp data in 30-day window
    # useful for flagging fish with incomplete windows
    temp_30day_n = {
      if(is.na(benicia_date)) 0L else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 30) &
            temp_flow_clean$date <   benicia_date]
        sum(!is.na(window))
      }
    }
    
  ) %>%
  ungroup()

cat("\n=== TEMPERATURE COVARIATE SUMMARY ===\n")
cat("Missing temp_occ2:", sum(is.na(fish_temp$temp_occ2_raw)), "\n")
cat("Missing temp_occ3:", sum(is.na(fish_temp$temp_occ3_raw)), "\n")
cat("Missing temp_30day:", sum(is.na(fish_temp$temp_30day_raw)), "\n")

cat("\n30-day window completeness:\n")
print(table(fish_temp$temp_30day_n))
cat("Fish with full 30-day window (30 days):",
    sum(fish_temp$temp_30day_n == 30, na.rm = TRUE), "\n")
cat("Fish with partial window (<30 days):",
    sum(fish_temp$temp_30day_n > 0 &
          fish_temp$temp_30day_n < 30, na.rm = TRUE), "\n")

cat("\nRaw temperature summaries:\n")
cat("Occ2 same-day temp (C):\n")
print(summary(fish_temp$temp_occ2_raw))
cat("Occ3 same-day temp (C):\n")
print(summary(fish_temp$temp_occ3_raw))
cat("30-day mean temp (C):\n")
print(summary(fish_temp$temp_30day_raw))

# check WY2007 fish specifically
cat("\nWY2007 fish temperature coverage:\n")
fish_temp %>%
  filter(water_year == 2007) %>%
  dplyr::select(animal_id, benicia_date,
                temp_30day_raw, temp_30day_n) %>%
  print()

#==============================================================================
# SECTION 4: STANDARDIZE TEMPERATURE COVARIATES
#==============================================================================

temp_occ2_mean  <- mean(fish_temp$temp_occ2_raw,  na.rm = TRUE)
temp_occ2_sd    <- sd(fish_temp$temp_occ2_raw,    na.rm = TRUE)
temp_occ3_mean  <- mean(fish_temp$temp_occ3_raw,  na.rm = TRUE)
temp_occ3_sd    <- sd(fish_temp$temp_occ3_raw,    na.rm = TRUE)
temp_30day_mean <- mean(fish_temp$temp_30day_raw, na.rm = TRUE)
temp_30day_sd   <- sd(fish_temp$temp_30day_raw,   na.rm = TRUE)

fish_temp <- fish_temp %>%
  mutate(
    temp_occ2_std  = (temp_occ2_raw  - temp_occ2_mean)  / temp_occ2_sd,
    temp_occ3_std  = (temp_occ3_raw  - temp_occ3_mean)  / temp_occ3_sd,
    temp_30day_std = (temp_30day_raw - temp_30day_mean) / temp_30day_sd
  )

cat("\n=== STANDARDIZED TEMPERATURE SUMMARIES ===\n")
cat("Mean raw occ2 temp:", round(temp_occ2_mean, 2), "C\n")
cat("SD raw occ2 temp:",   round(temp_occ2_sd, 2),   "C\n")
cat("Mean raw occ3 temp:", round(temp_occ3_mean, 2), "C\n")
cat("SD raw occ3 temp:",   round(temp_occ3_sd, 2),   "C\n")
cat("Mean raw 30-day temp:", round(temp_30day_mean, 2), "C\n")
cat("SD raw 30-day temp:",   round(temp_30day_sd, 2),   "C\n")

print(summary(fish_temp$temp_occ2_std))
print(summary(fish_temp$temp_occ3_std))
print(summary(fish_temp$temp_30day_std))

#==============================================================================
# SECTION 5: CHECK CORRELATIONS WITH FLOW
# Critical before deciding whether to include temp alongside flow
#==============================================================================

cat("\n=== CORRELATION BETWEEN FLOW AND TEMPERATURE ===\n")

# load flow data if not in environment
if(!exists("fish_flow_09")){
  load("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_mcmc_flow_09_run1.RData")
}

# join flow and temp for correlation check
flow_temp_check <- fish_temp %>%
  dplyr::select(animal_id, water_year,
                temp_occ2_raw, temp_occ3_raw, temp_30day_raw) %>%
  left_join(
    fish_flow_09 %>%
      dplyr::select(animal_id, water_year,
                    flow_occ2_raw, flow_occ3_raw, flow_30day_raw),
    by = c("animal_id", "water_year")
  )

cat("Same-day flow vs same-day temp (occ2):\n")
cat("r =", round(cor(flow_temp_check$flow_occ2_raw,
                     flow_temp_check$temp_occ2_raw,
                     use = "complete.obs"), 3), "\n")

cat("Same-day flow vs same-day temp (occ3):\n")
cat("r =", round(cor(flow_temp_check$flow_occ3_raw,
                     flow_temp_check$temp_occ3_raw,
                     use = "complete.obs"), 3), "\n")

cat("30-day flow vs 30-day temp:\n")
cat("r =", round(cor(flow_temp_check$flow_30day_raw,
                     flow_temp_check$temp_30day_raw,
                     use = "complete.obs"), 3), "\n")

cat("Occ2 temp vs occ3 temp:\n")
cat("r =", round(cor(flow_temp_check$temp_occ2_raw,
                     flow_temp_check$temp_occ3_raw,
                     use = "complete.obs"), 3), "\n")

cat("Same-day temp vs 30-day temp:\n")
cat("r =", round(cor(flow_temp_check$temp_occ2_raw,
                     flow_temp_check$temp_30day_raw,
                     use = "complete.obs"), 3), "\n")

#==============================================================================
# SECTION 6: EXTRACT VECTORS FOR NIMBLE
#==============================================================================

# verify row order matches detection_history
identical(fish_temp$animal_id, detection_history$animal_id)
identical(as.numeric(fish_temp$water_year),
          as.numeric(detection_history$water_year))

# extract standardized vectors — impute 0 for missing
temp_occ2_vec  <- fish_temp$temp_occ2_std
temp_occ3_vec  <- fish_temp$temp_occ3_std
temp_30day_vec <- fish_temp$temp_30day_std

temp_occ2_vec[is.na(temp_occ2_vec)]   <- 0
temp_occ3_vec[is.na(temp_occ3_vec)]   <- 0
temp_30day_vec[is.na(temp_30day_vec)] <- 0

cat("\nFinal temperature vectors:\n")
cat("temp_occ2 range:",  round(range(temp_occ2_vec),  2), "\n")
cat("temp_occ3 range:",  round(range(temp_occ3_vec),  2), "\n")
cat("temp_30day range:", round(range(temp_30day_vec), 2), "\n")
cat("NAs imputed to 0 (occ2):",
    sum(is.na(fish_temp$temp_occ2_std)), "\n")
cat("NAs imputed to 0 (occ3):",
    sum(is.na(fish_temp$temp_occ3_std)), "\n")
cat("NAs imputed to 0 (30day):",
    sum(is.na(fish_temp$temp_30day_std)), "\n")

#==============================================================================
# SECTION 7: SAVE
#==============================================================================

save(fish_temp,
     temp_occ2_vec, temp_occ3_vec, temp_30day_vec,
     temp_occ2_mean, temp_occ2_sd,
     temp_occ3_mean, temp_occ3_sd,
     temp_30day_mean, temp_30day_sd,
     flow_temp_check,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/temp_covariates.RData")

cat("\nTemperature covariates saved\n")
cat("Share the correlation output before proceeding to model 11\n")


#TESTING WHICH TEMP TIME SCALE IS BEST
##################################################################
# first build all temporal scale temperature covariates
# then screen with simple logistic regression before running NIMBLE

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)

# build temperature windows at multiple temporal scales
# anchored to Benicia passage date like the 30-day flow covariate

fish_temp_multiscale <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  left_join(occ2_dates,   by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,   by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # same-day temperature at occ2 arrival
    temp_1day_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(temp_flow_clean$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else temp_flow_clean$temp_c[idx]
      }
    },
    
    # 7-day mean before Benicia
    temp_7day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 7) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 14-day mean before Benicia
    temp_14day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 14) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 30-day mean before Benicia
    temp_30day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 30) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 45-day mean before Benicia
    temp_45day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 45) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 60-day mean before Benicia
    temp_60day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 60) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    }
    
  ) %>%
  ungroup() %>%
  # join route outcomes for screening
  left_join(
    detection_history %>%
      dplyr::select(animal_id, water_year, occ_2, occ_3, status),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    took_geo    = as.integer(occ_2 == 2),
    took_ss     = as.integer(occ_3 == 4),
    failed      = as.integer(status %in% c("up_incomplete", "incomplete_dead"))
  )

cat("Multiscale temperature data built\n")
cat("Rows:", nrow(fish_temp_multiscale), "\n")

#==============================================================================
# SCREEN: SIMPLE LOGISTIC REGRESSION AT EACH TEMPORAL SCALE
# Fit logistic regression of each outcome on each temperature covariate
# Compare pseudo-R2 and significance as a quick screen before NIMBLE
#==============================================================================

library(pROC)

# temperature scales to test
temp_vars <- c("temp_1day_raw", "temp_7day_raw", "temp_14day_raw",
               "temp_30day_raw", "temp_45day_raw", "temp_60day_raw")
temp_labels <- c("Same-day", "7-day", "14-day", "30-day", "45-day", "60-day")

# outcomes to test
outcomes <- c("took_geo", "took_ss", "failed")
outcome_labels <- c("Georgiana routing", "SS routing", "Migration failure")

# storage
screen_results <- data.frame()

for(o in seq_along(outcomes)){
  outcome_var <- outcomes[o]
  
  for(t in seq_along(temp_vars)){
    temp_var <- temp_vars[t]
    
    # subset to complete cases for this combination
    df_sub <- fish_temp_multiscale %>%
      dplyr::select(all_of(c(outcome_var, temp_var))) %>%
      filter(complete.cases(.))
    
    # for SS routing exclude Georgiana fish
    if(outcome_var == "took_ss"){
      df_sub <- fish_temp_multiscale %>%
        filter(occ_2 != 2) %>%
        dplyr::select(all_of(c(outcome_var, temp_var))) %>%
        filter(complete.cases(.))
    }
    
    if(nrow(df_sub) < 10) next
    
    # fit logistic regression
    formula_str <- paste(outcome_var, "~", temp_var)
    fit <- glm(as.formula(formula_str),
               data   = df_sub,
               family = binomial(link = "logit"))
    
    s <- summary(fit)
    
    # McFadden pseudo R2
    null_ll <- logLik(glm(as.formula(paste(outcome_var, "~ 1")),
                          data = df_sub, family = binomial))
    full_ll <- logLik(fit)
    pseudo_r2 <- 1 - as.numeric(full_ll) / as.numeric(null_ll)
    
    # AUC
    auc_val <- tryCatch({
      as.numeric(auc(roc(df_sub[[outcome_var]],
                         fitted(fit), quiet = TRUE)))
    }, error = function(e) NA)
    
    screen_results <- rbind(screen_results, data.frame(
      outcome      = outcome_labels[o],
      temp_scale   = temp_labels[t],
      n            = nrow(df_sub),
      beta         = round(coef(fit)[2], 3),
      se           = round(s$coefficients[2, 2], 3),
      z            = round(s$coefficients[2, 3], 3),
      p_value      = round(s$coefficients[2, 4], 4),
      pseudo_r2    = round(pseudo_r2, 4),
      auc          = round(auc_val, 3)
    ))
  }
}

cat("\n=== TEMPORAL SCALE SCREENING RESULTS ===\n")
print(screen_results, row.names = FALSE)

# plot pseudo R2 by temporal scale and outcome
p_screen <- ggplot(screen_results,
                   aes(x = temp_scale, y = pseudo_r2,
                       color = outcome, group = outcome)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Georgiana routing" = "#56B4E9",
    "SS routing"        = "#009E73",
    "Migration failure" = "#D55E00"
  )) +
  scale_x_discrete(limits = temp_labels) +
  labs(
    x     = "Temperature averaging window",
    y     = "McFadden pseudo R²",
    color = "Outcome",
    title = "Temporal scale screening — temperature covariates",
    subtitle = "Higher pseudo R² indicates better explanatory power at that temporal scale"
  ) +
  theme_bw(base_size = 11, base_family = "Times New Roman") +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

print(p_screen)

ggsave(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/figures/temp_scale_screening.png",
  plot = p_screen, width = 8, height = 5, dpi = 300)

cat("\nScreening plot saved\n")
cat("Use pseudo R2 and AUC to select temporal scale before running NIMBLE\n")






#==============================================================================
# TURBIDITY COVARIATE EXTRACTION AND SCREENING
# Data: RioVista_Turbidity_Hourly.xlsx
# Hourly NTU at Rio Vista — aggregate to daily, check correlations with flow
#==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(zoo)
library(lubridate)
library(ggplot2)

#==============================================================================
# SECTION 1: LOAD AND CLEAN TURBIDITY DATA
#==============================================================================

turb_raw <- read_xlsx(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/RioVista_Turbidity_Hourly.xlsx"
)

# check structure
cat("Columns:", paste(names(turb_raw), collapse = ", "), "\n")
cat("Rows:", nrow(turb_raw), "\n")
cat("Class of each column:\n")
print(sapply(turb_raw, class))
cat("\nFirst few rows:\n")
print(head(turb_raw, 10))

# check DATA_FLAG values
cat("\nDATA_FLAG values:\n")
print(table(turb_raw$DATA_FLAG, useNA = "always"))

# check VALUE range
cat("\nVALUE summary:\n")
print(summary(turb_raw$VALUE))

# clean turbidity — convert VALUE to numeric and aggregate to daily
turb_daily <- turb_raw %>%
  mutate(
    datetime = as.POSIXct(`DATE TIME`, format = "%Y-%m-%d %H:%M",
                          tz = "America/Los_Angeles"),
    date  = as.Date(datetime),
    turb_ntu = as.numeric(VALUE)
  ) %>%
  filter(!is.na(turb_ntu),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  group_by(date) %>%
  dplyr::summarise(
    turb_daily_mean = mean(turb_ntu, na.rm = TRUE),
    turb_daily_n    = n(),   # hours of data per day
    turb_daily_max  = max(turb_ntu, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date)

cat("Daily turbidity records:", nrow(turb_daily), "\n")
cat("Date range:", as.character(range(turb_daily$date)), "\n")
cat("Days with <12 hourly obs:", sum(turb_daily$turb_daily_n < 12), "\n")
cat("\nTurbidity summary (NTU):\n")
print(summary(turb_daily$turb_daily_mean))
print(summary(turb_daily$turb_daily_max))

# check for extreme outliers — turbidity can spike during storms
hist(turb_daily$turb_daily_mean, breaks = 50,
     main = "Daily mean turbidity distribution",
     xlab = "NTU")

# what are the highest values?
turb_daily %>%
  arrange(desc(turb_daily_mean)) %>%
  head(20) %>%
  print()

# log transform turbidity
turb_daily <- turb_daily %>%
  mutate(log_turb = log(turb_daily_mean))

cat("Log turbidity summary:\n")
print(summary(turb_daily$log_turb))

# complete to full daily record with NAs for missing dates
turb_complete <- tidyr::complete(
  turb_daily,
  date = seq(as.Date("2006-09-01"),
             as.Date("2017-09-30"),
             by = "day")
) %>%
  # interpolate short gaps only — don't interpolate across the long
  # WY2007/early WY2008 gap since there's no surrounding data
  mutate(log_turb = zoo::na.approx(log_turb, maxgap = 7, na.rm = FALSE))

cat("\nTotal days:", nrow(turb_complete), "\n")
cat("Missing log_turb:", sum(is.na(turb_complete$log_turb)), "\n")

# check which fish are affected by missing turbidity
# by checking benicia dates against data availability
benicia_dates %>%
  mutate(
    has_turb_30day = benicia_date >= as.Date("2008-03-02"),  # Feb 1 + 30 days
    has_turb_7day  = benicia_date >= as.Date("2008-02-08"),  # Feb 1 + 7 days
    wy_group = case_when(
      water_year == 2007 ~ "WY2007 (no data)",
      water_year == 2008 ~ "WY2008 (partial)",
      TRUE ~ "WY2009+ (full data)"
    )
  ) %>%
  group_by(wy_group, has_turb_30day) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  print()

# also check WY2008 fish specifically
benicia_dates %>%
  filter(water_year == 2008) %>%
  arrange(benicia_date) %>%
  mutate(
    has_30day_window = benicia_date >= as.Date("2008-03-02")
  ) %>%
  print()

# build multiscale turbidity covariates for screening
fish_turb_multiscale <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  left_join(occ2_dates,   by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,   by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    # same-day log turbidity at occ2 arrival
    turb_1day_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(turb_complete$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else turb_complete$log_turb[idx]
      }
    },
    
    # 7-day mean log turbidity before Benicia
    turb_7day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- turb_complete$log_turb[
          turb_complete$date >= (benicia_date - 7) &
            turb_complete$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 14-day mean log turbidity before Benicia
    turb_14day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- turb_complete$log_turb[
          turb_complete$date >= (benicia_date - 14) &
            turb_complete$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 30-day mean log turbidity before Benicia
    turb_30day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- turb_complete$log_turb[
          turb_complete$date >= (benicia_date - 30) &
            turb_complete$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 45-day mean log turbidity before Benicia
    turb_45day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- turb_complete$log_turb[
          turb_complete$date >= (benicia_date - 45) &
            turb_complete$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    },
    
    # 60-day mean log turbidity before Benicia
    turb_60day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- turb_complete$log_turb[
          turb_complete$date >= (benicia_date - 60) &
            turb_complete$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    }
    
  ) %>%
  ungroup() %>%
  left_join(
    detection_history %>%
      dplyr::select(animal_id, water_year, occ_2, occ_3, status),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    took_geo = as.integer(occ_2 == 2),
    took_ss  = as.integer(occ_3 == 4),
    failed   = as.integer(status %in% c("up_incomplete", "incomplete_dead"))
  )

cat("Multiscale turbidity data built\n")
cat("Missing turb_1day:", sum(is.na(fish_turb_multiscale$turb_1day_raw)), "\n")
cat("Missing turb_7day:", sum(is.na(fish_turb_multiscale$turb_7day_raw)), "\n")
cat("Missing turb_30day:", sum(is.na(fish_turb_multiscale$turb_30day_raw)), "\n")

# check correlations with flow before screening
flow_turb_check <- fish_turb_multiscale %>%
  dplyr::select(animal_id, water_year,
                turb_1day_raw, turb_7day_raw,
                turb_30day_raw) %>%
  left_join(
    fish_flow_09 %>%
      dplyr::select(animal_id, water_year,
                    flow_occ2_raw, flow_30day_raw),
    by = c("animal_id", "water_year")
  )

cat("\n=== TURBIDITY-FLOW CORRELATIONS ===\n")
cat("Same-day log turb vs same-day flow:\n")
cat("r =", round(cor(flow_turb_check$flow_occ2_raw,
                     flow_turb_check$turb_1day_raw,
                     use = "complete.obs"), 3), "\n")
cat("30-day log turb vs 30-day flow:\n")
cat("r =", round(cor(flow_turb_check$flow_30day_raw,
                     flow_turb_check$turb_30day_raw,
                     use = "complete.obs"), 3), "\n")
cat("7-day log turb vs 30-day flow:\n")
cat("r =", round(cor(flow_turb_check$flow_30day_raw,
                     flow_turb_check$turb_7day_raw,
                     use = "complete.obs"), 3), "\n")

# run multiscale screening
turb_vars   <- c("turb_1day_raw", "turb_7day_raw", "turb_14day_raw",
                 "turb_30day_raw", "turb_45day_raw", "turb_60day_raw")
turb_labels <- c("Same-day", "7-day", "14-day", "30-day", "45-day", "60-day")
outcomes      <- c("took_geo", "took_ss", "failed")
outcome_labels <- c("Georgiana routing", "SS routing", "Migration failure")

turb_screen <- data.frame()

for(o in seq_along(outcomes)){
  outcome_var <- outcomes[o]
  for(t in seq_along(turb_vars)){
    turb_var <- turb_vars[t]
    
    df_sub <- fish_turb_multiscale %>%
      dplyr::select(all_of(c(outcome_var, turb_var, "occ_2"))) %>%
      filter(complete.cases(dplyr::select(., all_of(c(outcome_var, turb_var)))))
    
    if(outcome_var == "took_ss"){
      df_sub <- df_sub %>% filter(occ_2 != 2)
    }
    
    if(nrow(df_sub) < 10) next
    
    fit <- glm(as.formula(paste(outcome_var, "~", turb_var)),
               data = df_sub, family = binomial(link = "logit"))
    s   <- summary(fit)
    
    null_ll   <- logLik(glm(as.formula(paste(outcome_var, "~ 1")),
                            data = df_sub, family = binomial))
    pseudo_r2 <- 1 - as.numeric(logLik(fit)) / as.numeric(null_ll)
    
    auc_val <- tryCatch({
      as.numeric(pROC::auc(pROC::roc(df_sub[[outcome_var]],
                                     fitted(fit), quiet = TRUE)))
    }, error = function(e) NA)
    
    turb_screen <- rbind(turb_screen, data.frame(
      outcome    = outcome_labels[o],
      turb_scale = turb_labels[t],
      n          = nrow(df_sub),
      beta       = round(coef(fit)[2], 3),
      se         = round(s$coefficients[2, 2], 3),
      p_value    = round(s$coefficients[2, 4], 4),
      pseudo_r2  = round(pseudo_r2, 4),
      auc        = round(auc_val, 3)
    ))
  }
}

cat("\n=== TURBIDITY TEMPORAL SCALE SCREENING ===\n")
print(turb_screen, row.names = FALSE)

# compare failed vs completed fish turbidity
fish_turb_multiscale %>%
  group_by(failed) %>%
  dplyr::summarise(
    n = n(),
    mean_turb_7day_ntu  = round(exp(mean(turb_7day_raw,  na.rm = TRUE)), 1),
    mean_turb_30day_ntu = round(exp(mean(turb_30day_raw, na.rm = TRUE)), 1),
    .groups = "drop"
  ) %>%
  print()

save(fish_turb_multiscale, turb_daily, turb_complete,
     turb_screen, flow_turb_check,
     file = "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/turb_covariates.RData")

cat("Turbidity covariates saved — not pursued further due to:\n")
cat("1. Lower explanatory power than temperature (pseudo R2 0.05 vs 0.44)\n")
cat("2. Borderline collinearity with 30-day flow (r = 0.727)\n")
cat("3. Biologically less direct cue than temperature or flow\n")








rv_flow <- read.csv(
  "C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv"
)

rv_flow_complete <- rv_flow %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_cfs = value) %>%
  filter(!is.na(flow_cfs),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_cfs = zoo::na.approx(flow_cfs, na.rm = FALSE))

geo_flow <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/RioVista_Confluence_Flows.csv"
)

ges_flow_complete <- geo_flow %>%
  mutate(date = as.Date(dateTime, format = "%m/%d/%Y %H:%M")) %>%
  group_by(date) %>%
  dplyr::summarise(flow_GES = mean(GES, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(date),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_GES = zoo::na.approx(flow_GES, na.rm = FALSE))

temp_flow_raw <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/daily_tidalfilter_riovista_temp_flow.csv"
)

temp_flow_clean <- temp_flow_raw %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  dplyr::select(date, temp_c = daily_temp_c) %>%
  filter(!is.na(date),
         date >= as.Date("2006-09-01"),
         date <= as.Date("2017-09-30")) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(as.Date("2006-09-01"),
                             as.Date("2017-09-30"),
                             by = "day")) %>%
  mutate(temp_c = zoo::na.approx(temp_c, na.rm = FALSE))

benicia_dates <- model1_events %>%
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(benicia_date = as.Date(min(first_detection)),
                   .groups = "drop")

occ2_dates <- model1_events %>%
  filter(occasion == 2) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(arr_date_occ2 = as.Date(min(first_detection)),
                   .groups = "drop")

occ3_dates <- model1_events %>%
  filter(occasion == 3) %>%
  group_by(animal_id, water_year) %>%
  dplyr::summarise(arr_date_occ3 = as.Date(min(first_detection)),
                   .groups = "drop")

fish_cov_12b <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(occ2_dates,    by = c("animal_id", "water_year")) %>%
  left_join(occ3_dates,    by = c("animal_id", "water_year")) %>%
  left_join(benicia_dates, by = c("animal_id", "water_year")) %>%
  rowwise() %>%
  mutate(
    
    flow_occ2_raw = {
      if(is.na(arr_date_occ2)) NA_real_ else {
        idx <- which(rv_flow_complete$date == arr_date_occ2)
        if(length(idx) == 0) NA_real_ else rv_flow_complete$flow_cfs[idx]
      }
    },
    
    flow_occ3_raw = {
      if(is.na(arr_date_occ3)) NA_real_ else {
        idx <- which(ges_flow_complete$date == arr_date_occ3)
        if(length(idx) == 0) NA_real_ else ges_flow_complete$flow_GES[idx]
      }
    },
    
    temp_7day_raw = {
      if(is.na(benicia_date)) NA_real_ else {
        window <- temp_flow_clean$temp_c[
          temp_flow_clean$date >= (benicia_date - 7) &
            temp_flow_clean$date <   benicia_date]
        window <- window[!is.na(window)]
        if(length(window) == 0) NA_real_ else mean(window, na.rm = TRUE)
      }
    }
    
  ) %>%
  ungroup()

cat("Missing flow_occ2:", sum(is.na(fish_cov_12b$flow_occ2_raw)), "\n")
cat("Missing flow_occ3:", sum(is.na(fish_cov_12b$flow_occ3_raw)), "\n")
cat("Missing temp_7day:", sum(is.na(fish_cov_12b$temp_7day_raw)), "\n")

flow_occ2_mean_12b <- mean(fish_cov_12b$flow_occ2_raw, na.rm = TRUE)
flow_occ2_sd_12b   <- sd(fish_cov_12b$flow_occ2_raw,   na.rm = TRUE)
flow_occ3_mean_12b <- mean(fish_cov_12b$flow_occ3_raw, na.rm = TRUE)
flow_occ3_sd_12b   <- sd(fish_cov_12b$flow_occ3_raw,   na.rm = TRUE)
temp_7day_mean_12b <- mean(fish_cov_12b$temp_7day_raw, na.rm = TRUE)
temp_7day_sd_12b   <- sd(fish_cov_12b$temp_7day_raw,   na.rm = TRUE)

fish_cov_12b <- fish_cov_12b %>%
  mutate(
    flow_occ2_std = (flow_occ2_raw - flow_occ2_mean_12b) / flow_occ2_sd_12b,
    flow_occ3_std = (flow_occ3_raw - flow_occ3_mean_12b) / flow_occ3_sd_12b,
    temp_7day_std = (temp_7day_raw - temp_7day_mean_12b) / temp_7day_sd_12b
  )

cat("\nflow_occ2 mean:", round(flow_occ2_mean_12b, 0), "cfs SD:",
    round(flow_occ2_sd_12b, 0), "\n")
cat("flow_occ3 mean:", round(flow_occ3_mean_12b, 0), "cfs SD:",
    round(flow_occ3_sd_12b, 0), "\n")
cat("temp_7day mean:", round(temp_7day_mean_12b, 2), "C SD:",
    round(temp_7day_sd_12b, 2), "\n")

identical(fish_cov_12b$animal_id, detection_history$animal_id)
identical(as.numeric(fish_cov_12b$water_year),
          as.numeric(detection_history$water_year))

flow_occ2_12b <- fish_cov_12b$flow_occ2_std
flow_occ3_12b <- fish_cov_12b$flow_occ3_std
temp_7day_12b <- fish_cov_12b$temp_7day_std

flow_occ2_12b[is.na(flow_occ2_12b)] <- 0
flow_occ3_12b[is.na(flow_occ3_12b)] <- 0
temp_7day_12b[is.na(temp_7day_12b)] <- 0

cat("flow_occ2 range:", round(range(flow_occ2_12b), 2), "\n")
cat("flow_occ3 range:", round(range(flow_occ3_12b), 2), "\n")
cat("temp_7day range:", round(range(temp_7day_12b), 2), "\n")

# save the covariate file and standardization parameters
# run this after your existing covariate extraction code

library(writexl)

# fish level covariates
fish_cov_save <- fish_cov_12b %>%
  dplyr::select(animal_id, water_year,
                flow_occ2_raw, flow_occ2_std,
                flow_occ3_raw, flow_occ3_std,
                temp_7day_raw, temp_7day_std)

write.csv(fish_cov_save,
          "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_model12b_covariates.csv",
          row.names = FALSE)

# standardization parameters for back-transforming later
std_params <- data.frame(
  covariate = c("flow_occ2", "flow_occ3", "temp_7day"),
  mean      = c(flow_occ2_mean_12b, flow_occ3_mean_12b, temp_7day_mean_12b),
  sd        = c(flow_occ2_sd_12b,   flow_occ3_sd_12b,   temp_7day_sd_12b)
)

write.csv(std_params,
          "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_model12b_std_params.csv",
          row.names = FALSE)

cat("Covariate files saved\n")
cat("flow_occ2 mean:", round(flow_occ2_mean_12b, 0), "SD:", round(flow_occ2_sd_12b, 0), "\n")
cat("flow_occ3 mean:", round(flow_occ3_mean_12b, 0), "SD:", round(flow_occ3_sd_12b, 0), "\n")
cat("temp_7day mean:", round(temp_7day_mean_12b, 2), "SD:", round(temp_7day_sd_12b, 2), "\n")
