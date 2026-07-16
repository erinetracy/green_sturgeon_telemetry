#==============================================================================
# 00_BUILD_ENVIRONMENTAL_DATA.R
# Builds a single clean daily environmental file from raw sources
# Run once — saves gs_daily_environmental_data.csv
# Covers: Rio Vista flow, GES flow, Rio Vista temperature
#builds off of temp_flow_turbidity_covariate_extraction code
#==============================================================================

library(dplyr)
library(tidyr)
library(zoo)
library(readxl)

start_date <- as.Date("2006-09-01")
end_date   <- as.Date("2017-09-30")

# --- Rio Vista tidally filtered flow ---
rv_flow <- read.csv(
  "C:/Users/eetracy/Desktop/Post_doc_GS/daily_tidalfilter_riovista.csv"
) %>%
  mutate(date = as.Date(time)) %>%
  dplyr::select(date, flow_rv_cfs = value) %>%
  filter(!is.na(date), date >= start_date, date <= end_date) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(start_date, end_date, by = "day")) %>%
  mutate(flow_rv_cfs = zoo::na.approx(flow_rv_cfs, na.rm = FALSE))

cat("RV flow: ", nrow(rv_flow), "days |",
    "Missing:", sum(is.na(rv_flow$flow_rv_cfs)), "\n")

# --- GES flow (Sacramento above Steamboat/Sutter junction) ---
ges_flow <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/RioVista_Confluence_Flows.csv"
) %>%
  mutate(date = as.Date(dateTime, format = "%m/%d/%Y %H:%M")) %>%
  group_by(date) %>%
  dplyr::summarise(flow_ges_cfs = mean(GES, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(date), date >= start_date, date <= end_date) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%
  mutate(flow_ges_cfs = zoo::na.approx(flow_ges_cfs, na.rm = FALSE))

cat("GES flow:", nrow(ges_flow), "days |",
    "Missing:", sum(is.na(ges_flow$flow_ges_cfs)), "\n")

# --- Rio Vista temperature ---
temp_data <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/flow/daily_tidalfilter_riovista_temp_flow.csv"
) %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  dplyr::select(date, temp_c = daily_temp_c) %>%
  filter(!is.na(date), date >= start_date, date <= end_date) %>%
  arrange(date) %>%
  tidyr::complete(date = seq(start_date, end_date, by = "day")) %>%
  mutate(temp_c = zoo::na.approx(temp_c, na.rm = FALSE))

cat("Temperature:", nrow(temp_data), "days |",
    "Missing:", sum(is.na(temp_data$temp_c)), "\n")
cat("Temperature starts:", as.character(min(temp_data$date[!is.na(temp_data$temp_c)])), "\n")

# --- Merge into single daily file ---
daily_env <- rv_flow %>%
  left_join(ges_flow,  by = "date") %>%
  left_join(temp_data, by = "date")

cat("\nMerged daily environmental file:\n")
cat("Rows:", nrow(daily_env), "\n")
cat("Missing flow_rv_cfs:",  sum(is.na(daily_env$flow_rv_cfs)), "\n")
cat("Missing flow_ges_cfs:", sum(is.na(daily_env$flow_ges_cfs)), "\n")
cat("Missing temp_c:",       sum(is.na(daily_env$temp_c)), "\n")
cat("Date range:", as.character(range(daily_env$date)), "\n")

# --- Save ---
write.csv(daily_env,
          "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/gs_daily_environmental_data.csv",
          row.names = FALSE)

cat("\nSaved: gs_daily_environmental_data.csv\n")
