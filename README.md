# Green Sturgeon Upstream Migration Telemetry Analysis

## Overview
This repository contains R scripts for analyzing acoustic telemetry data from adult green sturgeon (*Acipenser medirostris*) making upstream spawning migrations in the Sacramento River system. The analysis estimates route-specific survival, detection, and routing probabilities using a discrete hidden Markov multistate model (dDHMMo) following Perry et al. (2018). A secondary analysis examines the influence of flow and temperature on migration initiation timing.

## Data
Acoustic telemetry detection data were downloaded from the Pacific Ocean Tracking Network (OTN) PATH database. Receiver metadata were compiled from ArcGIS spatial analysis. Flow data are tidally-filtered daily discharge from the USGS Rio Vista gauge (station 11455420). Temperature data are daily values from the Rio Vista monitoring station.

- **Model period:** Water years 2007–2017
- **Fish:** 228 confirmed upstream migrants (140 up_complete, 83 up_incomplete, 5 incomplete_dead)
- **Receiver arrays:** Benicia/Carquinez, Rio Vista, SR_MOUTH, SR_BLWSTEAM, SR_KK345R/SR_FREEPORT, upper Sacramento, spawning ground

---

## Script Pipeline

Scripts should be run in order. Each script saves outputs required by the next.

---

### `receiver_cleaning_gsflow.R`
**Purpose:** Receiver metadata cleaning, spatial grouping, and mapping.

Loads raw receiver deployment records from the OTN PATH database and assigns spatial group labels (bay, benicia, carquinez, sacramento, georgiana, mok_deltacross, steamboat_sutter, yolo_bypass, spawning_ground) based on ArcGIS spatial analysis. Produces interactive maps for QC and exports cleaned receiver metadata.

---

### `detection_cleaning_gsflow.R`
**Purpose:** Detection data download, false detection filtering, and event reduction.

Downloads raw green sturgeon detections from the OTN PATH database (requires UCD VPN and credentials), runs the GLATOS false detection filter (tf = 3600 seconds, ~0.8% flagged), reduces raw detections to detection events, removes double-tagged fish, assigns water years, and joins receiver group labels from `arc_receivers_update.csv`.

**Inputs:**
- PATH database (or `green_sturgeon_detection_OTN_012926.csv` if already downloaded)
- `arc_receivers_update.csv`

**Outputs:**
- `events_with_receivergroups_032026.csv` — cleaned detection events with receiver group labels and water year

### `04_multistate_data_prep.R`
**Purpose:** Prepare detection history matrix for the multistate model.

Filters events to model years (2007–2017) and upstream migration statuses using an explicit inclusion list built from migration_status to avoid join artifacts. Assigns occasions (1–7) and states (1–4) to each detection event. Builds the 228 × 7 detection history matrix using last(state) per occasion to capture the final committed route rather than exploratory forays. Applies fixes for impossible detection sequences caused by staging behavior. Saves all objects needed for the model script.
Uses detection history relative to key receiver arrays (golden_gate, bay, benicia/carquinez, sacramento, spawning_ground) to classify fish as up_complete, up_incomplete, down_complete, down_incomplete, incomplete_dead, or bad. Applies manual corrections for two fish with anomalous detection histories and labels confirmed mortality events for specific water years only (not all years for multi-year tagged fish).

**Status categories:**
- `up_complete` — confirmed full upstream spawning migration (ocean → Benicia → spawning ground)
- `up_incomplete` — started upstream migration but did not reach spawning ground
- `incomplete_dead` — confirmed or likely mortality during upstream migration
- `down_complete` — confirmed full downstream migration (spawning ground → ocean)
- `down_incomplete` — started downstream migration but did not reach ocean
- `bad` — ambiguous movement, delta resident, or data issues

**Key decisions documented in script:**
- Decker Island excluded from occasion 2 (ambiguous route indicator)
- DCC restricted to SR_DCC and SR_DCC2 only (confirmed channel entry)
- last(state) used for detection history (committed route, not exploration)
- Steamboat/Sutter pass-through at occ4, rejoin at occ5 (SR_KK345R above both rejoin points)
- Exploration behavior preserved in `exploration_summary` for IBM analysis

**Inputs:**
- `events_with_receivergroups_032026.csv`
- `migratory_status_03162026.csv` (migration_status object)
- `arc_receivers_update.csv`

**Outputs:**
- `gs_multistate_data.RData` — all model objects:
  - `ch_mat_nimble` — 228 × 7 matrix formatted for dDHMMo (0 recoded to nstate+1)
  - `ch_mat` — 228 × 7 matrix (0 = not detected)
  - `detection_history` — fish-level detection history with status and exploration variables
  - `fish_info` — fish-level metadata
  - `exploration_summary` — fish-level exploration behavior variables
  - `model1_events` — filtered detection events with occasions and states
  - `nstate` — number of states (6)

---

### `05_multistate_model.R`
**Purpose:** Build and run the discrete hidden Markov multistate model (dDHMMo).

Implements a multistate mark-recapture model following Perry et al. (2018) using nimbleEcology. Estimates route-specific survival (S), detection (p), routing probabilities (psi), migration failure probability (phi_fail), and spawning ground arrival probability (lambda). Uses uninformative flat priors dbeta(1,1) on all parameters and slice samplers for MCMC following Perry et al.

**Model structure:**
- States: 1=Sacramento, 2=Georgiana, 3=DCC, 4=Steamboat/Sutter, 5=Dead (absorbing), 6=Failed migration (absorbing)
- Occasions: 7 (Benicia → Rio Vista → SR_MOUTH → SR_BLWSTEAM → SR_KK345R/SR_FREEPORT → upper Sacramento → spawning ground)
- MCMC: 3 chains × 50,000 iterations, 10,000 burnin, thin = 10

**Inputs:**
- `gs_multistate_data.RData`

**Outputs:**
- `gs_mcmc_full_run2_228fish.RData` — MCMC posterior samples
- `gs_mcmc_traces_run2.pdf` — trace plots for key parameters

---

### `06_flow_multistate.R`
**Purpose:** Flow and temperature covariate development and migration timing analysis.

Explores the relationship between Sacramento River hydrology (tidally-filtered daily discharge at USGS Rio Vista gauge) and green sturgeon migration behavior at two scales: (1) migration initiation — what flow and temperature conditions trigger fish to pass Benicia/Carquinez; (2) route selection — what flow conditions at each junction predict alternative route use.

**Analyses:**
- Hydrographs centered on fish passage dates (relative day 0)
- Flow at multiple temporal scales (day of, 3-day, 7-day, 15-day, 30-day, 45-day, 60-day mean prior to passage)
- ANOVA of flow temporal scales against route selection (F-statistic comparison)
- Temperature distribution at migration start across temporal scales
- Flow rate of change at migration start
- Route selection proportions by flow quartile

**Key findings:**
- 30-day mean flow before Benicia passage best separates route choices (F = 2.13, p = 0.097)
- Steamboat/Sutter use decreases with increasing flow; DCC use increases dramatically at high flow
- Temperature at migration start is more normally distributed and consistent across temporal scales than flow
- Most fish migrate during stable flow conditions (rate of change near zero)
- Flow signal is driven primarily by wet years (2016–2017 in model period)

**Inputs:**
- `gs_multistate_data.RData`
- `daily_tidalfilter_riovista.csv` — tidally-filtered daily flow at Rio Vista (USGS 11455420)
- `Rio_Vista_daily_temp.csv` — daily temperature at Rio Vista

**Outputs:** Exploratory plots (not saved by default)

---

### `GSSMB_*.R` scripts and code_perry_2018
**Purpose:** Perry et al. (2018, 2026) simulation and verification code.

These scripts contain the original Perry et al. (2018) JAGS-based multistate model code adapted for green sturgeon, used for model structure verification and comparison and 2026 nimble code Not part of the main analysis pipeline.

---

## Dependencies

```r
# Core
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

# Telemetry
library(glatos)        # install from OTN r-universe
library(RPostgres)     # PATH database connection

# Spatial
library(leaflet)

# Bayesian model
library(nimble)
library(nimbleEcology)
library(abind)
library(MCMCvis)
```

Install glatos from OTN r-universe:
```r
install.packages('glatos', repos = c('https://ocean-tracking-network.r-universe.dev',
                                      'https://cloud.r-project.org'))
```

---

## Reference
Perry, R.W., et al. (2018). Flow-mediated effects on travel time, survival, and routing of juvenile Chinook salmon in a spatially complex, tidally forced river delta. *Transactions of the American Fisheries Society.*
