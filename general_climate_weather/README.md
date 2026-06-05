# General climate and weather analyses

This directory contains climate and weather analyses used in the primary manuscript text for **“Stable states in an unstable landscape: microbial resistance at the front line of climate change.”**

The workflow places the 2011–2017 temporal metagenomic study period in the context of long-term air temperature and precipitation records from Abisko/Stordalen. It generates figures describing annual and seasonal climate context, including temperature trajectories, precipitation patterns, and comparisons between the temporal study years and longer observational records.

## Directory contents

```text
general_climate_weather/
└── 01_weather_analyses.R
```

The script creates local output directories when run:

```text
general_climate_weather/
├── figures/
└── outputs/
```

## Script overview

| Script                  | Purpose                                                                                                                                                                                 |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `01_weather_analyses.R` | Reads weather and climate data from Abisko and field site, cleans and combines temperature/precipitation records, defines plotting aesthetics, and generates long-term climate-context figures for the manuscript. |

## Required R packages

| Script                  | Required R packages                                                                                                               |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `01_weather_analyses.R` | `tidyverse`, `viridis`, `cowplot`, `lubridate`, `broom`, `broom.mixed`, `MuMIn`, `here`, `readxl`, `RColorBrewer`, `zoo`, `knitr` |

## Expected inputs

The script sources the top-level project setup file:

```r
source(here("setup.R"))
```

It also reads publicly available weather data that it expects to find in the repository-level `data/` directory, including:

```text
data/WeatherData/Abisko_1913-2016.xls
data/smhi-opendata_2_188800_20230613_170724.csv
data/smhi-opendata_1_188800_20230613_173232.csv
data/smhi-opendata_5_188800_20230613_165112.csv
```

Some large or restricted input files may not be committed directly to GitHub. See the top-level repository README for instructions on preparing the `data/` directory.

## Workflow

Run from the repository root:

```bash
Rscript general_climate_weather/01_weather_analyses.R
```

