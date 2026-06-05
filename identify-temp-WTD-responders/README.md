# Temperature and water-table-depth responders

This directory contains analyses used in the primary manuscript text to identify microbial community features associated with temperature and water-table-depth variation.

The workflow evaluates relationships between environmental gradients and community composition, including Mantel tests, response models, and plots of temperature- or WTD-associated MAGs/pathways. 

The directory is organized as two analysis scripts: one focused on temperature response (used in manuscript) and one focused on water-table-depth response (not used in mansucript).

## Directory contents

```text
identify-temp-WTD-responders/
├── 00_Temp_response_analyses.R
└── 01_WTD_response_analyses.R
```

The scripts create local output directories when run:

```text
identify-temp-WTD-responders/
├── figures/
└── outputs/
```

## Script overview

| Script                        | Purpose                                                                                                                                                                    |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `00_Temp_response_analyses.R` | Identifies and visualizes temperature-associated community patterns, including air-temperature and soil-temperature Mantel analyses, response models, and responder plots. |
| `01_WTD_response_analyses.R`  | Identifies and visualizes relationships between water-table-depth metrics and community composition, including WTD Mantel tests and methane/WTD exploratory plots.         |

## Required R packages

| Script                        | Required R packages                                                                                                                                                       |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `00_Temp_response_analyses.R` | `tidyverse`, `broom`, `broom.mixed`, `vegan`, `viridis`, `GUniFrac`, `cowplot`, `kableExtra`, `multcompView`, `propr`, `here`, `RColorBrewer`, `phyloseq`, `ape`, `knitr` |
| `01_WTD_response_analyses.R`  | `tidyverse`, `vegan`, `viridis`, `GUniFrac`, `cowplot`, `multcompView`, `here`, `RColorBrewer`, `knitr`                                                                   |

Some functions are called through package-qualified syntax rather than `library()`, including `RColorBrewer::brewer.pal()`, `ape::keep.tip()`, and `phyloseq::UniFrac()`.

## Expected inputs

Both scripts source the top-level setup file:

```r
source(here("setup.R"))
```

Expected objects loaded or created by the project setup include:

| Object               | Description                                                                                    |
| -------------------- | ---------------------------------------------------------------------------------------------- |
| `input_ra`           | Relative-abundance input object used as the main community input.                              |
| `sample_metadata`    | Sample metadata containing habitat, depth, year, temperature, and water-table-depth variables. |
| `trimmed_mean`       | MAG relative-abundance table used for community dissimilarity calculations.                    |
| `metapathway_groups` | Functional pathway grouping table used for plotting colors and labels.                         |
| `input$tree`         | Phylogenetic tree used for UniFrac-based analyses where relevant.                              |

The temperature-response script also reads, generated from the WGCNA analysis:

```text
data/mags_vip_temperature.csv
```

## Outputs

Outputs are written to:

```text
identify-temp-WTD-responders/figures/
identify-temp-WTD-responders/outputs/
```
