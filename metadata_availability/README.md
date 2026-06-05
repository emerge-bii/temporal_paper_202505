# Metadata availability

This directory contains metadata-availability analyses it was not used in the primary manuscript, but aided in decision-making around what analyses were possible.

The workflow summarizes which environmental metadata variables are available for sequenced samples and visualizes data completeness across year, habitat, core, and depth. It also generates a cleaned metadata table used for plotting.

## Directory contents

```text
metadata_availability/
├── R/
│   └── 01_plot_available_data_in_metadata_sheet.R
└── outputs/
    └── metadata_used_for_plots.csv
```

The script creates a local figure directory when run:

```text
metadata_availability/
└── figures/
    └── data_completeness_plots/
```

## Script overview

| Script                                         | Purpose                                                                                                                          |
| ---------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `R/01_plot_available_data_in_metadata_sheet.R` | Plots availability of environmental metadata variables across sequenced samples and writes the metadata table used for plotting. |

## Required R packages

| Script                                         | Required R packages                                                                  |
| ---------------------------------------------- | ------------------------------------------------------------------------------------ |
| `R/01_plot_available_data_in_metadata_sheet.R` | `tidyverse`, `vegan`, `viridis`, `cowplot`, `ggh4x`, `here`, `RColorBrewer`, `knitr` |

## Expected inputs

The script sources the top-level setup file:

```r
source(here("setup.R"))
```

Expected objects loaded by setup include:

| Object                     | Description                                                                            |
| -------------------------- | -------------------------------------------------------------------------------------- |
| `input_ra`                 | Relative-abundance object containing sample metadata.                                  |
| `input_ra$sample_metadata` | Metadata used to evaluate availability of environmental variables.                     |
| `sample_metadata`          | Sample metadata used to plot sample depths, water-table depth, and active-layer depth. |

Expected metadata columns include:

```text
temporal_sample_id
Year__
Habitat__
Core__
DepthAvg__
DepthMin__
DepthMax__
DepthLumping
WTD
ALD
ALD_detectable
```

## Outputs

Generated outputs include:

```text
metadata_availability/outputs/metadata_used_for_plots.csv
metadata_availability/figures/data_completeness_plots/*.png
```
The `.png` files show whether each metadata variable is present or missing for each sequenced sample.
