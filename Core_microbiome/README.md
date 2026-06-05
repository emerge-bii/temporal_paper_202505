# Core microbiome

This directory contains exploratory core-microbiome analyses for the EMERGE temporal manuscript repository.

**Note:** The analyses in this folder contributed to the development of ideas and interpretation for the manuscript, but they did not make it into the final manuscript text. They are retained here for transparency, provenance, and potential reuse in future analyses.

## Directory contents

```text id="j96u3y"
Core_microbiome/
├── 00_core_microbiome.R
├── sncm.fit_metagenome.R
└── outputs/
```

The main script also creates a local `figures/` directory when run:

```text id="cyxmsx"
Core_microbiome/
├── figures/
└── outputs/
```

## Overview

The core-microbiome workflow was designed to identify persistent and influential MAGs/OTUs across EMERGE habitats over time. The analysis combines occupancy, relative abundance, Bray-Curtis contribution, and neutral-model approaches to rank taxa and identify candidate “core” members of the community.

The workflow is based conceptually on the core-microbiome framework described by Shade and Stopnisek 2019, with additional filtering and visualization for habitat- and depth-specific EMERGE metagenomic data.

## Scripts

| Script                  | Purpose                                                                                                                                                                                                                                                            |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `00_core_microbiome.R`  | Main exploratory analysis script. Loads project data, defines helper functions, ranks OTUs/MAGs by occupancy and Bray-Curtis contribution, applies core thresholds, fits neutral-model comparisons, and writes habitat/depth-specific core-OTU tables and figures. |
| `sncm.fit_metagenome.R` | Helper functions for fitting Sloan neutral community models to metagenomic abundance data. This script is sourced by `00_core_microbiome.R`.                                                                                                                       |

## Setup

Run scripts from the top level of the repository so that `here()` resolves paths correctly.

```bash id="ev9ko3"
cd temporal_paper_202505
```

Create and activate the project conda environment using the top-level repository instructions:

```bash id="h6mb8z"
./install_dependencies.sh
conda activate temporal_paper_vX
```

Replace `temporal_paper_vX` with the current environment name used by the repository.

The main script sources the top-level setup file:

```r id="w3pkf5"
source(here("setup.R"))
```

Because of this, the workflow assumes that shared project objects are available after `setup.R` is sourced.

## Required R packages

| Script                  | Required R packages                                                                                               |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `00_core_microbiome.R`  | `tidyverse`, `vegan`, `viridis`, `here`, `RColorBrewer`, `knitr`                                                  |
| `sncm.fit_metagenome.R` | Base R/statistical functions used by the Sloan neutral-model fitting workflow; sourced by `00_core_microbiome.R`. |

## Expected input objects

The main script assumes that several shared objects are available after sourcing `setup.R`.

| Object                         | Description                                                                                                                                 |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `input_ra`                     | Relative-abundance input object. The script initially renames this object to `input`.                                                       |
| `input_counts`                 | Count-based input object containing at least an OTU/MAG table and sample metadata.                                                          |
| `input_counts$otu_table`       | OTU/MAG abundance table. Rows are genomes or OTUs, and columns are samples.                                                                 |
| `input_counts$sample_metadata` | Sample metadata used to subset samples by habitat and depth. Expected fields include `temporal_sample_id`, `Habitat__`, and `DepthLumping`. |
| `metapathway_groups`           | Functional/pathway grouping table used for setting pathway color palettes.                                                                  |

These objects are generated elsewhere in the repository workflow and are not created from scratch in this directory.

## Outputs

Generated outputs are written to:

```text id="f5eehx"
Core_microbiome/outputs/
Core_microbiome/figures/
```

## Manuscript citation

If using or adapting these scripts, please cite the associated manuscript:

Cronin, D. R., Holland-Moritz, H., Smith, D. A., Aroney, S. T. N., Hodgkins, S. B., Borton, M., Li, Y.-F., Zayed, A., Healy, K., Persson, A., IsoGenie Field & Analytic Teams 2010–2017, EMERGE Institute Coordinators, Tfaily, M. M., Crill, P., McCalley, C. K., Wrighton, K., Varner, R. K., Tyson, G. W., Woodcroft, B. J., Bagby, S. C., Ernakovich, J., and Rich, V. I. **Stable states in an unstable landscape: microbial resistance at the front line of climate change.** bioRxiv 2025.02.07.636677.
