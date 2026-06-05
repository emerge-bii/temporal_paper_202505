# Quantifying stability with time

This directory contains analyses and figure-generation scripts used in the primary manuscript text to quantify microbial-community stability through time.

The workflow combines multiple “lenses” of stability, including organismal diversity, functional or metabolic summaries, community assembly, network properties, and beta-diversity through time. The outputs support the manuscript figure and associated interpretation of temporal stability across habitats and depth intervals.

## Directory contents

```text
quantify_stability_with_time_figure/
├── Figure1_beta_diversity_over_time/
│   ├── beta_diversity_over_time.Rmd
│   ├── mean_beta_diversity.csv
│   ├── mean_beta_diversity_within_ref.csv
│   ├── minmax_beta_diversity.csv
│   ├── minmax_beta_diversity_within_ref.csv
│   ├── ref_diversity_long.csv
│   ├── ref_diversity_long_against_ref.csv
│   └── ref_diversity_long_within_ref.csv
├── outputs/
├── Figure1_beta_diversity_over_time.zip
├── Stability_fig.R
└── agu_plots.R
```

## Script overview

| Script or subdirectory                                          | Purpose                                                                                                                                                           |
| --------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Stability_fig.R`                                               | Main manuscript stability-figure script. Combines alpha diversity, granularity, assembly, network, and beta-diversity outputs into stability summaries and plots. |
| `agu_plots.R`                                                   | Generates additional ordination and presentation-style plots, including NMDS/PCoA visualizations of community dissimilarity across habitat, depth, and year.      |
| `Figure1_beta_diversity_over_time/beta_diversity_over_time.Rmd` | RMarkdown workflow for calculating and visualizing beta diversity through time.                                                                                   |
| `Figure1_beta_diversity_over_time/*.csv`                        | Precomputed beta-diversity summary tables used by the stability figure workflow.                                                                                  |

## Required R packages

| Script                                                          | Required R packages                                                                                                          |
| --------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `Stability_fig.R`                                               | `tidyverse`, `ggnewscale`, `here`, `RColorBrewer`, `data.table`, `furrr`, `future`, `vegan`, `broom`, `purrr`, `knitr`       |
| `agu_plots.R`                                                   | `tidyverse`, `here`, `ggnewscale`, `RColorBrewer`, `vegan`, `ape`, `grid`, `broom`, `mgcv`, `lme4`                           |
| `Figure1_beta_diversity_over_time/beta_diversity_over_time.Rmd` | RMarkdown beta-diversity workflow; likely requires `tidyverse`, `vegan`, `here`, and plotting packages used in the notebook. |

Some packages are called through namespace syntax or through functions loaded by `tidyverse`, including `purrr::map()`, `broom::tidy()`, `broom::glance()`, and `furrr::future_map()`.

## Expected inputs

The scripts source the top-level setup file:

```r
source(here("setup.R"))
```

Expected objects loaded by setup include:

| Object               | Description                                                                            |
| -------------------- | -------------------------------------------------------------------------------------- |
| `sample_metadata`    | Sample metadata with habitat, depth, year, and temporal sample identifiers.            |
| `input_ra`           | Relative-abundance object used for sample metadata and abundance analyses.             |
| `trimmed_mean`       | MAG relative-abundance table used for ordination and community dissimilarity analyses. |
| `network_stats`      | Network statistics used to evaluate network stability over time.                       |
| `singlem_otu_tables` | SingleM OTU/gene table used to calculate strain-level alpha diversity.                 |

Additional input files read by `Stability_fig.R` include:

```text
data/20221109-cv-mags.csv
data/granularity-BH - granularity-corrected-10092023.csv
Assembly-analysis/outputs/mantel_with_time.txt
data/20240403_network_comparisons.txt
data/singlem_alpha.txt
```

If `data/singlem_alpha.txt` does not exist, the script attempts to generate it from `singlem_otu_tables`.

## Workflow
The beta-diversity RMarkdown workflow can be rendered from the repository root or from within its subdirectory. It should be run first.

```bash
Rscript -e "rmarkdown::render('quantify_stability_with_time_figure/Figure1_beta_diversity_over_time/beta_diversity_over_time.Rmd')"
```

Main figure generation, & statsitical analysis:

```bash
Rscript quantify_stability_with_time_figure/Stability_fig.R
```

Optional/supporting plots:

```bash
Rscript quantify_stability_with_time_figure/agu_plots.R
```

## Outputs

Outputs are written to:

```text
quantify_stability_with_time_figure/outputs/
quantify_stability_with_time_figure/figures/
quantify_stability_with_time_figure/Figure1_beta_diversity_over_time/
```
