# CN versatility

This directory contains scripts for visualizing carbon and nitrogen metabolic versatility among EMERGE MAGs, with a focus on MAGs identified as important members of WGCNA modules or environmental-variable-associated VIP lists.

The analyses combine pathway calls, MAG metadata, specialization assignments, VIP membership tables, and environmental correlations to generate figures describing overlap between carbon- and nitrogen-cycling potential.

## Directory contents

```text
CN_versatility/
├── CN_Versatility.R
├── CN_network_vis.R
└── CN_versatility_type2b.R
```

The scripts create local output folders when run:

```text
CN_versatility/
├── outputs/
└── figures/
```

These folders are used for generated tables, diagnostic plots, and manuscript or exploratory figures.

## Scripts

| Script                    | Purpose                                                                                                                                                                                                 | Main outputs                                                                                                                       |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `CN_Versatility.R`        | Defines carbon and nitrogen versatility categories from pathway annotations and generates early Type 1 and Type 2 VIP-versatility plots.                                                                | PNG figures showing VIP MAGs by carbon and nitrogen pathway membership, habitat, environmental variable, and specialization class. |
| `CN_versatility_type2b.R` | Generates a revised Type 2b version of the carbon/nitrogen versatility plot. This version uses distribution-style dot plots and patchwork panels to show overlap between carbon and nitrogen functions. | `type2b_*_VIP_specials.png` figures in `CN_versatility/figures/`.                                                                  |
| `CN_network_vis.R`        | Generates network- and correlation-focused visualizations linking VIP MAGs, environmental variables, specialization classes, carbon/nitrogen versatility, and network centrality.                       | Hive plots, MAG/environment correlation tables, correlation-check PDFs, combined VIP correlation tables, and CN capability plots.  |

## Setup

Run these scripts from the top level of the repository so that `here()` resolves paths correctly.

```bash
cd temporal_paper_202505
```

Create and activate the project conda environment using the top-level repository instructions:

```bash
./install_dependencies.sh
conda activate temporal_paper_vX
```

Replace `temporal_paper_vX` with the current environment name used by the repository.

Each script sources the top-level setup file:

```r
source(here("setup.R"))
```

Because of this, most shared objects are expected to be loaded or defined by `setup.R` and the repository-level data-preparation workflow.

## Required R packages

The recommended installation route is the top-level conda environment. The table below summarizes package requirements by script.

| Script                    | Required R packages                                                                                                                            |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `CN_Versatility.R`        | `tidyverse`, `here`, `RColorBrewer`                                                                                                            |
| `CN_versatility_type2b.R` | `tidyverse`, `here`, `RColorBrewer`, `ggdist`, `patchwork`, `purrr`, `ggpubr`                                                                  |
| `CN_network_vis.R`        | `tidyverse`, `here`, `RColorBrewer`, `ggraph`, `igraph`, `ggplot2`, `tidygraph`, `dplyr`, `HiveR`, `grid`, `cowplot`, `broom`, `fs`, `ggrepel` |

## Expected input objects and files from setup.R

The scripts assume that several shared objects are already available after sourcing `setup.R`. These likely include:

| Object or file                             | Description                                                                                            |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------ |
| `product_refined`                          | Refined pathway/function table for MAGs, including pathway/subpathway calls.                           |
| `genome_clusters`                          | MAG clustering or representative-genome table.                                                         |
| `metapathway_groups`                       | Mapping of subpathways to broader metapathway groups.                                                  |
| `vip_members`                              | Table of VIP MAGs, environmental variables, habitat/data-origin labels, module colors, and VIP scores. |
| `specialisations`                          | MAG specialization assignments, including carbon and nitrogen specialization categories.               |
| `sample_metadata`                          | Sample metadata, including habitat and environmental variables.                                        |
| `trimmed_mean`                             | MAG relative-abundance table used for environmental correlation analyses.                              |
| `input_ra` / `input_counts`                | MAG abundance and taxonomy objects used for plotting and labeling.                                     |
| `network_stats`                            | Network centrality/statistics table used in network visualizations.                                    |
| `genome_sizes`                             | Genome-size metadata used in CN/network plots.                                                         |
| `data/vip_members_pathways-08-29-2024.txt` | VIP/pathway membership table read by the network-visualization script.                                 |
| `data/all_organisms_for_naming_v6.csv`     | Organism naming table used to create an updated naming table.                                          |


## Outputs

Generated outputs are written to:

```text
CN_versatility/figures/
CN_versatility/outputs/
```

Common output types include:

* `.png` figures for CN versatility plots,
* `.pdf` files of correlation-check plots,
* `.tsv` correlation tables,
* updated organism/pathway naming files.

Examples of expected output naming patterns include:

```text
type1_*_VIP_specials.png
type2_*_VIP_specials.png
type2b_*_VIP_specials.png
*_VIP_hive_plot.png
*_VIP_corr_table.txt
all_vip_env_var_corr.tsv
```
## Manuscript citation

If using or adapting these scripts, please cite the associated manuscript:

Cronin, D. R., Holland-Moritz, H., Smith, D. A., Aroney, S. T. N., Hodgkins, S. B., Borton, M., Li, Y.-F., Zayed, A., Healy, K., Persson, A., IsoGenie Field & Analytic Teams 2010–2017, EMERGE Institute Coordinators, Tfaily, M. M., Crill, P., McCalley, C. K., Wrighton, K., Varner, R. K., Tyson, G. W., Woodcroft, B. J., Bagby, S. C., Ernakovich, J., and Rich, V. I. **Stable states in an unstable landscape: microbial resistance at the front line of climate change.** bioRxiv 2025.02.07.636677.
