# Assembly analysis

This directory contains scripts and selected outputs for the microbial community assembly analyses used in the manuscript **“Stable states in an unstable landscape: microbial resistance at the front line of climate change.”** These analyses evaluate patterns of phylogenetic and functional assembly across EMERGE metagenomic datasets, including analyses based on MAGs, SingleM OTUs, and functional modules.

## Directory structure

```text
Assembly-analysis/
├── R/
│   ├── 00_all_MAGs_assembly/
│   ├── 01_SingleM_assembly/
│   ├── 02_Module_assembly/
│   ├── 03_Misc_figures/
│   └── helper_scripts/
├── bash/
├── slurm/
└── outputs/
```

### `R/`

R scripts for preparing inputs, running assembly analyses, summarizing results, and generating figures.

* `00_all_MAGs_assembly/`: **Primary analysis in manuscript** analyses based on the full MAG dataset, including phylogenetic signal tests, betaNTI, RCBC, random forest analyses (not shown in paper), and supporting ordination/distance-based analyses.
* `01_SingleM_assembly/`: analyses based on SingleM OTUs, including tree construction, input preparation, and betaNTI calculations.
* `02_Module_assembly/`: analyses based on functional modules, including preparation of module inputs, betaNTI, RCBC, and combined result summaries.
* `03_Misc_figures/`: supporting scripts for additional figure panels and exploratory analyses, including specialization, phylofactorization, and tree-based plotting.
* `helper_scripts/`: helper functions used by the assembly-analysis workflow.

### `bash/`

Shell scripts used for command-line or high-performance-computing steps that support the R analyses, including phylogenetic tree construction for MAGs and SingleM genes.

### `slurm/`

SLURM job scripts for running computationally intensive analyses on an HPC cluster. These include jobs for phylogenetic signal testing, GTDB tree construction, betaNTI, RCBC, SingleM protein-tree workflows, module-based assembly analyses, and random forest analyses.

Before running these scripts on a new cluster, check and edit:

* account/partition settings,
* requested memory and wall time,
* module or conda activation commands,
* input and output paths,
* job-array settings, if applicable.

### `outputs/`

Selected analysis outputs used by downstream scripts or manuscript figures. These include betaNTI and RCBC outputs, SingleM OTU relative-abundance and taxonomy tables, representative-sequence FASTA files, alignments, trees, and module-specific output folders.

Some large or restricted files may not be included in the public repository. See the top-level repository README and `.gitignore` for required external data inputs.

## Setup

From the top level of the repository, create and activate the project conda environment:

```bash
./install_dependencies.sh
conda activate temporal_paper_vX
```

Replace `temporal_paper_vX` with the current environment name/version created by `install_dependencies.sh`.

Most R scripts assume that they are run from the repository root or that paths are resolved using the shared project setup. In new analysis scripts, source the top-level setup file before loading inputs:

```r
source("setup.R")
```

### R package requirements by analysis subdirectory

The table below summarizes R packages used by scripts in each subdirectory of `Assembly-analysis/R/`, in general each sub directory represents an independent analysis path.

| Subdirectory              | Main analysis purpose                                                                                                                                                                       | Required R packages                                                                                                                                                                                                                                                                                                       |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `R/00_all_MAGs_assembly/` | MAG-level phylogenetic signal tests, betaNTI, RCBC, pairwise predictor preparation, random forest analyses, random forest plotting, dbRDA/MRM analyses, and assembly-result interpretation. | `tidyverse`, `here`, `knitr`, `ape`, `vegan`, `analogue`, `parallel`, `picante`, `tidytree`, `ggtree`, `ggnewscale`, `ggstance`, `GGally`, `viridis`, `cowplot`, `multcompView`, `data.table`, `iterators`, `Rmpi`, `doMPI`, `doParallel`, `ranger`, `fastshap`, `caret`, `janitor`, `RColorBrewer`, `ecodist`, `shapviz` |
| `R/01_SingleM_assembly/`  | SingleM OTU table preparation, SingleM tree handling, SingleM taxonomy preparation, and SingleM betaNTI analyses.                                                                           | `tidyverse`, `here`, `knitr`, `vegan`, `viridis`, `cowplot`, `ape`, `phytools`, `picante`, `Rmpi`, `doMPI`, `doParallel`, `iCAMP`                                                                                                                                                                                         |
| `R/02_Module_assembly/`   | Functional/module-level input preparation, module-specific tree plotting, module betaNTI, module RCBC, and combined module assembly summaries.                                              | `tidyverse`, `here`, `knitr`, `ggrepel`, `ape`, `picante`, `vegan`, `tidytree`, `ggtree`, `ggnewscale`, `ggstance`, `data.table`, `iterators`, `Rmpi`, `doMPI`, `doParallel`, `viridis`                                                                                                                                   |
| `R/03_Misc_figures/`      | Supporting figure generation, network centrality/specialization plots, phylofactorization, tree plots, and basic specialization analyses.                                                   | `tidyverse`, `here`, `gganimate`, `ggrepel`, `ape`, `tidytree`, `ggtree`, `ggnewscale`, `ggstance`, `viridis`, `phylofactor`, `vegan`, `GUniFrac`, `cowplot`, `multcompView`, `RColorBrewer`, `rstatix`, `av`                                                                                                             |
| `R/helper_scripts/`       | Helper functions for phylogenetic trait-conservatism calculations.                                                                                                                          | `data.table`, `adephylo`, `ape`                                                                                                                                                                                                                                                                                           |


## Input data

The assembly analyses depend on processed EMERGE data products prepared outside this directory. Required public or shareable inputs should be linked into the top-level `data/` directory as described in the repository README. Restricted or large data products are not committed directly to GitHub.

Typical inputs include:

* EMERGE MAG data,
* SingleM OTU tables,
* functional annotation/distillate files,
* metadata tables,
* processed metatranscriptomic or abundance data,
* phylogenetic trees and representative-sequence files.

## Suggested workflow

The exact workflow depends on which analysis is being regenerated, but the scripts are generally numbered in the order they should be run.

### 1. MAG-based assembly analyses

Use scripts in, run in order of numbers preceeding script:

```text
R/00_all_MAGs_assembly/
```

These scripts prepare pairwise data, calculate betaNTI and RCBC metrics, summarize null-model results, run random forest analyses on Assembly output (not shown in paper), and generate plots from MAG-based assembly outputs.

### 2. SingleM-based assembly analyses

Use scripts in:

```text
bash/
R/01_SingleM_assembly/
slurm/
```

The bash and SLURM scripts support tree creation and alignment steps. The R scripts prepare SingleM inputs and calculate SingleM-based betaNTI.

### 3. Module-based assembly analyses

Use scripts in:

```text
R/02_Module_assembly/
slurm/
```

These scripts prepare module-level abundance/function inputs, calculate betaNTI and RCBC metrics for modules, and combine results for downstream interpretation.

### 4. Supporting figures and exploratory analyses

Use scripts in:

```text
R/03_Misc_figures/
```

These scripts generate additional tree-based, phylogenetic, and specialization-related analyses used for supporting figures or interpretation.

## Running analyses

For smaller analyses, scripts can generally be run interactively from R or with `Rscript` from the repository root:

```bash
Rscript Assembly-analysis/R/00_all_MAGs_assembly/03_prepare_pairwise_data.R
```

For larger analyses, submit the corresponding SLURM script:

```bash
sbatch Assembly-analysis/slurm/02_run_betaNTI.slurm
```

Check each script before running to confirm that paths, environment activation, and cluster-specific settings are correct.

## Outputs

Outputs generated by these workflows should be written to analysis-specific output directories. Where possible, intermediate files should be named descriptively and grouped by analysis type, for example:

```text
outputs/
├── all_module_otus/
├── bog_module_otus/
├── fen_module_otus/
├── palsa_module_otus/
├── module_assembly/
└── consentrait/
```

## Reproducibility notes

* Run scripts from the repository root unless otherwise specified.
* Source `setup.R` to load shared paths and project settings.
* Verify that required data symlinks exist in `data/`.
* Confirm that cluster-specific SLURM options match the system being used.

## Manuscript citation

If using or adapting these scripts, please cite the associated manuscript:

Cronin, D. R., Holland-Moritz, H., Smith, D. A., Aroney, S. T. N., Hodgkins, S. B., Borton, M., Li, Y.-F., Zayed, A., Healy, K., Persson, A., IsoGenie Field & Analytic Teams 2010–2017, EMERGE Institute Coordinators, Tfaily, M. M., Crill, P., McCalley, C. K., Wrighton, K., Varner, R. K., Tyson, G. W., Woodcroft, B. J., Bagby, S. C., Ernakovich, J., and Rich, V. I. **Stable states in an unstable landscape: microbial resistance at the front line of climate change.** bioRxiv 2025.02.07.636677.
