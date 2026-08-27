# M. stor (PLGY01) Network and VIP Analysis

This directory contains analyses of multilayer network and environmental correlations for *Methanoflorens stordalenirensis* (*M. stor* / actor ID `PLGY01`) in the Stordalen Mire Bog habitat. 

The workflow integrates multilayer network data (years 2011–2017) with Variable Importance in Projection (VIP) scores and Spearman correlations for environmental parameters to characterize the niche and key associations of *M. stor*'s primary neighbors.

## Directory contents

```text
bog-mtn-associations-mstor/
├── mstor-plots.Rmd       # Core analysis and plotting script
├── Input/                # VIP scores, correlation tables, and ID references
└── Output/               # Directory for generated PDF/SVG plots
```

## Script overview

| Script | Purpose |
| :--- | :--- |
| `mstor-plots.Rmd` | Processes the Bog multilayer network, extracts the neighbors of `PLGY01`, merges network statistics with VIP and correlation values, generates publication-quality plots (PDF and SVG) for key neighbors, and exports Cytoscape edge attributes. |

## Expected inputs

The analysis script expects the following input files:

### Network data
* `../visualization/Input/multinet_network_input_bog.txt` (the multilayer network representation of the bog)

### VIP values and correlation tables
* `Input/Bog.08.VIP_values_with_CO2.mM__blue.txt`
* `Input/Bog.08.VIP_values_with_N.percent__blue.txt`
* `Input/Bog.08.VIP_values_with_alphaCblue.txt`
* `Input/Bog.08.VIP_values_with_d13C_CO2__blue.txt`
* `Input/all_vip_env_var_corr_DSedits.xlsx`
* `Input/heptapus_ids.txt`

## Outputs

Running the script creates/populates:
* **`Output/`**: A directory containing individual SVG and PDF plots (e.g. `plot_PLGY01.svg`, `plot_3300037175_6.svg`) showing VIP values and Spearman's rho correlations for environmental variables (CO_2 mM, %N, alphaC) for *M. stor* and its close network neighbors.
* **`edge_attributes.csv`**: A CSV file containing Cytoscape-formatted edge properties (`Source`, `Target`, `edge_line`, `num_layers`) to visualize network associations.

## Workflow

To run the analysis and generate the output files and plots:

```bash
# Render the R Markdown document
Rscript -e "rmarkdown::render('mstor-plots.Rmd')"
```
