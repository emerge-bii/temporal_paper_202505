# *M. stor* (PLGY01) Network & VIP Analysis

Analysis of multilayer network persistence and environmental correlations for *Methanoflorens stordalenirensis* (*M. stor* / `PLGY01`) in the Stordalen Mire Bog dataset.

This pipeline combines 2011–2017 multilayer network connectivity with PLS Variable Importance in Projection (VIP) scores and Spearman correlations to analyze environmental drivers for *M. stor* and its primary co-occurring neighbors.

## Directory Structure

```text
bog-mtn-associations-mstor/
├── mstor-plots.Rmd       # Main R Markdown analysis & figure generation script
├── Input/                # Input data (VIP values, correlation tables, ID mappings)
└── Output/               # Output directory for exported SVG/PDF plots
```

## Inputs

- **Network file**: `../visualization/Input/multinet_network_input_bog.txt`
- **VIP & Correlation files**:
  - `Input/Bog.08.VIP_values_with_CO2.mM__blue.txt`
  - `Input/Bog.08.VIP_values_with_N.percent__blue.txt`
  - `Input/Bog.08.VIP_values_with_alphaCblue.txt`
  - `Input/Bog.08.VIP_values_with_d13C_CO2__blue.txt`
  - `Input/all_vip_env_var_corr_DSedits.xlsx`

## Outputs

- **`Output/`**: Individual PDF and SVG plots (e.g., `plot_PLGY01.pdf`) showing VIP scores and environmental correlations (CO2, %N, alphaC) for *M. stor* and key neighbors.
- **`edge_attributes.csv`**: Cytoscape edge attribute file formatted for network visualization.

## Running the Analysis

You can run the script interactively in RStudio or render it from the terminal:

```bash
Rscript -e "rmarkdown::render('mstor-plots.Rmd')"
```
