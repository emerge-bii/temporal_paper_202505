# WGCNA & VIP Analysis

Weighted Gene Co-expression Network Analysis (WGCNA) and Partial Least Squares (PLS) regression for MAG abundance profiles across Stordalen Mire habitats (Palsa, Bog, Fen, and All Combined).

The pipeline constructs signed co-abundance networks, calculates module eigengenes, evaluates module-trait correlations across environmental metadata, and identifies high-VIP driver MAGs associated with key parameters (e.g., CO2, alphaC, soil temperature).

## Directory Structure

```text
wgcna/
├── 01_wgcna_refactored.Rmd   # Main WGCNA & sPLS regression analysis script
├── README.md                 # Workflow overview and execution instructions
└── output/                   # Directory for exported figures, tables, and TOM matrices
```

## Inputs

Expects data cached via `setup.R` (`data/setup_cache.RData`) containing:
- MAG count matrix (`otu_table`)
- Sample environmental metadata (`sample_metadata`)
- GTDB taxonomy table (`taxonomy`)

## Outputs

Generated in `output/`:
- **Soft-thresholding plots**: `Palsa_Thresholding.pdf`, `Bog_Thresholding.pdf`, etc.
- **Module membership & Eigengenes**: `[Habitat]_Module_membership.csv`, `[Habitat].04.MAGs_Module_eigen_values.csv`
- **Module-Trait Heatmaps & Tables**: `[Habitat].05.Module-trait relationships.pdf`, `[Habitat]_Combined_ModuleTrait_Table.csv`, `[Habitat]_LongForm_ModuleTrait_Table.csv`
- **sPLS VIP Driver Tables & Plots**: `[Habitat].08.VIP_values_with_[Parameter][Module].txt`, `[Habitat].11.measured_vs_predicted_[Module]-vs-[Parameter].pdf`

## Running the Analysis

Render the document in RStudio or run:

```bash
Rscript -e "rmarkdown::render('01_wgcna_refactored.Rmd')"
```
