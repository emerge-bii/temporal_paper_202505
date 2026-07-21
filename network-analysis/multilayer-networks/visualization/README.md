# Multilayer Network Visualization

This directory contains scripts and data for visualizing multilayer networks (years 2011–2017): Palsa, Bog, and Fen. 

The workflow constructs multiplex network visualizations using `muxViz` and calculates year-over-year network degree correlation statistics to output a combined, multi-panel publication figure.

## Directory contents

```text
visualization/
├── muxviz.Rmd        # Core visualization and correlation script
├── Input/            # Layer-specific edge lists and multinet input files
└── Output/           # Directory for generated publication figures
```

## Script overview

| `muxviz.Rmd` | Processes year-specific edge lists for Palsa, Bog, and Fen; generates multiplex network layout plots; computes pairwise layer degree correlations (Spearman's rho); and compiles these plots into a single multi-panel figure. |


## Expected inputs

The script expects the following files in the `Input/` directory:

### Edge lists per layer (years 2011–2017)
* `Input/layer1_palsa_unique.csv` to `Input/layer7_palsa_unique.csv` (Palsa edges)
* `Input/layer1_bog_unique.csv` to `Input/layer7_bog_unique.csv` (Bog edges)
* `Input/layer1_fen_unique.csv` to `Input/layer7_fen_unique.csv` (Fen edges)

### Multinet networks
* `Input/multinet_network_input_palsa.txt`
* `Input/multinet_network_input_bog.txt`
* `Input/multinet_network_input_fen.txt`

### Infomap node lists (for community mapping)
* `Input/infomap_nodes_palsa.txt`
* `Input/infomap_nodes_bog.txt`
* `Input/infomap_nodes_fen.txt`

## Outputs

Running the script generates the following publication-quality figures:
* **`Output/network.png`**: High-resolution (600 DPI) multi-panel PNG.
* **`Output/network.svg`**: Multi-panel SVG.

Each figure integrates:
1. 2D multiplex network layouts representing co-occurrence across the 7 study years for each habitat.
2. Boxplots of year-over-year degree correlations (rho) to illustrate network structural consistency over time.

## Workflow

To execute the script and update the figures:

```bash
# Render the R Markdown document
Rscript -e "rmarkdown::render('muxviz.Rmd')"
```
