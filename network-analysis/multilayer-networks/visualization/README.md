# Multilayer Network Visualization

Visualizes multilayer co-occurrence networks (2011–2017) across the three Stordalen Mire habitats: Palsa, Bog, and Fen.

The pipeline processes annual layer edge lists with `muxViz` to generate multiplex network layout plots, calculates year-over-year degree correlations (Spearman's rho), and compiles them into a publication figure (`Output/network.png` and `network.svg`).

## Directory Structure

```text
visualization/
├── muxviz.Rmd        # Main R Markdown visualization script
├── Input/            # Annual layer edge lists, multinet files, and infomap node files
└── Output/           # Exported figures (network.png, network.svg)
```

## Inputs

- **Layer edge lists**: `Input/layer[1-7]_[palsa|bog|fen]_unique.csv`
- **Multinet network files**: `Input/multinet_network_input_[palsa|bog|fen].txt`
- **Infomap node files**: `Input/infomap_nodes_[palsa|bog|fen].txt`

## Outputs

- **`Output/network.png`** / **`network.svg`**: Combined multi-panel publication figure containing:
  1. 2D multiplex network layouts representing co-occurrence across the 7 study years for each habitat.
  2. Boxplots of year-over-year degree correlations (rho) showing network structural consistency over time.

## Running the Analysis

Render the document in RStudio or run:

```bash
Rscript -e "rmarkdown::render('muxviz.Rmd')"
```
