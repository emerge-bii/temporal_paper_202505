#' Creating a stability with time figure

#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(ggnewscale)
library(here)


# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("quantify_stability_with_time_figure", "outputs")
figures.fp <- here("quantify_stability_with_time_figure", "figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}
#### ====================================================================== ####
#### Read in data for each part of figure
#### ====================================================================== ####

# Lens 1 - Organisms
# Phylum
# Species
Species_CV <- read_csv(here("data/20221109-cv-mags.csv"))

# Dylan's granularity BH
granularity.raw <- read.csv(here("data/granularity-BH - granularity-corrected-10092023.csv"), skip = 1)

# 
# Lens 2 - Functions

# Lens 3 - Organization

# Assembly data 
assembly_stab <- read_delim(file = here("Assembly-analysis/outputs/mantel_with_time.txt"))

# Network
network_stats %>% names()

# network comparisons from Dylan
network_comparisons <- read_tsv(here("data", "20240403_network_comparisons.txt"))

# Basic alpha diversity, beta diversity
# Use SingleM data


#### ====================================================================== ####

# Setup plotting parameters
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
shape_habitat  <- c(15, 16, 17)
fill_habitat   <- colour_habitat
depth_levels <- c("0-9", "10-19", "20-29", "30-39")
depth_labels <- c("0", "10", "20", "30")
fill_depth   <- RColorBrewer::brewer.pal(5, "YlOrBr")[-1]
colour_depth <- fill_depth
year_levels <- c("2011", "2012", "2013", "2014", "2015", "2016", "2017")
year_labels <- year_levels
fill_year   <- RColorBrewer::brewer.pal(7, "OrRd")
colour_year <- fill_year
assembly_levels <- c("Homogenous selection", "Heterogenous selection",
                     "Homogenizing dispersal", "Dispersal limitation and drift",
                     "Drift")
assembly_labels <- assembly_levels
#colour_assembly <- c("#133253", "#738CA6","#1C4A00","#80BA5D", "#B1DF95") # blue/green
#colour_assembly <- c("#3F002E", "#BE7FAD","#871200","#DA3015", "#FF816D") # purple/red
#colour_assembly <- c("#6C358D", "#AB81C4","#085A4F","#2D8478", "#79BDB4") # purple/teal
colour_assembly <- c("#521168", "#8e318f","#BD4D0C","#FF8000", "#FADBAC") # purple/orange
fill_assembly <- colour_assembly

#### ====================================================================== ####

#### Helper Functions
#### ====================================================================== ####
# Makes triangles to overlay on heatmap
# make_triangles <- function(x, y, point = "up") {
#   x <- as.integer(as.factor((x)))
#   y <- as.integer(as.factor((y)))
#   
#   if (point == "up") {
#     newx <- sapply(x, function(x) {
#       c(x - 0.5, x - 0.5, x + 0.5)
#     }, simplify = FALSE)
#     newy <- sapply(y, function(y) {
#       c(y - 0.5, y + 0.5, y + 0.5)
#     }, simplify = FALSE)
#   } else if (point == "down") {
#     newx <- sapply(x, function(x) {
#       c(x - 0.5, x + 0.5, x + 0.5)
#     }, simplify = FALSE)
#     newy <- sapply(y, function(y) {
#       c(y - 0.5, y - 0.5, y + 0.5)
#     }, simplify = FALSE)
#   }
#   data.frame(x = unlist(newx), y = unlist(newy))
# }

# make placeholders for functions data and organizations data

#### ====================================================================== ####

# Combine and reorder data
#### ====================================================================== ####
# Alpha diversity
if(!file.exists(here("data","singlem_alpha.txt"))) {
  plan(multisession, workers = 3)
  # Genus-level representation table
  singlem_genus <- singlem_otu_tables %>%
    select(temporal_sample_id, coverage, all.genus.binned) %>%
    group_by(temporal_sample_id, all.genus.binned) %>%
    summarise(coverage = sum(coverage)) %>%
    pivot_wider(names_from = all.genus.binned, values_from = coverage) %>%
    mutate(prop_genus = `1` / (`0` + `1`)) %>%
    select(temporal_sample_id, prop_genus) %>%
    left_join(sample_metadata %>% select(temporal_sample_id, Habitat__, DepthLumping, Year__)) %>%
    filter(!is.na(Habitat__))
  
  both_domain_genes <- c(
    "S3.1.ribosomal_protein_L2_rplB",
    "S3.2.ribosomal_protein_L3_rplC",
    "S3.3.ribosomal_protein_L14b_L23e_rplN",
    "S3.4.ribosomal_protein_L16_L10E_rplP",
    "S3.5.ribosomal_protein_S2_rpsB",
    "S3.6.ribosomal_protein_S5",
    "S3.7.ribosomal_protein_S7",
    "S3.8.ribosomal_protein_S12_S23",
    "S3.9.ribosomal_protein_S15P_S13e",
    "S3.10.ribosomal_protein_S19_rpsS",
    "S3.11.pheS",
    "S3.12.ribosomal_L1",
    "S3.13.ribosomal_S9"
  )
  
  singlem_alpha <- singlem_otu_tables %>%
    filter(gene %in% both_domain_genes) %>%
    select(temporal_sample_id, gene, sequence, coverage) %>%
    group_by(temporal_sample_id, gene) %>%
    nest() %>%
    mutate(nseq = map_dbl(data, nrow)) %>%
    filter(nseq >= 100) %>% as.data.table() %>%
    mutate(
      sample = furrr::future_map(data,
                                 ~ pivot_wider(., names_from = sequence, values_from = coverage)),
      richness = map_dbl(sample, ~ specnumber(.)),
      shannon = map_dbl(sample, diversity, index = "shannon"),
      pielou = shannon / log(richness),
      simpson = map_dbl(sample, diversity, index = "simpson"),
      invsimpson = map_dbl(sample, diversity, index = "invsimpson"),
    ) %>%
    group_by(temporal_sample_id) %>%
    summarise(
      richness = mean(richness),
      shannon = mean(shannon),
      pielou = mean(pielou),
      simpson = mean(simpson),
      invsimpson = mean(invsimpson)
    )
  
  write_tsv(x = singlem_alpha, here("data","singlem_alpha.txt"))
} else {
  singlem_alpha <- read_tsv(file = here("data","singlem_alpha.txt"))
}

# Calculate significance of alpha diversity regression over time
singlem_alpha_directional <- singlem_alpha %>%
  left_join(sample_metadata %>% select(temporal_sample_id, Habitat__, DepthLumping, Year__)) %>%
  filter(!is.na(Habitat__)) %>%
  select(Habitat__, DepthLumping, Year__, temporal_sample_id, richness, shannon, pielou, simpson, invsimpson) %>%
  pivot_longer(richness:invsimpson, names_to = "DiversityIndex", values_to = "Diversity") %>%
  group_by(Habitat__, DepthLumping, `DiversityIndex`) %>%
  nest() %>%
  mutate(lm.model = purrr::map(data, ~lm(Diversity~Year__, data = .)),
         lm.tidied = purrr::map(lm.model, broom::tidy),
         lm.summaried = purrr::map(lm.model, summary),
         lm.model.rsqr = map_dbl(lm.model, ~broom::glance(.)$r.squared),
         lm.model.pval = map_dbl(lm.model, ~broom::glance(.)$p.value)) %>%
  mutate(lm.padj = p.adjust(lm.model.pval, method = "BH")) %>%
  mutate(Test_direction = "Directional") %>%
  mutate(`DiversityIndex` = gsub("_", " ", str_to_title(`DiversityIndex`)))

# select only alpha diversity mentioned in text: richness, evenness, shannon's div.
# match table to style of granularity table
singlem_alpha_directional_gran <- singlem_alpha_directional %>%
  filter(DiversityIndex %in% c("Richness", "Pielou", "Shannon")) %>%
  mutate(DiversityIndex = gsub("Pielou", replacement = "Pielou's evenness", DiversityIndex),
         DiversityIndex = gsub("$", " - strain", DiversityIndex),
         DiversityIndex = gsub("^", "Diversity - ", DiversityIndex)) %>%
  select(Habitat__, DepthLumping, DiversityIndex, lm.model.rsqr, lm.padj, Test_direction) %>%
  rename(granularity = DiversityIndex,
         Habitat = Habitat__, Depth.Interval = DepthLumping,
         R2 = lm.model.rsqr, corrected_.pvalue = lm.padj) %>%
  # add significance marker
  mutate(significant = if_else(corrected_.pvalue < 0.05, "*", NA)) %>%
  rename(corrected_.pvalue_sigflag = significant)


# filter(DiversityIndex == "shannon") %>%
#   ggplot(aes(x = Year__, y = Diversity)) +
#   geom_point() +
#   stat_smooth(method = "lm", aes(group = DiversityIndex)) +
#   facet_grid(Habitat__ ~ DepthLumping)

# network stability data (from Setup.R)
network_stab_test <- network_stats %>%
  group_by(Habitat, layer) %>%
  summarize(across(matches(c("degree_centrality", "betweeness_centrality", "eigenvector_centrality")),
               ~mean(., na.rm = T))) %>%
  ungroup() %>%
  pivot_longer(contains("centrality"), names_to = "centrality metric", values_to = "AvgCentrality") %>%
  group_by(Habitat, `centrality metric`) %>%
  nest() %>%
  mutate(lm.model = purrr::map(data, ~lm(AvgCentrality~layer, data = .)),
         lm.tidied = purrr::map(lm.model, broom::tidy),
         lm.summaried = purrr::map(lm.model, summary),
         lm.model.rsqr = map_dbl(lm.model, ~broom::glance(.)$r.squared),
         lm.model.pval = map_dbl(lm.model, ~broom::glance(.)$p.value)) %>%
  mutate(lm.padj = p.adjust(lm.model.pval, method = "BH")) %>%
  mutate(Test_direction = "Directional") %>%
  mutate(`centrality metric` = gsub("_", " ", str_to_title(`centrality metric`)))

# More granular network stats
test_time_correlation_net <- function(network_data){
  metadata_filt <- input_ra$sample_metadata %>% 
    select(Year__) %>%
    na.omit() %>%
    # Add in a Time since 2011 column
    mutate(TimeSince2011 = Year__-2011) %>%
    distinct()
  
  # Variants on year to test
  # euc_Year: Simple euclidean distance in years
  euc_Year <- dist(metadata_filt %>% select(Year__), method = "euclidean")
  # bin_Year: Whether two samples are from the same (0) or different years (1), no continuous distance 
  bin_Year <- as.matrix(dist(model.matrix(~0 + as.character(metadata_filt$Year__)), method = "binary"))
  rownames(bin_Year) <- rownames(metadata_filt)
  colnames(bin_Year) <- rownames(metadata_filt)
  # euc_2011_baseline: Whether 2 samples are the same distance from the year 2011 or not; Samples that are a similar distance from
  # 2011 are closer to 0, samples that are dissimilarly distant from 2011 (eg 2012 and 2017) are considered far apart. These results
  # should be functionally similar to euc_Year except they use an arbitrary 3rd year as the baseline
  euc_2011_baseline <- dist(metadata_filt %>% select(TimeSince2011), method = "euclidean")
  
  # Add network stats
  network_stat_mat <- network_data %>%
    pivot_wider(names_from = "Year2", values_from = "value") %>%
    select(Year1, `2011`, everything()) %>%
    column_to_rownames(var = "Year1") %>% as.matrix()

  # Test if Year distances and betaNTI are correlated
  euc_Year_mantel <- ecodist::mantel(as.dist(network_stat_mat) ~ euc_Year, mrank = T)
  
  df_network_stat_return <- data.frame(MantelX = c("EucYear"),
                               MantelY = "NetworkStat",
                               bind_rows(euc_Year_mantel))
  return(df_network_stat_return)
}

Network_compare_granularity <- network_comparisons %>%
  group_by(Habitat, Metric) %>%
  nest() %>%
  mutate(mantel_results = purrr::map(data, ~test_time_correlation_net(network_data = .x))) %>%
  unnest(mantel_results) %>%
  mutate(padj = p.adjust(pval3, method = "fdr")) %>%
  select(Habitat, Metric, mantelr, padj) %>%
  mutate(Test_direction = "Directional",
         Habitat = str_to_title(Habitat),
         mantelr = abs(mantelr)) %>%
  left_join(granularity.raw %>% select(Habitat, Depth.Interval) %>% distinct(),
            by = "Habitat") %>%
  rename(granularity = Metric,
         R2 = mantelr,
         corrected_.pvalue = padj) %>%
  mutate(significant = ifelse(Depth.Interval == "0-9" & corrected_.pvalue < 0.05, "*", NA),
         granularity = gsub("dissimilarity.degree", "Networks - Dissimilarity of degree distribution - coarse", granularity),
         granularity = gsub("rho.degree", "Networks - Rho correlation of degree structure - medium", granularity),
         granularity = gsub("jaccard.triangles", "Networks - Overlap of tripartite relationships - fine", granularity)) %>%
  rename(corrected_.pvalue_sigflag = significant) %>%
  filter(grepl("-", granularity))

# granularity: via coarse-, medium- or fine-resolution metrics of topology: interannual dissimilarity of degree distributions (agnostic to specific edge-node relationships), correlations of degree structure (tracking specific nodes across years), and overlap of tripartite relationships among nodes, respectively.”

network_stab <- network_stab_test %>%
  select(Habitat, `centrality metric`, lm.model.rsqr, lm.padj, Test_direction) %>%
  rename(granularity = `centrality metric`,
  R2 = lm.model.rsqr, corrected_.pvalue = lm.padj) %>%
  left_join(granularity.raw %>% select(Habitat, Depth.Interval) %>% distinct(),
            by = "Habitat") %>%
  # add only one significant marker for all four depths since this is not a depth-resolved metric
  mutate(significant = ifelse(Depth.Interval == "0-9" & corrected_.pvalue < 0.05, "*", NA)) %>%
  rename(corrected_.pvalue_sigflag = significant) %>%
  # Filter out unnecessary metrics
  filter(granularity %in% "Eigenvector centrality")

# Assembly Stability Stats
assembly_stats <- assembly_stab %>% 
#  filter(MantelY == "BetaNTI") %>% 
  select(Habitat__, DepthLumping, MantelY, MantelX, mantelr, padj) %>%
  #rename(mantelr = mantelr) %>%
  mutate(mantelr = abs(mantelr),
         significant = ifelse(padj < 0.05, "*", NA)) %>%
  # Scale metrics
  #mutate_at(c("Metric1"), ~((.-min(.))/(max(.)-min(.)))) %>%
  #pivot_longer(starts_with("Metric"), values_to = "Metric", names_to = "type") %>%
  mutate(Test_direction = ifelse(MantelX == "BinYear", "Non-directional", "Directional"),
         granularity = ifelse(MantelY == "BetaNTI", "Community assembly - deterministic - MAGs", "Community assembly - stochastic - MAGs")) %>%
  rename(Habitat = Habitat__, Depth.Interval = DepthLumping, corrected_.pvalue_sigflag = significant,
         R2 = mantelr, corrected_.pvalue = padj)


## Combine assembly + network + raw granularity data into stats
granularity <- granularity.raw %>%
  select(Habitat, Depth.Interval, granularity, corrected_.pvalue_.numeric, R2_.numeric) %>%
  mutate(Test_direction = "Directional") %>%
  rename(corrected_.pvalue = corrected_.pvalue_.numeric, R2 = R2_.numeric) %>%
  # Hacky pivot longer - put the categorical data as rows below the numeric
  # numeric is directional, categorical is non-directional
  bind_rows(granularity.raw %>%
              select(Habitat, Depth.Interval, granularity, corrected_.pvalue_.categorical, R2_.categorical) %>%
              mutate(Test_direction = "Non-directional") %>%
              rename(corrected_.pvalue = corrected_.pvalue_.categorical, R2 = R2_.categorical)) %>%
  # fix order of factors
  mutate(Habitat = factor(Habitat, levels = habitat_levels),
         Depth.Interval = factor(Depth.Interval, levels = depth_levels),
         Depth.Interval = fct_rev(Depth.Interval),
         granularity = gsub("singlem-", "Composition - Marker genes - ", granularity),
         granularity = gsub("MAGs", "Composition - MAGs", granularity),
         granularity = ifelse(granularity == "genes", "Composition - Protein Clusters", granularity),
         granularity = ifelse(granularity == "Metabolism", "Composition - Pathways", granularity),
         granularity = factor(granularity, levels = c("Composition - Marker genes - strain", "Composition - MAGs", 
                                                      "Composition - Marker genes - genus", 
                                                      "Composition - Marker genes - phylum", "Composition - Protein Clusters", "Composition - Pathways"))) %>%
  # add significance symbols
  mutate(corrected_.pvalue_sigflag = ifelse(corrected_.pvalue < 0.05, "*", NA)) %>%
  bind_rows(assembly_stats) %>%
#  bind_rows(network_stab) %>%
  bind_rows(Network_compare_granularity) %>%
  bind_rows(singlem_alpha_directional_gran) %>%
  # fix the order and names of things
  mutate(granularity = factor(granularity, levels = c(# organisations
                                                      "Community assembly - deterministic - MAGs", 
                                                      "Community assembly - stochastic - MAGs",
                                                      "Networks - Dissimilarity of degree distribution - coarse",
                                                      "Networks - Rho correlation of degree structure - medium",
                                                      "Networks - Overlap of tripartite relationships - fine",
                                                      # organisms
                                                      "Diversity - Shannon - strain",
                                                      "Diversity - Pielou's evenness - strain",
                                                      "Diversity - Richness - strain",
                                                      "Composition - Marker genes - strain", "Composition - MAGs", 
                                                      "Composition - Marker genes - genus", "Composition - Marker genes - phylum",
                                                      # functions
                                                      "Composition - Protein Clusters",
                                                      "Composition - Pathways"))) %>%
  mutate(Habitat = factor(Habitat, levels = habitat_levels),
         Depth.Interval = factor(Depth.Interval, levels = depth_levels),
         Depth.Interval = fct_rev(Depth.Interval),
         Test_direction = ifelse(Test_direction == "Non-directional", "Fluctuation over time", "Change with time")) %>%
  # Add column designating lenses (not used in plotting, but helpful for those looking at table)
  mutate(lens = ifelse(as.numeric(granularity) <= 5, "Organisations",
                       ifelse(as.numeric(granularity) <= 12 & as.numeric(granularity) > 5, "Organisms",
                       "Functions")))


#### ====================================================================== ####

#### ====================================================================== ####
# Export plot data for Kieran
write.table(granularity, file = paste0(outputs.fp, "/granularity_fig_table.txt"),
            sep = "\t", row.names = F)

#### ====================================================================== ####
# 
# # Revised Granularity Figure after 11/15/2023 discussion
# #### ====================================================================== ####
# 
# # Base layout of trustus plot:
# trustus_layers <- list(
#     geom_tile(aes(fill = Habitat, alpha = R2)),
#     scale_fill_manual(values = colour_habitat, breaks = habitat_levels, na.value = "grey50"),
#     scale_alpha_continuous(name = expression(atop(Metric,
#                                                   (Mantel~r~or~R^{2}))),
#                            limits=c(0, 1), breaks=seq(from = 0, to = 1, by = 0.2), na.value = 0),
#     geom_text(aes(label = corrected_.pvalue_sigflag), size = rel(10), vjust = 0.8),
#     facet_grid(Habitat~`Test_direction`, space = "free_x"),
#     xlab(""), 
#     ylab("Depth (cm)"),
#     theme_bw(),
#     scale_x_discrete(expand = c(0,0)),
#     scale_y_discrete(expand = c(0,0)),
#     guides(fill = guide_legend(nrow = 1, byrow = TRUE,
#                                direction = "horizontal", 
#                                title.position = "left", title.hjust = 0, title.vjust = 1,
#                                title.theme = element_text(size = rel(15)), 
#                                label.position = "bottom",
#                                label.theme = element_text(angle = 45, size = rel(15), hjust = 0, vjust = 0),
#                                text = element_text(size = rel(12))),
#            alpha = guide_legend(nrow = 1, byrow = TRUE,
#                                direction = "horizontal", 
#                                title.position = "left", title.hjust = 0, title.vjust = 1,
#                                title.theme = element_text(size = rel(15)), 
#                                label.position = "bottom",
#                                label.theme = element_text(size = rel(15), angle = 45, hjust = 0, vjust = 0),
#                                text = element_text(size = rel(12)))),
#     theme(text = element_text(size = rel(4)),
#           axis.text.x = element_text(size = rel(5), angle = 45, hjust = 1, face = "bold"),
#           axis.text.y = element_text(size = rel(5), face = "bold"),
#           axis.title = element_text(size = rel(1), color = "black", face = "bold"),
#           panel.grid = element_blank(),
#           strip.background.y = element_blank(),
#           strip.background.x = element_rect(color = "white", fill = "white"),
#           strip.text.y = element_blank(),
#           strip.text.x = element_text(face = "bold", color = "black"),
#           legend.position = "bottom",
#           legend.justification = "center",
#           legend.key.size = unit(rel(1.5), "lines"),
#           plot.margin = margin(10, 10, 10, 100))
# )
# 
# label_size <- rel(5)
# trustus_plot <- ggplot(data = granularity, aes(y = Depth.Interval, x = granularity)) +
#   trustus_layers +
#   # add line and labels to delinieate lenses (can't be in trustus layers b/c changes with plot subset)
#   geom_vline(xintercept = 5.5, linewidth = rel(2)) +
#   geom_vline(xintercept = 9.5, linewidth = rel(2)) #+
# #  annotate("text", label = expression(bold(Organisms)), x = 7.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Organisation)), x = 2.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Functions)), x = 10.5, y = 4.3, color = "grey30", size = label_size)
# trustus_plot
# 
# ggsave(trustus_plot, 
#        filename = paste0(figures.fp, "/trustus.png"),
#        width = 10, height = 10, dpi = 400)
# 
# 
# # Subset so to only include linear change with time
# trustus_plot_cwt <- granularity %>% filter(Test_direction == "Change with time") %>%
#   ggplot(aes(y = Depth.Interval, x = granularity)) +
#   trustus_layers +
#   # add line and labels to delinieate lenses
#   geom_vline(xintercept = 5.5, linewidth = rel(2)) +
#   geom_vline(xintercept = 9.5, linewidth = rel(2)) #+
# #  annotate("text", label = expression(bold(Organisms)), x = 7.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Organisation)), x = 2.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Functions)), x = 10.5, y = 4.3, color = "grey30", size = label_size)
# trustus_plot_cwt
# 
# ggsave(trustus_plot_cwt, 
#        filename = paste0(figures.fp, "/trustus_Directional_only.png"),
#        width = 10, height = 10, dpi = 400)
# 
# # Subset to only include fluctuation over time
# trustus_plot_fot <- granularity %>% filter(Test_direction == "Fluctuation over time") %>%
#   ggplot(aes(y = Depth.Interval, x = granularity)) +
#   trustus_layers +
#   # add line and labels to delinieate lenses
#   geom_vline(xintercept = 2.5, linewidth = rel(2)) +
#   geom_vline(xintercept = 6.5, linewidth = rel(2)) #+
# #  annotate("text", label = expression(bold(Organisms)), x = 4.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Organisation)), x = 1.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Functions)), x = 7.5, y = 4.3, color = "grey30", size = label_size)
# trustus_plot_fot
# 
# ggsave(trustus_plot_fot, 
#        filename = paste0(figures.fp, "/trustus_fluctuation_only.png"),
#        width = 10, height = 10, dpi = 400)
# #### ====================================================================== ####


# Revised Granularity Figure after 03/2024 discussions
#### ====================================================================== ####

# Base layout of trustus plot:

trustus_layers <- list(
  geom_tile(aes(fill = Habitat, alpha = R2)),
  scale_fill_manual(values = colour_habitat, breaks = habitat_levels, na.value = "white"),
  scale_alpha_continuous(name = expression(atop(Metric,
                                                (Mantel~r~or~R^{2}))),
                         limits=c(0, 1), 
                         #breaks=seq(from = 0, to = 1, by = 0.1),
                         na.value = 0),
  geom_text(aes(label = corrected_.pvalue_sigflag), color = "white", size = rel(17), fontface = "bold", vjust = 0.8),
  facet_grid(lens~Habitat, scales = "free_y", space = "free_y"),
  xlab("Depth (cm)"), 
  ylab(""),
  theme_bw(),
  scale_x_discrete(expand = c(0,0)),
  scale_y_discrete(expand = c(0,0)),
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             direction = "horizontal", 
                             title.position = "left", title.hjust = 0, title.vjust = 1,
                             title.theme = element_text(size = rel(5)), 
                             label.position = "bottom",
                             label.theme = element_text(angle = 45, size = rel(5), hjust = 0, vjust = 0),
                             text = element_text(size = rel(3))),
         alpha = guide_legend(nrow = 1, byrow = TRUE,
                              direction = "horizontal", 
                              title.position = "left", title.hjust = 0, title.vjust = 1,
                              title.theme = element_text(size = rel(5)), 
                              label.position = "bottom",
                              label.theme = element_text(size = rel(5), angle = 45, hjust = 0, vjust = 0),
                              text = element_text(size = rel(3)))),
  theme(text = element_text(size = rel(4)),
        axis.text.x = element_text(size = rel(5), angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = rel(5), face = "bold"),
        axis.title = element_text(size = rel(1), color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background.y = element_blank(),
        strip.background.x = element_rect(color = "white", fill = "white"),
        strip.text.y = element_text(size = rel(5),face = "bold", color = "black"),
        strip.text.x = element_text(size = rel(5), face = "bold", color = "black"),
        legend.position = "bottom",
        legend.justification = "center",
        legend.key.size = unit(rel(1), "lines"),
        plot.margin = margin(10, 10, 10, 100))
)

label_size <- rel(5) #rel(5)
# trustus_plot <- ggplot(data = granularity, aes(y = Depth.Interval, x = granularity)) +
#   trustus_layers +
#   # add line and labels to delinieate lenses (can't be in trustus layers b/c changes with plot subset)
#   geom_vline(xintercept = 5.5, linewidth = rel(2)) +
#   geom_vline(xintercept = 9.5, linewidth = rel(2)) #+
# #  annotate("text", label = expression(bold(Organisms)), x = 7.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Organisation)), x = 2.5, y = 4.3, color = "grey30", size = label_size) +
# #  annotate("text", label = expression(bold(Functions)), x = 10.5, y = 4.3, color = "grey30", size = label_size)
# trustus_plot
# 
# ggsave(trustus_plot, 
#        filename = paste0(figures.fp, "/trustus.png"),
#        width = 10, height = 10, dpi = 400)


# Subset so to only include linear change with time
trustus_plot_cwt <- granularity %>% filter(Test_direction == "Change with time") %>%
  filter(!grepl("Networks - [RO]", granularity)) %>%
  mutate(lens = factor(lens, levels = c("Organisms", "Functions", "Organisations")),
         Depth.Interval = factor(Depth.Interval, levels = depth_levels)) %>%
  ggplot(aes(x = Depth.Interval, y = granularity)) +
  trustus_layers +
  theme(legend.position = "none")
trustus_plot_cwt

ggsave(trustus_plot_cwt, 
       filename = paste0(figures.fp, "/trustus_Directional_only_no_networkfine_med.png"),
       width = 18, height = 10, dpi = 400)

trustus_legend <- trustus_plot_cwt + 
  theme(legend.position = "bottom",
        plot.margin = margin(100, 100, 100, 100))


ggsave(trustus_legend, 
       filename = paste0(figures.fp, "/trustus_Directional_only_legend.png"),
       width = 20, height = 10, dpi = 400)

# Subset to only include fluctuation over time
trustus_plot_fot <- granularity %>% filter(Test_direction == "Fluctuation over time") %>%
  mutate(lens = factor(lens, levels = c("Organisms", "Functions", "Organisations")),
         Depth.Interval = factor(Depth.Interval, levels = depth_levels)) %>%
  ggplot(aes(x = Depth.Interval, y = granularity)) +
  trustus_layers + 
  theme(legend.position = "none")
  # trustus_layers +
  # # add line and labels to delinieate lenses
  # geom_vline(xintercept = 2.5, linewidth = rel(2)) +
  # geom_vline(xintercept = 6.5, linewidth = rel(2)) #+
#  annotate("text", label = expression(bold(Organisms)), x = 4.5, y = 4.3, color = "grey30", size = label_size) +
#  annotate("text", label = expression(bold(Organisation)), x = 1.5, y = 4.3, color = "grey30", size = label_size) +
#  annotate("text", label = expression(bold(Functions)), x = 7.5, y = 4.3, color = "grey30", size = label_size)
trustus_plot_fot

ggsave(trustus_plot_fot, 
       filename = paste0(figures.fp, "/trustus_fluctuation_only.png"),
       width = 16, height = 10, dpi = 400)
#### ====================================================================== ####


#### ====================================================================== ####
# MAG time response
#### ====================================================================== ####
library(propr)

# Mag list from slide 161 in https://docs.google.com/presentation/d/18bqCpz-xO4J3HrTOoeVfJ4ua9rSXZT_qf3CD-a-qYuw/edit#slide=id.g25655bf0596_0_5
mag_list <- c("20120700_E2D_6",
              "20110800_E1D_6",
              "20120700_S1X_2",
              "20110800_E3D_26",
              "20100900_E1D_4",
              "3300036929_5",
              "20120600_E1D_53",
              "20120800_E3X_27",
              "20120600_E2D_15",
              "PMBS01")
mags_of_interest <- mag_list
abund_TM_RA_clr <- trimmed_mean$rel_abund %>%
  # translate dataframe to work with propr expectation
  column_to_rownames(var = "genome") %>%
  t() %>% as_tibble(rownames = "temporal_sample_id") %>%
  column_to_rownames(var = "temporal_sample_id") %>% 
  # Filter out MAGs with >90% zero abundance
  #dplyr::select(where(~ sum(. != 0) > length(.) / 10)) %>% 
  # apply zeros transformation corrections
  mutate(across(everything(), transform_perc)) %>%
  # apply clr transformation
  propr(metric = "rho") %>%
  pluck("logratio") %>%
  rownames_to_column(var = "temporal_sample_id") %>%
  select(temporal_sample_id, any_of(mags_of_interest)) %>%
  pivot_longer(-temporal_sample_id, names_to = "genome", values_to = "CLR_cumu_relabund")

# Combine with metadata
# Assemble all relevant data about VIPs
otu_otu_resp <- sample_metadata %>%
  dplyr::select(temporal_sample_id, Year__, Habitat__, T_soil.deg_C, DepthAvg__, DepthLumping ) %>%
  # only keep samples from Fen 30
  dplyr::filter(DepthLumping == "30-39" & Habitat__ == "Fen") %>% 
  # Filter abundance table for samples in the sample metadata and with soil temp data
  left_join(abund_TM_RA_clr, by = c("temporal_sample_id")) %>% 
  left_join(input_ra$taxonomy, by = "genome") %>%
  # Fix taxonomy names
  mutate(across(Domain:Species, ~gsub("^[pcofgs]__", "", .))) %>%
  # Change the name of the abundance column for convenience
  rename(coverage = CLR_cumu_relabund) %>%
  # Fix habitat levels
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>% 
  # Use mutate to standardize our explantory variables of interest
  mutate(across(all_of(c("T_soil.deg_C", "DepthAvg__")), ~scale(.)[,1], .names = "{.col}_scl")) %>%
  # Calculate log10 of coverage for models where residuals are not normally distributed
  mutate(coverage_log10 = log10(coverage + 1 - min(coverage)),  # Because we're using clr transformed data, we will translate the data to have no negatives first by adding the minimum value to each group
         coverage_cuberoot = sign(coverage) * (abs(coverage))^(1/3))


# Run correlations
correlations_otus_time <- otu_otu_resp %>%
  group_by(genome) %>%
  nest() %>%
  mutate(cor.model = purrr::map(data, ~cor.test(~coverage + Year__, data = ., method = "spearman")),
         cor.tidied = purrr::map(cor.model, broom::tidy)) %>%
  unnest(cor.tidied) %>% ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

# note: genome = "20120600_E2D_15" cannot compute p-value b/c ties
correlations_otus_time

# Plot the data
overview_plot <- otu_otu_resp %>%
  ggplot(aes(y = coverage, x = Year__)) + 
  geom_point(aes(fill = DepthLumping ), shape = 21, size = 3, alpha = 0.5) + 
  facet_wrap(Habitat__~genome, scales = "free_y", nrow = 2)

if(!dir.exists(paste0(figures.fp, "/", dataType))) {dir.create(paste0(figures.fp,"/", dataType))}

ggsave(overview_plot, 
       filename = paste0(figures.fp, "/", "10_MAG_time_correlations", "_overview_plot.png"),
       device = "png", height = 14, width = 20)


#### ====================================================================== ####
# GAMs
#### ====================================================================== ####
########################################
### Bray-Curtis dissimilarity matrix ###
########################################
get_dis_matrix <- function(rel_abund) {
  dis_matrix <- rel_abund %>%
    t() %>%
    vegdist(method = "bray")
  
  return(dis_matrix)
}

dis_matrix <- trimmed_mean$rel_abund %>%
  select(all_of(sample_metadata$temporal_sample_id)) %>%
  get_dis_matrix()

dis_matrix_palsa <- trimmed_mean$rel_abund %>%
  select(all_of(
    sample_metadata %>%
      filter(Habitat__ == "Palsa") %>%
      `$`(temporal_sample_id)
  )) %>%
  get_dis_matrix()

dis_matrix_bog <- trimmed_mean$rel_abund %>%
  select(all_of(
    sample_metadata %>%
      filter(Habitat__ == "Bog") %>%
      `$`(temporal_sample_id)
  )) %>%
  get_dis_matrix()

dis_matrix_fen <- trimmed_mean$rel_abund %>%
  select(all_of(
    sample_metadata %>%
      filter(Habitat__ == "Fen") %>%
      `$`(temporal_sample_id)
  )) %>%
  get_dis_matrix()



# PCOA Plot Function
plot_pcoa_habitat <- function(dis_matrix, habitat, legible = TRUE, output_dir = figures.fp) {
  abund_pcoa_raw <- ape::pcoa(dis_matrix) # this is from the package ape; 
  abund_pcoa.pct_ex <- round((abund_pcoa_raw$values$Relative_eig) * 100, 1)
  
  abund_pcoa <- abund_pcoa_raw$vectors[,1:2] %>% data.frame() %>%
    rownames_to_column(var = "temporal_sample_id") %>%
    left_join(sample_metadata, by = c("temporal_sample_id")) %>%
    filter(!is.na(Habitat__))
  
  ################
  ### Plotting ###
  ################
  pcoa_plotting_layers <- list(
    stat_ellipse(type = "t", geom = "polygon", alpha = 1/10, colour = NA),
    stat_ellipse(type = "norm", linetype = 2, colour = "black"),
    scale_shape_manual(name = "Habitat", values = shape_habitat, breaks = habitat_levels),
    scale_colour_manual(name = "Habitat", values = colour_habitat, breaks = habitat_levels),
    scale_fill_manual(name = "Habitat", values = fill_habitat, breaks = habitat_levels),
    xlab(paste0("PCoA 1 (", abund_pcoa.pct_ex[1],"%)")), # extract the percent for the first axis
    ylab(paste0("PCoA 2 (", abund_pcoa.pct_ex[2],"%)")), # extract the percent for the second axis
    theme_bw(), 
    theme(legend.text = element_text(size = rel(5)), 
          text = element_text(size = rel(5)))
  )
  
  # Point settings
  point_size <- rel(1)
  point_alpha <- 1
  
  if(legible) {
    nmds_plotting_layers <- list(
      pcoa_plotting_layers,
      theme(legend.text = element_text(size = rel(5)), 
            text = element_text(size = rel(5)))
    )
    # Point settings
    point_size <- rel(5)
    point_alpha <- 0.7
  }
  
  
  abund_pcoa %>%
    ggplot(aes(Axis.1, Axis.2, shape = Habitat__, colour = Habitat__, fill = Habitat__)) +
    pcoa_plotting_layers +
    geom_point()
  ggsave(str_c("PCOA", habitat, "12.png", sep = "_"), 
         height = 7, width = 7,
         path = output_dir, dpi = 900)
  
  abund_pcoa %>%
    ggplot(aes(Axis.1, Axis.2, shape = Habitat__, fill = Habitat__)) +
    pcoa_plotting_layers +
    new_scale_fill() +
    geom_point(aes(shape = Habitat__, colour = NA, fill = DepthLumping), size = point_size, alpha = point_alpha) +
    scale_fill_manual(name = "Depth (cm)", values = colour_depth, breaks = depth_levels) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  ggsave(str_c("PCOA", habitat, "depth.png", sep = "_"), 
         height = 7, width = 8.5,
         path = output_dir, dpi = 900)
  
  abund_pcoa %>%
    ggplot(aes(Axis.1, Axis.2, shape = Habitat__, fill = Habitat__)) +
    pcoa_plotting_layers +
    new_scale_fill() +
    geom_point(aes(shape = Habitat__, colour = NA, fill = factor(Year__)), size = point_size, alpha = point_alpha) +
    scale_fill_manual(name = "Year", values = colour_year, breaks = year_levels) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  ggsave(str_c("PCOA", habitat, "year.png", sep = "_"), 
         height = 7, width = 8,
         path = output_dir, dpi = 900)
  
  return(abund_pcoa)
}

abund_pcoa <- plot_pcoa_habitat(dis_matrix, "")
abund_pcoa_palsa <- plot_pcoa_habitat(dis_matrix_palsa, "palsa")
abund_pcoa_bog <- plot_pcoa_habitat(dis_matrix_bog, "bog")
abund_pcoa_fen <- plot_pcoa_habitat(dis_matrix_fen, "fen")


temp_depth_gam <- gam(T_soil.deg_C ~ s(Axis.1, Axis.2, by = Habitat__) + DepthLumping + Habitat__,
                      data = abund_pcoa, method = "REML")
gam.check(temp_depth_gam)
summary(temp_depth_gam)



abund_pcoa <- mutate(abund_pcoa, Habitat__ = factor(Habitat__, levels = habitat_levels),
                       DepthLumping = factor(DepthLumping, levels = depth_levels))
  
  
  
  temp_only_gam <- gam(T_soil.deg_C ~s(Axis.1, Axis.2),
                       data = abund_pcoa, method = "REML") # adjs. r^2 = 0.342
  gam.check(temp_only_gam)
  summary(temp_only_gam)
  
  temp_gam <- gam(T_soil.deg_C ~Habitat__*DepthLumping + s(Axis.1, Axis.2),
                  data = abund_pcoa, method = "REML") # adjs. r^2 = 0.617
  gam.check(temp_gam)
  summary(temp_gam)
  # Interpretation: once habitat and depth are included in the model, they take up the bulk of the variation
  # in temperature. and the ordination axes are non significant
  vis.gam(x = temp_gam,                # GAM object
          view = c("Axis.1", "Axis.2"),   # variables
          plot.type = "persp", se = T, theta = 250, phi = 0)
  
  time_only_gam <- gam(Year__ ~ s(Axis.1, Axis.2),
                       data = abund_pcoa, method = "REML") # adjs. r^2 = 0.000651
  gam.check(time_only_gam)
  summary(time_only_gam)
  
  time_gam <- gam(Year__ ~Habitat__*DepthLumping + s(Axis.1, Axis.2),
                  data = abund_pcoa, method = "REML") # adjs. r^2 = 0.056
  gam.check(time_gam)
  summary(time_gam)
  
  vis.gam(x = temp_gam,                # GAM object
          view = c("Axis.1", "Axis.2"),   # variables
          plot.type = "persp", se = T, theta = 250, phi = 0)
  
  
  vis.gam(x = temp_gam,                # GAM object
          view = c("Axis.1", "Habitat__"),   # variables
          plot.type = "persp", se = T, theta = 250, phi = 0)
  
  
  time_gam <- gam(Year__ ~ s(Axis.1, Axis.2, by = DepthLumping) + Habitat__ + DepthLumping,
                  data = abund_pcoa, method = "REML") # adjs. r^2 = 0.102
  gam.check(time_gam)
  summary(time_gam)
  
  
  
  temp_gam <- gam(T_soil.deg_C ~ te(Axis.1, Axis.2, DepthLumping, bs = "fs") +
                    s(Axis.1, Axis.2),
                  data = abund_pcoa, method = "REML")
  summary(temp_gam)
  plot(temp_gam, pages = 1)
  
  coef(temp_gam)
  gam.check(temp_gam, cex = 3)
  
  vis.gam(x = temp_gam,                # GAM object
          view = c("Axis.1", "Axis.2"),   # variables
          plot.type = "persp", se = T, theta = 200, phi = 0)
  vis.gam(x = temp_gam,                # GAM object
          view = c("Axis.1", "Axis.2"),   # variables
          plot.type = "persp", se = T, theta = 250, phi = 0)
  
  
  plot(temp_gam, pages = 1, cex = 5)
  
  plot(temp_gam, scheme = 1, cex = 5)
  curvi
  
  summary(temp_gam)
  
  
  # Axis gam
  axis_1_gam <- gam(Axis.1 ~ s(T_soil.deg_C, Habitat__, bs = "fs") + DepthLumping,
                    data = abund_pcoa, method = "REML")
  summary(axis_1_gam)
  plot(axis_1_gam)
  
  
  ggplot(abund_pcoa, aes(x = Axis.1, y = Year__)) +
    geom_point(aes(color = Habitat__)) +
    facet_wrap(DepthLumping ~ Habitat__)
  
  ggplot(abund_pcoa, aes(x = Axis.1, y = T_soil.deg_C)) +
    geom_point(aes(color = Habitat__)) +
    facet_wrap(DepthLumping ~ Habitat__, scales = "free_y")
  
  
  ggplot(abund_pcoa, aes(x = Axis.2, y = Year__)) +
    geom_point(aes(color = Habitat__)) +
    facet_wrap(DepthLumping ~ Habitat__)
  
  ggplot(abund_pcoa, aes(x = Axis.2, y = T_soil.deg_C)) +
    geom_point(aes(color = Habitat__)) +
    facet_wrap(DepthLumping ~ Habitat__, scales = "free_y")
  
  
  plot(temp_gam, pages = 1)
  
  
  
  
  extract.xyz <- function(obj) {
    xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
    xyz <- cbind(xy, c(obj$grid$z))
    names(xyz) <- c("x", "y", "z")
    return(xyz)
  }
  
  contour.vals <- extract.xyz(obj = ordi_temp)
  
}
