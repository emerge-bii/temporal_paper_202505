#' ## Calculate Random Forest Assembly
#' This step runs random forest models on assembly process

#+ include=FALSE
# some setup options for outputing markdown files; feel free to ignore these
knitr::opts_chunk$set(eval = FALSE, 
                      include = TRUE, 
                      warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      dpi = 300,
                      fig.dim = c(9, 9),
                      out.width = '98%',
                      out.height = '98%')
 
#'
#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(cowplot)
library(fastshap) # for shapley values
library(shapviz)
library(caret)
library(viridis)
library(janitor)
library(here)

# Load required data
#source(here("setup.R"))
source(here("Assembly-analysis", "R", "00_all_MAGs_assembly","03_prepare_pairwise_data.R")) # note this also sources setup.R

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
inputs.fp <- here("Assembly-analysis", "outputs")
outputs_mags.fp <- here("Assembly-analysis", "outputs", "mag_rf_model")
outputs_func.fp <- here("Assembly-analysis", "outputs", "funct_rf_model")
figures.fp <- here("Assembly-analysis", "figures")

# If makes sense, why things might happen, how compares to microbial data

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

assembly_levels <- c("Homogenous.selection", "Heterogenous.selection",
                     "Homogenizing.dispersal", "Dispersal.limitation.and.drift",
                     "Drift")

#### ====================================================================== ####

#### Read in data 
#### ====================================================================== ####
# Read in OTU table and env data
input <- input_ra
input$otu_table <- input$otu_table[-1]

# Read in differential data
betanull.lf.diff <- read_csv(file = paste0(inputs.fp, "/betanull.lf.diff.csv"))

# Read in RF outputs mags
all_hardcode_model <- readRDS(file = paste0(outputs_mags.fp, "/all_hardcode_model.RDS"))
palsa_hardcode_model <- readRDS(paste0(outputs_mags.fp, "/palsa_hardcode_model.RDS"))
bog_hardcode_model <- readRDS(paste0(outputs_mags.fp, "/bog_hardcode_model.RDS"))
fen_hardcode_model <- readRDS(paste0(outputs_mags.fp, "/fen_hardcode_model.RDS"))

# Read in RF outputs functions
all_hardcode_model_func <- readRDS(file = paste0(outputs_func.fp, "/all_hardcode_model.RDS"))
palsa_hardcode_model_func <- readRDS(paste0(outputs_func.fp, "/palsa_hardcode_model.RDS"))
bog_hardcode_model_func <- readRDS(paste0(outputs_func.fp, "/bog_hardcode_model.RDS"))
fen_hardcode_model_func <- readRDS(paste0(outputs_func.fp, "/fen_hardcode_model.RDS"))
#### ====================================================================== ####

# Setup plotting parameters
#### ====================================================================== ####

colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
shape_habitat  <- c(15, 16, 17)
fill_habitat   <- colour_habitat
habitat_comp_fill <- fill_habitat
names(habitat_comp_fill) <- c("Palsa_Palsa", "Bog_Bog", "Fen_Fen")
depth_levels <- c("0-9", "10-19", "20-29", "30-39")
depth_labels <- c("0", "10", "20", "30")
fill_depth   <- RColorBrewer::brewer.pal(5, "YlOrBr")[-1]
colour_depth <- fill_depth
year_levels <- c("2011", "2012", "2013", "2014", "2015", "2016", "2017")
year_labels <- year_levels
fill_year   <- RColorBrewer::brewer.pal(7, "OrRd")
colour_year <- fill_year
phylum_levels <- c(
  # Bacteria
  "p__Acidobacteriota",
  "p__Actinobacteriota",
  "p__Bacteroidota",
  "p__Chloroflexota",
  "p__Desulfobacterota",
  "p__Eremiobacterota",
  "p__Proteobacteria",
  "p__Verrucomicrobiota",
  # Archaea
  "p__Halobacteriota",
  "p__Methanobacteriota",
  "p__Thermoplasmatota",
  "p__Thermoproteota",
  "other",
  "Environmental Predictor"
)
phylum_labels <- c(
  # Bacteria
  "Acidobacteriota",
  "Actinobacteriota",
  "Bacteroidota",
  "Chloroflexota",
  "Desulfobacterota",
  "Eremiobacterota",
  "Proteobacteria",
  "Verrucomicrobiota",
  # Archaea
  "Halobacteriota",
  "Methanobacteriota",
  "Thermoplasmatota",
  "Thermoproteota",
  "other",
  "Environmental Predictor"
)
colour_phylum <- c(RColorBrewer::brewer.pal(12, "Set3"), "white", "grey30")
fill_phylum <- colour_phylum
metapathway_levels <- metapathway_groups %>% arrange(metapath) %>% `$`(metapath) %>% unique()
colour_metapathway <- RColorBrewer::brewer.pal(10, "Set3")
fill_metapathway <- colour_metapathway

special_levels <- c("generalist", "methanogen", "fermenter", "macromolecule_degrader", "methanotroph", "homoacetogen", "monomer_degrader", NA)
special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-8]
special_fill <- special_colour
n_special_levels <- c("nitric_oxide_oxidiser", "nitrogen_fixer", "nitrate_reducer", "denitrifier", "nitrite_oxidiser", NA)
n_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-6:-8]
n_special_fill <- n_special_colour
s_special_levels <- c("sulfur_oxidiser", "sulfate_reducer", NA)
s_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-3:-8]
s_special_fill <- s_special_colour

all_special_levels <- c(na.exclude(special_levels), na.exclude(n_special_levels), na.exclude(s_special_levels), NA)
all_special_labels <- gsub("_", " ", all_special_levels)
all_special_colour <- c(RColorBrewer::brewer.pal(7, "OrRd")[c(7:1)],RColorBrewer::brewer.pal(5, "Greens")[c(5,2,4,1,3)],
                        RColorBrewer::brewer.pal(9, "Purples")[c(5,8)])
all_special_fill <- all_special_colour

#### ====================================================================== ####

# Useful functions
#### ====================================================================== ####

prepare_plot_vip_importance <- function(hardcode_model = all_hardcode_model, 
                                        plotN = NULL,
                                        RFStats = T,
                                        RFStatsVerb = T) {
  # hardcode_model = output from running rf model
  # define phylum plotting characteristics
  phylum_levels <- c(
    # Bacteria
    "p__Acidobacteriota",
    "p__Actinobacteriota",
    "p__Bacteroidota",
    "p__Chloroflexota",
    "p__Desulfobacterota",
    "p__Eremiobacterota",
    "p__Proteobacteria",
    "p__Verrucomicrobiota",
    # Archaea
    "p__Halobacteriota",
    "p__Methanobacteriota",
    "p__Thermoplasmatota",
    "p__Thermoproteota",
    "other",
    "Environmental Predictor"
  )
  phylum_labels <- c(
    # Bacteria
    "Acidobacteriota",
    "Actinobacteriota",
    "Bacteroidota",
    "Chloroflexota",
    "Desulfobacterota",
    "Eremiobacterota",
    "Proteobacteria",
    "Verrucomicrobiota",
    # Archaea
    "Halobacteriota",
    "Methanobacteriota",
    "Thermoplasmatota",
    "Thermoproteota",
    "other",
    "Environmental Predictor"
  )
  colour_phylum <- c(RColorBrewer::brewer.pal(12, "Set3"), "white", "grey30")
  fill_phylum <- colour_phylum
  
  special_levels <- c("generalist", "methanogen", "fermenter", "macromolecule_degrader", "methanotroph", "homoacetogen", "monomer_degrader", NA)
  special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-8]
  special_fill <- special_colour
  n_special_levels <- c("nitric_oxide_oxidiser", "nitrogen_fixer", "nitrate_reducer", "denitrifier", "nitrite_oxidiser", NA)
  n_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-6:-8]
  n_special_fill <- n_special_colour
  s_special_levels <- c("sulfur_oxidiser", "sulfate_reducer", NA)
  s_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-3:-8]
  s_special_fill <- s_special_colour
  
  all_special_levels <- c(na.exclude(special_levels), na.exclude(n_special_levels), na.exclude(s_special_levels), NA)
  all_special_labels <- gsub("_", " ", all_special_levels)
  all_special_colour <- c(RColorBrewer::brewer.pal(7, "OrRd")[c(7:1)],RColorBrewer::brewer.pal(5, "Greens")[c(5,2,4,1,3)],
                          RColorBrewer::brewer.pal(9, "Purples")[c(5,8)])
  all_special_fill <- all_special_colour
  
  # Get RF statistics
  ntrees <- hardcode_model$rf_permutation$num.trees
  node.size <-  hardcode_model$rf_permutation$min.node.size
  mtry <- hardcode_model$rf_permutation$mtry
  splitrule <- hardcode_model$rf_permutation$splitrule
  accur <- round(hardcode_model$conf_mat_permute$overall["Accuracy"], digits = 2)
  accurlwr <- round(hardcode_model$conf_mat_permute$overall["AccuracyLower"], digits = 2)
  accurupr <- round(hardcode_model$conf_mat_permute$overall["AccuracyUpper"], digits = 2)
  accurpval <- round(hardcode_model$conf_mat_permute$overall["AccuracyPValue"], digits = 2)
  RF_stats_det <- paste0("Ntrees = ", ntrees, "\n",
                         "node.size = ", node.size, "\n",
                         "mtry = ", mtry, "\n",
                         "splitrule = ", splitrule, "\n") 
  
  RF_stats_few <- paste0("Accuracy = ", accur, ",", " CI: (", accurlwr, ",", accurupr, ")\n",
                         "Accuracy Pval = ", accurpval)
  
  
  RF_stats_verb <- paste0(RF_stats_det, "\n", RF_stats_few)
  
  # Specializations
  specialisations_lf <- specialisations %>% 
    rename(carbon_specialisation = specialisation) %>%
    pivot_longer(contains("specialisation"), names_to = "specialisation_type",
                 values_to = "specialisation")
  
  # Get plot df
  plot_df <- hardcode_model$vip_permutation %>% 
    filter(pvalue < 0.05) %>%
    arrange(desc(Importance)) %>%
    left_join(names_translation_data.rf, by = c("Variable" = "new")) %>%
    mutate(Variable = original) %>%
    mutate(Variable = gsub("_mag_diff_abund", "", Variable)) %>%
    left_join(input$taxonomy, by = c("Variable" = "genome")) %>%
    mutate(PhylumGlom = ifelse(!is.na(Phylum) & !(Phylum %in% phylum_levels), "other", Phylum)) %>% 
    mutate(PhylumLabel = ifelse(!is.na(PhylumGlom), PhylumGlom, "Environmental Predictor")) %>%
    left_join(specialisations_lf, by = c("Variable" = "genome")) %>%
    mutate(Variable = fct_reorder(Variable, Importance, .desc = T))
  
  if(!is.null(plotN)) {
    plot_df <- plot_df %>%
      slice_max(order_by = Importance, n = plotN) 
  }
  
  
  bar_layers <- list(
    geom_bar(stat = "identity", aes(fill = PhylumLabel), color = "grey20"),
    scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels),
    scale_y_continuous(expand = expansion(mult = c(0, 0))),
    ylab("Importance"),
    theme_cowplot(),
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white"))
  )
  
  specialisation_bars <- list(
    geom_bar(stat = "identity", aes(fill = specialisation), color = "grey20"),
    scale_fill_manual("Specialisation", values = all_special_colour, breaks = all_special_levels, labels = all_special_labels),
    scale_y_continuous(expand = expansion(mult = c(0, 0))),
    facet_wrap(~specialisation_type),
    ylab("Importance"),
    theme_cowplot(),
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white"))
  )
  
  
  if(RFStats && RFStatsVerb) {
    bar_layers <- append(bar_layers, list(annotate("text", label = RF_stats_verb, x = Inf, y = Inf, 
                                                   hjust = 1, vjust = 1)))
  } else {
    if(RFStats) {
      bar_layers <- append(bar_layers, list(annotate("text", label = RF_stats_few, x = Inf, y = Inf, 
                                                     hjust = 1, vjust = 1)))
    }}
  
  return_list <- list(plot_df = plot_df, 
                      bar_layers = bar_layers,
                      specialisation_bars = specialisation_bars)
  
  return(return_list)
}

prepare_plot_funct_importance <- function(hardcode_model = all_hardcode_model_func, 
                                        plotN = NULL,
                                        RFStats = T,
                                        RFStatsVerb = T) {
  # hardcode_model = output from running rf model
  # define pathway plotting characteristics
  metapathway_levels <- c(metapathway_groups %>% arrange(metapath) %>% `$`(metapath) %>% unique(),
                          "Environmental Predictor")
  colour_metapathway <- c(RColorBrewer::brewer.pal(10, "Set3"),  "grey30")
  fill_metapathway <- colour_metapathway
  
  # Get RF statistics
  ntrees <- hardcode_model$rf_permutation$num.trees
  node.size <-  hardcode_model$rf_permutation$min.node.size
  mtry <- hardcode_model$rf_permutation$mtry
  splitrule <- hardcode_model$rf_permutation$splitrule
  accur <- round(hardcode_model$conf_mat_permute$overall["Accuracy"], digits = 2)
  accurlwr <- round(hardcode_model$conf_mat_permute$overall["AccuracyLower"], digits = 2)
  accurupr <- round(hardcode_model$conf_mat_permute$overall["AccuracyUpper"], digits = 2)
  accurpval <- round(hardcode_model$conf_mat_permute$overall["AccuracyPValue"], digits = 2)
  RF_stats_det <- paste0("Ntrees = ", ntrees, "\n",
                         "node.size = ", node.size, "\n",
                         "mtry = ", mtry, "\n",
                         "splitrule = ", splitrule, "\n") 
  
  RF_stats_few <- paste0("Accuracy = ", accur, ",", " CI: (", accurlwr, ",", accurupr, ")\n",
                         "Accuracy Pval = ", accurpval)
  
  
  RF_stats_verb <- paste0(RF_stats_det, "\n", RF_stats_few)
  
  # Get plot df
  plot_df <- hardcode_model$vip_permutation %>% 
    filter(pvalue < 0.05) %>%
    arrange(desc(Importance)) %>%
    left_join(names_translation_data.rf, by = c("Variable" = "new")) %>%
    mutate(Variable = original) %>%
    mutate(Variable = gsub("_path_diff_abund", "", Variable)) %>%
    left_join(metapathway_groups, by = c("Variable" = "subpathway")) %>%
    mutate(PathLabel = ifelse(!is.na(metapath), metapath, "Environmental Predictor")) %>%
    mutate(Variable = fct_reorder(Variable, Importance, .desc = T))
  
  if(!is.null(plotN)) {
    plot_df <- plot_df %>%
      slice_max(order_by = Importance, n = plotN) 
  }
  
  
  bar_layers <- list(
    geom_bar(stat = "identity", aes(fill = PathLabel), color = "grey20"),
    scale_fill_manual("Pathway", values = colour_metapathway, breaks = metapathway_levels, labels = metapathway_levels),
    scale_y_continuous(expand = expansion(mult = c(0, 0))),
    ylab("Importance"),
    theme_cowplot(),
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white"),
          plot.margin = margin(10, 10, 10, 100))
  )
  
  if(RFStats && RFStatsVerb) {
    bar_layers <- append(bar_layers, list(annotate("text", label = RF_stats_verb, x = Inf, y = Inf, 
                                                   hjust = 1, vjust = 1)))
  } else {
    if(RFStats) {
      bar_layers <- append(bar_layers, list(annotate("text", label = RF_stats_few, x = Inf, y = Inf, 
                                                     hjust = 1, vjust = 1)))
    }}
  
  return_list <- list(plot_df = plot_df, 
                      bar_layers = bar_layers)
  
  return(return_list)
}

prepare_shap_object <- function(shap_val_rf = palsa_hardcode_model_func$shap_vals_perm,
                                shapdat = palsa_hardcode_model_func$train_dat) {
  return_list <- list()
  
  HomSel <- shap_val_rf %>%
    as.data.frame() %>%
    select(matches("_Homogenous.selection")) %>%
    rename_all(~gsub("_Homogenous.selection", "", .))
  
  if(ncol(HomSel) > 0) {
    class(HomSel) <- c("tbl_df", "tbl", "data.frame", "explain")
    HomSelshap <- shapviz(object = HomSel, X = shapdat[-1])
    return_list <- append(return_list, list(HomSel = HomSelshap))
  }
  
  HetSel <- shap_val_rf %>%
    as.data.frame() %>%
    select(matches("_Heterogenous.selection")) %>%
    rename_all(~gsub("_Heterogenous.selection", "", .))
  
  
  if(ncol(HetSel) > 0) {
    class(HetSel) <- c("tbl_df", "tbl", "data.frame", "explain")
    HetSelshap <- shapviz(object = HetSel, X = shapdat[-1])
    return_list <- append(return_list, list(HetSel = HetSelshap))
  }
  
  
  HomDisp <- shap_val_rf %>%
    as.data.frame() %>%
    select(matches("_Homogenizing.dispersal")) %>%
    rename_all(~gsub("_Homogenizing.dispersal", "", .))
  
  if(ncol(HomDisp) > 0) {
    class(HomDisp) <- c("tbl_df", "tbl", "data.frame", "explain")
    HomDispshap <- shapviz(object = HomDisp, X = shapdat[-1])
    return_list <- append(return_list, list(HomDisp = HomDispshap))
  }
  
  DispLim <- shap_val_rf %>%
    as.data.frame() %>%
    select(matches("_Dispersal.limitation.and.drift")) %>%
    rename_all(~gsub("_Dispersal.limitation.and.drift", "", .))
  
  
  if(ncol(DispLim) > 0) {
    class(DispLim) <- c("tbl_df", "tbl", "data.frame", "explain")
    DispLimshap <- shapviz(object = DispLim, X = shapdat[-1])
    return_list <- append(return_list, list(DispLim = DispLimshap))
  }
  
  
  Drift <- shap_val_rf %>%
    as.data.frame() %>%
    select(matches("_Drift")) %>%
    rename_all(~gsub("_Drift", "", .))
  
  if(ncol(Drift) > 0) {
    class(Drift) <- c("tbl_df", "tbl", "data.frame", "explain")
    Driftshap <- shapviz(object = Drift, X = shapdat[-1])
    return_list <- append(return_list, list(Drift = Driftshap))
  }
  
  
  return(return_list)
}
#### ====================================================================== ####

# Run random forest models with Biotic Side (aka, taxa)
#### ====================================================================== ####
# Remove na values and assign training and testing data (with proportional groups)
data.rf <- betanull.lf.diff %>% 
  select(Assembly_Process, Habitat_comp, Depth_diff, T_air_7d_diff, 
         pct_time_below_WTD_21d_diff, Year_diff, 
         contains("mag_diff_abund"), 
         Site1, Site2) %>% 
  select(!contains("scl_")) %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>% 
  mutate(Assembly_Process = gsub(" ", ".", Assembly_Process)) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = assembly_levels))  %>% 
  filter(Habitat_comp %in% c("Fen_Fen", "Bog_Bog", "Palsa_Palsa")) %>%
  mutate(Habitat_comp = droplevels(Habitat_comp)) %>%
  filter(Depth_diff < 10) # Only use depth differences in the same depthLumping 
names_translation_data.rf <- data.frame(original = names(data.rf)) 
data.rf <- data.rf %>%
  janitor::clean_names(case = "none")

names_translation_data.rf$new <- names(data.rf)



  
all_imp_plot <- prepare_plot_vip_importance(hardcode_model = all_hardcode_model, plotN = 50,
                                            RFStats = T, RFStatsVerb = F)

all_pred_plot <- ggplot(data = all_imp_plot$plot_df %>% select(!contains("specialisation")) %>%
                          distinct(), aes(x = Variable, y = Importance)) +
  all_imp_plot$bar_layers

all_pred_plot  
ggsave(all_pred_plot, filename = paste0(figures.fp, "/all_mag_predictors_classification_rf.png"),
       width = 14, height = 10, dpi = 300)


all_pred_plot_spec <- ggplot(data = all_imp_plot$plot_df, aes(x = Variable, y = Importance)) +
  all_imp_plot$specialisation_bars
ggsave(all_pred_plot_spec, filename = paste0(figures.fp, "/all_mag_predictors_classification_rf_spec.png"),
       width = 18, height = 10, dpi = 300)
#### ====================================================================== ####

#### ====================================================================== ####
# palsa
palsa_data <- data.rf %>% filter(Habitat_comp == "Palsa_Palsa") %>%
  select(-Habitat_comp) %>%
  mutate(Assembly_Process = fct_drop(Assembly_Process))


palsa_hardcode_model$conf_mat_permute

palsa_imp_plot <- prepare_plot_vip_importance(hardcode_model = palsa_hardcode_model, plotN = 50,
                                            RFStats = T, RFStatsVerb = F)

palsa_pred_plot <- ggplot(data = palsa_imp_plot$plot_df %>% select(!contains("specialisation")) %>%
                            distinct(), aes(x = Variable, y = Importance)) +
  palsa_imp_plot$bar_layers

palsa_pred_plot

ggsave(palsa_pred_plot, filename = paste0(figures.fp, "/palsa_mag_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)

palsa_imp_plot$plot_df %>%
  left_join(network_stats %>% filter(Habitat == "Palsa"), 
            by = c("Variable" = "genome_id")) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  ggplot(aes(y = eigenvector_centrality, x = fct_reorder(Variable, Importance, .desc = T))) +
  geom_bar(aes(fill = PhylumLabel, group = Year), stat = "identity", position = "dodge") +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("Eigenvector Centrality") + xlab("Taxa Ordered by Assembly Importance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




palsa_pred_plot_spec <- ggplot(data = palsa_imp_plot$plot_df %>% filter(!is.na(specialisation_type)), 
                               aes(x = Variable, y = Importance)) +
  palsa_imp_plot$specialisation_bars
ggsave(palsa_pred_plot_spec, filename = paste0(figures.fp, "/palsa_mag_predictors_classification_rf_spec.png"),
       width = 18, height = 10, dpi = 300)


bog_imp_plot <- prepare_plot_vip_importance(hardcode_model = bog_hardcode_model, plotN = 50,
                                              RFStats = T, RFStatsVerb = F)

bog_pred_plot <- ggplot(data = bog_imp_plot$plot_df%>% select(!contains("specialisation")) %>%
                          distinct(), aes(x = Variable, y = Importance)) +
  bog_imp_plot$bar_layers

ggsave(bog_pred_plot, filename = paste0(figures.fp, "/bog_mag_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)

bog_pred_plot_spec <- ggplot(data = bog_imp_plot$plot_df %>% filter(!is.na(specialisation_type)),
                             aes(x = Variable, y = Importance)) +
  bog_imp_plot$specialisation_bars
ggsave(bog_pred_plot_spec, filename = paste0(figures.fp, "/bog_mag_predictors_classification_rf_spec.png"),
       width = 18, height = 10, dpi = 300)

bog_imp_plot$plot_df %>%
  left_join(network_stats %>% filter(Habitat == "Bog"), 
            by = c("Variable" = "genome_id")) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  ggplot(aes(y = eigenvector_centrality, x = fct_reorder(Variable, Importance, .desc = T))) +
  geom_bar(aes(fill = PhylumLabel, group = Year), stat = "identity", position = "dodge") +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("Eigenvector Centrality") + xlab("Taxa Ordered by Assembly Importance") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


bog_imp_plot$plot_df %>%
  left_join(network_stats %>% filter(Habitat == "Bog"), 
            by = c("Variable" = "genome_id")) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  ggplot(aes(x = eigenvector_centrality, y = Importance)) +
  geom_point(aes(fill = PhylumLabel), shape = 21) +
  # geom_text(aes(label = Year), hjust = 1) +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("Assembly Importance") + xlab("Network Centrality (eigenvector centrality)") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

fen_imp_plot <- prepare_plot_vip_importance(hardcode_model = fen_hardcode_model, plotN = 50,
                                              RFStats = T, RFStatsVerb = F)

fen_pred_plot <- ggplot(data = fen_imp_plot$plot_df %>% select(!contains("specialisation")) %>%
                          distinct(), aes(x = Variable, y = Importance)) +
  fen_imp_plot$bar_layers

ggsave(fen_pred_plot, filename = paste0(figures.fp, "/fen_mag_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)


fen_imp_plot$plot_df %>%
  left_join(network_stats %>% filter(Habitat == "Fen"), 
            by = c("Variable" = "genome_id")) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  ggplot(aes(x = eigenvector_centrality, y = Importance)) +
  geom_point(aes(fill = PhylumLabel), shape = 21) +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("Assembly Importance") + xlab("Network Centrality (eigenvector centrality)") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


fen_imp_plot$plot_df %>%
  left_join(network_stats %>% filter(Habitat == "Fen"), 
            by = c("Variable" = "genome_id")) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  ggplot(aes(y = eigenvector_centrality, x = fct_reorder(Variable, Importance, .desc = T))) +
  geom_bar(aes(fill = PhylumLabel, group = Year), stat = "identity", position = "dodge") +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("Eigenvector Centrality") + xlab("Taxa Ordered by Assembly Importance") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


fen_pred_plot_spec <- ggplot(data = fen_imp_plot$plot_df %>% filter(!is.na(specialisation_type)), 
                             aes(x = Variable, y = Importance)) +
  fen_imp_plot$specialisation_bars
ggsave(fen_pred_plot_spec, filename = paste0(figures.fp, "/fen_mag_predictors_classification_rf_spec.png"),
       width = 18, height = 10, dpi = 300)



# All habitats Relationship between centrality and Important Assembly Predictors
bind_rows(
  palsa_imp_plot$plot_df %>%
    mutate(Habitat = "Palsa"),
  bog_imp_plot$plot_df %>%
  mutate(Habitat = "Bog"),
  fen_imp_plot$plot_df %>%
    mutate(Habitat = "Fen")) %>%
  left_join(network_stats %>%
            group_by(Habitat, genome_id) %>%
              summarize(across(contains("centrality"), list(mean = ~mean(., na.rm = T), 
                                                            sd = ~sd(., na.rm = T),
                                                            se = ~sd(.,na.rm = T)/sqrt(n())))), 
            by = c("Variable" = "genome_id", 
                   "Habitat" = "Habitat")) %>%
  mutate(lower = eigenvector_centrality_mean - eigenvector_centrality_se,
         upper = eigenvector_centrality_mean + eigenvector_centrality_se,
         Habitat = factor(Habitat, levels = habitat_levels)) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  filter(pvalue < 0.05) %>%
  ggplot(aes(x = eigenvector_centrality_mean, y = Importance)) +
  geom_linerange(aes(xmin = lower, 
                      xmax = upper)) +
  geom_point(aes(fill = PhylumLabel), shape = 21) +
#  geom_text(aes(label = Year), hjust = 1) +
  facet_wrap(~Habitat, scales = "free") +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("MAG Assembly Predictor Importance") + xlab("Network Centrality (eigenvector centrality over 7 years)") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  
all_imp_plot$plot_df %>%
  left_join(network_stats %>%
              pivot_longer(contains("centrality"), names_to = "centrality_metric", values_to = "centrality") %>%
              group_by(Habitat, centrality_metric, genome_id) %>%
              summarize(across(all_of("centrality"), list(mean = ~mean(., na.rm = T), 
                                                            sd = ~sd(., na.rm = T),
                                                            se = ~sd(.,na.rm = T)/sqrt(n())))),
            by = c("Variable" = "genome_id")) %>%
  mutate(lower = centrality_mean - centrality_se,
         upper = centrality_mean + centrality_se,
         Habitat = factor(Habitat, levels = habitat_levels)) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  filter(pvalue < 0.05) %>% filter(!is.na(centrality_mean)) %>%
  ggplot(aes(x = centrality_mean, y = Importance)) +
  geom_linerange(aes(xmin = lower, 
                     xmax = upper)) +
  geom_point(aes(fill = PhylumLabel), shape = 21) +
  #  geom_text(aes(label = Year), hjust = 1) +
  facet_wrap(centrality_metric~Habitat, scales = "free") +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("MAG Assembly Predictor Importance") + xlab("Network Centrality (eigenvector centrality over 7 years)") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


bind_rows(
  palsa_imp_plot$plot_df %>%
    mutate(Habitat = "Palsa"),
  bog_imp_plot$plot_df %>%
    mutate(Habitat = "Bog"),
  fen_imp_plot$plot_df %>%
    mutate(Habitat = "Fen")) %>%
  left_join(network_stats %>%
              pivot_longer(contains("centrality"), names_to = "centrality_metric", values_to = "centrality") %>%
              group_by(Habitat, centrality_metric, genome_id) %>%
              summarize(across(all_of("centrality"), list(mean = ~mean(., na.rm = T), 
                                                          sd = ~sd(., na.rm = T),
                                                          se = ~sd(.,na.rm = T)/sqrt(n())))),
            by = c("Variable" = "genome_id",
                   "Habitat" = "Habitat")) %>%
  mutate(lower = centrality_mean - centrality_se,
         upper = centrality_mean + centrality_se,
         Habitat = factor(Habitat, levels = habitat_levels)) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  filter(pvalue < 0.05) %>% filter(!is.na(centrality_mean)) %>%
  ggplot(aes(x = centrality_mean, y = Importance)) +
  geom_linerange(aes(xmin = lower, 
                     xmax = upper)) +
  geom_point(aes(fill = PhylumLabel), shape = 21) +
  facet_wrap(centrality_metric~Habitat, scales = "free") +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("MAG Assembly Predictor Importance") + xlab("Network Centrality (eigenvector centrality over 7 years)") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


bog_imp_plot$plot_df %>%
  left_join(network_stats %>% filter(Habitat == "Bog"), 
            by = c("Variable" = "genome_id")) %>%
  select(!contains("specialisation")) %>% distinct() %>%
  ggplot(aes(x = eigenvector_centrality, y = Importance)) +
  geom_point(aes(fill = PhylumLabel), shape = 21) +
  geom_text(aes(label = Year), hjust = 1) +
  scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
  ylab("Assembly Importance") + xlab("Network Centrality (eigenvector centrality)") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

#### ====================================================================== ####
# Shapley
#### ====================================================================== ####

# Importance plot
palsa_mag_shap <- prepare_shap_object(shap_val_rf = palsa_hardcode_model$shap_vals_perm, 
                                        shapdat = palsa_hardcode_model$train_dat)
sv_importance(palsa_mag_shap$HomSel, kind = "both") +
  ggtitle("Explaining Homogenous Selection") 

sv_importance(palsa_mag_shap$HomDisp, kind = "both") +
  ggtitle("Explaining Homogenizing Dispersal") 

sv_importance(palsa_mag_shap$DispLim, kind = "both") +
  ggtitle("Explaining Dispersal Limitation") 

sv_importance(palsa_mag_shap$Drift, kind = "both") +
  ggtitle("Explaining Drift") 


sv_dependence(palsa_mag_shap$HomSel, v = "X20120700_P2M_1_mag_diff_abund") + ggtitle("Homogenizing Selection")
#sv_dependence(bog_funct_shap$HetSel, v = "X20120700_P2M_1_mag_diff_abund") + ggtitle("Heterogeneous Selection")
sv_dependence(palsa_mag_shap$HomDisp, v = "X20120700_P2M_1_mag_diff_abund") + ggtitle("Homogenizing Dispersal")
sv_dependence(palsa_mag_shap$DispLim, v = "X20120700_P2M_1_mag_diff_abund") + ggtitle("Dispersal Limitation")
sv_dependence(palsa_mag_shap$Drift, v = "X20120700_P2M_1_mag_diff_abund") + ggtitle("Drift")

product_refined %>%
  filter(genome %in%c("20120700_P2M_1")) %>%
  filter(call == TRUE)


sv_dependence(palsa_mag_shap$HomSel, v = "X20120700_P2M_5_mag_diff_abund") + ggtitle("Homogenizing Selection")
#sv_dependence(bog_funct_shap$HetSel, v = "X20120700_P2M_1_mag_diff_abund") + ggtitle("Heterogeneous Selection")
sv_dependence(palsa_mag_shap$HomDisp, v = "X20120700_P2M_5_mag_diff_abund") + ggtitle("Homogenizing Dispersal")
sv_dependence(palsa_mag_shap$DispLim, v = "X20120700_P2M_5_mag_diff_abund") + ggtitle("Dispersal Limitation")
sv_dependence(palsa_mag_shap$Drift, v = "X20120700_P2M_5_mag_diff_abund") + ggtitle("Drift")



# Methanogens and Friends

bog_mag_shap <- prepare_shap_object(shap_val_rf = bog_hardcode_model$shap_vals_perm, 
                                      shapdat = bog_hardcode_model$train_dat)

sv_dependence(bog_mag_shap$HomSel, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Homogenizing Selection")
sv_dependence(bog_mag_shap$HetSel, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Heterogeneous Selection")
sv_dependence(bog_mag_shap$HomDisp, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(bog_mag_shap$DispLim, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(bog_mag_shap$Drift, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Drift")

sv_dependence(bog_mag_shap$HomSel, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Homogenizing Selection")
sv_dependence(bog_mag_shap$HetSel, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Heterogeneous Selection")
sv_dependence(bog_mag_shap$HomDisp, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(bog_mag_shap$DispLim, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(bog_mag_shap$Drift, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Drift")


sv_dependence(bog_mag_shap$HomSel, v = "X20170700_S25_9_mag_diff_abund", "auto") + ggtitle("Homogenizing Selection")
# Het sel. more ikely as differential abundance of S1X_2 increases
sv_dependence(bog_mag_shap$HetSel, v = "X20170700_S25_9_mag_diff_abund", "X20170700_S25_26_mag_diff_abund") + ggtitle("Heterogeneous Selection")
sv_dependence(bog_mag_shap$HomDisp, v = "X20170700_S25_9_mag_diff_abund", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(bog_mag_shap$DispLim, v = "X20170700_S25_9_mag_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(bog_mag_shap$Drift, v = "X20170700_S25_9_mag_diff_abund", "auto") + ggtitle("Drift")



fen_mag_shap <- prepare_shap_object(shap_val_rf = fen_hardcode_model$shap_vals_perm, 
                                    shapdat = fen_hardcode_model$train_dat)

# A.stor and M.Stor friend
sv_dependence(fen_mag_shap$HomSel, v = "X20170700_S25_9_mag_diff_abund", "PLGY01_mag_diff_abund") + ggtitle("Homogenizing Selection")
# Het sel. more likely as differential abundance of S1X_2 increases
sv_dependence(fen_mag_shap$HetSel, v = "X20170700_S25_9_mag_diff_abund", "") + ggtitle("Heterogeneous Selection")
sv_dependence(fen_mag_shap$HomDisp, v = "X20170700_S25_9_mag_diff_abund", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(fen_mag_shap$DispLim, v = "X20170700_S25_9_mag_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(fen_mag_shap$Drift, v = "X20170700_S25_9_mag_diff_abund", "auto") + ggtitle("Drift")


sv_dependence(fen_mag_shap$HomSel, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Homogenizing Selection")
# Het sel. more likely as differential abundance of S1X_2 increases
sv_dependence(fen_mag_shap$HetSel, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Heterogeneous Selection")
sv_dependence(fen_mag_shap$HomDisp, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(fen_mag_shap$DispLim, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(fen_mag_shap$Drift, v = "X20120700_S1X_2_mag_diff_abund", "auto") + ggtitle("Drift")


sv_dependence(fen_mag_shap$HomSel, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Homogenizing Selection")
sv_dependence(fen_mag_shap$HetSel, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Heterogeneous Selection")
sv_dependence(fen_mag_shap$HomDisp, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(fen_mag_shap$DispLim, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(fen_mag_shap$Drift, v = "X20110800_S1D_5_mag_diff_abund", "auto") + ggtitle("Drift")



sv_importance(fen_mag_shap$HomDisp, kind = "both")
sv_importance(fen_mag_shap$DispLim, kind = "both")
sv_importance(fen_mag_shap$Drift, kind = "both")
sv_importance(fen_mag_shap$HomSel, kind = "both")
sv_importance(fen_mag_shap$HetSel, kind = "both")





product_refined %>%
  filter(genome %in%c("20120700_S1X_15")) %>%
  filter(call == TRUE) %>%
  left_join(input$taxonomy, by = "genome") %>% View()

sv_dependence(fen_mag_shap$HomSel, v = "X20120700_S1X_15_mag_diff_abund", "auto") + ggtitle("Homogenizing Selection")
sv_dependence(fen_mag_shap$HetSel, v = "X20120700_S1X_15_mag_diff_abund", "auto") + ggtitle("Heterogeneous Selection")
sv_dependence(fen_mag_shap$HomDisp, v = "X20120700_S1X_15_mag_diff_abund", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(fen_mag_shap$DispLim, v = "X20120700_S1X_15_mag_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(fen_mag_shap$Drift, v = "X20120700_S1X_15_mag_diff_abund", "auto") + ggtitle("Drift")

#### ====================================================================== ####

# Functional data
#### ====================================================================== ####

# Run random forest models with functional data
#### ====================================================================== ####
# Remove na values and assign training and testing data (with proportional groups)
data.rf <- betanull.lf.diff %>% 
  select(Assembly_Process, Habitat_comp, Depth_diff, T_air_7d_diff, 
         pct_time_below_WTD_21d_diff, Year_diff, 
         contains("path_diff_abund"), 
         Site1, Site2) %>% 
  select(!contains("scl_")) %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>% 
  mutate(Assembly_Process = gsub(" ", ".", Assembly_Process)) %>%
  mutate(Assembly_Process = factor(Assembly_Process, levels = assembly_levels))  %>% 
  filter(Habitat_comp %in% c("Fen_Fen", "Bog_Bog", "Palsa_Palsa")) %>%
  mutate(Habitat_comp = droplevels(Habitat_comp)) %>%
  filter(Depth_diff < 10) # Only use depth differences in the same depthLumping 
names_translation_data.rf <- data.frame(original = names(data.rf)) 
data.rf <- data.rf %>%
  janitor::clean_names(case = "none")

names_translation_data.rf$new <- names(data.rf)




all_imp_plot_func <- prepare_plot_funct_importance(hardcode_model = all_hardcode_model_func, plotN = 50,
                                            RFStats = T, RFStatsVerb = F)

all_pred_plot_func <- ggplot(data = all_imp_plot_func$plot_df, aes(x = Variable, y = Importance)) +
  all_imp_plot_func$bar_layers

all_pred_plot_func  
ggsave(all_pred_plot_func, filename = paste0(figures.fp, "/all_func_predictors_classification_rf.png"),
       width = 14, height = 10, dpi = 300)

#### ====================================================================== ####

#### ====================================================================== ####
# palsa
palsa_data <- data.rf %>% filter(Habitat_comp == "Palsa_Palsa") %>%
  select(-Habitat_comp) %>%
  mutate(Assembly_Process = fct_drop(Assembly_Process))


palsa_imp_plot_func <- prepare_plot_funct_importance(hardcode_model = palsa_hardcode_model_func, plotN = 50,
                                              RFStats = T, RFStatsVerb = F)

palsa_pred_plot_func <- ggplot(data = palsa_imp_plot_func$plot_df, aes(x = Variable, y = Importance)) +
  palsa_imp_plot_func$bar_layers

palsa_pred_plot_func
ggsave(palsa_pred_plot_func, filename = paste0(figures.fp, "/palsa_func_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)

bog_imp_plot_func <- prepare_plot_funct_importance(hardcode_model = bog_hardcode_model_func, plotN = 50,
                                                     RFStats = T, RFStatsVerb = F)

bog_pred_plot_func <- ggplot(data = bog_imp_plot_func$plot_df, aes(x = Variable, y = Importance)) +
  bog_imp_plot_func$bar_layers

bog_pred_plot_func
ggsave(bog_pred_plot_func, filename = paste0(figures.fp, "/bog_func_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)

fen_imp_plot_func <- prepare_plot_funct_importance(hardcode_model = fen_hardcode_model_func, plotN = 50,
                                                     RFStats = T, RFStatsVerb = F)

fen_pred_plot_func <- ggplot(data = fen_imp_plot_func$plot_df, aes(x = Variable, y = Importance)) +
  fen_imp_plot_func$bar_layers

fen_pred_plot_func
ggsave(fen_pred_plot_func, filename = paste0(figures.fp, "/fen_func_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)
#### ====================================================================== ####


# Exploring specific functional pathways
#### ====================================================================== ####
# Palsa

# Glycerol degradation


# Importance plot
palsa_funct_shap <- prepare_shap_object(shap_val_rf = palsa_hardcode_model_func$shap_vals_perm, 
                            shapdat = palsa_hardcode_model_func$train_dat)
sv_importance(palsa_funct_shap$HomSel, kind = "both") +
  ggtitle("Explaining Homogenous Selection") 

sv_importance(palsa_funct_shap$HomDisp, kind = "both") +
  ggtitle("Explaining Homogenizing Dispersal") 

sv_importance(palsa_funct_shap$DispLim, kind = "both") +
  ggtitle("Explaining Dispersal Limitation") 

sv_importance(palsa_funct_shap$Drift, kind = "both") +
  ggtitle("Explaining Drift") 

sv_dependence(test$HomSel, v = "ethanol_fermentation_all_path_diff_abund")

# Higher likelihood of dispersal limitation with glycerol degradation; may be some interaction with xylose pathway
a <- sv_dependence(palsa_funct_shap$HomSel, v = "glycerol_degradation_glycerol_path_diff_abund") + ggtitle("Homogenizing Selection")
sv_dependence(palsa_funct_shap$HomDisp, v = "glycerol_degradation_glycerol_path_diff_abund") + ggtitle("Homogenizing Dispersal")
b <- sv_dependence(palsa_funct_shap$DispLim, v = "glycerol_degradation_glycerol_path_diff_abund") + ggtitle("Dispersal Limitation")
sv_dependence(palsa_funct_shap$Drift, v = "glycerol_degradation_glycerol_path_diff_abund") + ggtitle("Drift")
test <- plot_grid(a, b)
save_plot(filename = "~/Downloads/test.png", plot = test, base_height = 7,base_width = 14, dpi = 300)


sv_importance(palsa_funct_shap$HomDisp, kind = "both") +
  ggtitle("Explaining Homogenizing Dispersal") 

# When depth difference is greater, Homogenizing dispersal is less likely
sv_dependence(palsa_funct_shap$HomDisp, v = "Depth_diff") +
  xlab("Depth Difference (cm)") + ylab("Importance of Depth difference in explaining homogenizing dispersal") +
  geom_smooth(color = "black")


sv_dependence(palsa_funct_shap$HomDisp, v = "T_air_7d_diff") +
  xlab("Temp Difference (DegC/7days)") + ylab("Importance of 7-day Air temperature difference in explaining homogenizing dispersal") +
  geom_smooth(color = "black")


sv_dependence(palsa_funct_shap$HomDisp, v = "CAZy_Xyloglucan_path_diff_abund")
sv_dependence(palsa_funct_shap$HomSel, v = "CAZy_Xyloglucan_path_diff_abund")

# Bog
bog_funct_shap <- prepare_shap_object(shap_val_rf = bog_hardcode_model_func$shap_vals_perm, 
                                      shapdat = bog_hardcode_model_func$train_dat)
sv_importance(bog_funct_shap$HomSel, kind = "both") +
  ggtitle("Explaining Homogenous Selection") 

sv_importance(bog_funct_shap$HetSel, kind = "both") +
  ggtitle("Explaining Heterogeneous Selection") 
sv_importance(bog_funct_shap$HomDisp, kind = "both") +
  ggtitle("Explaining Homogenizing Dispersal") 
sv_importance(bog_funct_shap$DispLim, kind = "both") +
  ggtitle("Explaining Dispersal Limitation") 
sv_importance(bog_funct_shap$Drift, kind = "both") +
  ggtitle("Explaining Drift") 

# Urea in Bog
a <- sv_dependence(bog_funct_shap$HomSel, v = "urea_degradation_all_path_diff_abund") + ggtitle("Homogenizing Selection") +
  theme(legend.position = "none")
sv_dependence(bog_funct_shap$HetSel, v = "urea_degradation_all_path_diff_abund", "auto") + ggtitle("Heterogeneous Selection")

b <- sv_dependence(bog_funct_shap$HomDisp, v = "urea_degradation_all_path_diff_abund") + ggtitle("Homogenizing Dispersal")
bleg <- get_legend(b)
b <- b + theme(legend.position = "none")
sv_dependence(bog_funct_shap$DispLim, v = "urea_degradation_all_path_diff_abund", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(bog_funct_shap$Drift, v = "urea_degradation_all_path_diff_abund", "auto") + ggtitle("Drift")


test <- plot_grid(a, b)
save_plot(filename = "~/Downloads/test.png", plot = test, base_height = 7,base_width = 14, dpi = 300)

# Depth/WTD in Bog
sv_dependence(bog_funct_shap$HomSel, v = "Depth_diff") + ggtitle("Homogenizing Selection")
sv_dependence(bog_funct_shap$HetSel, v = "Depth_diff") + ggtitle("Heterogeneous Selection")
a <- sv_dependence(bog_funct_shap$HomDisp, v = "Depth_diff") + ggtitle("Homogenizing Dispersal")
sv_dependence(bog_funct_shap$DispLim, v = "Depth_diff") + ggtitle("Dispersal Limitation")
sv_dependence(bog_funct_shap$Drift, v = "Depth_diff") + ggtitle("Drift")

ggsave("~/Downloads/test.png", plot = a, dpi = 300)


sv_dependence(bog_funct_shap$HomSel, v = "pct_time_below_WTD_21d_diff") + ggtitle("Homogenizing Selection")
sv_dependence(bog_funct_shap$HetSel, v = "pct_time_below_WTD_21d_diff") + ggtitle("Heterogeneous Selection")
sv_dependence(bog_funct_shap$HomDisp, v = "pct_time_below_WTD_21d_diff") + ggtitle("Homogenizing Dispersal")
sv_dependence(bog_funct_shap$DispLim, v = "pct_time_below_WTD_21d_diff") + ggtitle("Dispersal Limitation")
a <- sv_dependence(bog_funct_shap$Drift, v = "pct_time_below_WTD_21d_diff") + ggtitle("Drift")

ggsave("~/Downloads/test.png", plot = a, dpi = 300)




all_funct_shap <- prepare_shap_object(shap_val_rf = all_hardcode_model_func$shap_vals_perm, 
                                        shapdat = all_hardcode_model_func$train_dat)

# Habitat
sv_dependence(all_funct_shap$HomSel, v = "Habitat_comp", "auto") + ggtitle("Homogenous Selection")
sv_dependence(all_funct_shap$HetSel, v = "Habitat_comp", "auto") + ggtitle("Heterogeneous Selection")
sv_dependence(all_funct_shap$HomDisp, v = "Habitat_comp", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(all_funct_shap$DispLim, v = "Habitat_comp", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(all_funct_shap$Drift, v = "Habitat_comp", "auto") + ggtitle("Drift")


# Habitat
sv_dependence(all_funct_shap$HomSel, v = "Habitat_comp", "hydrogenotrophic_methanogenesis_all_path_diff_abund") + ggtitle("Homogenous Selection")
sv_dependence(all_funct_shap$HetSel, v = "Habitat_comp", "acetoclastic_methanogenesis_all_path_diff_abund") + ggtitle("Heterogeneous Selection")
sv_dependence(all_funct_shap$HomDisp, v = "Habitat_comp", "auto") + ggtitle("Homogenizing Dispersal")
sv_dependence(all_funct_shap$DispLim, v = "Habitat_comp", "auto") + ggtitle("Dispersal Limitation")
sv_dependence(all_funct_shap$Drift, v = "Habitat_comp", "auto") + ggtitle("Drift")

# Chitin
sv_dependence(all_funct_shap$HomSel, v = "CAZy_Chitin_path_diff_abund", "Habitat_comp") + ggtitle("Homogenizing Selection")

sv_dependence(all_funct_shap$Drift, v = "CAZy_Chitin_path_diff_abund", "auto") + ggtitle("Drift") +
  scale_color_manual("Habitat", values = habitat_comp_fill, labels = c("Palsa", "Bog", "Fen"))

sv_dependence(all_funct_shap$HomSel, v = "CAZy_Chitin_path_diff_abund", "auto") + ggtitle("Homogenizing Selection") +
  scale_color_manual("Habitat", values = habitat_comp_fill, labels = c("Palsa", "Bog", "Fen"))

sv_dependence(all_funct_shap$HetSel, v = "CAZy_Chitin_path_diff_abund", "Habitat_comp") + ggtitle("Heterogenous Selection")
sv_dependence(all_funct_shap$HomDisp, v = "CAZy_Chitin_path_diff_abund", "Habitat_comp") + ggtitle("Homogenizing Dispersal")
sv_dependence(all_funct_shap$DispLim, v = "CAZy_Chitin_path_diff_abund", "Habitat_comp") + ggtitle("Dispersal Limitation")
sv_dependence(all_funct_shap$Drift, v = "CAZy_Chitin_path_diff_abund", "Habitat_comp") + ggtitle("Drift")


# all Galacturonic Acid
sv_dependence(all_funct_shap$HomSel, v = "galacturonic_acid_degradation_all_path_diff_abund", "Habitat_comp") + ggtitle("Homogeneous Selection") +
  scale_color_manual("Habitat", values = habitat_comp_fill, labels = c("Palsa", "Bog", "Fen"))
sv_dependence(all_funct_shap$HetSel, v = "galacturonic_acid_degradation_all_path_diff_abund", "Habitat_comp") + ggtitle("Heterogenous Selection") +
  scale_color_manual("Habitat", values = habitat_comp_fill, labels = c("Palsa", "Bog", "Fen"))
sv_dependence(all_funct_shap$HomDisp, v = "galacturonic_acid_degradation_all_path_diff_abund", "Habitat_comp") + ggtitle("Homogenizing Dispersal") +
  scale_color_manual("Habitat", values = habitat_comp_fill, labels = c("Palsa", "Bog", "Fen"))
sv_dependence(all_funct_shap$DispLim, v = "galacturonic_acid_degradation_all_path_diff_abund", "Habitat_comp") + ggtitle("Dispersal Limitation") +
  scale_color_manual("Habitat", values = habitat_comp_fill, labels = c("Palsa", "Bog", "Fen"))
sv_dependence(all_funct_shap$Drift, v = "galacturonic_acid_degradation_all_path_diff_abund", "Habitat_comp") + ggtitle("Drift") +
  scale_color_manual("Habitat", values = habitat_comp_fill, labels = c("Palsa", "Bog", "Fen"))

# Fen Galacturonic Acid

fen_funct_shap <- prepare_shap_object(shap_val_rf = fen_hardcode_model_func$shap_vals_perm, 
                                        shapdat = fen_hardcode_model_func$train_dat)
sv_importance(fen_funct_shap$HomSel, kind = "both") +
  ggtitle("Explaining Homogenous Selection") 

sv_importance(fen_funct_shap$HetSel, kind = "both") +
  ggtitle("Explaining Heterogeneous Selection") 
sv_importance(fen_funct_shap$HomDisp, kind = "both") +
  ggtitle("Explaining Homogenizing Dispersal") 
sv_importance(fen_funct_shap$DispLim, kind = "both") +
  ggtitle("Explaining Dispersal Limitation") 
sv_importance(fen_funct_shap$Drift, kind = "both") +
  ggtitle("Explaining Drift") 






sv_dependence(fen_funct_shap$HomSel, v = "galacturonic_acid_degradation_all_path_diff_abund", "auto") + ggtitle("Homogeneous Selection")

sv_dependence(fen_funct_shap$HetSel, v = "galacturonic_acid_degradation_all_path_diff_abund", "auto") + ggtitle("Heterogenous Selection")
a <- sv_dependence(fen_funct_shap$HomDisp, v = "galacturonic_acid_degradation_all_path_diff_abund") + ggtitle("Homogenizing Dispersal") +
  geom_hline(yintercept = 0, linetype = "dashed")
a
sv_dependence(fen_funct_shap$DispLim, v = "galacturonic_acid_degradation_all_path_diff_abund", "auto") + ggtitle("Dispersal Limitation")
b <- sv_dependence(fen_funct_shap$Drift, v = "galacturonic_acid_degradation_all_path_diff_abund") + ggtitle("Drift") +
  geom_hline(yintercept = 0, linetype = "dashed")
b
test <- plot_grid(a, b)
save_plot(filename = "~/Downloads/test.png", plot = test, base_height = 5,base_width = 10, dpi = 300)



sv_dependence(fen_funct_shap$HomSel, v = "CAZy_Alpha_mannan_path_diff_abund", "auto") + ggtitle("Homogeneous Selection")
sv_dependence(fen_funct_shap$HetSel, v = "CAZy_Alpha_mannan_path_diff_abund", "auto") + ggtitle("Heterogenous Selection")
sv_dependence(fen_funct_shap$Drift, v = "CAZy_Alpha_mannan_path_diff_abund", "auto") + ggtitle("Drift")

#### ====================================================================== ####
# Relationship between importance and VIP score (or not, as it turns out)
#### ====================================================================== ####
vip_imp_plot <- function(imp_plot_df = all_imp_plot$plot_df,
                         Habitat = c("All")) {
  vip_imp <- imp_plot_df %>%
    left_join(vip_members,
              by = c("Variable" = "MAG")) %>%
    select(!contains("specialisation")) %>%
    distinct()
  
  ggplot(vip_imp, aes(x = `VIP Score`, y = Importance)) +
    geom_point(aes(fill = PhylumGlom), shape = 21) +
    facet_wrap(~`Data Origin`) +
    scale_fill_manual("Phylum", values = colour_phylum, breaks = phylum_levels, labels = phylum_labels) +
    ylab(paste0("Assembly Importance ", Habitat))
  
}

vip_imp_plot(imp_plot_df = all_imp_plot$plot_df)
vip_imp_plot(imp_plot_df = palsa_imp_plot$plot_df, Habitat = "Palsa")
vip_imp_plot(imp_plot_df = bog_imp_plot$plot_df, Habitat = "Bog")
vip_imp_plot(imp_plot_df = fen_imp_plot$plot_df, Habitat = "Fen")


#### ====================================================================== ####