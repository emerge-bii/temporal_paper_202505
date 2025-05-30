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
outputs.fp <- here("Assembly-analysis", "outputs")
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
    mutate(Variable = gsub("_mag_diff_abund", "", Variable)) %>%
    left_join(input$taxonomy, by = c("Variable" = "genome")) %>%
    mutate(PhylumGlom = ifelse(!is.na(Phylum) & !(Phylum %in% phylum_levels), "other", Phylum)) %>% 
    mutate(PhylumLabel = ifelse(!is.na(PhylumGlom), PhylumGlom, "Environmental Predictor")) %>%
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

prepare_plot_funct_importance <- function(hardcode_model = all_hardcode_model, 
                                        plotN = NULL,
                                        RFStats = T,
                                        RFStatsVerb = T) {
  # hardcode_model = output from running rf model
  # define pathway plotting characteristics
  metapathway_levels <- metapathway_groups %>% arrange(metapath) %>% `$`(metapath) %>% unique()
  colour_metapathway <- RColorBrewer::brewer.pal(10, "Set3")
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
    mutate(Variable = gsub("_func_diff_abund", "", Variable)) %>%
    left_join(input$taxonomy, by = c("Variable" = "genome")) %>%
    mutate(PhylumGlom = ifelse(!is.na(Phylum) & !(Phylum %in% phylum_levels), "other", Phylum)) %>% 
    mutate(PhylumLabel = ifelse(!is.na(PhylumGlom), PhylumGlom, "Environmental Predictor")) %>%
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

all_pred_plot <- ggplot(data = all_imp_plot$plot_df, aes(x = Variable, y = Importance)) +
  all_imp_plot$bar_layers
  
ggsave(all_pred_plot, filename = paste0(figures.fp, "/all_mag_predictors_classification_rf.png"),
       width = 14, height = 10, dpi = 300)

#### ====================================================================== ####

#### ====================================================================== ####
# palsa
palsa_data <- data.rf %>% filter(Habitat_comp == "Palsa_Palsa") %>%
  select(-Habitat_comp) %>%
  mutate(Assembly_Process = fct_drop(Assembly_Process))


palsa_hardcode_model$conf_mat_permute

palsa_imp_plot <- prepare_plot_vip_importance(hardcode_model = palsa_hardcode_model, plotN = 50,
                                            RFStats = T, RFStatsVerb = F)

palsa_pred_plot <- ggplot(data = palsa_imp_plot$plot_df, aes(x = Variable, y = Importance)) +
  palsa_imp_plot$bar_layers

ggsave(palsa_pred_plot, filename = paste0(figures.fp, "/palsa_mag_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)

bog_imp_plot <- prepare_plot_vip_importance(hardcode_model = bog_hardcode_model, plotN = 50,
                                              RFStats = T, RFStatsVerb = F)

bog_pred_plot <- ggplot(data = bog_imp_plot$plot_df, aes(x = Variable, y = Importance)) +
  bog_imp_plot$bar_layers

ggsave(bog_pred_plot, filename = paste0(figures.fp, "/bog_mag_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)

fen_imp_plot <- prepare_plot_vip_importance(hardcode_model = fen_hardcode_model, plotN = 50,
                                              RFStats = T, RFStatsVerb = F)

fen_pred_plot <- ggplot(data = fen_imp_plot$plot_df, aes(x = Variable, y = Importance)) +
  fen_imp_plot$bar_layers

ggsave(fen_pred_plot, filename = paste0(figures.fp, "/fen_mag_predictors_classification_rf.png"),
       width = 18, height = 10, dpi = 300)

#### ====================================================================== ####

# Exploring 20170700_S25_9 - a close associate of M.Stor and A.Stor
#### ====================================================================== ####

fen_hardcode_model$shap_vals_perm$X20170700_S25_9_mag_diff_abund_Homogenous.selection

fen_hardcode_model$rf_permutation$
fen_hardcode_model$shap_vals_perm %>% shapviz::sv_dependence()

#### ====================================================================== ####

# Functional data
#### ====================================================================== ####


#### ====================================================================== ####

