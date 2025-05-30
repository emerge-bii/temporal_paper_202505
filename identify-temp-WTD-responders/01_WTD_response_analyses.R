#' ## Interpret Assembly Analysis results
#' This step takes the output from 
#' 03_prepare_and_combine_assembly_analysis_results.R and runs initial figures 
#' and analyses on them.

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
library(vegan); packageVersion("vegan") # for ecological applications
library(viridis)
library(GUniFrac)
library(cowplot)
library(multcompView)
library(here)

# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("identify-temp-WTD-responders", "outputs")
figures.fp <- here("identify-temp-WTD-responders", "figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#### ====================================================================== ####

#### Read in data 
#### ====================================================================== ####
# Rename input
input <- input_ra

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

#### ====================================================================== ####

# Define some helper functions
#### ====================================================================== ####
# mantel plotting settings
mant_plot <- list(
  geom_point(aes(fill = sig), shape = 21, size = 3),
  scale_fill_manual(values = c("black", "white"), breaks = c("yes", "no"),
                    guide = "none"),
  xlab(""),
  theme_bw()
)

get_wtd_mantels <- function(metadata = input$sample_metadata) {
  
  metadata <- metadata %>%
    select(temporal_sample_id, Year__, Habitat__, WTD, contains("pct_time_below_WTD")) %>%
    rename_all(~gsub("pct_time_below_WTD", "", .)) %>%
    rename_if(is.numeric, ~gsub("_", "", .))
  
  # mrm_result <- ecodist::MRM(as.dist(dm) ~ dist(Year) + dist(Tair.degC),
  #     data =metadat_filt)
  mant_df <- data.frame()
  for(i in c("WTD", "21d", "28d", "growing",
             "allgrowing")) {
    metadat_filt <- metadata %>%
      select(all_of(c("temporal_sample_id", i))) %>%
      filter_all(all_vars(!is.na(.)))
    
    dm <- trimmed_mean$rel_abund %>% 
      select(all_of(metadat_filt$temporal_sample_id)) %>%
      t() %>%
      vegdist(method = "bray")
    
    mant_result <- mantel(as.dist(dm), dist(metadat_filt[,i]), method = "spearman")
    
    mant_df_new <- data.frame(Temp_metric = i, 
                              significance = mant_result$signif, 
                              rho = mant_result$statistic)
    
    mant_df <- bind_rows(mant_df, mant_df_new)
  }
  mant_df <- mant_df %>%
    mutate(sig = ifelse(significance < 0.05, "yes", "no")) %>%
    mutate(Temp_metric = factor(Temp_metric, 
                                levels = c("WTD", "21d", "28d", "growing",
                                           "allgrowing")))
  return(mant_df)
}
#### ====================================================================== ####

# What is the variation of time spent below the WT in our samples?
#### ====================================================================== ####
WTD_plot <- input$sample_metadata %>%
  select(temporal_sample_id, Habitat__, Year__, Samp_dist_WT,  
         contains("pct_time_below_WTD")) %>%
  mutate(pct_time_below_WTD_samplingdate = ifelse(Samp_dist_WT > 0, 1, 0)) %>% # Pos = sample below WT, neg = sample above WT
  select(-Samp_dist_WT) %>%
  pivot_longer(cols = contains("WTD"), 
               names_to = "WTD_metric",
               values_to = "pct_time_wtd") %>% 
  mutate(WTD_metric = gsub("pct_time_below_WTD_", "", WTD_metric),
         WTD_metric = gsub("_", "", WTD_metric),
         WTD_metric = factor(WTD_metric, levels = c("samplingdate", "21d", "28d", "growing", "allgrowing"))) %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
  ggplot(aes(x = Year__, y = pct_time_wtd)) + 
  geom_point( size = rel(2)) + 
  facet_grid(WTD_metric~Habitat__, switch = "y") +
  xlab("Year") + ylab("Percent Time Below WT") +
  theme_bw() + 
  theme(axis.title = element_text(size = rel(2)), 
        strip.text = element_text(size = rel(1)))
WTD_plot
ggsave(paste0(figures.fp, "/wet_stordalen.png"), plot = WTD_plot, device = "png", dpi = 300,
       height = 11, width = 10)
#### ====================================================================== ####

# WTD and community dissimilarity
#### ====================================================================== ####
mant_wtd_bogfen <- get_wtd_mantels(metadata = input$sample_metadata %>%
                                  filter(Habitat__ %in% c("Bog", "Fen")))

mant_wtd_bog <- get_wtd_mantels(metadata = input$sample_metadata %>% 
                                   filter(Habitat__ == "Bog"))
mant_wtd_fen <- get_wtd_mantels(metadata = input$sample_metadata %>% 
                                   filter(Habitat__ == "Fen"))

# plot mantel results

mant_temp_bogfen.plot <- mant_wtd_bogfen %>%
  ggplot(aes(x = Temp_metric, y = rho)) +
  mant_plot +
  ggtitle("Mantel correlations: Bray-Curtis with WTD - Bog/Fen")
ggsave(paste0(figure.fp, "/BC_wtd_mant_bogfen.png"), plot = mant_temp_bogfen.plot, device = "png", dpi = 300,
       height = 5, width = 10)

mant_temp_bog.plot <- mant_wtd_bog %>%
  ggplot(aes(x = Temp_metric, y = rho)) +
  mant_plot +
  ggtitle("Mantel correlations: Bray-Curtis with WTD - Bog")
ggsave(paste0(figure.fp, "/BC_wtd_mant_bog.png"), plot = mant_temp_bog.plot, device = "png", dpi = 300,
       height = 5, width = 10)

mant_temp_fen.plot <- mant_wtd_fen %>%
  ggplot(aes(x = Temp_metric, y = rho)) +
  mant_plot +
  ggtitle("Mantel correlations: Bray-Curtis with WTD - Fen")
ggsave(paste0(figure.fp, "/BC_wtd_mant_fen.png"), plot = mant_temp_fen.plot, device = "png", dpi = 300,
       height = 5, width = 10)

#### ====================================================================== ####

# Methane Flux and WTD
#### ====================================================================== ####
methanogenesis_wtd <- sample_metadata %>% 
  filter(Habitat__ %in% c("Bog", "Fen")) %>%
  select(temporal_sample_id, Habitat__, Samp_dist_WT, contains("CH4_Flux_AVE"), 
         contains("pct_time_below_WTD")) %>%
  mutate(pct_time_below_WTD_samplingdate = ifelse(Samp_dist_WT > 0, 1, 0)) %>% # Pos = sample below WT, neg = sample above WT
  select(-Samp_dist_WT) %>%
  pivot_longer(cols = contains("WTD"), 
               names_to = "WTD_metric",
               values_to = "pct_time_wtd") %>% 
  mutate(WTD_metric = gsub("pct_time_below_WTD_", "", WTD_metric),
         WTD_metric = gsub("_", "", WTD_metric),
         WTD_metric = factor(WTD_metric, levels = c("samplingdate", "21d", "28d", "growing", "allgrowing"))) %>% 
  pivot_longer(cols = contains("CH4_Flux_AVE"), 
               names_to = "CH4_Flux_metric",
               values_to = "CH4_Flux") %>%
  mutate(CH4_Flux_metric = gsub("CH4_Flux_AVE_", "", CH4_Flux_metric),
         CH4_Flux_metric = gsub("_before_coring", "", CH4_Flux_metric),
         CH4_Flux_metric = fct_rev(CH4_Flux_metric)) %>%
  ggplot(aes(y = CH4_Flux, x = pct_time_wtd)) +
  geom_point(aes(color = Habitat__)) + 
  scale_color_manual(name = "Habitat", values = colour_habitat, 
                     breaks = habitat_levels, 
                     labels = habitat_levels) +
  facet_grid(CH4_Flux_metric ~ WTD_metric, scales = "free", switch = "y") +
  xlab("Percent of time below water table") +
  theme_bw() +
  theme(strip.placement = "outside")
methanogenesis_wtd

ggsave(paste0(figures.fp, "/methaneflux_wtd.png"), 
       plot = methanogenesis_wtd, device = "png", dpi = 300,
       height = 5, width = 12)
#### ====================================================================== ####
