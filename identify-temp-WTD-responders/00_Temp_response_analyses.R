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
library(broom) # for easy extraction of model results
library(broom.mixed)
library(vegan); packageVersion("vegan") # for ecological applications
library(viridis)
library(GUniFrac)
library(cowplot)
library(kableExtra)
library(multcompView)
library(propr)
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

# Read in Temperature Responder VIPs from Dylan
temp_respVIP <- read_csv(here("data", "mags_vip_temperature.csv")) 
#### ====================================================================== ####

# Setup plotting parameters
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), 
                                 c("#737373", "#FFFFFF", "#000000")), 
                          c("blue", "darkblue", "green", "darkgreen", "red", 
                            "darkred", "orange", "darkorange", "purple", 
                            "darkpurple", "yellow", "brown", "grey", "white", 
                            "black"))
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

# Define some helper functions
#### ====================================================================== ####
# Air Temperature Mantel cacluations
get_temp_mantels <- function(metadata = input$sample_metadata) {
  
  metadata <- metadata %>%
    select(temporal_sample_id, Year__, Habitat__, T_air.deg_C, contains("mean_AirTemperature")) %>%
    rename_all(~gsub("mean_AirTemperature", "", .)) %>%
    rename_if(is.numeric, ~gsub("_", "", .))
  
  # mrm_result <- ecodist::MRM(as.dist(dm) ~ dist(Year) + dist(Tair.degC),
  #     data =metadat_filt)
  mant_df <- data.frame()
  for(i in c("Tair.degC", "samplingdate", "7d", "14d", "21d", "28d", "growing",
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
                                levels = c("Tair.degC", "samplingdate", "7d", 
                                           "14d", "21d", "28d", "growing",
                                           "allgrowing")))
  return(mant_df)
}

# temperature mantel plotting settings
mant_plot <- list(
  geom_point(aes(fill = sig), shape = 21, size = 3),
  scale_fill_manual(values = c("black", "white"), breaks = c("yes", "no"),
                    guide = "none"),
  xlab(""),
  theme_bw()
)


# Soil Temperature Mantel Calculations
get_soiltemp_mantels <- function(metadata = input$sample_metadata,
                                 distmet = "bray") {
  
  metadata <- metadata %>%
    dplyr::select(temporal_sample_id, Year__, Habitat__, T_soil.deg_C) %>%
    rename_if(is.numeric, ~gsub("_", "", .))
  
  metadat_filt <- metadata %>%
    dplyr::select(all_of(c("temporal_sample_id", "Tsoil.degC"))) %>%
    filter_all(all_vars(!is.na(.)))
  
  if(distmet == "bray") {
    dm <- trimmed_mean$rel_abund %>% 
      dplyr::select(all_of(metadat_filt$temporal_sample_id)) %>%
      t() %>%
      vegdist(method = "bray")
  } else if(distmet == "unifrac") { # unifrac
    genomes_left <- trimmed_mean$rel_abund %>% 
      dplyr::select(genome, all_of(metadat_filt$temporal_sample_id)) %>%
      mutate(Sum = rowSums(.[2:nrow(metadat_filt)])) %>%
      dplyr::select(genome, Sum) %>%
      filter(Sum > 0)
    
    filt_tree <- ape::keep.tip(input$tree, genomes_left$genome)
    
    otu_tab <- phyloseq::otu_table(trimmed_mean$rel_abund %>% 
                                     dplyr::select(genome, all_of(metadat_filt$temporal_sample_id)) %>%
                                     mutate(Sum = rowSums(.[2:nrow(metadat_filt)])) %>%
                                     filter(Sum > 0) %>%
                                     dplyr::select(-Sum) %>%
                                     column_to_rownames(var = "genome") %>%
                                     t(), taxa_are_rows = F)
    
    phyclas <- phyloseq::phyloseq(otu_tab, filt_tree)
    unifr <- phyloseq::UniFrac(phyclas, weighted = T)
    dm <- unifr
  }
  
  
  mant_result <- mantel(as.dist(dm), dist(metadat_filt[,"Tsoil.degC"]), method = "spearman")
  
  mant_df <- data.frame(Temp_metric = "Tsoil.degC", 
                        significance = mant_result$signif, 
                        rho = mant_result$statistic)
  
  mant_df <- mant_df %>%
    mutate(sig = ifelse(significance < 0.05, "yes", "no"))
  return(mant_df)
}

# Runs soil temperature partial mantels correcting for depth 
get_soiltemp_part.mantels <- function(metadata = input$sample_metadata,
                                      distmet = "bray") {
  
  metadata <- metadata %>%
    dplyr::select(temporal_sample_id, Year__, Habitat__, T_soil.deg_C, DepthAvg__) %>%
    rename_if(is.numeric, ~gsub("_", "", .))
  
  
  metadat_filt <- metadata %>%
    dplyr::select(all_of(c("temporal_sample_id", "Tsoil.degC", "DepthAvg"))) %>%
    filter_all(all_vars(!is.na(.)))
  
  
  if(distmet == "bray") {
    dm <- trimmed_mean$rel_abund %>% 
      dplyr::select(all_of(metadat_filt$temporal_sample_id)) %>%
      t() %>%
      vegdist(method = "bray")
  } else if(distmet == "unifrac") { # unifrac
    genomes_left <- trimmed_mean$rel_abund %>% 
      dplyr::select(genome, all_of(metadat_filt$temporal_sample_id)) %>%
      mutate(Sum = rowSums(.[2:nrow(metadat_filt)])) %>%
      dplyr::select(genome, Sum) %>%
      filter(Sum > 0)
    
    filt_tree <- ape::keep.tip(input$tree, genomes_left$genome)
    
    otu_tab <- phyloseq::otu_table(trimmed_mean$rel_abund %>% 
                                     dplyr::select(genome, all_of(metadat_filt$temporal_sample_id)) %>%
                                     mutate(Sum = rowSums(.[2:nrow(metadat_filt)])) %>%
                                     filter(Sum > 0) %>%
                                     dplyr::select(-Sum) %>%
                                     column_to_rownames(var = "genome") %>%
                                     t(), taxa_are_rows = F)
    
    phyclas <- phyloseq::phyloseq(otu_tab, filt_tree)
    unifr <- phyloseq::UniFrac(phyclas, weighted = T)
    dm <- unifr
  }
  
  mant_result <- mantel.partial(as.dist(dm), dist(metadat_filt[,"Tsoil.degC"]),
                                dist(metadat_filt[,"DepthAvg"]), method = "spearman")
  
  mant_df <- data.frame(Temp_metric = "Tsoil.degC", 
                        significance = mant_result$signif, 
                        rho = mant_result$statistic)
  
  mant_df <- mant_df %>%
    mutate(sig = ifelse(significance < 0.05, "yes", "no"))
  return(mant_df)
}

# Filter and plot metadata over time with temperature
metadata_time_plot <- function(metadata = input$sample_metadata, variable = "DepthAvg__",
                               Temp_var = "T_soil.deg_C") {
  plot_df <- metadata %>%
    filter(!is.na(!!as.name(variable))) %>%
    filter(!is.na(!!as.name(Temp_var))) %>%
    mutate(Habitat__ = factor(Habitat__, levels = habitat_levels),
           Year__ = factor(Year__)) 
  if(nrow(plot_df) < 1) {
    return(NULL)
  } else {
    xlabel <- ifelse(Temp_var == "T_soil.deg_C", "Ground Temperature (°C)", Temp_var)
    var.plot <- ggplot(data = plot_df, aes(y = !!as.name(variable), x = !!as.name(Temp_var), 
                                           group = Year__, fill = Year__, color = Year__)) + 
      geom_point() +
      facet_grid(~Habitat__) +
      geom_smooth(method = "lm") +
      xlab(xlabel) +
      scale_color_manual(values = RColorBrewer::brewer.pal(7, "YlOrBr"), name = "Year") +
      scale_fill_manual(values = RColorBrewer::brewer.pal(7, "YlOrBr"), name = "Year") +
      theme_bw()
    
    return(var.plot)
  }
}

# Temperature responders Functions
# Plotting models individualy with trend lines
responder_plot_funct <- function(dat, scale_vars = T) {
  # dat must contain columns: 
  # plottingmetric, T_soil.deg_C, plottingcoverage, lm.pval, lm.rsquared
  dat <- dat
  
  # Plotting axis options
  xvar <- "T_soil.deg_C"
  xlab <- "Ground Temperature (°C)"
  
  writeLines(paste0("Xvar is ", xvar))
  
  ylab <- unique(dat$coveragemetric)
  ylab <- ifelse(ylab == "coverage", "Cumu Relative Abundance (CLR transformed)", 
                 "Log10(Cumu Relative Abundance) (CLR transformed)")
  plotname <- unique(dat$plottitle)
  
  # Estimate response holding depth at 3 quantiles
  if(scale_vars) {
    pred_dat <- with(dat, 
                     data.frame(T_soil.deg_C = rep(seq(min(T_soil.deg_C), 
                                                       max(T_soil.deg_C), length.out = 100), 3),
                                DepthAvg__ = rep(quantile(DepthAvg__, probs = c(0.25, 0.5, 0.75)), each = 100)))
    
    pred_dat$Predictions <- predict(dat$lm.model[[1]], pred_dat)
    
    
    model_est <- pred_dat %>% 
      mutate(DepthCategory = as.character(as.numeric(as.factor(DepthAvg__))),
             DepthCategory = ifelse(DepthCategory == "1", "0.25 Quartile",
                                    ifelse(DepthCategory == "2", "Median",
                                           "0.75 Quartile"))) %>%
      select(Predictions, T_soil.deg_C, DepthAvg__, DepthCategory) %>%
      rename(plottingcoverage = Predictions)
  } else {
    pred_dat <- with(dat, 
                     data.frame(T_soil.deg_C = rep(seq(min(T_soil.deg_C), 
                                                       max(T_soil.deg_C), length.out = 100), 3),
                                DepthAvg__ = rep(quantile(DepthAvg__, probs = c(0.25, 0.5, 0.75)), each = 100)))
    
    pred_dat$Predictions <- predict(dat$lm.model[[1]], pred_dat)
    
    
    model_est <- pred_dat %>% 
      mutate(DepthCategory = as.character(as.numeric(as.factor(DepthAvg__))),
             DepthCategory = ifelse(DepthCategory == "1", "0.25 Quartile",
                                    ifelse(DepthCategory == "2", "Median",
                                           "0.75 Quartile"))) %>%
      select(Predictions, T_soil.deg_C, DepthAvg__, DepthCategory) %>%
      rename(plottingcoverage = Predictions)
  }
  model_est %>% select(contains("DepthAvg"), DepthCategory) %>% 
    distinct() %>% kableExtra::kable(format = "pipe",
                                     caption = plotname) %>% print()
  
  mod_anot <- dat %>%
    select(lm.rsquared, lm.pval) %>%
    distinct() %>%
    mutate(modlab = paste0(" pval = ", round(lm.pval, digits = 3),
                           "\nadj.r2 = ", round(lm.rsquared, digits = 2)))
  
  responder_status = dat %>% 
    mutate(T_responder = ifelse(T_responder, "T", ""),
           Depth_responder = ifelse(Depth_responder, "D", ""),
           TDepth_responder = ifelse(TDepth_responder, "TD", ""),
           responses = paste0("Significant terms: ", T_responder, ", ",
                              Depth_responder, ", ", TDepth_responder),
           responses = gsub(" , ", "", responses),
           responses = gsub(",$", "", responses))
  
  p <- ggplot(dat, aes(x = !!as.name(xvar), y = plottingcoverage)) +
    geom_point(aes(fill = DepthLumping), shape = 21) +
    #geom_point(aes(color = as.factor(Year__))) +
    annotate(geom = "text", x = Inf, y = Inf, 
             vjust = 1.3, hjust = 1.3,
             label = mod_anot$modlab[1]) +
    annotate(geom = "text", x = -Inf, y = Inf, 
             vjust = 1.3, hjust = -0.1, color = "orange",
             label = responder_status$responses[1]) +
    ylab(ylab) + xlab(xlab) +
    ggtitle(plotname) +
    theme_bw()
  
  plot <- p
  
  # Overwrite plot with lines if pval < 0.05, VIF is okay, and there is a significant temperature response
  if(unique(dat$lm.pval) < 0.05 & unique(dat$VIFOK) & unique(dat$T_responder)) {# & unique(dat$ShapiroOK)) {
    plot <- p + 
      geom_line(data = model_est, 
                aes(x = !!as.name(xvar), y = plottingcoverage,
                    color = as.factor(DepthAvg__))) +
      scale_color_manual(values = RColorBrewer::brewer.pal(3,"Oranges"),
                         name = "Depth (cm)")
  }
  
  return(plot)
}

compare_and_plot_slopes <- function(resp_model, slope = "Temperature",
                                    plot_dir, dataType = "metaG",
                                    plot_height = 7, plot_width = 14) {
  
  # Figure out what slope to filter on
  if(slope == "Temperature") {
    slopeterm <- "SoilTemp"
    xlabel <- "Temperature"
    plottitle <- "Slope of response to ground temperature"
    fname <- paste0(dataType, "_temperature_slope_comparison.png")
    writeLines("Plotting temperature slope terms")
  } else if(slope == "Depth") {
    slopeterm <- "Depth"
    xlabel <- "Depth"
    plottitle <- "Slope of response to depth"
    fname <- paste0(dataType, "_depth_slope_comparison.png")
    writeLines("Plotting depth slope terms")
  } else if (slope == "Interaction") {
    slopeterm <- "SoilTemp:Depth"
    xlabel <- "Temperature and Depth Interaction"
    plottitle <- "Slope of temperature and depth interaction"
    fname <- paste0(dataType, "_interaction_slope_comparison.png")
    writeLines("Plotting interaction slope terms")
  } else {
    writeLines("Please correctly specify the slope term, it can be either: Temperature, Depth, or Interaction")
  }
  
  # Different settings for genomes vs. pathways
  # Pathway/Genome rename
  path_vip_rename <- c(responder = "Pathways", responder = "genome")
  
  # Add taxonomy if we're working with genomes
  if("genome" %in% colnames(resp_model)) {
    resp_model <- left_join(resp_model,
                            input$taxonomy, by = "genome") %>%
      # Fix taxonomy names
      mutate(across(Domain:Species, ~gsub("^[pcofgs]__", "", .))) %>%
      # change responder name to include family designation
      mutate(responder = paste0(genome, " (", Family, ")"))
  } else if ("Pathways" %in% colnames(resp_model)) {
    resp_model <- resp_model %>%
      rename(responder = Pathways)
  } else {
    writeLines("Neither Pathways nor genome found, please check input data")
  }
  
  
  # Extract the slope term information
  temp_estimates <- resp_model %>% #filter(T_responder|TDepth_responder) %>%
    mutate(plottedcoverage = ifelse(lm.resid.shapiro < 0.05, "coverage_log10", "coverage")) %>%
    select(Habitat__, responder, lm.tidied, plottedcoverage, ShapiroOK,
           VIFOK, T_responder, T_padj, Depth_padj, TDepth_padj, lm.pval, lm.rsquared,
           any_of(c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "genome"))) %>%
    unnest(lm.tidied) %>% 
    filter(term != "(Intercept)") %>%
    # Make clearer term names
    mutate(term = case_when(
      term == "scale(T_soil.deg_C):scale(DepthAvg__)" ~ "SoilTemp:Depth",
      term == "scale(T_soil.deg_C)" ~ "SoilTemp",
      term == "scale(DepthAvg__)" ~ "Depth",
      TRUE ~ NA_character_),
      term = factor(term, levels = c("SoilTemp:Depth", "Depth", "SoilTemp"))) %>%
    # Identify when a term is significant
    mutate(significant = case_when(
      T_padj < 0.05 & term == "SoilTemp" ~ "sig",
      Depth_padj < 0.05 & term == "Depth" ~ "sig",
      TDepth_padj < 0.05 & term == "SoilTemp:Depth" ~ "sig",
      TRUE ~ "not sig"),
      # Specify VIFOK as "cannot determine" rather than significant/not significant
      significant = ifelse(VIFOK, significant, "cannot determine"),
      # Specify ShapiroOK as "not sig"
      #significant = ifelse(ShapiroOK, significant, "not sig"),
      TempEst = ifelse(term == "SoilTemp", estimate, NA)) %>% 
    ungroup() %>%
    arrange(VIFOK, desc(T_responder), TempEst) %>% 
    mutate(order = ifelse(term == "SoilTemp" & T_padj < 0.05 & VIFOK, row_number(), NA)) %>%
    group_by(Habitat__, responder) %>% fill(order) %>% ungroup()
  
  # Plot slope estimates
  slope_est_plot <- temp_estimates %>%
    # Filter based on slope terms
    filter(term == slopeterm) %>%
    ggplot(aes(x = estimate, y = Habitat__)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high,
                      group = responder, linetype = significant, color = Habitat__), 
                  width = 0.3) +
    geom_point(aes(group = responder, shape = significant, 
                   color = Habitat__)) +
    scale_shape_manual("Significant", values = c(19, 1, 8), 
                       breaks = c("sig", "not sig", "cannot determine"), 
                       labels = c("Significant", "Not Significant", "Temperature and depth\ntoo correlated\nto determine")) +
    scale_linetype_manual("Significant", values = c("solid", "dashed", "dotted"), 
                          breaks = c("sig", "not sig", "cannot determine"), 
                          labels = c("Significant", "Not Significant", "Temperature and depth\ntoo correlated\nto determine")) +
    scale_color_manual("Habitat", values = colour_habitat, 
                       breaks = habitat_levels, labels = habitat_levels) +
    facet_wrap(.~responder, ncol = 1, strip.position = "left") +
    xlab(paste0(xlabel, " slope estimate (+/- 95% confidence interval)")) + ylab("") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0),
          panel.spacing.y = unit(0, "lines")) +
    ggtitle(plottitle)
  
  # Save plot
  ggsave(path = plot_dir, filename = fname, plot = slope_est_plot,
         device = "png", height = plot_height, width = plot_width)
  
  return(temp_estimates)
}

# Temperature MAG responders
calc_temp_resp_mags <- function(data = trimmed_mean$rel_abund,
                                mag_list = temp_respVIP,
                                dataType = "VIP",
                                habitat = c("Palsa", "Bog", "Fen"),
                                scale_vars = T, skip_resids = T) {
  
  mags_of_interest <- mag_list$mag %>% unique()
  
  abund_TM_RA_clr <- data %>%
    # translate dataframe to work with propr expectation
    column_to_rownames(var = "genome") %>%
    t() %>% as_tibble(rownames = "temporal_sample_id") %>%
    column_to_rownames(var = "temporal_sample_id") %>% 
    # Filter out MAGs with >90% zero abundance
    dplyr::select(where(~ sum(. != 0) > length(.) / 10)) %>% 
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
  temp_otu_resp <- input$sample_metadata %>%
    dplyr::select(temporal_sample_id, Year__, Habitat__, T_soil.deg_C, DepthAvg__, DepthLumping ) %>%
    # Remove samples that have no soil temperature data
    dplyr::filter(!is.na(T_soil.deg_C)) %>% 
    # Filter abundance table for samples in the sample metadata and with soil temp data
    left_join(abund_TM_RA_clr, by = c("temporal_sample_id")) %>% 
    # Gather data about module membership, and taxonomy
    left_join(mag_list, by = c("genome" = "mag")) %>%
    left_join(input$taxonomy, by = "genome") %>%
    # Match the habitat of the module to the habitat of the samples, filter out those samples where the module habitat doesn't match the sample habitat; 
    filter(Habitat__ == habitat) %>%
    # Fix taxonomy names
    mutate(across(Domain:Species, ~gsub("^[pcofgs]__", "", .))) %>%
    # Change the name of the abundance column for convenience
    rename(coverage = CLR_cumu_relabund) %>%
    # Fix habitat levels
    mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>% 
    group_by(Habitat__, genome) %>%
    add_count() %>% 
    # Only use genomes who are present >4 times in a habitat; 
    # this is true for all VIPs at the moment
    filter(n > 4) %>% # lowest is 43, we're probably good
    dplyr::select(-n) %>%
    # Use mutate to standardize our explantory variables of interest
    mutate(across(all_of(c("T_soil.deg_C", "DepthAvg__")), ~scale(.)[,1], .names = "{.col}_scl")) %>%
    # Calculate log10 of coverage for models where residuals are not normally distributed
    mutate(coverage_log10 = log10(coverage + 1 - min(coverage)),  # Because we're using clr transformed data, we will translate the data to have no negatives first by adding the minimum value to each group
           coverage_cuberoot = sign(coverage) * (abs(coverage))^(1/3))
  # Plot the data
  overview_plot <- temp_otu_resp %>%
    ggplot(aes(y = coverage, x = T_soil.deg_C)) + 
    geom_point(aes(fill = DepthLumping ), shape = 21, size = 3, alpha = 0.5) + 
    facet_wrap(Habitat__~genome, scales = "free_y")
  
  if(!dir.exists(paste0(figures.fp, "/", dataType))) {dir.create(paste0(figures.fp,"/", dataType))}
  
  ggsave(overview_plot, 
         filename = paste0(figures.fp, "/", dataType, "/", dataType, "_overview_plot.png"),
         device = "png", height = 14, width = 20)
  
  # Models to use
  # If scale_vars = T, use the scaled variables; else use unscaled data
  if(scale_vars) {
    cov_mod <- as.formula("coverage ~ scale(T_soil.deg_C) + scale(DepthAvg__) + scale(T_soil.deg_C)*scale(DepthAvg__)")
    cov_log10_mod <- as.formula("coverage_log10 ~ scale(T_soil.deg_C) + scale(DepthAvg__) + scale(T_soil.deg_C)*scale(DepthAvg__)")
    t_var <- "scale(T_soil.deg_C)"
    d_var <- "scale(DepthAvg__)"
  } else {
    cov_mod <- as.formula("coverage ~ T_soil.deg_C + DepthAvg__ + T_soil.deg_C*DepthAvg__")
    cov_log10_mod <- as.formula("coverage_log10 ~ T_soil.deg_C + DepthAvg__ + T_soil.deg_C*DepthAvg__")
    t_var <- "T_soil.deg_C"
    d_var <- "DepthAvg__"
  }
  
  writeLines("Running Models")
  temp_otu_resp_models <- temp_otu_resp %>%
    nest() %>% 
    mutate(lm.model = purrr::map(data, ~lm(cov_mod, data = .)),
           # significant value indicates residuals not normally distributed
           lm.resid.shapiro = purrr::map_dbl(lm.model, ~ pluck(.x,"residuals") %>% 
                                               shapiro.test(.) %>% pluck("p.value")),
           # If the shapiro test fails, try to redo with log10 coverage instead
           lm.model = ifelse(lm.resid.shapiro < 0.05, 
                             purrr::map(data, ~lm(cov_log10_mod, data = .)),
                             lm.model),
           # coverage metric used after initial shapiro test
           coveragemetric = ifelse(lm.resid.shapiro < 0.05, 
                                   "coverage_log10",
                                   "coverage"),
           # redo shapiro test
           lm.resid.shapiro = purrr::map_dbl(lm.model, ~ pluck(.x,"residuals") %>% 
                                               shapiro.test(.) %>% pluck("p.value")),
           # Plot residual normality
           lm.resid.plot = purrr::map(lm.model, ~ pluck(.x,"residuals") %>% 
                                        data.frame(resids = .) %>%
                                        mutate(idealresidsx = seq(min(resids), max(resids), length.out = n()),
                                               idealresids = dnorm(x = idealresidsx, 
                                                                   mean = unique(mean(resids)),
                                                                   sd = unique(sd(resids)))) %>%
                                        ggplot() + 
                                        geom_density(mapping = aes(x = resids),
                                                     color = "blue") +
                                        geom_line(mapping = aes(x = idealresidsx, y = idealresids),
                                                  color = "red")),
           # add the shapiro annotations onto the residual plots
           lm.resid.plot = purrr::map2(lm.resid.plot, 
                                       lm.resid.shapiro, ~ .x + annotate(geom = "text", 
                                                                         x = Inf, y = Inf,
                                                                         hjust = 1, vjust = 1.5,
                                                                         label = paste0("Shapiro = ",
                                                                                        round(.y, digits = 3)))),
           residplotname = paste0(Habitat__, "_", genome, "_resid_plot.png"),
           # Tidy up lm results
           lm.tidied = purrr::map(lm.model, ~broom::tidy(., conf.int = T, conf.level = 0.95)),
           lm.rsquared = purrr::map_dbl(lm.model, ~broom::glance(.)$adj.r.squared),
           lm.pval = purrr::map_dbl(lm.model, ~broom::glance(.)$p.value),
           # Check variance inflation factor and shapiro
           lm.vif = purrr::map(lm.model, ~car::vif(.)), # we don't care about the interaction term vif
           ShapiroOK = ifelse(lm.resid.shapiro > 0.05, TRUE, FALSE),
           VIFOK = purrr::map_lgl(lm.vif, ~(pluck(., t_var) < 5 & pluck(., d_var) < 5)),
           T_effect = purrr::map_dbl(lm.tidied, ~pluck(., "estimate", 2)),
           T_pval = purrr::map_dbl(lm.tidied, ~pluck(., "p.value", 2)),
           Depth_effect = purrr::map_dbl(lm.tidied, ~pluck(., "estimate", 3)),
           Depth_pval = purrr::map_dbl(lm.tidied, ~pluck(., "p.value", 3)),
           TDepth_effect = purrr::map_dbl(lm.tidied, ~pluck(., "estimate", 4)),
           TDepth_pval = purrr::map_dbl(lm.tidied, ~pluck(., "p.value", 4))) %>%
    ungroup() %>%
    group_by(Habitat__) %>%
    mutate(T_padj = p.adjust(T_pval, method = "fdr"),
           Depth_padj = p.adjust(Depth_pval, method = "fdr"),
           TDepth_padj = p.adjust(TDepth_pval, method = "fdr"),
           T_responder = ifelse(T_padj < 0.05, T, F),
           Depth_responder = ifelse(Depth_padj < 0.05, T, F),
           TDepth_responder = ifelse(TDepth_padj < 0.05, T, F))
  
  # Plot residual assuptions plots
  if(!skip_resids) {
    writeLines("Saving Residual Plots")
    if(!dir.exists(paste0(figures.fp, "/", dataType, "/path_resid_plots"))) {dir.create(paste0(figures.fp,"/", dataType, "/path_resid_plots"))}
    temp_otu_resp_models %>%
      ungroup() %>%
      arrange(desc(lm.resid.shapiro)) %>%
      select(residplotname, lm.resid.plot) %>%
      purrr::pmap(., ~ggsave(path = paste0(figures.fp, "/", dataType, "/path_resid_plots"), filename = .x, plot = .y,
                             device = "png", height = 7, width = 7))
  }
  
  
  # Plotting models 1 at a time
  writeLines("Saving Model Plots")
  response_plots_dir <- here(figures.fp,"/", dataType ,"T_responder_vip_plots")
  if(!dir.exists(response_plots_dir)){dir.create(response_plots_dir)}
  
  temp_otu_resp_models %>%
    unnest(c(data)) %>%
    #filter(T_responder) %>%
    mutate(plottingcoverage = ifelse(coveragemetric == "coverage_log10", coverage_log10, coverage),
           plottitle = paste0("Temperature response for ", Habitat__, " ", module, " module: ", genome, " (", Family, ")")) %>%
    select(Habitat__, genome, T_soil.deg_C, DepthAvg__, DepthLumping, lm.model,
           T_soil.deg_C_scl, DepthAvg___scl, plottingcoverage, coveragemetric, ShapiroOK, VIFOK, 
           lm.pval, lm.rsquared, plottitle, T_responder, Depth_responder, TDepth_responder) %>%
    ungroup() %>% group_by(Habitat__, genome) %>% nest() %>% #pluck("data",3)
    # Add back in the model for line plotting
    # Add back in module infomation for plotting
    mutate(responder_plot = purrr::map(data, ~responder_plot_funct(., scale_vars = scale_vars)),
           responderplotname = paste0(Habitat__, "_", genome, "_responder_plot.png")) %>%
    #pluck("responder_plot", 3)
    ungroup() %>% # removes grouping variables for pmap
    select(responderplotname, responder_plot) %>%
    purrr::pmap(., ~ggsave(path = response_plots_dir, filename = .x, plot = .y,
                           device = "png", height = 7, width = 10))
  
  return(temp_otu_resp_models)
  
}
# Pathway responders
calc_temp_resp_pathways <- function(data = pathway_abundance_cumu,
                                    dataType = "metaG",
                                    habitat = c("Palsa", "Bog", "Fen"),
                                    scale_vars = T, skip_resids = T) {
  
  pathways_of_interest <- c("acetoclastic_methanogenesis-all",
                            "acetogenesis-simple", #(NO METATs)
                            "Wood_Ljugndahl-acetogen", 
                            "hydrogenotrophic_methanogenesis-all",
                            "carbon_redox-assimilatory_methanotrophy",
                            "carbon_redox-dissimilatory_methanotrophy", 
                            "nitrogen_redox-nitrogen_fixation",
                            "sulfur_redox-dissimilatory_sulfate_reduction", 
                            "sulfur_redox-dissimilatory_sulfur_oxidation",
                            "glycerol_degradation-glycerol",
                            "glycerol_degradation-glycerone", 
                            "ethanol_fermentation-all",
                            "lactate_fermentation-all", 
                            "propionate_fermentation-succinate_pathway",
                            "propionate_fermentation-acrylate_pathway",
                            "butanoate_fermentation-acetyl_coa_pathway",
                            "CAZy-Xylans", 
                            "CAZy-Chitin", # check with sam and derek that this is the correct annotation for chitin
                            "CAZy-Amorphous_Cellulose",
                            "CAZy-Xyloglucan")
  if(dataType == "metaG") {
    data <- rename(data, temporal_sample_id = sample)
  }
  
  abund_TM_RA_clr <- data %>%
    left_join(sample_metadata %>%
                dplyr::select(temporal_sample_id, Habitat__), 
              by = c("temporal_sample_id")) %>%
    #filter(Habitat__ %in% habitat) %>%
    column_to_rownames(var = "temporal_sample_id") %>%
    select(-c(Habitat__)) %>% 
    # Filter out pathways with >90% zero abundance
    dplyr::select(where(~ sum(. != 0) > length(.) / 10)) %>%
    mutate(across(everything(), transform_perc)) %>%
    propr(metric = "rho") %>%
    pluck("logratio") %>%
    rownames_to_column(var = "temporal_sample_id") %>%
    select(temporal_sample_id, any_of(pathways_of_interest)) %>%
    pivot_longer(-temporal_sample_id, names_to = "Pathways", values_to = "CLR_cumu_relabund")
  
  
  # Combine with metadata
  # Assemble all relevant data about VIPs
  temp_path_resp <- input$sample_metadata %>%
    dplyr::select(temporal_sample_id, Year__, Habitat__, T_soil.deg_C, DepthAvg__, DepthLumping ) %>%
    # Remove samples that have no soil temperature data
    dplyr::filter(!is.na(T_soil.deg_C)) %>% 
    # Remove samples that have no metaT/G data 
    dplyr::filter(temporal_sample_id %in% unique(abund_TM_RA_clr$temporal_sample_id)) %>%
    # Gather data about VIP abundance, module membership, and taxonomy
    left_join(abund_TM_RA_clr, by = c("temporal_sample_id")) %>%
    # Change the name of the abundance column for convenience
    rename(coverage = CLR_cumu_relabund) %>%
    # Fix habitat levels
    mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>% 
    group_by(Habitat__, Pathways) %>%
    add_count() %>% 
    # Only use pathways who are present >4 times in a habitat; 
    # this is true for all VIPs at the moment
    filter(n > 4) %>% # lowest is 43 for fen, we're probably good
    dplyr::select(-n) %>%
    # Use mutate to standardize our explantory variables of interest
    mutate(across(all_of(c("T_soil.deg_C", "DepthAvg__")), ~scale(.)[,1], .names = "{.col}_scl")) %>%
    # Calculate log10 of coverage for models where residuals are not normally distributed
    mutate(coverage_log10 = log10(coverage + 1 - min(coverage)),  # Because we're using clr transformed data, we will translate the data to have no negatives first by adding the minimum value to each group
           coverage_cuberoot = sign(coverage) * (abs(coverage))^(1/3))
  # Plot the data
  overview_plot <- temp_path_resp %>%
    ggplot(aes(y = coverage, x = T_soil.deg_C)) + 
    geom_point(aes(fill = DepthLumping ), shape = 21, size = 3, alpha = 0.5) + 
    facet_wrap(Habitat__~Pathways, scales = "free_y")
  
  if(!dir.exists(paste0(figures.fp, "/", dataType))) {dir.create(paste0(figures.fp,"/", dataType))}
  
  ggsave(overview_plot, 
         filename = paste0(figures.fp, "/", dataType, "/", dataType, "_overview_plot.png"),
         device = "png", height = 14, width = 20)
  
  # Models to use
  # If scale_vars = T, use the scaled variables; else use unscaled data
  if(scale_vars) {
    cov_mod <- as.formula("coverage ~ scale(T_soil.deg_C) + scale(DepthAvg__) + scale(T_soil.deg_C)*scale(DepthAvg__)")
    cov_log10_mod <- as.formula("coverage_log10 ~ scale(T_soil.deg_C) + scale(DepthAvg__) + scale(T_soil.deg_C)*scale(DepthAvg__)")
    t_var <- "scale(T_soil.deg_C)"
    d_var <- "scale(DepthAvg__)"
  } else {
    cov_mod <- as.formula("coverage ~ T_soil.deg_C + DepthAvg__ + T_soil.deg_C*DepthAvg__")
    cov_log10_mod <- as.formula("coverage_log10 ~ T_soil.deg_C + DepthAvg__ + T_soil.deg_C*DepthAvg__")
    t_var <- "T_soil.deg_C"
    d_var <- "DepthAvg__"
  }
  
  writeLines("Running Models")
  temp_path_resp_models <- temp_path_resp %>%
    nest() %>% 
    mutate(lm.model = purrr::map(data, ~lm(cov_mod, data = .)),
           # significant value indicates residuals not normally distributed
           lm.resid.shapiro = purrr::map_dbl(lm.model, ~ pluck(.x,"residuals") %>% 
                                               shapiro.test(.) %>% pluck("p.value")),
           # If the shapiro test fails, try to redo with log10 coverage instead
           lm.model = ifelse(lm.resid.shapiro < 0.05, 
                             purrr::map(data, ~lm(cov_log10_mod, data = .)),
                             lm.model),
           # coverage metric used after initial shapiro test
           coveragemetric = ifelse(lm.resid.shapiro < 0.05, 
                                   "coverage_log10",
                                   "coverage"),
           # redo shapiro test
           lm.resid.shapiro = purrr::map_dbl(lm.model, ~ pluck(.x,"residuals") %>% 
                                               shapiro.test(.) %>% pluck("p.value")),
           # Plot residual normality
           lm.resid.plot = purrr::map(lm.model, ~ pluck(.x,"residuals") %>% 
                                        data.frame(resids = .) %>%
                                        mutate(idealresidsx = seq(min(resids), max(resids), length.out = n()),
                                               idealresids = dnorm(x = idealresidsx, 
                                                                   mean = unique(mean(resids)),
                                                                   sd = unique(sd(resids)))) %>%
                                        ggplot() + 
                                        geom_density(mapping = aes(x = resids),
                                                     color = "blue") +
                                        geom_line(mapping = aes(x = idealresidsx, y = idealresids),
                                                  color = "red")),
           # add the shapiro annotations onto the residual plots
           lm.resid.plot = purrr::map2(lm.resid.plot, 
                                       lm.resid.shapiro, ~ .x + annotate(geom = "text", 
                                                                         x = Inf, y = Inf,
                                                                         hjust = 1, vjust = 1.5,
                                                                         label = paste0("Shapiro = ",
                                                                                        round(.y, digits = 3)))),
           residplotname = paste0(Habitat__, "_", Pathways, "_resid_plot.png"),
           # Tidy up lm results
           lm.tidied = purrr::map(lm.model, ~broom::tidy(., conf.int = T, conf.level = 0.95)),
           lm.rsquared = purrr::map_dbl(lm.model, ~broom::glance(.)$adj.r.squared),
           lm.pval = purrr::map_dbl(lm.model, ~broom::glance(.)$p.value),
           # Check variance inflation factor
           lm.vif = purrr::map(lm.model, ~car::vif(.)), # we don't care about the interaction term vif
           ShapiroOK = ifelse(lm.resid.shapiro > 0.05, TRUE, FALSE),
           VIFOK = purrr::map_lgl(lm.vif, ~(pluck(., t_var) < 5 & pluck(., d_var) < 5)),
           T_effect = purrr::map_dbl(lm.tidied, ~pluck(., "estimate", 2)),
           T_pval = purrr::map_dbl(lm.tidied, ~pluck(., "p.value", 2)),
           Depth_effect = purrr::map_dbl(lm.tidied, ~pluck(., "estimate", 3)),
           Depth_pval = purrr::map_dbl(lm.tidied, ~pluck(., "p.value", 3)),
           TDepth_effect = purrr::map_dbl(lm.tidied, ~pluck(., "estimate", 4)),
           TDepth_pval = purrr::map_dbl(lm.tidied, ~pluck(., "p.value", 4))) %>%
    ungroup() %>%
    group_by(Habitat__) %>%
    mutate(T_padj = p.adjust(T_pval, method = "fdr"),
           Depth_padj = p.adjust(Depth_pval, method = "fdr"),
           TDepth_padj = p.adjust(TDepth_pval, method = "fdr"),
           T_responder = ifelse(T_padj < 0.05, T, F),
           Depth_responder = ifelse(Depth_padj < 0.05, T, F),
           TDepth_responder = ifelse(TDepth_padj < 0.05, T, F))
  
  # Plot residual assuptions plots
  if(!skip_resids) {
    writeLines("Saving Residual Plots")
    if(!dir.exists(paste0(figures.fp, "/", dataType, "/path_resid_plots"))) {dir.create(paste0(figures.fp,"/", dataType, "/path_resid_plots"))}
    temp_path_resp_models %>%
      ungroup() %>%
      arrange(desc(lm.resid.shapiro)) %>%
      select(residplotname, lm.resid.plot) %>%
      purrr::pmap(., ~ggsave(path = paste0(figures.fp, "/", dataType, "/path_resid_plots"), filename = .x, plot = .y,
                             device = "png", height = 7, width = 7))
  }
  
  
  # Plotting models 1 at a time
  writeLines("Saving Model Plots")
  response_plots_dir <- here(figures.fp,"/", dataType ,"T_responder_path_plots")
  if(!dir.exists(response_plots_dir)){dir.create(response_plots_dir)}
  
  temp_path_resp_models %>%
    unnest(c(data)) %>% #View()
    #filter(T_responder) %>%
    mutate(plottingcoverage = ifelse(coveragemetric == "coverage_log10", coverage_log10, coverage),
           plottitle = paste0("Temperature response for ", Habitat__, " ", Pathways)) %>%
    select(Habitat__, Pathways, T_soil.deg_C, DepthAvg__, DepthLumping, lm.model,
           T_soil.deg_C_scl, DepthAvg___scl, plottingcoverage, coveragemetric, ShapiroOK, VIFOK, 
           lm.pval, lm.rsquared, plottitle, T_responder, Depth_responder, TDepth_responder) %>%
    ungroup() %>% group_by(Habitat__, Pathways) %>% nest() %>% #pluck("data",3)
    # Add back in the model for line plotting
    # Add back in module infomation for plotting
    mutate(responder_plot = purrr::map(data, ~responder_plot_funct(., scale_vars = scale_vars)),
           responderplotname = paste0(Habitat__, "_", Pathways, "_responder_plot.png")) %>%
    #pluck("responder_plot", 3)
    ungroup() %>% # removes grouping variables for pmap
    select(responderplotname, responder_plot) %>%
    purrr::pmap(., ~ggsave(path = response_plots_dir, filename = .x, plot = .y,
                           device = "png", height = 7, width = 10))
  
  return(temp_path_resp_models)
  
}


#### ====================================================================== ####

# Air Temperature and community dissimilarity
#### ====================================================================== ####
mant_temp_all <- get_temp_mantels(metadata = input$sample_metadata %>%
                                    filter(DepthAvg__ < 5))

mant_temp_palsa <- get_temp_mantels(metadata = input$sample_metadata %>% 
                                   filter(Habitat__ == "Palsa")  %>%
                                     filter(DepthAvg__ < 5))
mant_temp_bog <- get_temp_mantels(metadata = input$sample_metadata %>% 
                                     filter(Habitat__ == "Bog") %>%
                                    filter(DepthAvg__ < 5))
mant_temp_fen <- get_temp_mantels(metadata = input$sample_metadata %>% 
                                   filter(Habitat__ == "Fen") %>%
                                    filter(DepthAvg__ < 5))

# plot mantel results
mant_temp_all.plot <- mant_temp_all %>%
  ggplot(aes(x = Temp_metric, y = rho)) +
  mant_plot +
  ggtitle("Mantel correlations: Bray-Curtis with air temperature")
mant_temp_all.plot
ggsave(paste0(figures.fp, "/BC_temp_mant_all.png"), plot = mant_temp_all.plot, device = "png", dpi = 300,
       height = 5, width = 10)


mant_temp_palsa.plot <- mant_temp_palsa %>%
  ggplot(aes(x = Temp_metric, y = rho)) +
  mant_plot +
  ggtitle("Mantel correlations: Bray-Curtis with temperature - Palsa")
ggsave(paste0(figures.fp, "/BC_temp_mant_palsa.png"), plot = mant_temp_palsa.plot, device = "png", dpi = 300,
       height = 5, width = 10)

mant_temp_bog.plot <- mant_temp_bog %>%
  ggplot(aes(x = Temp_metric, y = rho)) +
  mant_plot +
  ggtitle("Mantel correlations: Bray-Curtis with air temperature - Bog")
ggsave(paste0(figures.fp, "/BC_temp_mant_bog.png"), plot = mant_temp_bog.plot, device = "png", dpi = 300,
       height = 5, width = 10)

mant_temp_fen.plot <- mant_temp_fen %>%
  ggplot(aes(x = Temp_metric, y = rho)) +
  mant_plot +
  ggtitle("Mantel correlations: Bray-Curtis with air temperature - Fen")
ggsave(paste0(figures.fp, "/BC_temp_mant_fen.png"), plot = mant_temp_fen.plot, device = "png", dpi = 300,
       height = 5, width = 10)


#### ====================================================================== ####

# Ground Temperature and community dissimilarity
#### ====================================================================== ####
# With Bray-curtis Distance
mant_soiltemp_all <- get_soiltemp_mantels(metadata = input_ra$sample_metadata)

mant_soiltemp_palsa <- get_soiltemp_mantels(metadata = input_ra$sample_metadata %>% 
                                      filter(Habitat__ == "Palsa")) 
mant_soiltemp_bog <- get_soiltemp_mantels(metadata = input_ra$sample_metadata %>% 
                                    filter(Habitat__ == "Bog"))
mant_soiltemp_fen <- get_soiltemp_mantels(metadata = input_ra$sample_metadata %>% 
                                    filter(Habitat__ == "Fen"))

mant_soiltemp_all
mant_soiltemp_palsa
mant_soiltemp_bog
mant_soiltemp_fen

# With weighted Unifrac distance

mantuf_soiltemp_all <- get_soiltemp_mantels(metadata = input_ra$sample_metadata, distmet = "unifrac")

mantuf_soiltemp_palsa <- get_soiltemp_mantels(metadata = input_ra$sample_metadata %>% 
                                              filter(Habitat__ == "Palsa"), distmet = "unifrac") 
mantuf_soiltemp_bog <- get_soiltemp_mantels(metadata = input_ra$sample_metadata %>% 
                                            filter(Habitat__ == "Bog"), distmet = "unifrac")
mantuf_soiltemp_fen <- get_soiltemp_mantels(metadata = input_ra$sample_metadata %>% 
                                            filter(Habitat__ == "Fen"), distmet = "unifrac")

mantuf_soiltemp_all
mantuf_soiltemp_palsa
mantuf_soiltemp_bog
mantuf_soiltemp_fen


#### ====================================================================== ####

# Ground Temperature and community dissimilarity (Depth corrected)
# use partial mantels to account for depth
#### ====================================================================== ####
# Bray-Curtis distance
mant_soiltemp_all.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata)

mant_soiltemp_palsa.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata %>% 
                                              filter(Habitat__ == "Palsa"))
mant_soiltemp_bog.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata %>% 
                                            filter(Habitat__ == "Bog"))
mant_soiltemp_fen.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata %>% 
                                            filter(Habitat__ == "Fen"))

mant_soiltemp_all.part
mant_soiltemp_palsa.part
mant_soiltemp_bog.part
mant_soiltemp_fen.part

# With weighted unifrac distance
mantuf_soiltemp_all.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata,
                                                      distmet = "unifrac")

mantuf_soiltemp_palsa.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata %>% 
                                                        filter(Habitat__ == "Palsa"),
                                                        distmet = "unifrac")
mantuf_soiltemp_bog.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata %>% 
                                                      filter(Habitat__ == "Bog"),
                                                      distmet = "unifrac")
mantuf_soiltemp_fen.part <- get_soiltemp_part.mantels(metadata = input$sample_metadata %>% 
                                                      filter(Habitat__ == "Fen"),
                                                      distmet = "unifrac")

mantuf_soiltemp_all.part
mantuf_soiltemp_palsa.part
mantuf_soiltemp_bog.part
mantuf_soiltemp_fen.part


#### ====================================================================== ####

# Do the relationships of temperature and metadata values vary over time?
#### ====================================================================== ####
# Air temperature and other metadata over time
metadata_airtemp_plots <- lapply(names(input$sample_metadata %>% select_if(is.numeric)), function(x) {
  # If column name is a forbidden name, skip it and return NULL
  if(x %in% c("Habitat__", "Year__", "Month__", "Core__",
              "UpdateDate__.Biogeochemistry","UpdateDate__.DepthInfo")) {
    return(NULL)
  } else {
    writeLines(x)
    metadat <- input$sample_metadata %>%
      select(Habitat__, Year__, DepthLumping, T_air.deg_C, !!as.name(x))
    
    temp_plot <- metadata_time_plot(metadata = metadat, variable = x, Temp_var = "T_air.deg_C")
    if(!dir.exists(paste0(figures.fp, "/TempRelated"))) {dir.create(paste0(figures.fp, "/TempRelated"))}
    plot_name <- paste0(figures.fp, "/TempRelated/airtemp_", x, "_plot.png")
    ggsave(plot_name, 
           plot = temp_plot, device = "png", dpi = 300,
           height = 11, width = 12)
    return(temp_plot)
  }
})

# Soil temperature and other metadata over time
metadata_soiltemp_plots <- lapply(names(input$sample_metadata %>% select_if(is.numeric)), 
                                  function(x) {
  # If column name is a forbidden name, skip it and return NULL
  if(x %in% c("Habitat__", "Year__", "Month__", "Core__",
              "UpdateDate__.Biogeochemistry","UpdateDate__.DepthInfo")) {
    return(NULL)
  } else {
    writeLines(x)
    metadat <- input$sample_metadata %>%
      select(Habitat__, Year__, DepthLumping, T_soil.deg_C, !!as.name(x))
    
    temp_plot <- metadata_time_plot(metadata = metadat, variable = x)
    plot_name <- paste0(figures.fp, "/TempRelated/soiltemp_", x, "_plot.png")
    if(!dir.exists(paste0(figures.fp, "/TempRelated"))) {dir.create(paste0(figures.fp, "/TempRelated"))}
    ggsave(plot_name, 
           plot = temp_plot, device = "png", dpi = 300,
                  height = 11, width = 12)
    return(temp_plot)
  }
})

#### ====================================================================== ####

# Methane Flux and Temperature
#### ====================================================================== ####
methanogenesis_temp.df <- sample_metadata %>% 
  #filter(Habitat__ %in% c("Bog")) %>%
  select(temporal_sample_id, Habitat__, contains("CH4_Flux_AVE"), 
         contains("mean_AirTemperature")) %>%
  pivot_longer(cols = contains("mean_AirTemperature"), 
               names_to = "Temp_metric",
               values_to = "mean_AirTemp") %>%
  mutate(Temp_metric = gsub("mean_AirTemperature", "", Temp_metric),
         Temp_metric = gsub("_", "", Temp_metric),
         Temp_metric = factor(Temp_metric, levels = c("Tair.degC", "samplingdate", "7d", "14d",
                                                                    "21d", "28d", "growing", "allgrowing"))) %>%
  pivot_longer(cols = contains("CH4_Flux_AVE"), 
               names_to = "CH4_Flux_metric",
               values_to = "CH4_Flux") %>%
  mutate(CH4_Flux_metric = gsub("CH4_Flux_AVE_", "", CH4_Flux_metric),
         CH4_Flux_metric = gsub("_before_coring", "", CH4_Flux_metric),
         CH4_Flux_metric = fct_rev(CH4_Flux_metric)) %>% 
  select(-temporal_sample_id) %>% distinct() %>% # removes misleading extra "observations" of CH4 flux that actually are just the same across samples
  group_by(Habitat__, Temp_metric, CH4_Flux_metric) %>%
  nest() %>%
  mutate(lm.model = map(data, ~lm(CH4_Flux~mean_AirTemp, data = .)),
         lm.tidied = map(lm.model, tidy),
         summaried = map(lm.model, summary), 
         lm.model.rsqr = map_dbl(lm.model, ~glance(.)$r.squared),
         lm.model.pval = map_dbl(lm.model, ~glance(.)$p.value)) %>%
  unnest(lm.tidied, names_sep = ".") %>%
  mutate(cor.model = map(data, ~cor.test(~CH4_Flux + mean_AirTemp, data = .)),
         cor.tidied = map(cor.model, tidy)) %>%
  unnest(cor.tidied, names_sep = ".") %>% 
  unnest(data) %>%
  mutate(sig = ifelse(cor.tidied.p.value < 0.05, TRUE, NA))
methanogenesis_temp.df


methanogenesis_temp <- methanogenesis_temp.df %>% 
  ggplot(aes(y = CH4_Flux, x = mean_AirTemp)) +
  geom_point(aes(color = Habitat__)) +
  geom_text(data = methanogenesis_temp.df %>%
              filter(sig == TRUE) %>%
              filter(lm.tidied.term == "mean_AirTemp") %>%
              select(CH4_Flux_metric, Temp_metric, Habitat__, lm.tidied.estimate) %>%
              mutate(lm.tidied.estimate = round(lm.tidied.estimate, digits = 2)) %>%
              distinct(),
            aes(label = paste0("slope = ", lm.tidied.estimate)), x = Inf, y = Inf, vjust = 1, hjust = 1) +
  geom_smooth(data = methanogenesis_temp.df %>%
                filter(sig == TRUE),
              aes(group = Habitat__, color = Habitat__), method = "lm") +
  scale_color_manual(name = "Habitat", values = colour_habitat, 
                     breaks = habitat_levels, 
                     labels = habitat_levels) +
  facet_grid(CH4_Flux_metric ~ Temp_metric, switch = "y", scales = "free_x") +
  theme_bw() +
  theme(strip.placement = "outside")
methanogenesis_temp

ggsave(paste0(figures.fp, "/methaneflux_temp.png"), plot = methanogenesis_temp, device = "png", dpi = 300,
       height = 5, width = 12)
#### ====================================================================== ####

# Ground Temperature responder VIPs
#### ====================================================================== ####

metaG_temp_vip_resp_models <- calc_temp_resp_mags(data = trimmed_mean$rel_abund, 
                    mag_list = temp_respVIP, 
                    dataType = "VIP",
                    scale_vars = T, skip_resids = T)


metaG_temp_vip_resp_models %>%
  filter(T_responder == T)
#### ====================================================================== #### 

# Compare slopes to each other
#### ====================================================================== #### 


VIP_temp_slope_est <- compare_and_plot_slopes(resp_model = metaG_temp_vip_resp_models,
                                           slope = "Temperature",
                                           plot_dir = here(figures.fp,"/", "VIP"),
                                           dataType = "VIP", plot_height = 20)
VIP_depth_slope_est <- compare_and_plot_slopes(resp_model = metaG_temp_vip_resp_models,
                                         slope = "Depth",
                                         plot_dir = here(figures.fp,"/", "VIP"),
                                         dataType = "VIP", plot_height = 20)

VIP_Tdepth_slope_est <- compare_and_plot_slopes(resp_model = metaG_temp_vip_resp_models,
                                               slope = "Interaction",
                                               plot_dir = here(figures.fp,"/", "VIP"),
                                               dataType = "VIP", plot_height = 20)


#### ====================================================================== #### 

# Temperature responder pathways 
#### ====================================================================== ####

metaG_temp_path_resp_models <- calc_temp_resp_pathways(data = pathway_abundance_cumu, 
                                                  dataType = "metaG", scale_vars = T)

# Calculate metaT cumulative abundance
pathway_expression_cumu <- metaT_exp %>%
    group_by(subpathway, temporal_sample_id) %>%
    summarise(abund = sum(tpm)) %>%
    pivot_wider(names_from = "subpathway", values_from = "abund")

metaT_temp_path_resp_models <- calc_temp_resp_pathways(data = pathway_expression_cumu, 
                        dataType = "metaT", scale_vars = T)
#### ====================================================================== ####

# Compare temperature slopes to each other
#### ====================================================================== ####
compare_and_plot_slopes <- function(resp_model, slope = "Temperature",
                                    plot_dir, dataType = "metaG") {
  
  # Figure out what slope to filter on
  if(slope == "Temperature") {
    slopeterm <- "SoilTemp"
    xlabel <- "Temperature"
    plottitle <- "Slope of response to ground temperature"
    fname <- paste0(dataType, "_temperature_slope_comparison.png")
    writeLines("Plotting temperature slope terms")
  } else if(slope == "Depth") {
    slopeterm <- "Depth"
    xlabel <- "Depth"
    plottitle <- "Slope of response to depth"
    fname <- paste0(dataType, "_depth_slope_comparison.png")
    writeLines("Plotting depth slope terms")
  } else if (slope == "Interaction") {
    slopeterm <- "SoilTemp:Depth"
    xlabel <- "Temperature and Depth Interaction"
    plottitle <- "Slope of temperature and depth interaction"
    fname <- paste0(dataType, "_interaction_slope_comparison.png")
    writeLines("Plotting interaction slope terms")
  } else {
    writeLines("Please correctly specify the slope term, it can be either: Temperature, Depth, or Interaction")
  }
  
  # Extract the slope term information
  temp_path_estimates <- resp_model %>% #filter(T_responder|TDepth_responder) %>%
    mutate(plottedcoverage = ifelse(lm.resid.shapiro < 0.05, "coverage_log10", "coverage")) %>%
    select(Habitat__, Pathways, lm.tidied, plottedcoverage, VIFOK, T_responder, T_padj, Depth_padj, TDepth_padj, lm.pval, lm.rsquared) %>%
    unnest(lm.tidied) %>% 
    filter(term != "(Intercept)") %>%
    # Make clearer term names
    mutate(term = case_when(
      term == "scale(T_soil.deg_C):scale(DepthAvg__)" ~ "SoilTemp:Depth",
      term == "scale(T_soil.deg_C)" ~ "SoilTemp",
      term == "scale(DepthAvg__)" ~ "Depth",
      TRUE ~ NA_character_),
      term = factor(term, levels = c("SoilTemp:Depth", "Depth", "SoilTemp"))) %>%
    # Identify when a term is significant
    mutate(significant = case_when(
      T_padj < 0.05 & term == "SoilTemp" ~ "sig",
      Depth_padj < 0.05 & term == "Depth" ~ "sig",
      TDepth_padj < 0.05 & term == "SoilTemp:Depth" ~ "sig",
      TRUE ~ "not sig"),
      # Specify VIFOK as "cannot determine" rather than significant/not significant
      significant = ifelse(VIFOK, significant, "cannot determine"),
      TempEst = ifelse(term == "SoilTemp", estimate, NA)) %>% 
    ungroup() %>%
    arrange(VIFOK, desc(T_responder), TempEst) %>% 
    mutate(order = ifelse(term == "SoilTemp" & T_padj < 0.05 & VIFOK, row_number(), NA)) %>%
    group_by(Habitat__, Pathways) %>% fill(order) %>% ungroup()
  
  # Plot slope estimates
  slope_est_plot <- temp_path_estimates %>%
    # Filter based on slope terms
    filter(term == slopeterm) %>%
    ggplot(aes(x = estimate, y = Habitat__)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high,
                      group = Pathways, linetype = significant, color = Habitat__), 
                  width = 0.3) +
    geom_point(aes(group = Pathways, shape = significant, 
                   color = Habitat__)) +
    scale_shape_manual("Significant", values = c(19, 1, 8), 
                       breaks = c("sig", "not sig", "cannot determine"), 
                       labels = c("Significant", "Not Significant", "Temperature and depth\ntoo correlated\nto determine")) +
    scale_linetype_manual("Significant", values = c("solid", "dashed", "dotted"), 
                          breaks = c("sig", "not sig", "cannot determine"), 
                          labels = c("Significant", "Not Significant", "Temperature and depth\ntoo correlated\nto determine")) +
    scale_color_manual("Habitat", values = colour_habitat, 
                       breaks = habitat_levels, labels = habitat_levels) +
    facet_wrap(.~Pathways, ncol = 1, strip.position = "left") +
    xlab(paste0(xlabel, " slope estimate (+/- 95% confidence interval)")) + ylab("") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0),
          panel.spacing.y = unit(0, "lines")) +
    ggtitle(plottitle)
  
  # Save plot
  ggsave(path = plot_dir, filename = fname, plot = slope_est_plot,
         device = "png", height = 7, width = 14)
  
  return(temp_path_estimates)
}

metaG_slope_est <- compare_and_plot_slopes(resp_model = metaG_temp_path_resp_models,
                                           slope = "Temperature",
                                           plot_dir = here(figures.fp,"/", "metaG"),
                                           dataType = "metaG")

metaG_slope_est <- compare_and_plot_slopes(resp_model = metaG_temp_path_resp_models,
                                           slope = "Depth",
                                           plot_dir = here(figures.fp,"/", "metaG"),
                                           dataType = "metaG")

metaT_slope_est <- compare_and_plot_slopes(resp_model = metaT_temp_path_resp_models,
                                           slope = "Temperature",
                                           plot_dir = here(figures.fp,"/", "metaT"),
                                           dataType = "metaT")
metaT_slope_est <- compare_and_plot_slopes(resp_model = metaT_temp_path_resp_models,
                                           slope = "Depth",
                                           plot_dir = here(figures.fp,"/", "metaT"),
                                           dataType = "metaT")
