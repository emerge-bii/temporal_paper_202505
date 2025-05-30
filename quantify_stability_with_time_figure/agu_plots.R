# AGU presentation nmds with time

#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(here)
library(ggnewscale)

# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("quantify_stability_with_time_figure", "outputs")
figures.fp <- here("quantify_stability_with_time_figure", "figures")

output_dir <- figures.fp
if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

# 
# library(ggnewscale)
# library(ggh4x)
# library(propr)
# library(cowplot)
# library(vegan)
# library(ape)
# library(tidyverse)
# library(here)
# source(here("setup.R"))
# 
# map <- purrr::map
# 
# main_dir <- here("Metabolic-analysis", "00_MAG_overview")
# dir.create(main_dir, recursive = TRUE)
#### ====================================================================== ####
# Set colors
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
shape_habitat  <- c(22, 21, 24)
fill_habitat   <- colour_habitat

depth_levels <- c("0-9", "10-19", "20-29", "30-39")
depth_labels <- c("0", "10", "20", "30")
fill_depth   <- RColorBrewer::brewer.pal(5, "YlOrBr")[-1]
colour_depth <- fill_depth

year_levels <- c(2011, 2012, 2013, 2014, 2015, 2016, 2017)
year_labels <- c(1, 2, 3, 4, 5, 6, 7)
fill_year   <- RColorBrewer::brewer.pal(7, "BuPu")
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
  "other"
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
  "other"
)
colour_phylum <- c(RColorBrewer::brewer.pal(12, "Set3"), "white")
fill_phylum <- colour_phylum

set.seed(42)

# euclidian distance for arrow calculations
euc.dist <- function(x1, x2, y1, y2){ 
  arrow_dist <- sqrt((x2-x1)^2 + (y2 - y1)^2)
  return(arrow_dist)
  }

#### ====================================================================== ####
#### Read in data for each part of figure
#### ====================================================================== ####

rel_abund <- trimmed_mean$rel_abund %>%
  select(genome, all_of(sample_metadata$temporal_sample_id)) %>%
  pivot_longer(-genome, names_to = "temporal_sample_id", values_to = "rel_abund")


## Ordination figures

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


#############
### Stats ###
#############
abund_stats <- adonis2(dis_matrix ~ Habitat__ * DepthLumping * Year__, data = sample_metadata, by = "terms")

# NMDS Plot Function
plot_nmds_habitat <- function(dis_matrix, habitat, legible = TRUE) {
  abund_nmds_raw <- metaMDS(dis_matrix)
  stress <- abund_nmds_raw$stress %>%
    round(2)
  abund_nmds <- as_tibble(abund_nmds_raw$points, rownames = "sample") %>%
    left_join(sample_metadata, by = c("sample" = "temporal_sample_id")) %>%
    filter(!is.na(Habitat__))
  
  ################
  ### Plotting ###
  ################
  nmds_plotting_layers <- list(
    stat_ellipse(type = "t", geom = "polygon", alpha = 1/10, colour = NA),
    stat_ellipse(type = "norm", linetype = 2, colour = "black"),
    scale_shape_manual(name = "Habitat", values = shape_habitat, breaks = habitat_levels),
    scale_colour_manual(name = "Habitat", values = colour_habitat, breaks = habitat_levels),
    scale_fill_manual(name = "Habitat", values = fill_habitat, breaks = habitat_levels),
    annotation_custom(grid::textGrob(label = paste("stress:", stress), x = unit(0.8, "npc"), y = unit(0.1, "npc"), hjust = 0)),
    theme_bw(), 
    theme(legend.text = element_text(size = rel(5)), 
          text = element_text(size = rel(5)))
  )
  
  # Point settings
  point_size <- rel(1)
  point_alpha <- 1
  
  if(legible) {
    nmds_plotting_layers <- list(
      nmds_plotting_layers,
      theme(legend.text = element_text(size = rel(5)), 
            text = element_text(size = rel(5)))
    )
    # Point settings
    point_size <- rel(5)
    point_alpha <- 0.7
  }
  
  
  abund_nmds %>%
    ggplot(aes(MDS1, MDS2, shape = Habitat__, colour = Habitat__, fill = Habitat__)) +
    nmds_plotting_layers +
    geom_point()
  ggsave(str_c("nmds", habitat, "12.png", sep = "_"), 
         height = 7, width = 7,
         path = output_dir, dpi = 900)
  
  abund_nmds %>%
    ggplot(aes(MDS1, MDS2, shape = Habitat__, fill = Habitat__)) +
    nmds_plotting_layers +
    new_scale_fill() +
    geom_point(aes(shape = Habitat__, colour = NA, fill = DepthLumping), size = point_size, alpha = point_alpha) +
    scale_fill_manual(name = "Depth (cm)", values = colour_depth, breaks = depth_levels) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  ggsave(str_c("nmds", habitat, "depth.png", sep = "_"), 
         height = 7, width = 8.5,
         path = output_dir, dpi = 900)
  
  abund_nmds %>%
    ggplot(aes(MDS1, MDS2, shape = Habitat__, fill = Habitat__)) +
    nmds_plotting_layers +
    new_scale_fill() +
    geom_point(aes(shape = Habitat__, colour = NA, fill = factor(Year__)), size = point_size, alpha = point_alpha) +
    scale_fill_manual(name = "Year", values = colour_year, breaks = year_levels) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  ggsave(str_c("nmds", habitat, "year.png", sep = "_"), 
         height = 7, width = 8,
         path = output_dir, dpi = 900)
}

plot_nmds_habitat(dis_matrix, "")
plot_nmds_habitat(dis_matrix_palsa, "palsa")
plot_nmds_habitat(dis_matrix_bog, "bog")
plot_nmds_habitat(dis_matrix_fen, "fen")


# PCOA Plot Function
plot_pcoa_habitat <- function(dis_matrix, habitat, legible = TRUE) {
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
}

plot_pcoa_habitat(dis_matrix, "")
plot_pcoa_habitat(dis_matrix_palsa, "palsa")
plot_pcoa_habitat(dis_matrix_bog, "bog")
plot_pcoa_habitat(dis_matrix_fen, "fen")






ggplot(iris, aes(x = Petal.Width, y = Sepal.Length)) + geom_point(aes(color = Species))


iris.gam <- gam(Sepal.Length ~ Petal.Width + s(Species, k = 3, bs = "re"), data = iris)
summary(iris.gam)
gam.check(iris.gam)

iris.lmer <- lmer(Sepal.Length ~ Petal.Width + (1|Species), data = iris)
summary(iris.lmer)



# Ordisurf

layer_temperature_ordisurf <- function(abund_pcoa) {
  ordi_temp <- ordisurf(abund_pcoa[,c(2,3)] ~ abund_pcoa$T_soil.deg_C)
  ordi_depth <- ordisurf(abund_pcoa[,c(2,3)] ~ abund_pcoa$DepthAvg__)
  ordi_Habitat <- ordisurf(abund_pcoa[,c(2,3)] ~ abund_pcoa$Habitat__)
  broom::glance(ordi_temp)
  broom::glance(ordi_depth)
  
  abund_pcoa <- mutate(abund_pcoa, Habitat__ = factor(Habitat__, levels = habitat_levels),
                       DepthLumping = factor(DepthLumping, levels = depth_levels))
  temp_depth_gam <- gam(T_soil.deg_C ~ s(Axis.1, Axis.2, by = Habitat__) + DepthLumping + Habitat__,
                  data = abund_pcoa, method = "REML")
  gam.check(temp_depth_gam)
  summary(temp_depth_gam)
  
  
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
  
  time_gam <- gam(Year__ ~ 1 + Habitat__*DepthLumping + s(Axis.1, Axis.2),
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

#### ====================================================================== ####
# Pathway nmdss
#### ====================================================================== ####

# Calculate function distance matrix
as_matrix <- function(x) {
  output <- as.matrix.data.frame(x[, -1])
  rownames(output) <- x[[1]]
  output
}

get_dis_matrix <- function(pathway_abundance_cumu) {
  dis_matrix <- pathway_abundance_cumu %>%
    as_matrix() %>%
    vegdist(method = "bray")
  
  return(dis_matrix)
}

dis_matrix_path <- get_dis_matrix(pathway_abundance_cumu)

plot_nmds_habitat(dis_matrix_path, "path")
plot_nmds_habitat(dis_matrix_path_palsa, "path_palsa")
plot_nmds_habitat(dis_matrix_path_bog, "path_bog")
plot_nmds_habitat(dis_matrix_path_fen, "path_fen")

#### ====================================================================== ####
# Plot arrows over time NMDS
#### ====================================================================== ####
# Define plotting function
plot_nmds_change_by_year <- function(dis_matrix, habitat, type, max_years = 2017, plot_points = T) {
  # type = "path" or "Mags" will name the output file accordingly
  abund_nmds_raw <- metaMDS(dis_matrix)
  stress <- abund_nmds_raw$stress %>%
    round(2)
  abund_nmds <- as_tibble(abund_nmds_raw$points, rownames = "sample") %>%
    left_join(sample_metadata, by = c("sample" = "temporal_sample_id")) %>%
    filter(!is.na(Habitat__))
  
  ################
  ### Plotting ###
  ################
  nmds_plotting_layers <- list(
    stat_ellipse(type = "t", geom = "polygon", alpha = 1/10, colour = NA),
    stat_ellipse(type = "norm", linetype = 2, colour = "black"),
    scale_shape_manual(name = "Habitat", values = shape_habitat, breaks = habitat_levels),
    scale_colour_manual(name = "Habitat", values = colour_habitat, breaks = habitat_levels),
    scale_fill_manual(name = "Habitat", values = fill_habitat, breaks = habitat_levels),
    annotation_custom(grid::textGrob(label = paste("stress:", stress), x = unit(0.8, "npc"), y = unit(0.1, "npc"), hjust = 0)),
    theme_bw()
  )
  
  # Filter for max years to plot to create annimation
  abund_nmds_filt <- abund_nmds %>%
    filter(Year__ <= max_years)
  
  # create arrows from average core locations each year (grouped by depth)
  abund_nmds_avg <- abund_nmds_filt %>%
    select(sample, MDS1, MDS2, Core__, Year__, Habitat__, DepthLumping) %>%
    mutate(Year__ = Year__ - 1) %>% # this way 2012 average value will be aligned with 2011 
    group_by(Year__, Habitat__, DepthLumping) %>%
    add_count(name = "grp_size") %>%
    summarize(AvgNextMDS1 = ifelse(all(grp_size) > 1, mean(MDS1, na.rm = T), MDS1),
              AvgNextMDS2 = ifelse(all(grp_size) > 1, mean(MDS2, na.rm = T), MDS2))
  
  segment_abund_nmds <- abund_nmds_filt %>%
    select(sample, MDS1, MDS2, Core__, Year__, Habitat__, DepthLumping) %>%
    group_by(Year__, Habitat__, DepthLumping) %>% 
    add_count(name = "grp_size") %>%
    summarize(AvgFirstMDS1 = ifelse(all(grp_size) > 1, mean(MDS1, na.rm = T), MDS1),
              AvgFirstMDS2 = ifelse(all(grp_size) > 1, mean(MDS2, na.rm = T), MDS2)) %>%
    left_join(abund_nmds_avg)
  
  
  # Fix NAs that are caused by only one observation by using just that observation
  
  
  hab_year_plot <- abund_nmds %>%
    select(sample, MDS1, MDS2, Year__, Habitat__, DepthLumping) %>%
    ggplot(aes(MDS1, MDS2, shape = Habitat__, fill = Habitat__)) +
    nmds_plotting_layers +
    new_scale_fill()
  
  if(plot_points) {
    hab_year_plot <- hab_year_plot + 
      geom_point(data = abund_nmds,
                 aes(shape = Habitat__, colour = NA, fill = factor(Year__)),
                 alpha = 1)
  } else {
    hab_year_plot <- hab_year_plot + 
      geom_point(data = abund_nmds,
                 aes(shape = Habitat__, colour = NA, fill = factor(Year__)),
                 alpha = 0)
  }
  
  # Set arrow color scale to purples for year
  colour_year <- rev(RColorBrewer::brewer.pal(n = 7, "Purples"))
  
  hab_year_plot <- hab_year_plot + 
    geom_segment(data = segment_abund_nmds, 
                 aes(x = AvgFirstMDS1, y = AvgFirstMDS2,
                     xend = AvgNextMDS1, yend = AvgNextMDS2,
                     color = factor(Year__)),
                 #color = "grey30",
                 lineend = "round", linejoin = "bevel",
                 linewidth = rel(1),
                 arrow = arrow(length = unit(0.02, "npc"))) +
    scale_fill_manual(name = "Year", values = colour_year, breaks = year_levels) +
    scale_color_manual(name = "Year", values = colour_year, breaks = year_levels) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  ggsave(str_c("nmds", type, habitat, max_years, "year.png", sep = "_"), 
         height = 7, width = 7,
         path = output_dir, dpi = 900)
  return(hab_year_plot)
}

# Plot pathways
for(i in 2011:2017) {
  plot_nmds_change_by_year(dis_matrix = dis_matrix, "", type = "path", max_years = i)
  plot_nmds_change_by_year(dis_matrix, "", type = "path_no_points", 
                           max_years = i, plot_points = F)
}



# Calculate MAGs distance matrix
get_dis_matrix <- function(rel_abund) {
  dis_matrix <- rel_abund %>%
    t() %>%
    vegdist(method = "bray")
  
  return(dis_matrix)
}

dis_matrix <- trimmed_mean$rel_abund %>%
  select(all_of(sample_metadata$temporal_sample_id)) %>%
  get_dis_matrix()

# Plot MAGs
for(i in 2011:2017) {
  plot_nmds_change_by_year(dis_matrix, "", type = "MAGs", max_years = i)
  plot_nmds_change_by_year(dis_matrix, "", type = "MAGs_no_points", 
                           max_years = i, plot_points = F)
}


#### ====================================================================== ####
#### Arrow lengths with time
#### ====================================================================== ####
# Define plotting function
plot_nmds_arrow_with_time <- function(dis_matrix = dis_matrix_palsa, habitat = "Palsa", 
                                      type = "Mags", max_years = 2017, plot_points = T) {
  # type = "path" or "Mags" will name the output file accordingly
  abund_nmds_raw <- metaMDS(dis_matrix)
  stress <- abund_nmds_raw$stress %>%
    round(2)
  abund_nmds <- as_tibble(abund_nmds_raw$points, rownames = "sample") %>%
    left_join(sample_metadata, by = c("sample" = "temporal_sample_id")) %>%
    filter(!is.na(Habitat__))
  
  ################
  ### Plotting ###
  ################
  nmds_plotting_layers <- list(
    stat_ellipse(type = "t", geom = "polygon", alpha = 1/10, colour = NA),
    stat_ellipse(type = "norm", linetype = 2, colour = "black"),
    scale_shape_manual(name = "Habitat", values = shape_habitat, breaks = habitat_levels),
    scale_colour_manual(name = "Habitat", values = colour_habitat, breaks = habitat_levels),
    scale_fill_manual(name = "Habitat", values = fill_habitat, breaks = habitat_levels),
    annotation_custom(grid::textGrob(label = paste("stress:", stress), x = unit(0.8, "npc"), y = unit(0.1, "npc"), hjust = 0)),
    theme_bw()
  )
  
  # Filter for max years to plot to create annimation
  abund_nmds_filt <- abund_nmds %>%
    filter(Year__ <= max_years)
  
  # create arrows from average core locations each year (grouped by depth)
  abund_nmds_avg <- abund_nmds_filt %>%
    select(sample, MDS1, MDS2, Core__, Year__, Habitat__, DepthLumping) %>%
    mutate(Year__ = Year__ - 1) %>% # this way 2012 average value will be aligned with 2011 
    group_by(Year__, Habitat__, DepthLumping) %>%
    add_count(name = "grp_size") %>%
    summarize(AvgNextMDS1 = ifelse(all(grp_size) > 1, mean(MDS1, na.rm = T), MDS1),
              AvgNextMDS2 = ifelse(all(grp_size) > 1, mean(MDS2, na.rm = T), MDS2))
  
  segment_abund_nmds <- abund_nmds_filt %>%
    select(sample, MDS1, MDS2, Core__, Year__, Habitat__, DepthLumping) %>%
    group_by(Year__, Habitat__, DepthLumping) %>% 
    add_count(name = "grp_size") %>%
    summarize(AvgFirstMDS1 = ifelse(all(grp_size) > 1, mean(MDS1, na.rm = T), MDS1),
              AvgFirstMDS2 = ifelse(all(grp_size) > 1, mean(MDS2, na.rm = T), MDS2)) %>%
    left_join(abund_nmds_avg) %>%
    # calculate euclidean distance between arrows
    mutate(arrow_distance = euc.dist(AvgFirstMDS1, AvgNextMDS1, AvgFirstMDS2, AvgNextMDS2)) %>%
    mutate(Habitat__ = factor(Habitat__, levels = habitat_levels))
  
  p <- ggplot(data = segment_abund_nmds, aes(x = paste0(Year__, " -\n ", Year__+ 1), y = arrow_distance)) +
    geom_point(aes(color = Habitat__), size = 3) + 
    geom_line(aes(group = Habitat__, color = Habitat__)) +
    scale_colour_manual(name = "Habitat", values = colour_habitat, breaks = habitat_levels) +
    facet_grid(DepthLumping~Habitat__, switch = "y") +
    xlab("Year Interval") + ylab("Arrow Length\n(Distance from center point of replicates at every depth each year)") +
    theme_bw() +
    theme(strip.placement = "outside")
  
  p
  ggsave("NMDS_arrows_over_time.png", plot = p, 
         height = 7, width = 12,
         path = output_dir, dpi = 300)
  
  # Fix NAs that are caused by only one observation by using just that observation
  
  
  hab_year_plot <- abund_nmds %>%
    select(sample, MDS1, MDS2, Year__, Habitat__, DepthLumping) %>%
    ggplot(aes(MDS1, MDS2, shape = Habitat__, fill = Habitat__)) +
    nmds_plotting_layers +
    new_scale_fill()
  
  if(plot_points) {
    hab_year_plot <- hab_year_plot + 
      geom_point(data = abund_nmds,
                 aes(shape = Habitat__, colour = NA, fill = factor(Year__)),
                 alpha = 1)
  } else {
    hab_year_plot <- hab_year_plot + 
      geom_point(data = abund_nmds,
                 aes(shape = Habitat__, colour = NA, fill = factor(Year__)),
                 alpha = 0)
  }
  
  # Set arrow color scale to purples for year
  colour_year <- rev(RColorBrewer::brewer.pal(n = 7, "Purples"))
  
  hab_year_plot <- hab_year_plot + 
    geom_segment(data = segment_abund_nmds, 
                 aes(x = AvgFirstMDS1, y = AvgFirstMDS2,
                     xend = AvgNextMDS1, yend = AvgNextMDS2,
                     color = factor(Year__)),
                 #color = "grey30",
                 lineend = "round", linejoin = "bevel",
                 linewidth = rel(1),
                 arrow = arrow(length = unit(0.02, "npc"))) +
    scale_fill_manual(name = "Year", values = colour_year, breaks = year_levels) +
    scale_color_manual(name = "Year", values = colour_year, breaks = year_levels) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))

  return(hab_year_plot)
}




#### ====================================================================== ####
#### Beta Diversity With Time
#### ====================================================================== ####

# Define plotting function
dis_matrix <- as.matrix(dis_matrix)
xy <- t(combn(colnames(dis_matrix), 2))
dis_mat_long <- data.frame(xy, dist = dis_matrix[xy])

short_sample_metadata <- sample_metadata %>%
  select(temporal_sample_id, Habitat__, Year__, DepthLumping)

dis_matrix_info <- left_join(dis_mat_long, short_sample_metadata, by = c("X1" = "temporal_sample_id")) %>%
  rename(Sample_A = X1) %>%
  left_join(short_sample_metadata, by = c("X2" = "temporal_sample_id"), suffix = c("_A", "_B")) %>%
  rename(Sample_B = X2)
  

dis_mat_2011_comp <- dis_matrix_info %>%
  filter(if_any(c(Year___A, Year___B), ~ . == "2011")) %>%
  filter(Habitat___A == Habitat___B) %>%
  filter(DepthLumping_A == DepthLumping_B)
  

dis_mat_2011_comp_avg <- dis_mat_2011_comp %>%
  group_by(Habitat___A, DepthLumping_A, Year___B) %>%
  summarize(Avg_BC = mean(dist, na.rm = T),
            sd_BC = sd(dist, na.rm = T))


# Temp_Time_depth.plot.df <- sample_metadata %>%
#   mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
#   dplyr::select(Habitat__, Year__, DepthAvg__, DepthLumping, T_soil.deg_C) %>%
#   filter(!is.na(T_soil.deg_C)) %>% # removes misleading extra "observations" of CH4 flux that actually are just the same across samples
#   dplyr::select(Habitat__, Year__, DepthLumping, T_soil.deg_C) %>%
#   group_by(Habitat__, DepthLumping) %>%
#   nest() %>%
#   mutate(lm.model = purrr::map(data, ~lm(T_soil.deg_C ~ Year__, data = .)),
#          lm.tidied = purrr::map(lm.model, tidy),
#          lm.model.rsqr = map_dbl(lm.model, ~glance(.)$adj.r.squared),
#          lm.model.pval = map_dbl(lm.model, ~glance(.)$p.value),
#          lm.model.yr_est = map_dbl(lm.tidied, ~.$estimate[2]),
#          lm.model.yr_int = map_dbl(lm.tidied, ~.$estimate[1])) %>%
#   mutate(Sig = ifelse(lm.model.pval < 0.05, T, F),
#          plotlab = paste0("\n r2 = ", round(lm.model.rsqr, 2),
#                           ", p = ",round(lm.model.pval, 2),
#                           "\n slope = ", round(lm.model.yr_est, 2))) %>%
#   unnest(data)

# Horizontal Plot;
BC_Time_depth.plot <- dis_mat_2011_comp %>%
  ggplot(aes(x = Year___B, y = dist)) +
  geom_point(aes(shape = DepthLumping_A), size = rel(5)) +
  # geom_segment(data = Temp_Time_depth.plot.df %>% 
  #                mutate(MinYear = min(Year__), MaxYear = max(Year__)) %>%
  #                select(-Year__, -T_soil.deg_C) %>% distinct() %>%
  #                mutate(lm.model.yr_est = ifelse(Sig == T, lm.model.yr_est, NA),
  #                       lm.model.yr_int = ifelse(Sig == T, lm.model.yr_int, NA)),
  #              aes(x = MinYear, xend = MaxYear,
  #                  y = lm.model.yr_int + lm.model.yr_est*MinYear,
  #                  yend = lm.model.yr_int + lm.model.yr_est*MaxYear,
  #                  group = DepthLumping, linetype = DepthLumping,
  #                  color = Habitat__), linewidth = rel(2)) +
  # geom_smooth(data = Temp_Time_depth.plot.df %>% filter(Sig == T),
  #             method = "lm", aes(group = DepthLumping, fill = Habitat__),
  #             alpha = 0.2, color = NA, linetype = "blank") +
  # remove model output text
  # geom_text(data = Temp_Time_depth.plot.df %>% filter(Sig == T) %>% 
  #             dplyr::select(DepthLumping, Habitat__, plotlab) %>% 
  #             distinct(),
  #           aes(label = plotlab), size = rel(2.5),
  #           x = -Inf, y = Inf, vjust = 0.8, hjust = 0) +
  geom_ribbon(data = dis_mat_2011_comp_avg,
            aes(y = Avg_BC, ymin = Avg_BC - sd_BC, ymax = Avg_BC + sd_BC,
            group = DepthLumping_A, fill = Habitat___A),
            alpha = 0.2, color = NA, linetype = "blank") +
  geom_line(data = dis_mat_2011_comp_avg, aes(y = Avg_BC, group = DepthLumping_A,
                                              linetype = DepthLumping_A, color = Habitat___A)) +
  facet_grid(~Habitat___A) +
  scale_x_continuous(limits = c(2011, 2017), breaks = seq(from = 2012, to = 2017, by = 2)) +
  coord_fixed(ratio = 6.9999) +
  scale_shape_manual(name = "Depth", breaks = c("0-9", "10-19", "20-29", "30-39"), 
                     values = c(1,2,3,4)) +
  scale_linetype_manual(name = "Depth", breaks = c("0-9", "10-19", "20-29", "30-39"),
                        values = c(1,2,3,4)) +
  scale_color_manual(name = "Habitat", breaks = habitat_levels, values = colour_habitat) +
  scale_fill_manual(name = "Habitat", breaks = habitat_levels, values = colour_habitat) +
  guides(linetype = guide_legend(keywidth = unit(5, "cm"), override.aes = list(fill = NA, color = "black", linewidth = rel(1)))) +
  theme_bw() +
  theme(axis.title = element_text(size = rel(2)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        strip.text = element_blank(),
        legend.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(2)),
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = NA)) +
  xlab("Year") + ylab("Bray-Curtis Distance\nfrom 2011")

BC_Time_depth.plot
ggsave(paste0(figures.fp, "BC_Dist_2011_horizontal.png"), plot = BC_Time_depth.plot, device = "png", dpi = 300,
       height = 5, width = 14)


