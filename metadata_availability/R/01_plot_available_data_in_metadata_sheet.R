#' ## 01 Plot available data in metadata sheet
#' This script is meant to go through the metadata sheet and generate plots of the data in that sheet.

#'
#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(vegan); packageVersion("vegan") # for ecological applications
library(viridis)
library(cowplot)
library(ggh4x)
library(here)

#+ clean environment, eval = FALSE
# clean environment to start
rm(list = ls())

# Load required data
source(here("setup.R"))

#### User Run Options
#### ====================================================================== ####
save_data_plots = TRUE # should the data completeness plots be created and saved?
#### ====================================================================== ####
# Plot aesthetics
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
#### ====================================================================== ####

#### Create a plotting function
#### ====================================================================== ####
# Write a function that plots the presence or absence of environmental data 
# variable on top of the data for which we have reads mapped to MAGS
plot_data_completeness <- function(env_data_table = env_data, 
                                   env_var, outdir = NULL) {
  # Record name of env variable
  var_name <- env_var
  # Create subsetted data set
  subset_data <- env_data %>%
    filter_at(vars(all_of(env_var)), all_vars(!is.na(.)))
  #print(subset_data[, env_var], n = 253)
  
  # create baseplot 
  baseplot <- ggplot(data = env_data_table, 
                     aes(x = Core__, y = DepthAvg__)) +
    geom_point(color = "grey40",
               fill = NA) +
    geom_errorbar(aes(ymin = DepthMin__, ymax = DepthMax__),
                  width = 0.05) +
    facet_nested(~Year__ + Habitat__)
  
  var_plot <- baseplot +
    geom_point(data = subset_data, 
               aes(x = Core__, y = DepthAvg__),
               color = "red") +
    scale_y_reverse() +
    ggtitle(paste0("Data availability for ", var_name, " : red = some value is",
    " present, black = value is NA"))
  
  # Save or printout plots; if outdir is NULL, plots are returned with function,
  # otherwise they are saved to the specified location
  if(is.null(outdir)) {
    return(var_plot)
  } else {
    # strip filepath name of odd characters
    var_name <- str_replace_all(var_name, "[\\s[:punct:]]", "\\_")
    fp <- paste0(outdir, "/", var_name, ".png")
    ggsave(plot = var_plot, filename = fp, dpi = 300, height = 8, width = 12)
    writeLines(paste0("Plot for ", var_name, " can be found here: ", fp))
    return(var_plot)
  }
  
}

#### Plot the timeline of samples and env data
#### ====================================================================== ####
# Extract environmental data from objects produced by setup.R
env_data <- input_ra$sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels),
         Core__ = as.factor(Core__))

# Create a base plot
baseplot <- ggplot(data = env_data, 
                   aes(x = Core__, y = DepthAvg__)) +
  geom_point(color = "grey40", fill = NA) +
  geom_errorbar(aes(ymin = DepthMin__, 
                    ymax = DepthMax__), width = 0.05) +
  facet_nested(~Year__ + Habitat__) +
  scale_y_reverse()
baseplot

# Print data availability plots to file
if(save_data_plots == TRUE) {
  # Set up folder for plots
  data_avail_plots.fp <- here("metadata_availability", "figures", 
                              "data_completeness_plots")
  if (!dir.exists(data_avail_plots.fp)) {dir.create(data_avail_plots.fp, 
                                                    recursive = TRUE)}
  
  # remove any plots already in the folder
  existing_files <- list.files(path = data_avail_plots.fp, pattern = "*.png", 
                               full.names = T)
  file.remove(existing_files)

  # Generate plots
  completeness_plots <- lapply(names(env_data), plot_data_completeness, 
                               env_data_table = env_data, 
                               outdir = data_avail_plots.fp)
}

#### ====================================================================== ####
# Plot samples and depth ranges over time
#### ====================================================================== ####
depth_lumping_rect <- data.frame(
  ymin = rep(c(0,10,20,30), times =3*3),
  ymax = rep(c(9,19,29,39), times =3*3),
  name = rep(depth_levels, times = 3*3),
  Core__ = as.factor(rep(c(1,2,3), each = 4)),
  Habitat__ = rep(c("Palsa", "Bog", "Fen"), each = 12),
  DepthAvg__ = rep(c(5,15,25,35), times =3*3)
) %>%
  left_join(distinct(sample_metadata %>% select(Habitat__, Year__)), 
            by = "Habitat__") %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels),
         xmin = as.numeric(Core__) -0.5,
         xmax = as.numeric(Core__) + 0.5)
baseplot <- ggplot(data = sample_metadata %>% 
                     mutate(Habitat__ = factor(Habitat__, levels = habitat_levels),
                            WTD_negated = -1 * WTD ), 
                   aes(x = Core__, y = DepthAvg__)) +
  geom_rect(data = depth_lumping_rect, aes(xmin = xmin, 
                                           ymin = ymin, 
                                           xmax = xmax, 
                                           ymax = ymax, 
                                           fill = name), 
            alpha = 0.5) +
  scale_fill_manual(name = "Depth Range", values = colour_depth) +
  # Water table
  geom_segment(aes(y = WTD_negated, yend = WTD_negated, x = Core__ - 0.5, xend = Core__ + 0.5), color = "blue") +
  # PF distance
  geom_segment(aes(y = ALD, yend = ALD, 
                   x = Core__ - 0.5, xend = Core__ + 0.5), color = "black") +
  geom_text(data = sample_metadata %>% 
              mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
              filter(ALD_detectable == "No"),
            label = "x", aes(y = ALD, x = Core__), color = "black") +
  geom_errorbar(aes(ymin = DepthMin__, 
                    ymax = DepthMax__), width = 0.05) +
  geom_point(color = "black", fill = "black") +
  facet_nested(~Year__ + Habitat__) +
  scale_y_reverse() +
  xlab("Core Number") + ylab("Depth") +
  theme_bw()
baseplot

#### ====================================================================== ####
# Write out the version of the metadata file that is used to generate the plots
outputs.fp <- here("metadata_availability", "outputs")
if (!dir.exists(outputs.fp)) {dir.create(outputs.fp, recursive = TRUE)}
write_csv(env_data, file = paste0(outputs.fp, "/metadata_used_for_plots.csv"))
