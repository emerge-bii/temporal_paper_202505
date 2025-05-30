#### ====================================================================== ####
#plot figure V2b (Gin's idea, Sarah mods)
#### ====================================================================== ####

# figure needs to have C abilities on the y axis and N abilities on the y
# each VIP will be represented by a small dot jittered so as not to overlap showing the metabolic potential of a vip for some VIP list
# goal is cluster of points representing all the vips with a function and their C/N overlap in abilities

library(ggdist)
library(patchwork)
library(purrr)
library(ggpubr)

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("CN_versatility", "outputs")
figures.fp <- here("CN_versatility", "figures")

output_dir <- figures.fp
if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

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

special_levels <- c("generalist", "methanogen", "fermenter", "macromolecule_degrader", "methanotroph", "homoacetogen", "monomer_degrader", NA)
special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-8]
special_fill <- special_colour

n_special_levels <- c("nitrite_oxidiser", "nitrogen_fixer", "nitrate_reducer", "denitrifier", NA)
n_special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-5:-8]
n_special_fill <- n_special_colour


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

#### ====================================================================== ####
#### Read in data
#### ====================================================================== ####
calculate_cn_versatility <- function(product_refined, genome_clusters) {
  carbon_metabolism <- function(pathways) {
    case_when(
      grepl("cazymes", pathways) ~ "carbon",
      grepl("fermentation", pathways) ~ "carbon",
      grepl("C fixation", pathways) ~ "carbon",
      grepl("methanotrophy", pathways) ~ "carbon",
      grepl("methanogenesis", pathways) ~ "carbon",
      grepl("degradation", pathways) ~ "carbon",
      TRUE ~ NA_character_
    )
  }
  
  nitrogen_metabolism <- function(pathways) {
    case_when(
      grepl("nitrogen", pathways) ~ "nitrogen",
      grepl("reductive_glycine", pathways) ~ "nitrogen",
      pathways == "CAZy-Chitin" ~ "nitrogen",
      pathways == "urea_degradation-all" ~ "nitrogen",
      TRUE ~ NA_character_
    )
  }
  
  
  CNVers <- metapathway_groups %>%
    mutate(CVersatility = carbon_metabolism(metapath),
           NVersatility = nitrogen_metabolism(subpathway))
  
  
  pathways <- product_refined %>%
    filter(genome %in% genome_clusters$representative) %>%
    left_join(CNVers, by = "subpathway") %>%
    filter(call)
  
  # Create CN paths for type 1 figure
  cpaths <- pathways %>% filter(CVersatility == "carbon") %>% select(-NVersatility) %>% rename(Versatility = CVersatility)
  npaths <- pathways %>% filter(NVersatility == "nitrogen") %>% select(-CVersatility) %>% rename(Versatility = NVersatility)
  
  cn_paths <- bind_rows(cpaths, npaths)
  
  # Create CN paths for type 2 figure
  cpaths <- pathways %>% filter(CVersatility == "carbon") %>% select(-NVersatility) %>% mutate(CVersatility = subpathway)
  npaths <- pathways %>% filter(NVersatility == "nitrogen") %>% select(-CVersatility) %>% mutate(NVersatility = subpathway)
  
  cn_paths_wide <- left_join(cpaths, npaths, by = "genome", suffix = c(".C", ".N")) %>%
    mutate(NVersatility = case_when(is.na(call.N) ~ "None Detected",
                                    TRUE ~ NVersatility))
  
  return(list(cn_paths_long = cn_paths, cn_paths_wide = cn_paths_wide))
}


pathways_CN <- calculate_cn_versatility(product_refined, genome_clusters)
#### ====================================================================== ####
#### ====================================================================== ####
# Normal faceting won't work because ggdist does strange things to the sizes of the dots across panels,
# so we make the panels individually, displaying x-axis labels and a fill legend only on one,
# and then we bodge them together with patchwork

# Function written with an eye to the alphaC version; not tested with other VIP variables
plot_vip_specials_typetwo_b_panel <- function(metapath, data) {
    metapath_df <- filter(data, metapath.C == metapath)
    xscale <- c("CAZy-Chitin", "reductive_glycine-serine_pathway",
                "reductive_glycine-acetyl_pathway","urea_degradation-all", "None Detected")
    
    p <- ggplot(data = metapath_df, aes(y = CVersatility, x = factor(NVersatility, levels = xscale))) +
      geom_dots(aes(fill = specialisation, group = NA), size = 1.2,
                binwidth = unit(0.015, "npc"),
                slab_color = "black", slab_linewidth = 0.3, 
                smooth = smooth_discrete(kernel = "ep", width = 0.8),
                side = "both", orientation = "horizontal") +
      scale_fill_manual(name = "Specialisation", breaks = special_levels, values = special_fill) +
      scale_x_discrete(limits = xscale) +
      scale_y_discrete(position = "right") +
      ggtitle(metapath) + 
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 8, face = "bold"),
            axis.text.y = element_text(size = 7),
            plot.margin = margin(t = 0, r = 0.1, b = 0.1, l = 0.1, unit = "lines"))

    if (metapath == "cazymes") {
        p <- p +
            guides(fill = guide_legend(override.aes = list(size = 2))) +
            labs(x = "\nN pathway") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                  legend.text = element_text(size = 7),
                  legend.title = element_text(size = 8),
                  axis.title.x = element_text(size = 10))
        p_legend <- ggpubr::get_legend(p)
        p <- p +
            guides(fill = "none") +
            inset_element(p_legend, align_to = "full", left = 0.7, bottom = 0.05, right = 0.95, top = 0.35)
    } else {
        p <- p + guides(fill = "none")
    }
    
    return(p)
}

plot_vip_specials_typetwo_b <- function(habitat, environmental_var) {
  vip_CN_plot_df <- vip_members %>% 
    filter(`Environmental Variable` == environmental_var) %>%
    filter(`Data Origin` == habitat)
  
  # number of MAGs
  MAG_num <- length(unique(vip_CN_plot_df$MAG))
  
  vip_wide <- pathways_CN$cn_paths_wide %>%
    filter(genome %in% vip_CN_plot_df$MAG) %>%
    left_join(specialisations, by = c("genome")) %>%
    group_by(CVersatility, NVersatility, specialisation) %>%
    add_tally() %>%
    rename(VIPCspec_count = n) %>%
    ungroup() %>% group_by(CVersatility, NVersatility, nitrogen_specialisation) %>%
    add_tally() %>%
    rename(VIPNspec_count = n) %>% ungroup()

  panel_metapaths <- unique(vip_wide$metapath.C)
  panels <- purrr::map(panel_metapaths, ~ plot_vip_specials_typetwo_b_panel(.x, data = vip_wide))
  panels <- panels[c(1:3, 5, 4)]
  panel_heights <- purrr::map_int(panel_metapaths,
                                  function(x) filter(vip_wide, metapath.C == x) %>% .$CVersatility %>% n_distinct())
  panel_heights <- panel_heights[c(1:3, 5, 4)]

  pprime <- patchwork::wrap_plots(panels, heights = panel_heights, ncol = 1, guides = "keep")
  
  ggsave(plot = pprime, str_c("type2b", environmental_var, habitat, "VIP_specials.png", 
                              sep = "_"), 
         height = 8, width = 6,
         path = output_dir, dpi = 300)
  return(panels)
}

plot_vip_specials_typetwo_b(habitat = "All", environmental_var = "alphaC")
