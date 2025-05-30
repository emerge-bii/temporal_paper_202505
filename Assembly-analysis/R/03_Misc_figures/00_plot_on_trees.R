#' ## Prepares modules for assembly analysis
#' This step subsets data by module membership so that BetaNTI and RCBC can be run on them

#+ include=TRUE
# Loading necessary packages and data
library(tidyverse); packageVersion("tidyverse") # for dataframe processing
library(here)
library(ggrepel)
library(ape)
library(tidytree) # For tree plotting
library(ggtree)
library(ggnewscale)
library(ggstance)
library(viridis)


# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("Assembly-analysis", "outputs")
figures.fp <- here("Assembly-analysis", "figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#### ====================================================================== ####

#### Define some helpful functions:
#### ====================================================================== ####
# This function filteres the input from setup.R by module membership 
filter_by_modules <- function(input = input, 
                  module_file = bog_modules, 
                  module_name = "yellow",
                  output_folder = "bog_module_otus",
                  output_prefix = "bog_module") {
  
  # NOTE: input should be the mapped count data not RA data
  module_taxa_list <- module_file %>%
    filter(color %in% module_name)
  
  otu_table <- input$otu_table
  
  module_otu_table <- otu_table %>%
    filter(genome %in% module_taxa_list$genome) %>%
    select(where( ~ is.numeric(.x) && sum(.x) != 0 | is.character(.x))) # Remove samples that the module is not present in

  # Drop tips that aren't from module 
  tips_to_drop <- input$tree$tip.label[!(input$tree$tip.label %in% module_otu_table$genome)]
  
  filt_tree <- ape::drop.tip(input$tree, tip = tips_to_drop)
  
  filt_taxonomy <- input$taxonomy %>%
    filter(genome %in% module_otu_table$genome)
  
  # Reorder to match tree order
  match.phylo.otu <- picante::match.phylo.data(filt_tree, 
                                               module_otu_table %>%
                                                 select(-genome))
  tree <- match.phylo.otu$phy
  otu_table <- match.phylo.otu$data %>%
    rownames_to_column(var = "genome") %>%
    mutate(rowname = genome) %>%
    column_to_rownames()
  
  # reorder taxonomy to match otu_table and tree tip order
  taxonomy <- filt_taxonomy[rownames(match.phylo.otu$data),]
  
  # filter mapping file to match otu table
  mapping_file <- input$sample_metadata %>%
    filter(temporal_sample_id %in% colnames(module_otu_table))
  
  input_filt <- input
  input_filt$otu_table <- otu_table
  input_filt$tree <- tree
  input_filt$taxonomy <- taxonomy
  input_filt$sample_metadata <- mapping_file

    
  # Write out files
  if (!dir.exists(here(outputs.fp, output_folder))) {dir.create(here(outputs.fp, output_folder))}
  filepath <- paste0(outputs.fp, "/", output_folder, "/", output_prefix, "_", module_name, "_otus.txt") 
  write.table(module_otu_table, file = filepath, row.names = FALSE)

  filepath <- paste0(outputs.fp, "/", output_folder, "/", output_prefix, "_", module_name, ".RDS") 
  saveRDS(input_filt, file = filepath)
  
}

# This function creates line plots of organisms in each module across the samples
plot_modules <- function(module_folder_name = "fen_module_otus",
                          special_otu_plotting = TRUE) {
  file_list <- list.files(here(outputs.fp, module_folder_name),
                          pattern = ".RDS", full.names = T)
  lapply(file_list,
       function(x) {
         input <- readRDS(x)
         title <- gsub(".RDS","", basename(x))
         print(title)
         notu <- nrow(input$otu_table)
         nsamples <- ncol(input$otu_table[-1])
         mod_otu_plot <- t(input$otu_table[-1]) %>%
           as.data.frame() %>%
           rownames_to_column() %>%
           pivot_longer(cols = !matches("rowname")) %>%
           left_join(input$sample_metadata, by = c("rowname" = "temporal_sample_id")) %>%
           arrange(Year__, DepthAvg__, Core__) %>%
           mutate(rowname = fct_inorder(rowname)) %>%
           group_by(name) %>%
           mutate(leadrowname = dplyr::lead(rowname)) %>%
           ungroup() %>%
           mutate(rowname_n = as.numeric(rowname),
                  leadrowname_n = as.numeric(leadrowname),
                  leadrowname_n = ifelse(is.na(leadrowname_n), max(leadrowname_n, na.rm = T) + 1, 
                                         leadrowname_n))
         sample_order <- mod_otu_plot$rowname %>% unique()
         mod_otu_plot$rowname <- factor(mod_otu_plot$rowname, levels = sample_order)
         
         if(special_otu_plotting) {
           special_otus <- c("PMEE01", "20100900_E1D_10", "PLGY01", "20170700_S25_26",
                             "20120600_S3D_1", "20120600_S3M_3", "20120700_S1X_2",
                             "20110800_S1D_5") 
         }
         
         p <- ggplot(mod_otu_plot, aes(x = rowname_n,
                                       y = value)) +
           geom_rect(data = mod_otu_plot %>%
                       select(rowname_n, leadrowname_n, DepthLumping) %>%
                       distinct() %>%
                       mutate(value = max(mod_otu_plot$value)*1.1),
                     aes(fill = DepthLumping, xmin = rowname_n-0.5,
                         xmax = leadrowname_n -0.5),
                     ymin = -Inf, ymax = Inf,
                     alpha = 0.3) +
           scale_x_continuous(breaks = c(1:length(sample_order)),
                              labels = sample_order,
                              limits = c(0.5, length(sample_order)+0.5),
                              expand = c(0,0)) +
           scale_fill_grey(start = 0.8, end = 0.2) +
           geom_point(aes(color = name)) +
           geom_label(data = mod_otu_plot %>%
                        filter(name %in% special_otus) %>%
                        group_by(name) %>%
                        filter(value == max(value)) %>%
                        select(name, rowname_n, value) %>%
                        mutate(hjust = ifelse(grepl("P", name), 1, 1)),
                      aes(label = name, color = name, hjust = hjust), y = Inf, 
                      vjust = 1) +
           geom_line(aes(group = name, color = name)) +
           geom_line(data = mod_otu_plot %>% 
                       filter(name %in% special_otus),
                     aes(group = name), color = "black", size = 2, alpha = 0.3) +
           xlab("SampleID") + ylab("Mapped Abundance") + 
           scale_y_continuous(limits = c(0, NA),
                              expand = c(0,0.5)) +
           guides(color = "none") +
           theme_classic() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
           ggtitle(paste0(title, ", ", notu," members, ", nsamples, " samples"))
         
         ggsave(paste0(outputs.fp, "/", module_folder_name, "/", title, ".png"), plot = p,
                width = 14, height = 7)
       })
}

# This function plots specialisations on a tree
plot_specialisations_on_tree <- function(treedata = specialisations_tre,
                                         specialisation_type = c("S", "N", "C"),
                                         hab_list = c("Palsa", "Bog", "Fen"),
                                         grouptips = NULL,
                                         groupname = NULL,
                                         filter_tree = TRUE) {
  
  # Step 1: Filter tree to drop tips not in list
  if(filter_tree == TRUE) {
    treedata <- treedata %>%
      filter(Habitat %in% hab_list)
    
    tree_filter_list <- input_ra$taxonomy %>%
      filter(genome %in% treedata$genome)
    
    tree.filt <- keep.tip(input_ra$tree, tree_filter_list$genome)
    
    writeLines(paste0("Filtered tree, ", Ntip(tree.filt), " tips remain."))
  } else {
    tree.filt <- input_ra$tree
    tree_filter_list <- input_ra$taxonomy
  }
  
  # Step 2: Prepare tree annotation data; Add taxonomy data to tree tibble and 
  # prepare module heatmap
  tib_tree <- as_tibble(tree.filt)
  
  
  tib_tree_data.circ <- left_join(tib_tree, 
                                  tree_filter_list,
                                  by = c("label" = "genome")) %>%
    as.treedata()
  
  # Setup module pallette
  # module_palette <- unique(module_list$color)
  # names(module_palette) <- module_palette
  
  # Step 3: Create heatmap for circular tree
  heatmap <- treedata %>%
    select(genome, Habitat, all_of(specialisation_type))
  
  
  # Step 4: Plot tree by layers
  # Get archaea node:
  archaea_tips <- treedata %>%
    dplyr::filter(Domain == "d__Archaea")
  
  Arc_node <- getMRCA(tree.filt,
                      tip = archaea_tips$genome)
  # Plot base tree
  p <- ggtree(tib_tree_data.circ,
              layout = "fan", open.angle = 0.1)
  
  # Plot base tree with groups if it exits 
  if(!is.null(grouptips)) {
    tib_tree_data.circ <- groupOTU(tib_tree_data.circ, 
                                   grouptips,
                                   group_name = groupname)
    VIPlab <- data.frame(label = tree.filt$tip.label) %>%
      mutate(VIPs =ifelse(label %in% grouptips, 'symbol("\\267")', '')) # bullet symbol
   # tib_tree_data.circ <- full_join(tib_tree_data.circ, VIPlab, by = "label")
    
    p <- ggtree(tib_tree_data.circ,
                aes(color = .data[[groupname]]),
                layout = "fan", open.angle = 0.1)  %<+% VIPlab +
      geom_tiplab(aes(label = VIPs), # vip tips
                    linetype = NULL, size = rel(1), 
                  offset = 0.25, 
                  parse = T, 
                  align = T) +
      scale_color_manual(values=c("black", "red"), 
                         labels=c("", groupname),
                         guide = "none")
  }
  
  
  # Add archaea annotations if an archaea and arachaea node exists
  if(nrow(archaea_tips) > 0) {
    if(!is.null(Arc_node)) {
      p <- p + geom_point2(aes(subset=(node==Arc_node)), # this is the archaea
                           shape=21, size=2, fill='red')
    } else { # only 1 archaea
      p <- p + 
        geom_tippoint(aes(subset=(label %in% archaea_tips$genome)), # Archaeal tips
                      shape = 21, size = 2, fill = "red", alpha = 0.5)
    }
  }
  
  # p <- p +
  #   geom_tippoint(aes(subset=(label %in% c("PLGY01", "20170700_S25_26"))), # a.stor/m.stor yellow module
  #                 shape = 21, size = 2, fill = "yellow", alpha = 0.5) +    # a.stor/m.crilli green module
  #   geom_tippoint(aes(subset=(label %in% c("PMEE01", "20100900_E1D_10"))),
  #                 shape = 21, size = 2, fill = "green", alpha = 0.5)
  # p
  
  # Set up clade_labels for phyla
  phylumlist <- treedata$Phylum %>%
    unique()
  classlist <- treedata %>%
    dplyr::filter(Phylum == "p__Proteobacteria") %>%
    dplyr::filter(Class != "NA") %>%
    select(Class) %>%
    unique()
  cladelist <- c(phylumlist[phylumlist!="p__Proteobacteria"], classlist$Class)
  
  # remove NAs from cladelist
  cladelist <- cladelist[!is.na(cladelist)]
  clade_labels <- rep(NA, length = length(cladelist))
  
  # IMPORTANT: This only labels phylum with more than 3 % of genomes in them. There are
  #  phyla that have fewer ASVs that are unlabelled (this helps with the spacing of the phylum labels on the tree).
  # For my own knowledge: "M.stor:PLGY01", "A.stor:20170700_S25_26"; M.crilli PMEE01, A.stor2: 20100900_E1D_10
  genome_length_cutoff <- max(floor(0.03*Ntip(tree.filt)), 1) # mark phyla with more than 3% of genomes; or 1 genome, whichever is higher
  
  for (i in 1:length(cladelist)) {
    # if proteobacteria
    if( grepl("proteobacteria", cladelist[i]) ) {
      tip_list <- treedata %>%
        dplyr::filter(Class == cladelist[i])
    } else { # not proteobacteria
      tip_list <- treedata %>%
        dplyr::filter(Phylum == cladelist[i])
    }
    print(paste0("Group ", cladelist[i]))
    print(length(tip_list$genome))
    clade_labels[i] <- ifelse(length(tip_list$genome) > genome_length_cutoff, getMRCA(tree.filt,
                                                                                      tip = tip_list$genome), NA)
    #clade_labels[i] <- getMRCA(sb.tree.rt.filt, tip = tip_list$taxa)
  }
  names(clade_labels) <- cladelist
  
  clade_labels <- na.omit(clade_labels)
  clade_labels <- sort(clade_labels)
  
  # make colors for phyla
  clade_labels_color <- rep(c("grey40", "grey20", "grey80"), 
                            length.out = length(clade_labels))
  
  # Annotate clades
  p1 <- p
  
  for (i in 1:length(clade_labels)) {
    p1 <- p1 + geom_cladelabel(node = clade_labels[i],
                               label = names(clade_labels)[i],
                               align = TRUE, barsize = 2,
                               offset = 3*length(hab_list)*0.2 + 0.4, # 0.4 for 0.2 buffer on each side, then 3 of room for the habitats,
                               offset.text = 0.3,
                               fontsize = 4,
                               hjust = 0,
                               #alpha = 0.5,
                               angle = "auto",
                               color = "grey20")
    #color = clade_labels_color[i])
  }
  # Check clade annotations
  # p1 +
  #   geom_nodelab(aes(label = node)) +
  #   geom_hilight(node = 80 )
  
  # Plot Heatmap
  p2 <- p1
  
  for(special in seq_along(specialisation_type)) {
    if(special == 1) {
      set_offset = 0.2
    } else {
      set_offset = set_offset + 0.2*length(hab_list) # for 3 specialisations
    }

    heatmap_plot <- heatmap %>%
      select(genome, Habitat, all_of(specialisation_type[special])) %>%
      pivot_longer(!all_of(c("genome", "Habitat")), 
                   names_to = "SpecialisationType", values_to = "Specialisation") %>%
      dplyr::filter(Habitat %in% hab_list)  %>%
      pivot_wider(id_cols = genome, names_from = Habitat, values_from = c(Specialisation)) %>%
      select(all_of(c("Palsa", "Bog", "Fen", "genome"))) %>% # reorder columns acording to thaw
      column_to_rownames(var = "genome")
    
    p2 <- gheatmap(p2, heatmap_plot, offset = set_offset,
                   width = .28, colnames = TRUE, hjust = 1,
                   colnames_offset_y = 4, 
                   font.size = rel(4),
                   colnames_angle = 85, color = NA) +
      scale_fill_manual(name = paste0("Specialisations"),
                        values = all_special_colour,
                        breaks = all_special_labels[1:length(all_special_labels)-1],
                        labels = all_special_labels[1:length(all_special_labels)-1],
                        na.value = "grey20") +
      xlim_tree(1.5) + ylim(NA, Ntip(tree.filt) + 0.1*Ntip(tree.filt))
  } 
  
  p2 <- p2 + theme(panel.background = element_rect(fill = NA, color = NA),
                   legend.background = element_rect(fill = NA, color = NA),
                   plot.background = element_rect(fill = NA, color = NA))
  #p2
  ggsave(plot = p2, device = "png",
         filename = paste0(figures.fp, "/", paste(specialisation_type, collapse = "-"), 
                           "_specials_tree.png"),
         height = 10, width = 10, dpi = 600, units = "in")
  ggsave(plot = p2, device = "svg",
         filename = paste0(figures.fp, "/", paste(specialisation_type, collapse = "-"), 
                           "_specials_tree.svg"),
         height = 10, width = 10, dpi = 600, units = "in")
}

# This function plots a single continuous valuegrouped by depth and habitat on a tree
plot_on_tree <- function(treedata = CV_dat, 
                         filter_tree = TRUE,
                         column_to_plot = "CV",
                         column_name = "Coefficient of Variation",
                         hab_list = c("Palsa", "Bog", "Fen"),
                         grouptips = NULL,
                         groupname = NULL,
                         special_tips = TRUE, # mutually exclusive with grouptips
                         plot_type = "png") {

  # Set up colors
  colour_habitat <- c("#703C1B", "#058000", "#0001FF")
  special_levels <- c("generalist", "methanogen", "fermenter", "macromolecule_degrader", "methanotroph", "homoacetogen", "monomer_degrader", NA)
  special_colour <- RColorBrewer::brewer.pal(9, "Set1")[-8]
  named_specials <- special_colour
  names(named_specials) <- special_levels
  special_fill <- special_colour
  names(colour_habitat) <- c("Palsa", "Bog", "Fen")
  
  # Set up tree asthetics
  font_size <- 3 # heatmap font size
  tree_angle <- 20
  
  # Get archaea:
  archaea_tips <- treedata %>%
    dplyr::filter(Domain == "d__Archaea")
  
  # Step 1: Filter tree to drop tips not in list
  if(filter_tree == TRUE) {
    treedata <- treedata %>%
      filter(Habitat %in% hab_list)
    
    tree_filter_list <- input_ra$taxonomy %>%
      filter(genome %in% treedata$genome)
      
    tree.filt <- keep.tip(input_ra$tree, tree_filter_list$genome)
    
    Arc_node <- getMRCA(tree.filt,
                        tip = archaea_tips$genome)
    
    tree.filt <- ape::root(tree.filt, node = Arc_node)
    
    # reorder tree data
    treedata <- treedata[match(tree.filt$tip.label, treedata$genome),]

    writeLines(paste0("Filtered tree, ", Ntip(tree.filt), " tips remain."))
  } else {
    tree.filt <- input_ra$tree
    
    Arc_node <- getMRCA(tree.filt,
                        tip = archaea_tips$genome)
    
    tree.filt <- ape::root(tree.filt, node = Arc_node)
    tree_filter_list <- input_ra$taxonomy
    
    # reorder tree data
    treedata <- treedata[match(tree.filt$tip.label, treedata$genome),]
  }
  
  # Step 2: Prepare tree annotation data; Add taxonomy data to tree tibble and 
  # prepare module heatmap
  tib_tree <- as_tibble(tree.filt)
  
  
  tib_tree_data.circ <- left_join(tib_tree, 
                                  tree_filter_list,
                                  by = c("label" = "genome")) %>%
    as.treedata()
  
  # Setup module pallette
  # module_palette <- unique(module_list$color)
  # names(module_palette) <- module_palette
  
  # Step 3: Create heatmap for circular tree
  heatmap <- treedata %>%
    filter(Habitat %in% hab_list) %>%
    select(all_of(c("genome", column_to_plot, "Habitat", "Depth")))
  
  specials_heatmap <- heatmap %>%
    left_join(specialisations, by = "genome") %>%
    select(genome, specialisation) %>%
    distinct()
  
  # Step 4: Plot tree by layers

  # Plot base tree

  p <- ggtree(tib_tree_data.circ,
              layout = "fan", open.angle = tree_angle,
              root.position = Arc_node) 
  # Plot base tree with groups if it exits 
  # Plot base tree with groups if it exits 
  if(!is.null(grouptips)) {
    tib_tree_data.circ <- groupOTU(tib_tree_data.circ, 
                                   grouptips,
                                   group_name = groupname)
    VIPlab <- data.frame(label = tree.filt$tip.label) %>%
      mutate(VIPs =ifelse(label %in% grouptips, 'symbol("\\267")', '')) # bullet symbol

    p <- ggtree(tib_tree_data.circ,
                aes(color = .data[[groupname]]),
                layout = "fan", open.angle = tree_angle)  %<+% VIPlab +
      geom_tiplab(aes(label = VIPs), # vip tips
                  linetype = NULL, size = rel(3), 
                  offset = 0.2, 
                  parse = T, 
                  align = T) +
      scale_color_manual(values=special_fill, 
                         labels=special_levels,
                         guide = "none")
  }
  if(!is.null(special_tips)) {
    
    speciallab <- specialisations %>% 
      select(genome, specialisation) %>%
      mutate(speciallab = ifelse(!is.na(specialisation), 'symbol("\\267")', '')) %>%
      rename(label = genome)
    tib_tree_data.circ <- left_join(tib_tree, 
                                    tree_filter_list,
                                    by = c("label" = "genome")) %>%
      left_join(speciallab, by = c("label")) %>%
      mutate(specialisation = factor(specialisation, levels = special_levels)) %>%
      as.treedata() 
    
    
    p <- ggtree(tib_tree_data.circ,
                #root.position = Arc_node, 
                aes(color = .data[["specialisation"]]),
                layout = "fan", open.angle = tree_angle, size = rel(0.25)) +
      scale_color_manual(values=special_fill, 
                         labels=special_levels,
                         guide = "none", na.value = "grey20")
  }
  

  
  # Add archaea annotations if an archaea and arachaea node exists
  if(nrow(archaea_tips > 0)) {
    if(!is.null(Arc_node)) {
      p <- p + geom_point2(aes(subset=(node==Arc_node)), # this is the archaea
                           shape=21, size=2, fill='red')
    } else { # only 1 archaea
      p <- p + 
        geom_tippoint(aes(subset=(label %in% archaea_tips$genome)), # Archaeal tips
                      shape = 21, size = 2, fill = "red", alpha = 0.5)
    }
  }
  
  # Plot Heatmap for habitats
  p1 <- p
  for(hab in seq_along(hab_list)) {
    if(hab == 1) {
      set_offset = 0.2
    } else {
      set_offset = set_offset + 0.3*4 # for offset = distance between four depths
    }
    heatmap_plot <- heatmap %>% 
      dplyr::filter(Habitat == hab_list[hab]) %>%
      pivot_wider(names_from = "Depth", values_from = column_to_plot)  %>%
      select(-Habitat) %>%
      column_to_rownames(var = "genome")
    
    p1 <- p1 + new_scale_fill()
    p1 <- gheatmap(p1, heatmap_plot, offset = set_offset,
                   width = .35, colnames = TRUE, hjust = 1,
                   colnames_offset_y = 0.2, font.size = font_size*0.9,
                   colnames_angle = 85, colnames_level = c("0-9", "10-19", "20-29", "30-39"),
                   color = NA) +
      # scale_fill_viridis_c(name = paste0(column_name," - ", hab_list[hab]),
      #                       #low = "white", high = colour_habitat[hab_list[hab]],
      #                       na.value = "grey30", option = "C") +
      scale_fill_continuous(name = paste0(column_name," - ", hab_list[hab]),
                            low = "grey70", high = colour_habitat[hab_list[hab]],
                            na.value = "grey30")
  } 
  p1
  
  # Add specialization heatmap
  spec_heatmap_plot <- specials_heatmap %>% 
    mutate(specialisation = factor(specialisation, levels = special_levels)) %>%
    column_to_rownames(var = "genome")
  p2 <- gheatmap(p1 + new_scale_fill(), spec_heatmap_plot, offset = set_offset + 0.15*5,
                 width = .15, colnames = TRUE, hjust = 1,
                 colnames_offset_y = 0.2, font.size = font_size,
                 colnames_angle = 85, 
                 color = NA) +
  scale_fill_manual(values=named_specials,
                     guide = "none", na.value = "grey30")
    
  set_offset <- set_offset + 0.15*7

  
  
  # Set up clade_labels for phyla
  phylumlist <- treedata$Phylum %>%
    unique()
  classlist <- treedata %>%
    dplyr::filter(Phylum == "p__Proteobacteria") %>%
    dplyr::filter(Class != "NA") %>%
    select(Class) %>%
    unique()
  cladelist <- c(phylumlist[phylumlist!="p__Proteobacteria"], classlist$Class)
  
  # remove NAs from cladelist
  cladelist <- cladelist[!is.na(cladelist)]
  clade_labels <- rep(NA, length = length(cladelist))
  
  # IMPORTANT: This only labels phylum with more than 3 % of genomes in them. There are
  #  phyla that have fewer ASVs that are unlabelled (this helps with the spacing of the phylum labels on the tree).
  # For my own knowledge: "M.stor:PLGY01", "A.stor:20170700_S25_26"; M.crilli PMEE01, A.stor2: 20100900_E1D_10
  genome_length_cutoff <- max(floor(0.03*Ntip(tree.filt)), 1) # mark phyla with more than 3% of genomes; or 1 genome, whichever is higher
  
  for (i in 1:length(cladelist)) {
    # if proteobacteria
    if( grepl("proteobacteria", cladelist[i]) ) {
      tip_list <- treedata %>%
        dplyr::filter(Class == cladelist[i])
    } else { # not proteobacteria
      tip_list <- treedata %>%
        dplyr::filter(Phylum == cladelist[i])
    }
    print(paste0("Group ", cladelist[i]))
    print(length(tip_list$genome))
    clade_labels[i] <- ifelse(length(tip_list$genome) > genome_length_cutoff, getMRCA(tree.filt,
                                                                                      tip = tip_list$genome), NA)
    #clade_labels[i] <- getMRCA(sb.tree.rt.filt, tip = tip_list$taxa)
  }
  names(clade_labels) <- gsub("p__|c__", "", cladelist) # remove the "p__ or c__" from clade names
  
  clade_labels <- na.omit(clade_labels)
  clade_labels <- sort(clade_labels)
  
  
  # make colors for phyla
  clade_labels_color <- rep(c("grey40", "grey20", "grey80"), 
                            length.out = length(clade_labels))
  
  # Annotate clades
  p3 <- p2
  
  for (i in 1:length(clade_labels)) {
    p3 <- p3 + 
      geom_cladelab(node = clade_labels[[i]],
                      label = names(clade_labels)[i],
                      align = TRUE, barsize = 1,
                      offset = set_offset + 0.35, #3*0.02*0.04,
                      #                    offset = 3*length(hab_list)*0.2 + 0.4, # 0.4 for 0.2 buffer on each side, then 3 of room for the habitats,
                      offset.text = 0.1,
                      fontsize = font_size,
                      hjust = 0,
                      #alpha = 0.5,
                      angle = "auto",
                      color = "grey20")
  }
  p3
  
  
  p4 <- p3 + theme(panel.background = element_rect(fill = NA, color = NA),
             legend.background = element_rect(fill = NA, color = NA),
             plot.background = element_rect(fill = NA, color = NA)) + 
    xlim_tree(7.5)
  p4
  
  ggsave(plot = p4, device = plot_type,
         filename = paste0(figures.fp, "/",
                           paste(hab_list, collapse = "-"),
                           "_", column_to_plot, "_tree.", plot_type),
         height = 10, width = 10, dpi = 1000, units = "in", bg = "white")
}



#### ====================================================================== ####
# Plotting colors and settings
#### ====================================================================== #### 
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
all_special_colour <- c("grey80", RColorBrewer::brewer.pal(7, "OrRd")[c(7:2)],
                        RColorBrewer::brewer.pal(5, "Greens")[c(5,2,4,1,3)],
                        RColorBrewer::brewer.pal(9, "Purples")[c(5,8)])
all_special_fill <- all_special_colour
#### ====================================================================== ####

#### ====================================================================== ####
#### Read in data 
#### ====================================================================== ####
# Read in input
input <- input_counts

# Read in the CV data
CV_dat <- read_csv(file = here("data/20221109-cv-mags.csv")) %>%
  left_join(input_ra$taxonomy, by = c("MAG" = "genome")) %>%
  rename(genome = MAG) %>%
  mutate(CV_recip = 1/CV,
         logCV_recip = log(CV_recip)) %>%
  left_join(specialisations, by = "genome")


#### ====================================================================== ####

# Plot Specialisations on tree

#### ====================================================================== ####
# Rearrange specialisations
specialisations_tre <- input_counts$otu_table %>%
  pivot_longer(-genome, names_to = "SampleID", values_to = "counts") %>%
  left_join(input_counts$sample_metadata %>% select(temporal_sample_id, Habitat__),
            by = c("SampleID" = "temporal_sample_id")) %>%
  group_by(genome, Habitat__) %>%
  summarise(HabCounts = sum(counts)) %>%
  filter(HabCounts > 0) %>%
  ungroup() %>%
  group_by(genome) %>%
#  summarize(Habitat = paste(Habitat__, collapse = ",")) %>%
  right_join(specialisations, by = "genome") %>%
  left_join(input_counts$taxonomy, by = "genome") %>%
  mutate(across(contains("specialisation"), ~gsub("_", " " ,.x))) %>%
  mutate(across(contains("specialisation"), ~factor(.x, levels = all_special_labels[1:length(all_special_labels)- 1]))) %>%
  rename(Habitat = Habitat__) %>%
  rename(C = specialisation, N = nitrogen_specialisation, S = sulfur_specialisation) %>%
  ungroup()


plot_specialisations_on_tree(treedata = specialisations_tre,
                             specialisation_type = c("C", "N", "S"), 
                             hab_list = c("Palsa", "Bog", "Fen"))
#### ====================================================================== ####

# Plot avg relative abundance over years on tree
#### ====================================================================== ####
all_vip <- vip_members %>%
  filter(`Data Origin` == "All") %>%
  select(MAG) %>%
  distinct()

plot_on_tree(CV_dat %>%
               mutate(LogAbundance = log(Average_Abundance)), 
             filter_tree = TRUE,
             hab_list = c("Palsa", "Bog", "Fen"),
             column_to_plot = "LogAbundance",
             column_name = "Log(Average Abundance)",
             grouptips = all_vip$MAG,
             groupname = "VIP")




#### ====================================================================== ####
# Plot coefficient of variation on tree
#### ====================================================================== ####
plot_on_tree(CV_dat, 
             filter_tree = TRUE,
             hab_list = c("Palsa", "Bog", "Fen"),
             column_to_plot = "CV",
             column_name = "Coefficient of Variation")

plot_on_tree(CV_dat, 
             filter_tree = TRUE,
             hab_list = c("Palsa", "Bog", "Fen"),
             column_to_plot = "CV",
             column_name = "Coefficient of Variation",
             plot_type = "svg")

plot_on_tree(CV_dat, filter_tree = TRUE,
             hab_list = c("Palsa"))

palsa_vip <- vip_members %>%
  filter(`Data Origin` == "Palsa") %>%
  select(MAG) %>%
  distinct()
plot_on_tree(CV_dat, filter_tree = TRUE,
             hab_list = c("Palsa"),
             grouptips = palsa_vip$MAG,
             groupname = "VIP")


plot_on_tree(CV_dat, filter_tree = TRUE,
             hab_list = c("Bog"))

plot_on_tree(CV_dat, filter_tree = TRUE,
             hab_list = c("Fen"))

all_vip <- vip_members %>%
  filter(`Data Origin` == "All") %>%
  select(MAG) %>%
  distinct()

plot_on_tree(CV_dat, filter_tree = TRUE,
             hab_list = c("Palsa", "Bog", "Fen"),
             column_to_plot = "logCV_recip",
             grouptips = all_vip$MAG,
             groupname = "VIP")

plot_on_tree(CV_dat, filter_tree = TRUE,
             hab_list = c("Bog", "Fen"))
plot_on_tree(CV_dat, filter_tree = TRUE,
             hab_list = c("Palsa", "Bog"))

library(phytools)

CV_dat_filt <- CV_dat %>%
  filter(Depth == "0-9") %>%
  filter(Habitat == "Palsa")

tree.filt <- keep.tip(input_ra$tree, tip = CV_dat_filt$genome)

CV_dat


# match order in tree
CV_dat_filt <- CV_dat_filt[match(tree.filt$tip.label, CV_dat_filt$genome),]


phylosig(tree = tree.filt, x = CV_dat_filt$CV, method = "K", test = TRUE,
         nsim = 10000)

phylosig(tree = tree.filt, x = CV_dat_filt$CV, method = "lambda", test = TRUE,
         nsim = 10000)

source(here("Assembly-analysis/R/helper_scripts/consentrait.R"))

traits <- CV_dat_filt %>%
  select(genome, CV) %>%
  as.data.frame()

calc_phylosig <- function(tree, data, phytest) {
  # get subtree
  tree.filt <- keep.tip(input_ra$tree, tip = data$genome)
  
  # match order in tree
  data_filt <- data[match(tree.filt$tip.label, data$genome),] %>%
    as.data.frame()
  rownames(data_filt) <- data_filt$genome
  
  # run phylogenetic test
  phy_result <- phylosig(tree = tree.filt, x = data_filt$CV_recip, method = phytest, 
                           test = TRUE,
                           nsim = 10000)
  
  return(phy_result)
}


# Reminder:
# https://static1.squarespace.com/static/5459da8ae4b042d9849b7a7b/t/57ea64eae58c62718aa34769/1474979059782/Nesin_Winternitz_Practical_1and2.pdf
# pagel's lamda - compares trait to brownian mostion. lambda = 0, no phylosgenetic signal, traits independent of phylogeny
#                 lambda = 1 brownian motion model - random genetic drift
#                 0 < lambda < 1 between no signal and genetic drift
#                 lambda > 1 selection (presumably)

# Blomberg's K - Do closely related species resemble each other more or less than expected by browninan motion model
#                K goes 0-Inf, K= 1 -> browninan motion, K> 1 speciesmore similar than expected under random drift, 
#                K <1 species less similar than expected under random drift

test <- CV_dat %>%
  group_by(Habitat, Depth) %>%
  nest() %>%
#  pluck(3, 4) %>% select(genome)
  mutate(blombergK = purrr::map(data, ~calc_phylosig(tree = input_ra$tree, data = ., phytest = "K")),
         pagellambda = purrr::map(data, ~calc_phylosig(tree = input_ra$tree, data = ., phytest = "lambda"))) %>%
  mutate(bK = map_dbl(blombergK, ~.$K),
         bP = map_dbl(blombergK, ~.$P),
         pglP = map_dbl(pagellambda, ~.$P)) %>%
  mutate(InterpretK = ifelse(bP >= 0.05, "Not significantly different from BM", 
                             ifelse(bK < 1, "Less similar than random", "more similar than random")))



# More exploration of CV
CV_dat %>% 
  mutate(Habitat = factor(Habitat, levels = c("Palsa", "Bog", "Fen"))) %>%
  group_by(Habitat) %>% 
#  nest() %>% 
#  pluck(2, 3) %>% 
#  pluck(2, 2) %>%
#  pluck(2, 1) %>% 
  ggplot(aes(x = Habitat, y = logCV_recip)) + 
  geom_boxplot(aes(group = Habitat, color = Habitat)) +
  facet_wrap(~Depth, ncol = 4) +
  # annotate(geom = "text", y = -Inf, x = -Inf, vjust = -0.1, hjust = 0, 
  #          label = "less stable MAGs", size = rel(2)) +
  # annotate(geom = "text", y = Inf, x = -Inf, vjust = 1, hjust = 0, 
  #          label = "more stable MAGs", size = rel(2)) +
  #coord_cartesian(ylim = c(0,NA)) +
  scale_color_manual(values = colour_habitat) +
  ylab("log Reciprocal Coefficient of Variation (higher = more stable)") +
  theme_bw() +
  theme(panel.spacing = unit(0,"lines"))


# Get ancestral states for branches
# Set up colors
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
names(colour_habitat) <- c("Palsa", "Bog", "Fen")
# run fastANC
test <- CV_dat %>%
  select(genome, logCV_recip, Habitat, Depth) %>%
  group_by(Habitat, genome) %>%
  summarize(meanlogCV = mean(logCV_recip)) %>%
  nest() %>% #pluck(2, 1) %>% pull(meanlogCV, genome) %>% length()
  mutate(cv_mat = purrr::map(data, ~pull(., meanlogCV, genome))) %>%
  mutate(filt_tree = purrr::map(data, ~keep.tip(input_ra$tree, .$genome))) %>%
  mutate(fit = purrr::map2(cv_mat, filt_tree, ~phytools::fastAnc(tree = .y, x = .x, vars = TRUE, CI = TRUE))) %>%
  mutate(td = purrr::map2(cv_mat, filt_tree, ~data.frame(node = nodeid(.y, names(.x)),
                                                         trait = .x))) %>% 
  mutate(nd = purrr::map(fit, ~data.frame(node = names(.x$ace), trait = .x$ace))) %>%
  mutate(d = purrr::map2(td, nd, ~rbind(.x,.y) %>% mutate(node = as.numeric(node)))) %>%
  mutate(plot_tree = purrr::map2(filt_tree, d, ~full_join(.x, .y, by = "node")),
         plot_tree = purrr::map(plot_tree, ~left_join(.x, input_ra$taxonomy, 
                                                       by = c("label" = "genome")))) %>%
  # remove unnecessary columns
  select(-td, -nd, -d, -filt_tree ) %>%
  mutate(plot = purrr::map2(plot_tree, Habitat, ~ggtree(.x, aes(color = trait), layout = "fan", 
                                              ladderize = FALSE, continuous = 'colour', 
                                              open.angle = 10) +
                             scale_color_gradient(low = "white",
                                                  name = "Log Reciprocal of CV",
                                                  high = colour_habitat[.y], na.value = "grey30") +
#                             scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
                             #geom_tiplab(hjust = -.1) + 
                             xlim(0, 1.2) + 
                             ggtitle(paste0(.y, " Tree RecipCV")) +
                             theme(legend.position = c(.03, .85)) ))

test %>% pluck(5,1) %>% as.treedata() %>% pull(Phylum, label)

test %>% pluck(6,1) 


pt <- test %>% pluck(5,2) %>%
  ggtree(aes(color = trait), layout = "fan", 
         ladderize = FALSE, continuous = 'colour', 
         open.angle = 10) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  #geom_tiplab(hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85)) 

pt

phylumlist <- treedata$Phylum %>%
  unique()
classlist <- treedata %>%
  dplyr::filter(Phylum == "p__Proteobacteria") %>%
  dplyr::filter(Class != "NA") %>%
  select(Class) %>%
  unique()
cladelist <- c(phylumlist[phylumlist!="p__Proteobacteria"], classlist$Class)

# remove NAs from cladelist
cladelist <- cladelist[!is.na(cladelist)]
clade_labels <- rep(NA, length = length(cladelist))

# IMPORTANT: This only labels phylum with more than 3 % of genomes in them. There are
#  phyla that have fewer ASVs that are unlabelled (this helps with the spacing of the phylum labels on the tree).
# For my own knowledge: "M.stor:PLGY01", "A.stor:20170700_S25_26"; M.crilli PMEE01, A.stor2: 20100900_E1D_10
genome_length_cutoff <- max(floor(0.03*Ntip(tree.filt)), 1) # mark phyla with more than 3% of genomes; or 1 genome, whichever is higher

for (i in 1:length(cladelist)) {
  # if proteobacteria
  if( grepl("proteobacteria", cladelist[i]) ) {
    tip_list <- treedata %>%
      dplyr::filter(Class == cladelist[i])
  } else { # not proteobacteria
    tip_list <- treedata %>%
      dplyr::filter(Phylum == cladelist[i])
  }
  print(paste0("Group ", cladelist[i]))
  print(length(tip_list$genome))
  clade_labels[i] <- ifelse(length(tip_list$genome) > genome_length_cutoff, getMRCA(tree.filt,
                                                                                    tip = tip_list$genome), NA)
  #clade_labels[i] <- getMRCA(sb.tree.rt.filt, tip = tip_list$taxa)
}
names(clade_labels) <- cladelist

clade_labels <- na.omit(clade_labels)
clade_labels <- sort(clade_labels)

# make colors for phyla
clade_labels_color <- rep(c("grey40", "grey20", "grey80"), 
                          length.out = length(clade_labels))



#tips : <- as_tibble() %>% filter(node < min(parent))
#### ====================================================================== ####
#### Filter OTU tables and input for modules
#### ====================================================================== ####
# All Modules
lapply(seq_along(unique(all_modules$color)),
       function(x) {
         input$sample_metadata <- input$sample_metadata
         input$otu_table <- input$otu_table %>%
           select(genome, matches(input$sample_metadata$temporal_sample_id))
         
         filter_by_modules(input = input,
                           module_file = all_modules,
                           module_name = unique(all_modules$color)[x],
                           output_folder = "all_module_otus",
                           output_prefix = "all_module")
       })

# Palsa Modules
lapply(seq_along(unique(palsa_modules$color)),
       function(x) {
         input$sample_metadata <- input$sample_metadata %>%
           filter(Habitat__ == "Palsa")
         input$otu_table <- input$otu_table %>%
           select(genome, matches(input$sample_metadata$temporal_sample_id))
         
         filter_by_modules(input = input,
                           module_file = palsa_modules,
                           module_name = unique(palsa_modules$color)[x],
                           output_folder = "palsa_module_otus",
                           output_prefix = "palsa_module")
       })

# Bog Modules
lapply(seq_along(unique(bog_modules$color)),
       function(x) {
         input$sample_metadata <- input$sample_metadata %>%
           filter(Habitat__ == "Bog")
         input$otu_table <- input$otu_table %>%
           select(genome, matches(input$sample_metadata$temporal_sample_id))
         
         filter_by_modules(input = input,
                           module_file = bog_modules,
                           module_name = unique(bog_modules$color)[x],
                           output_folder = "bog_module_otus",
                           output_prefix = "bog_module")
       })
# Fen Modules
lapply(seq_along(unique(fen_modules$color)),
       function(x) {
         input$sample_metadata <- input$sample_metadata %>%
           filter(Habitat__ == "Fen")
         input$otu_table <- input$otu_table %>%
           select(genome, matches(input$sample_metadata$temporal_sample_id))
         
         filter_by_modules(input = input,
                           module_file = fen_modules,
                           module_name = unique(fen_modules$color)[x],
                           output_folder = "fen_module_otus",
                           output_prefix = "fen_module")
       })


#### ====================================================================== ####

# Plot organisms in modules in each sample to see their co-occurance patterns
# across samples
#### ====================================================================== ####
plot_modules(module_folder_name = "all_module_otus")

plot_modules(module_folder_name = "palsa_module_otus")

plot_modules(module_folder_name = "bog_module_otus")

plot_modules(module_folder_name = "fen_module_otus")

#### ====================================================================== ####
# Phylofactor for temperature response
#### ====================================================================== ####
#### ====================================================================== ####

# Plot Modules on Trees
#### ====================================================================== ####
plot_module_on_tree(modules_to_plot = "All")

plot_module_on_tree(modules_to_plot = "Palsa", filter_tree = TRUE)

plot_module_on_tree(modules_to_plot = "Bog")

plot_module_on_tree(modules_to_plot = "Fen")

#### ====================================================================== ####