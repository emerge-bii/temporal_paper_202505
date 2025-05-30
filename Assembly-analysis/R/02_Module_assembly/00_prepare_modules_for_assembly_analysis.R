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

# This function plots modules on a tree
plot_module_on_tree <- function(modules_to_plot = c("Palsa"),
                                module_list = module_list_comb, 
                                filter_tree = TRUE) {
  # Get module lists
  module_list <- module_list %>%
    dplyr::filter(Module_type %in% modules_to_plot) %>%
    dplyr::filter(!grepl("grey", color)) # Remove grey modules; the "trash bin module"
  
  # Step 1: Filter tree to drop tips not in module list(s)
  if(filter_tree == TRUE) {
    tree_filter_list <- module_list %>%
      dplyr::filter(!grepl("grey", color)) %>% # remove grey taxa (these are the trash bin taxa)
      select(genome, Domain:Species) %>%
      distinct()
      
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
  module_palette <- unique(module_list$color)
  names(module_palette) <- module_palette
  
  # Step 3: Create heatmap for circular tree
  heatmap <- module_list %>%
    select(genome, color, Module_type) %>%
    nest(module_data = -Module_type) %>% 
    mutate(module_data = purrr::map(module_data, 
                             ~.x %>%
                               dplyr::filter(!grepl("grey", color)) %>%
                               mutate(presence = color) %>%
                               pivot_wider(id_cols = c("genome"),names_from = color, 
                                           values_from = presence, values_fill = NA)))
  
  
  # Step 4: Plot tree by layers
  # Get archaea node:
  archaea_tips <- module_list %>%
    dplyr::filter(!grepl("grey", color)) %>%
    dplyr::filter(Domain == "d__Archaea")

  Arc_node <- getMRCA(tree.filt,
          tip = archaea_tips$genome)
  # Plot base tree
  p <- ggtree(tib_tree_data.circ,
              layout = "fan", open.angle = 5)
  
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
  
  p <- p +
    geom_tippoint(aes(subset=(label %in% c("PLGY01", "20170700_S25_26"))), # a.stor/m.stor yellow module
                  shape = 21, size = 2, fill = "yellow", alpha = 0.5) +    # a.stor/m.crilli green module
    geom_tippoint(aes(subset=(label %in% c("PMEE01", "20100900_E1D_10"))),
                  shape = 21, size = 2, fill = "green", alpha = 0.5)
  p
  
  # Set up clade_labels for phyla
  phylumlist <- module_list$Phylum %>%
    unique()
  classlist <- module_list %>%
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
      tip_list <- module_list %>%
        dplyr::filter(Class == cladelist[i])
    } else { # not proteobacteria
      tip_list <- module_list %>%
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
                               offset = 1,
                               offset.text = 0.3,
                               fontsize = 3,
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
  
  # Plot Heatmap for modules
  p2 <- p1
  for(module_type in modules_to_plot) {
    heatmap_plot <- heatmap %>% 
      dplyr::filter(Module_type == module_type) %>% 
      pull(module_data) %>% pluck(1) %>% 
      column_to_rownames(var = "genome")
    p2 <- p2 + new_scale_fill()
    p2 <- gheatmap(p2, heatmap_plot, offset = .2,
                   width = .3, colnames = FALSE) +
      scale_fill_manual(values = module_palette, name = "Module membership", 
                        na.value = "white") +
      xlim_tree(1.5)
  }
  
  ggsave(plot = p2, filename = paste0(figures.fp, "/", modules_to_plot, "_modules_tree.png"),
         height = 10, width = 10, dpi = 600, units = "in")
}

#### ====================================================================== ####


#### ====================================================================== ####
#### Read in data 
#### ====================================================================== ####
# Read in input
input <- input_counts

# Read in the module membership data
# All modules
all_modules <- read_csv(file = here("data", "module_membership", 
                                    "All_Module_membership.csv"),
                        col_names = c("genome", "module", "color"),
                        skip = 1) %>%
  left_join(input_ra$taxonomy, by = "genome") %>%
  filter(color != "grey") # remove the grey module

# Palsa modules
palsa_modules <- read_csv(file = here("data", "module_membership", 
                                    "Palsa_Module_membership.csv"),
                          col_names = c("genome", "module", "color"),
                        skip = 1) %>%
  left_join(input_ra$taxonomy, by = "genome") %>%
  filter(color != "grey") # remove the grey module

# Bog modules
bog_modules <- read_csv(file = here("data", "module_membership", 
                                    "Bog_Module_membership.csv"),
                        col_names = c("genome", "module", "color"),
                        skip = 1) %>%
  left_join(input_ra$taxonomy, by = "genome") %>%
  filter(color != "grey") # remove the grey module
# Fen modules
fen_modules <- read_csv(file = here("data", "module_membership", 
                                    "Fen_Module_membership.csv"),
                        col_names = c("genome", "module", "color"),
                        skip = 1) %>%
  left_join(input_ra$taxonomy, by = "genome") %>%
  filter(color != "grey") # remove the grey module

# Create a combined list of modules

module_list_comb <- bind_rows(
  all_modules %>% mutate(Module_type = "All"),
  palsa_modules %>% mutate(Module_type = "Palsa"),
  bog_modules %>% mutate(Module_type = "Bog"),
  fen_modules %>% mutate(Module_type = "Fen")
)



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


# Plot Modules on Trees
#### ====================================================================== ####
plot_module_on_tree(modules_to_plot = "All")

plot_module_on_tree(modules_to_plot = "Palsa", filter_tree = TRUE)

plot_module_on_tree(modules_to_plot = "Bog")

plot_module_on_tree(modules_to_plot = "Fen")

#### ====================================================================== ####