# Install packages
here::i_am("setup.R")
library(here)
library(ape)
library(phytools)
library(picante)
library(jsonlite)
library(tidyverse)
library(data.table)

# Set V1 directory
v1_mags_directory <- here("data", "Emerge_MAGs_v1")
v4_singlem_directory <- here("data", "SingleM_otu_tables_v4")
v9_emerge_distillate_directory <- here("data", "EMERGE_distillate_v9")
v2_dram_distillate_directory <- here("data", "DRAM_distillate_v2")
v1_manual_methanogen_directory <- here("data", "Manual_methanogen_calls_v1")
v4_dram_annotations_directory <- here("data", "DRAM_annotations_v4")
v5_emerge_metaT_directory <- here("data", "Emerge_metaTs_v5")
v7_metaT_processed_directory <- here("data", "Emerge_metaTs_processed_v7")
v1_aa_frequencies_directory <- here("data", "AA_frequencies_v1")

# Memory saving? - should large intermediate objects be removed to save memory after they are no longer needed:
memory_saving <- TRUE # default is false

read_genome_info <- function() {
    d <- bind_rows(
        read_tsv(here(v1_mags_directory, "Field_bins", "Woodcroft_2018", "tax_complete_contam.txt")) %>% mutate(genome_set = "Woodcroft_2018"),
        read_tsv(here(v1_mags_directory, "Field_bins", "Cronin_2021", "tax_complete_contam.txt")) %>% mutate(genome_set = "Cronin_2021"),
        read_tsv(here(v1_mags_directory, "SIP_bins", "tax_complete_contam.txt")) %>% mutate(genome_set = "SIP_bins"),
        read_tsv(here(v1_mags_directory, "JGI_bins", "tax_complete_contam.txt")) %>% mutate(genome_set = "JGI_bins"),
        )
    return(d)
}

# Adds gtdb2.0 taxonomy to genome_info
read_gtdb2_taxonomy <- function(genome_info = genome_info) {
  gtdb_bac <- read_tsv(here("data", "gtdbtk.bac120.summary.tsv"))
  gtdb_arch <- read_tsv(here("data", "gtdbtk.ar53.summary.tsv"))
  
  gtdb2_taxa <- bind_rows(gtdb_bac, gtdb_arch)
  
  
  genome_info <- genome_info %>%
    left_join(gtdb2_taxa %>%
                select(user_genome, classification), 
              by = c("genome" = "user_genome")) %>%
    rename(gtdb2_taxonomy = classification)

  return(genome_info)
}

read_genome_clusters <- function() {
    d <- read_tsv(here("data", "95_genome_clusters.tsv"))

    return(d)
}

# Load the MAG abundances
read_trimmed_mean <- function() {
    d <- read_csv(here("data", "MAG_otu_Trimmed_Mean.txt")) %>%
        rename(genome = Genome) %>%
        pivot_longer(cols = !genome, names_to = "temporal_sample_id", values_to = "coverage") %>%
        filter(coverage > 0) %>%
        group_by(temporal_sample_id) %>%
        mutate(relabund_of_recovered = coverage / sum(coverage))
    
    rel_abund <- d %>%
        select(genome, temporal_sample_id, relabund_of_recovered) %>%
        pivot_wider(names_from = temporal_sample_id, values_from = relabund_of_recovered,
                    values_fill = 0)
    coverage <- d %>%
        select(genome, temporal_sample_id, coverage) %>%
        pivot_wider(names_from = temporal_sample_id, values_from = coverage,
                    values_fill = 0)
    
    read_trimmed_mean <- list(trimmed_mean = d, 
                              rel_abund = rel_abund, coverage = coverage)
    
    return(read_trimmed_mean)
}

# Use temporal_sample_id to match
read_sample_metadata <- function() {
    d <- read_csv(here("data", "coring_geochem_sequenced_samples_1.0.0.csv")) %>%
        mutate(
            sample = SampleID__,
            temporal_sample_id = metaG_JGI_NovaSeq__,
            temporal_sample_id = gsub(" ", "_", temporal_sample_id)
        )

    return(d)
}

# Get read counts for each sample
read_read_counts <- function() {
    cmr_read_name_conversion <- read_tsv(here("data", "cmr_sample_info.tsv"))
    d <- read_csv(here("data", "cmr_read_counts.csv"), col_names = c("read_basename", "read_count")) %>%
        left_join(cmr_read_name_conversion, by = "read_basename") %>%
        filter(`included?`) %>%
        select(-`included?`, -sequencer, -read_basename)
}

# Load all SingleM otu tables
read_singlem_otu_tables <- function() {
    cmr_read_name_conversion <- read_tsv(here("data", "cmr_sample_info.tsv"))
    # Read in collated sample otu table
    d <- read_tsv(here(v4_singlem_directory, "collated_appraised.tsv")) %>%
        left_join(cmr_read_name_conversion, by = c("sample" = "read_basename")) %>%
        filter(`included?`) %>%
        select(-`included?`, -sequencer)
}

# Load the SingleM appraise summary table
read_singlem_appraise_summary <- function() {
    cmr_read_name_conversion <- read_tsv(here("data", "cmr_sample_info.tsv"))
    d <- read_csv(here("data", "singleM_appraise_summary.csv")) %>%
        left_join(cmr_read_name_conversion, by = c("sample" = "read_basename")) %>%
        filter(`included?`) %>%
        select(-`included?`, -sequencer, -sample)
}


read_genome_tree <- function(otu_table, visualize = FALSE) {
    require(ape)
    require(phytools)
    
    arch.tree.raw <- phytools::read.newick(here("data", "gtdbtk.ar53.decorated.tree.nwk"))
    bac.tree.raw <- phytools::read.newick(here("data", "gtdbtk.bac120.decorated.tree.nwk"))
    
    # Artificially merge the two trees based on the rooting proposed in paper out of 
    # Rob Knight's lab: Zhu et al. 2019 https://www.nature.com/articles/s41467-019-13443-4: 
    # "Phylogenomics of 10,575 genomes reveals evolutionary proximity between domains Bacteria and Archaea"
    bac.arch.tree.raw <- ape::bind.tree(bac.tree.raw, arch.tree.raw)
    
    # check tree merge
    if(visualize) {
        writeLines("Plotting initial tree")
        plot(bac.arch.tree.raw)
    }
    
    # Fix the names by removing the extra quotes
    bac.arch.tree.raw$tip.label <- gsub("'", "", bac.arch.tree.raw$tip.label)
    
    # Drop tips that aren't from EMERGE MAGs 
    #(since we used gtdb-tk to make the tree these should be tips that are in gtdb, but not our MAGs)
    gtdb_tips <- bac.arch.tree.raw$tip.label[!(bac.arch.tree.raw$tip.label %in% otu_table$genome)]
    
    bac.arch.tree.nogtdb <- ape::drop.tip(bac.arch.tree.raw, tip = gtdb_tips)
    bac.tree <- ape::drop.tip(bac.tree.raw, tip = gtdb_tips)
    arch.tree <- ape::drop.tip(arch.tree.raw, tip = gtdb_tips)
    
    if(visualize) {
        writeLines("Plotting filtered tree")
        # plot resulting tree
        plot(bac.arch.tree.nogtdb)
    }
    
    # Root the tree
    ba.tree.rt <- phytools::midpoint.root(bac.arch.tree.nogtdb)
    
    # check tree
    if(visualize) {
        writeLines("Plotting filtered tree")
        # plot resulting tree
        plot(ba.tree.rt)
    }
    return(list(bacteria_tree = bac.tree, 
                archaea_tree = arch.tree,
                bac_arch_tree = ba.tree.rt))
}



# Check for matching MAGs between datasets and get taxonomy
get_taxonomy <- function(otu_table, genome_info, taxonomy_version = "gtdb2.0") {
    # If gtdb2.0 taxonomy is selected, rename variable accordingly
    if(taxonomy_version == "gtdb2.0") {
      genome_info <- genome_info %>%
        mutate(gtdb1.5_taxonomy = taxonomy) %>%
        mutate(taxonomy = gtdb2_taxonomy)
    }
    
    # Find genomes missing from OTU table (if any)
    missing_genomes_data <- otu_table$genome[!(otu_table$genome %in% 
                                               genome_info$genome)]
    writeLines(paste("Genomes missing from the otu_table that are present in",
                     "genome_info:", paste(missing_genomes_data, collapse = ", ")))
 
    
    # Prepare taxonomy table
    taxonomy_table <- genome_info %>%
        filter(genome %in% otu_table$genome) %>% # select only genomes in current data
        select(genome, taxonomy) %>%
        separate(taxonomy, 
                 into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                 sep = ";") %>%
        select(genome, Domain:Species)
    
    return(taxonomy_table)
}

get_methane_flux <- function() {
  ch4_flx_raw <- readxl::read_excel(here("data", 
                                     "MethaneDataforDylan.xlsx"))
  
  ch4_flx <- ch4_flx_raw %>%
    rename(`d13CH4_AVE_1wk_before_coring`= 3,
           `d13CH4_SE_1wk_before_coring` = 4,	
           `CH4_Flux_AVE_1wk_before_coring`	= 5,
           `CH4_Flux_SE_1wk_before_coring` = 6,
           `d13CH4_AVE_2wk_before_coring`= 7,
           `d13CH4_SE_2wk_before_coring` = 8,	
           `CH4_Flux_AVE_2wk_before_coring`	= 9,
           `CH4_Flux_SE_2wk_before_coring` = 10,
           `d13CH4_AVE_3wk_before_coring`= 11,
           `d13CH4_SE_3wk_before_coring` = 12,	
           `CH4_Flux_AVE_3wk_before_coring`	= 13,
           `CH4_Flux_SE_3wk_before_coring` = 14,
           `d13CH4_AVE_4wk_before_coring`= 15,
           `d13CH4_SE_4wk_before_coring` = 16,	
           `CH4_Flux_AVE_4wk_before_coring`	= 17,
           `CH4_Flux_SE_4wk_before_coring` = 18,
    ) %>%
    filter(!grepl("d13CH4", d13CH4_AVE_1wk_before_coring)) %>%
    mutate_at(vars(contains("CH4")), as.numeric) %>%
    mutate(Year__ = lubridate::year(`Coring Date`),
           Month__ = lubridate::month(`Coring Date`),
           Date__ = lubridate::date(`Coring Date`)) %>%
    rename(Habitat__ = Site) %>%
    select(-`Coring Date`) %>%
    select(Date__, Habitat__, Year__, Month__, everything())
  return(ch4_flx)

}

clean_sample_metadata <- function(sample_metadata) {
    # Bring in Suzanne's WTD and temperature calculations
    t_wtd_sum <- read_csv(here("identifying-outlier-years", "results", 
                               "t_wtd_summaries_July2011-2017samplings.csv"))
    
    # pivot wider so there's a column for every time interval-summary-statistic
    t_wtd_sum_wide <- t_wtd_sum %>% 
      select(Site__, Year__, Interval, contains("AirTemperature")) %>%
      mutate(SiteYear = paste(Site__, Year__)) %>%
      pivot_wider(id_cols = c(Site__, Year__), 
                  names_from = Interval, values_from = contains("AirTemperature")) %>%
      select(Site__, Year__, samplingdate_min_AirTemperature_7d, # select only one variant of all of these identical columns
             samplingdate_max_AirTemperature_7d,
             samplingdate_mean_AirTemperature_7d,
             n_AirTemperature_7d:sd_AirTemperature_all_growing) %>%
      rename(samplingdate_max_AirTemperature = samplingdate_max_AirTemperature_7d,
             samplingdate_min_AirTemperature = samplingdate_min_AirTemperature_7d,
             samplingdate_mean_AirTemperature = samplingdate_mean_AirTemperature_7d)
  
    wtd_perc_sum <- read_csv(here("identifying-outlier-years", "results", 
                               "wtd_summaries_July2011-2017samples.csv"))
    
    wtd_perc_sum_simple <- wtd_perc_sum %>%
      select(SampleID__, contains("pct_time"))

    
    # Read in and organize methane data from Carrie
    ch4_flx <- get_methane_flux()
    
    # Read in list of outliers (decided 05/10/2022)
    # Powerpoint of reasonings here: https://docs.google.com/presentation/d/1Xon7PO6dLge-NxfSh3UEj6ABCjDnX-pgS2xS1IHc1w0/edit#slide=id.g1277bcffdd5_0_83
    # outlier list here: https://docs.google.com/spreadsheets/d/1f0xB5rd1QY1qxPPKVSXt0itwFomuIpFs6T4ZhLvZCUg/edit#gid=0
    outlier_collection <- readxl::read_xlsx(here("data", "Sample Outlier Collection.xlsx"))
    
        
    sample_metadata_mod <- sample_metadata %>%
        # Filter out sites that aren't the autochamber sites
        filter(CoreGroup__ == "MainAutochamber") %>%
        # Filter out samples that don't have a temporal sample id
        # Noticed on 1/17/2023 there's a discrepancy between samples that are MainAutochamber
        # and temporal_sample_id samples - I think this is because the temporal_sample_ids are samples for which
        # there were JGI NovaSequenced MetaGs available at the start of this analysis. There are some samples in the
        # updated sample metadatasheet that don't have jgi metaGs. 
        filter(!is.na(temporal_sample_id)) %>%
        # Add in an outlier column for temporal paper's outlier designations
        mutate(Outlier_status = ifelse(temporal_sample_id %in% outlier_collection$`Sample ID`, "temporal.outlier", NA)) %>%
        # Bring in Suzanne's temperature and WTD calculations
        left_join(t_wtd_sum_wide, by = c("Site__", "Year__")) %>%  
        left_join(wtd_perc_sum_simple, by = c("SampleID__")) %>%
        # Bring in Carrie's Methane Flux data
        left_join(ch4_flx,
                by = c("Date__", "Habitat__", "Year__", "Month__")) %>%
        # alphaC value calculated as per Suzanne's instructions
        # METHANOGENESIS ONLY; Rule of thumb:  acetoclastic > 1, hydrogenotrophic >> 1
        mutate(alphaC = (d13C_CO2__ + 1000)/(d13C_CH4__ + 1000)) %>% 
        # Adding in the methane co2 ratio as described in the W&S paper and also as Derek said is typically done
        mutate(CH4_CO2_porewater_percent = CH4.mM__/(CH4.mM__ + CO2.mM__),
             CH4_CO2_porewater_ratio = CO2.mM__/CH4.mM__) %>%
        # Make a column for depth based on Virginia Rich's basecamp suggestion: 
        # https://3.basecamp.com/4869263/buckets/18955755/messages/4572257670
        # (lump 0-9, 10-19, 20-29 cm etc. based on DepthAvg__ column)
        # Currently keeping the DepthCodeMod column for backwards compatibility,
        # but all code should move towards using DepthLumping instead.
        mutate(DepthLumping =  ifelse(DepthAvg__ >= 1 & DepthAvg__ < 9, "0-9", 
                               ifelse(DepthAvg__ >= 10 & DepthAvg__ < 20, "10-19", 
                               ifelse(DepthAvg__ >= 20 & DepthAvg__ < 30, "20-29",
                               ifelse(DepthAvg__ >= 30 & DepthAvg__ < 40, "30-39",
                               ifelse(DepthAvg__ >= 40 & DepthAvg__ < 50, "40-49",
                               ifelse(DepthAvg__ >= 50 & DepthAvg__ < 60, "50-59",
                               ifelse(DepthAvg__ >= 60 & DepthAvg__ < 70, "60-69",
                               ifelse(DepthAvg__ >= 70 & DepthAvg__ < 80, "70-79",
                               NA))))))))) %>%
       # ALD and temp fill in NAs for samples that have values in 
       # other columns. This fixes a quirk with how the data was entered
       # in the original field data sheets in one year that meant only
       # the first sample for a core had a recorded ALD and Air temp.
       # NOTE: this may no longer be necessary after the new release of the metadata sheet (05/10/2022)
       group_by(Year__, Month__, Habitat__, Site__, Core__) %>%
       fill(ALD.cm__) %>%
       fill(`T_air.deg_C`) %>%
       ungroup() %>%
       # Try to solve the WTD/ALD detection issue (namely that when WTD
       # is below detection, we want to record this b/c it's particularly
       # relevant biologically but that the limit of detection changes
       # from year to year as the instruments changed)
       # Set samples "below detection" to -101 (lower than all observed values)
       # for WTD and to 101 for ALD (lower than all observed values and lower 
       # than the longest probe used); All these saved in new variables (e.g. 
       # "ALD" rather than "ALD.cm__") so they don't change any code that is 
       # dependent on the old versions.
       mutate(WTD = ifelse(grepl("Detect", WTD.cm_neg_is_below_sfc__), "-101", 
                           WTD.cm_neg_is_below_sfc__),
              WTD = ifelse(grepl("notes|Difficult|Pending|pending|Not Measured", WTD.cm_neg_is_below_sfc__), 
                           NA, WTD),
              WTD = ifelse(Habitat__ == "Palsa", NA, WTD), # as per discussion, it's nonsensical to assign WTD in Palsa
              WTD = as.numeric(WTD)) %>%
       mutate(ALD = ifelse(grepl("Not Measured", ALD.cm__), NA, 
                           ALD.cm__),
              ALD = ifelse(grepl("Detect", ALD), "101", 
                           ALD),
              ALD = as.numeric(ALD)) %>% 
       # As per group discussions, reframing ALD int terms of distance from sample.
       # Min_PF_dist_from_sample represents the exact or minimum (in the case of 
       # Below detection samples) distance the frozen soil is from the average depth
       # of the sample. 
       mutate(Min_ALD = ifelse(grepl("Not Measured", ALD.cm__), NA, 
                                                ALD.cm__),
              Min_ALD = ifelse(grepl("Detect", Min_ALD),
                               gsub("(Below Detection \\(>)(.{1,3})(\\))", "\\2", Min_ALD), 
                               Min_ALD),
              Min_ALD = gsub("=", "", Min_ALD), # fix the cases where there is a >= sign
              Min_ALD = ifelse(grepl("Detect", Min_ALD), NA, Min_ALD), # fix "Below detection" with no measurement of ALD probe (yr = 2013, PalsaHole, Hodgkins Fen samples)
              Min_ALD = as.numeric(Min_ALD),
              ALD_detectable = ifelse(grepl("Detect", ALD.cm__), "No", NA),
              Min_PF_dist_from_sample = Min_ALD - DepthAvg__) %>%
       # Doing something similar for WTD
       mutate(DepthAvg__neg = -1*DepthAvg__, # convert depthavg to be on same scale as WTD; neg = below surface
              Samp_dist_WT = WTD - DepthAvg__neg)  %>% # Pos = WT above sample, neg = WT below sample
       # Set the measurements we're unsure of being less than or equal to 0 to NA
       mutate(`T_soil.deg_C` = ifelse(grepl("<=0", `T_soil.deg_C`), NA, `T_soil.deg_C`),
              `T_soil.deg_C` = as.numeric(`T_soil.deg_C`))
    
    return(sample_metadata_mod)
}


# Check for matching names between datasets
check_sample_metadata <- function(otu_table, sample_metadata, 
                                  sample_id_column = "temporal_sample_id",
                                  taxonomy = NULL,
                                  tree = NULL) {
    # Checking OTU table and sample metadata
    # Find the missing samples for otu_table and sample metadata
    missing_genomic_data <- na.omit(sample_metadata$temporal_sample_id[!(sample_metadata$temporal_sample_id %in% names(otu_table[-1]))]) # Samples missing from OTU table
    writeLines(paste("Samples missing from the OTU table that are present in the",
                     "metadata:", paste(missing_genomic_data, collapse = ", ")))
    
    missing_env_data <- names(otu_table[-1])[!(names(otu_table[-1]) %in% sample_metadata$temporal_sample_id )] # Samples missing from env data
    writeLines(paste("Samples missing from the metadata that are present in the",
                     "OTU table:", paste(missing_env_data, collapse = ", ")))
    writeLines(paste0("Now filtering out missing sample(s)..."))
    
    # Filter out the samples that don't occur in both
    otu_table <- otu_table %>% select(genome, !all_of(missing_env_data))
    
    sample_metadata <- sample_metadata %>% 
        filter((!!as.name(sample_id_column) %in% names(otu_table[-1])))
    
    # After filtering samples, remove any columns which are entirely NAs (or 0s for MAGs)
    sample_metadata <- sample_metadata[,colSums(is.na(sample_metadata))<nrow(sample_metadata)]
    
    # Identify env_data columns that are all NA
    all_na_cols <- names(sample_metadata[,colSums(is.na(sample_metadata))==nrow(sample_metadata)])
    # filter out metadata variables that are all NAs
    sample_metadata <- sample_metadata[,colSums(is.na(sample_metadata))<nrow(sample_metadata)]

    
    # Reorder otu table columns to match order in sample metadata
    otu_table <- otu_table[,c("genome", sample_metadata$temporal_sample_id)]

    writeLines(paste(nrow(sample_metadata), "samples in sample_metadata"))
    writeLines(paste(ncol(otu_table[-1]), "samples in otu_table"))
    

    # Checking taxonomy, tree, and otu table (if taxonomy and tree are present)
    # and reordering to match tree order
    if(!is.null(tree) & !is.null(taxonomy)) {
        # reorder otu table to match order of tips in tree
        match.phylo.otu <- picante::match.phylo.data(tree, 
                                                     otu_table %>% as_tibble() %>%
                                                         column_to_rownames(var = "genome"))
        tree <- match.phylo.otu$phy
        otu_table <- match.phylo.otu$data %>%
            rownames_to_column(var = "genome") %>%
            mutate(rowname = genome) %>%
            column_to_rownames()
        
        # Check that otu table genomes are all in the tree and taxonomy file
        genomes_missing_from_taxonomy <- otu_table$genome[!otu_table$genome %in% taxonomy$genome]
        writeLines(paste("Genomes missing from the taxonomy that are present in the",
                         "OTU table:", paste(genomes_missing_from_taxonomy, collapse = ", ")))
        
        genomes_missing_from_tree <- otu_table$genome[!otu_table$genome %in% tree$tip.label]
        writeLines(paste("Genomes missing from the tree that are present in the",
                         "OTU table:", paste(genomes_missing_from_tree, collapse = ", ")))
        
        writeLines(paste0(nrow(otu_table), " taxa in otu table"))
        writeLines(paste0(length(tree$tip.label), " tips in tree"))
        writeLines(paste0(nrow(taxonomy), " taxa in taxonomy"))
                                                   
        # reorder taxonomy to match otu_table and tree tip order
        taxonomy <- taxonomy %>% as_tibble() %>% column_to_rownames("genome")
        taxonomy <- taxonomy[rownames(match.phylo.otu$data),] %>%
            rownames_to_column(var = "genome") %>%
            mutate(rowname = genome) %>%
            column_to_rownames()

    }

    
    # Create list of outputs to return:
    return_list <- list(sample_metadata = sample_metadata, otu_table = otu_table)
    
    if(!is.null(taxonomy)) {
        return_list$taxonomy <- taxonomy
    }
    
    if(!is.null(tree)) {
        return_list$tree <- tree
    }
    
    return(return_list)
}

# Output from TranscriptM rawcount_tpm_table.txt concatenated with extra column SampleID__
# Prefiltered to remove genes with no mapped transcripts
read_dram_annotations <- function(filename = "all_annotations.tsv") {
  annotations <- read_tsv(here(v4_dram_annotations_directory, filename)) %>%
    rename(gene_id = 1, genome = fasta) %>%
    select(genome, gene_id, ko_id, cazy_id)
}

read_metaT_data <- function(metadata, annotations = read_dram_annotations()) {
  d <- read_tsv(here(v5_emerge_metaT_directory, "compiled_tpm_tables.tsv")) %>%
    select(gene_id, contig, start, end, strand, length, gene, product, raw_count, tpm, SampleID__) %>%
    inner_join(metadata %>% select(SampleID__, temporal_sample_id)) %>%
    # Fix gene_id to match DRAM annotations gene_id
    mutate(gene_id = map_chr(gene_id, str_replace_all,
           pattern = "-", replacement = "_")) %>%
    inner_join(annotations)

  return(d)
}


read_reaction_definitions <- function() {
  d <- read_tsv(
      here("data", "DRAM_product_refined_resources", "custom_input_modules", "EMERGE_pathways_module.tsv")
      ) %>%
    mutate(
      reaction = map_chr(module_name, str_extract, "(?<=number:)[:digit:]+"),
      ko_id = str_split(definition, "\\(|\\)|\\+|,")
      ) %>%
    unnest(ko_id) %>%
    filter(ko_id != "") %>%
    select(pathway = complex, reaction, ko_id)

  return(d)
}

read_pathway_definitions <- function(reaction_definitions) {
  d <- read_tsv(here("data", "DRAM_product_refined_resources", "pathways_refined.tsv")) %>%
    mutate(
      reaction = str_split(reaction, "\\+|,")
    ) %>%
    unnest(reaction) %>%
    select(pathway, subpathway, reaction) %>%
    unique() %>%
    inner_join(reaction_definitions, by = c("pathway", "reaction"), multiple = "all") %>%
    mutate(subpathway = map_chr(subpathway, str_replace, "dependent$", "dep"))

  return(d)
}

read_cazy_definitions <- function() {
  d <- read_tsv(here("data", "DRAM_CAZy_definitions.tsv")) %>%
    mutate(
      subpathway = map_chr(function_name, str_replace_all,
        pattern = "-|[:space:]|\\(|\\)", replacement = "_"
        ),
      reaction = map2_chr(function_name, long_function_name,
        ~ .y %>%
            str_replace(pattern = .x, replacement = "") %>%
            str_trim() %>%
            str_replace_all(pattern = "-|[:space:]|\\(|\\)", replacement = "_")
        )
      ) %>%
    select(pathway = category, subpathway, reaction, cazy_id = function_ids) %>%
    separate_rows(cazy_id, sep = ", ") %>%
    distinct()

  return(d)
}

apply_pathway_labels_by_ko <- function(df, pathway_definitions, cazy_definitions) {
  labels <- bind_rows(
      pathway_definitions %>% rename(label_id = ko_id),
      cazy_definitions %>% rename(label_id = cazy_id)
    ) %>%
    unique()

  d <- df %>%
    mutate(
      label_id = map2_chr(ko_id, cazy_id,
        ~ c(
            str_split(.x, ",")[[1]],
            str_split(.y, "; ")[[1]] %>% str_remove("_.*") %>% unique()
            ) %>%
            str_replace_na(replacement = "") %>%
            str_c(sep = ",", collapse = ",")
      )) %>%
    separate_rows(label_id, sep = ",") %>%
    filter(label_id != "") %>%
    left_join(labels, by = "label_id") %>%
    filter(!is.na(pathway))

  return(d)
}

read_metaT_genomes <- function() {
  d <- read_tsv(here(v7_metaT_processed_directory, "metaT_genomes.tsv"))

  return(d)
}

read_metaT_pathways <- function() {
  d <- read_tsv(here(v7_metaT_processed_directory, "metaT_pathways.tsv"))

  return(d)
}

read_metaT_exp <- function() {
  d <- read_csv(here(v7_metaT_processed_directory, "metaT_pathway_exp.csv")) %>%
    select(-pathway) %>%
    complete(SampleID__, genome, subpathway, fill = list(tpm = 0))

  return(d)
}

read_metaT_react <- function() {
  d <- read_csv(here(v7_metaT_processed_directory, "metaT_reaction_exp.csv")) %>%
    complete(SampleID__, genome, pathway, reaction, fill = list(tpm = 0))

  return(d)
}

read_emerge_distillate <- function() {
    d <- read_tsv(here(v9_emerge_distillate_directory, "EMERGE_20220905_distillate.tsv")) %>%
        pivot_longer(!c(gene_id, gene_description, module, header, subheader), names_to = "genome", values_to = "call") %>%
        rename(reaction = module)

    return(d)
}

read_emerge_product <- function(memory_option = memory_saving) {
    if(memory_option) {
      # use data.table functions instead
      dr <- fread(file = here(v9_emerge_distillate_directory, "EMERGE_20220905_product.tsv"),
                  sep = "\t") %>%
        melt(., id.vars = c("genome"), measure.vars = 2:226, variable.name = "complex_info",
             value.name = "call") %>%
        setorder(., cols = "genome", "complex_info")
      dr <- dr[, call := ifelse(call >= 0.6, T, F)]
      dr <- dr[, pathway := gsub(x = complex_info, pattern = "(Complex )(.+)(: .+)", replacement = "\\2")]
      dr <- dr[, reaction := gsub(x = complex_info, pattern = "(Complex )(\\w+\\S+: )(.+)( number:)(.+)", replacement = "\\5")]
      d <- dr %>% as_tibble() # convert back to tibble
    } else { # original pattern extraction
      d <- read_tsv(here(v9_emerge_distillate_directory, "EMERGE_20220905_product.tsv")) %>%
        pivot_longer(-genome, names_to = "complex_info", values_to = "call") %>%
        mutate(
            call = map_lgl(call, ~ . >= 0.6),
            pathway = map_chr(complex_info, str_extract, pattern = "(?<=Complex ).+(?=: )"),
            reaction = map_chr(complex_info, str_extract, pattern = "(?<=number:)[:digit:]+"),
            )
    }

    return(d)
}

read_emerge_product_refined <- function() {
    d <- read_tsv(here(v9_emerge_distillate_directory, "EMERGE_20220905_product_refined.tsv")) %>%
        pivot_longer(-genome, names_to = "subpathway", values_to = "call")

    return(d)
}

read_dram_distillate <- function() {
    d <- bind_rows(
      readxl::read_excel(here(v2_dram_distillate_directory, "metabolism_summary.xlsx"), sheet = "MISC") %>% mutate(sheet = "MISC"),
      readxl::read_excel(here(v2_dram_distillate_directory, "metabolism_summary.xlsx"), sheet = "carbon utilization") %>% mutate(sheet = "carbon utilization"),
      readxl::read_excel(here(v2_dram_distillate_directory, "metabolism_summary.xlsx"), sheet = "Transporters") %>% mutate(sheet = "Transporters"),
      readxl::read_excel(here(v2_dram_distillate_directory, "metabolism_summary.xlsx"), sheet = "Energy") %>% mutate(sheet = "Energy"),
      readxl::read_excel(here(v2_dram_distillate_directory, "metabolism_summary.xlsx"), sheet = "Organic Nitrogen") %>% mutate(sheet = "Organic Nitrogen"),
      readxl::read_excel(here(v2_dram_distillate_directory, "metabolism_summary.xlsx"), sheet = "carbon utilization (Woodcroft)") %>% mutate(sheet = "carbon utilization (Woodcroft)")
      ) %>%
      pivot_longer(!c(gene_id, gene_description, module, header, subheader, sheet), names_to = "genome", values_to = "call") %>%
      rename(reaction = module) %>%
      filter(genome != "fasta")

    return(d)
}

read_dram_product <- function(definitions, distillate, memory_option = memory_saving) {
  if(memory_option){
    # use data.table functions instead
    suppressWarnings(dr <- fread(file = here(v2_dram_distillate_directory, "product.tsv"),
                sep = "\t") %>%
      melt(., id.vars = c("genome"), measure.vars = 2:99, variable.name = "complex_info",
           value.name = "call") %>% # note CAZymes are set to TRUE/FALSE and are coerced to 1/0 in this operation; throws a warning 
      setorder(., cols = "genome", "complex_info")) # match order with what pivot_longer would provide
    dr <- dr[, call := ifelse(call >= 0.6, T, F)]
    dr <- dr[, c("pathway", "reaction") := tstrsplit(complex_info, ": ", fixed = FALSE)] # fastest
    dr <- dr[, pathway := ifelse(is.na(reaction), NA, pathway)] # fix the NAs from tstrsplit
    d <- dr %>% as_tibble() # convert back to tibble
  } else {
    # original version
    d <- read_tsv(here(v2_dram_distillate_directory, "product.tsv")) %>%
      pivot_longer(-genome, names_to = "complex_info", values_to = "call") %>%
      mutate(
        call = map_lgl(call, ~ . >= 0.6),
        pathway = map_chr(complex_info, str_extract, pattern = "^.+(?=: )"),
        reaction = map_chr(complex_info, str_extract, pattern = "(?<=: ).+$"),
      )
  }

    # Remove positive calls containing no unique (signature) cazy_ids to the subpathway
    # Conservative approach to remove high levels of pathway overlap for some cazy_ids
    # Also mark families with <50% matching characterised proteins as non-signature
    # Add CE4 as a heterogenous family with chitin, peptidoglycan and unknown
    low_matching_families <- c("GH63", "GH49", "GH127", "CE2", "GH23", "PL30", "PL8", "PL33", "GH88", "GH67", "GH98", "CE4")
    cazy_signature <- definitions %>%
      select(subpathway, cazy_id) %>%
      distinct() %>%
      group_by(cazy_id) %>%
      count() %>%
      filter(n == 1, !cazy_id %in% low_matching_families)
    
    if(memory_option) {
      # Modified to use data.table()
      dr1 <- dr[pathway == "CAZy",] # select only CAZy pathways
      dr1 <- dr1[call == TRUE,] # select only called pathways
      dr1 <- dr1[, call := NULL] # remove call
      dr1 <- dr1[, subpathway := gsub("-|[[:space:]]|\\(|\\)", "_", reaction)] # replace -, spaces, and parentheses with _
      # Merge cazyme pathway definitions with table
      definitions_dtbl <- definitions %>% select(pathway, subpathway, cazy_id) %>% as.data.table()
      setkey(definitions_dtbl, pathway, subpathway)
      setkey(dr1, pathway, subpathway)
      dr1 <- dr1[definitions_dtbl, allow.cartesian = T]
      # add a column for signature genes
      dr1 <- dr1[, signature := ifelse(cazy_id %in% cazy_signature$cazy_id, T, F)]
      # Merge distillate
      distillate_dtbl <- distillate %>% select(genome, cazy_id = gene_id, call) %>% as.data.table()
      dr1 <- merge(dr1, distillate_dtbl, by = c("genome", "cazy_id"), all.x = T)
      dr1 <- dr1[call > 0,]
      # Summarie cazy calls by presence of signature
      dr1 <- dr1[, .(call_cazy = any(signature)), by = list(genome, reaction)]
      overlap_correction <- dr1 %>% as_tibble() # convert back to tibble
      
      # Clean up memory
      rm(definitions_dtbl)
      rm(distillate_dtbl)
      rm(dr1) # remove dr1 to save space
      rm(dr) # remove dr to save space
      gc() # call gc() to mop up
    } else {
      # original version
      overlap_correction <- d %>%
        filter(call, pathway == "CAZy") %>%
        select(-call) %>% 
        mutate(
          subpathway = map_chr(reaction, str_replace_all, pattern = "-|[:space:]|\\(|\\)", replacement = "_")
        ) %>%
        left_join(definitions %>% select(pathway, subpathway, cazy_id), multiple = "all") %>%
        mutate(signature = map_lgl(cazy_id, ~ . %in% cazy_signature$cazy_id)) %>%
        left_join(distillate %>% select(genome, cazy_id = gene_id, call)) %>%
        filter(call > 0) %>% 
        group_by(genome, reaction) %>%
        summarise(call_cazy = any(signature))
    }

    d <- d %>%
      left_join(overlap_correction) %>%
      mutate(
        call = coalesce(call_cazy, call)
      ) %>%
      select(-call_cazy)

    return(d)
}

read_manual_methanogen_calls <- function() {
  d <- readxl::read_excel(here(v1_manual_methanogen_directory, "19May22_DRAM_and_metat_substrate_analysis_of_mgens_with_Bathy.xlsx"), sheet = "import") %>%
    rename(genome = Bin, taxonomy = Taxonomy, verbose_call = call) %>%
    mutate(
      hydrogenotroph = map_lgl(verbose_call, str_detect, pattern = "Hydrogenotroph"),
      acetoclast = map_lgl(verbose_call, str_detect, pattern = "Acetoclastic"),
      methylotroph = map_lgl(verbose_call, str_detect, pattern = "Methylotroph"),
      hydrogen_dep = map2_lgl(methylotroph, hydrogenotroph, ~ .x & !.y),
      # Note: unclear if Methanobacterium genomes are H2-dependent (literature descriptions) or H2-independent (as indicated by genome calls)
      hydrogen_indep = map2_lgl(methylotroph, hydrogenotroph, ~ .x & .y),
    ) %>%
    select(genome, hydrogenotroph, acetoclast, hydrogen_dep, hydrogen_indep) %>%
    rename(
      `hydrogenotrophic_methanogenesis-all` = hydrogenotroph,
      `acetoclastic_methanogenesis-all` = acetoclast,
      `methylotrophic_methanogenesis-h2_dep` = hydrogen_dep,
      `methylotrophic_methanogenesis-h2_indep` = hydrogen_indep
    ) %>%
    pivot_longer(-genome, names_to = "subpathway", values_to = "call")

  return(d)
}

combine_distillates <- function(emerge, dram) {
  d <- bind_rows(
    emerge,
    dram %>% select(-sheet)
    )

  return(d)
}

combine_products <- function(emerge, dram, memory_option = memory_saving) {
  if(memory_option) {
    # modify dram distillate to merge with emerge
    dr <- as.data.table(dram)
    dr <- dr[header=="CAZY",]
    dr <- dr[, call := ifelse(call > 0, T, F)]
    # rename some columns 
    dr[, reaction := gene_id]
    dr[, pathway := header]
    dr <- dr[, complex_info := paste("Complex ", pathway, reaction, sep = ": ")]
    dr1 <- unique(dr[, c("gene_id", "gene_description", "header", "subheader", "sheet"):= NULL])
    dr1 <- as_tibble(dr1)
    
    d <- bind_rows(emerge, dr1)
    
    # clean up memory
    rm(dr)
    rm(dr1)
    gc()
  } else {
    # original version
    d <- bind_rows(
      emerge,
      dram %>%
        filter(header == "CAZY") %>%
        mutate(
          call = map_lgl(call, ~ . > 0),
          reaction = gene_id,
          pathway = header,
          complex_info = map2_chr(pathway, reaction, ~ str_c("Complex ", .x, .y, sep = ": "))
        ) %>%
        distinct(genome, complex_info, call, pathway, reaction)
    )
  }
  return(d)
}

combine_products_refined <- function(emerge, dram, manual_methanogen, clusters = genome_clusters) {
  d <- bind_rows(
    emerge %>%
      filter(!str_detect(subpathway, "methanogenesis"), subpathway != "acetogenesis-all"),
    dram %>%
      filter(pathway == "CAZy") %>%
      mutate(
        reaction = map_chr(reaction, str_replace_all, pattern = "-|[:space:]|\\(|\\)", replacement = "_"),
        subpathway = map2_chr(pathway, reaction, str_c, sep = "-")
      ) %>%
      select(genome, subpathway, call),
    manual_methanogen %>%
      filter(genome %in% clusters$representative)
  )

  return(d)
}

calculate_specialisations <- function(product_refined, genome_clusters) {
  specials <- function(pathways) {
    case_when(
      "hydrogenotrophic_methanogenesis-all" %in% pathways ~ "methanogen",
      "acetoclastic_methanogenesis-all" %in% pathways ~ "methanogen",
      "methylotrophic_methanogenesis-h2_dep" %in% pathways ~ "methanogen",
      "methylotrophic_methanogenesis-h2_indep" %in% pathways ~ "methanogen",
      "carbon_redox-dissimilatory_methanotrophy" %in% pathways ~ "methanotroph",
      "Wood_Ljungdahl-acetogen" %in% pathways ~ "homoacetogen",
      sum(startsWith(pathways, "CAZy")) >= 3 & sum(str_detect(pathways, "degradation")) >= 3 ~ "generalist",
      any(str_detect(pathways, "fermentation")) ~ "fermenter",
      any(startsWith(pathways, "CAZy")) ~ "macromolecule_degrader",
      any(str_detect(pathways, "degradation")) ~ "monomer_degrader",
      TRUE ~ NA_character_
    )
  }

  nitrogen_specials <- function(pathways) {
    case_when(
      "nitrogen_redox-comammox" %in% pathways ~ "comammox",
      "nitrogen_redox-ammonia_oxidation" %in% pathways ~ "ammonia_oxidiser",
      "nitrogen_redox-denitrification" %in% pathways ~ "denitrifier",
      "nitrogen_redox-nitrate_reduction" %in% pathways ~ "nitrate_reducer",
      "nitrogen_redox-nitrite_oxidation" %in% pathways ~ "nitrite_oxidiser",
      "nitrogen_redox-nitrogen_fixation" %in% pathways ~ "nitrogen_fixer",
      TRUE ~ "None detected"
    )
  }

  sulfur_specials <- function(pathways) {
    case_when(
      "sulfur_redox-dissimilatory_sulfate_reduction" %in% pathways ~ "sulfate_reducer",
      "sulfur_redox-dissimilatory_sulfur_oxidation" %in% pathways ~ "sulfur_oxidiser",
      TRUE ~ NA_character_
    )
  }

  pathways <- product_refined %>%
    filter(genome %in% genome_clusters$representative)

  pathways %>%
    filter(call) %>%
    group_by(genome) %>%
    summarise(
      specialisation = specials(subpathway),
      nitrogen_specialisation = nitrogen_specials(subpathway),
      sulfur_specialisation = sulfur_specials(subpathway)
      ) %>%
    full_join(pathways %>% distinct(genome) %>% select(genome))
}

read_kegg_brite <- function() {
    d <- fromJSON(
        here("data", "ko00001.json"),
        simplifyVector = FALSE
    ) %>%
    as_tibble() %>%
    unnest_wider(children, names_repair = "unique") %>%
    unnest(children) %>%
    unnest_wider(children, names_repair = "unique") %>%
    unnest(children) %>%
    unnest_wider(children, names_repair = "unique") %>%
    unnest(children) %>%
    unnest(children) %>%
    unnest(children) %>%
    select(-name...1) %>%
    rename(group_A = name...2, group_B = name...3, group_C = name...4) %>%
    separate(children, into = c("kegg_id", "kegg_description"), sep = "  ")

    return(d)
}

get_metapathway_groups <- function() {
  d <- tribble(
    ~"subpathway", ~"metapath", ~"subpathway_label",
    # CAZy
    "CAZy-Alpha_galactans",                 "cazymes", "a-galactans",
    "CAZy-Alpha_mannan",                    "cazymes", "a-mannan",
    "CAZy-Amorphous_Cellulose",             "cazymes", "amorphous cellulose",
    "CAZy-Arabinan",                        "cazymes", "arabinan",
    "CAZy-Arabinose_cleavage",              "cazymes", "arabinose cleavage",
    "CAZy-Beta_galactan__pectic_galactan_", "cazymes", "b-galactan",
    "CAZy-Beta_mannan",                     "cazymes", "b-mannan",
    "CAZy-Chitin",                          "cazymes", "chitin",
    "CAZy-Crystalline_Cellulose",           "cazymes", "crystalline cellulose",
    "CAZy-Fucose_Cleavage",                 "cazymes", "fucose cleavage",
    "CAZy-Mixed_Linkage_glucans",           "cazymes", "mixed",
    "CAZy-Mucin",                           "cazymes", "mucin",
    "CAZy-Pectin",                          "cazymes", "pectin",
    "CAZy-Polyphenolics",                   "cazymes", "polyphenolics",
    "CAZy-Rhamnose_cleavage",               "cazymes", "rhamnose cleavage",
    "CAZy-Starch",                          "cazymes", "starch",
    "CAZy-Sulf_Polysachharides",            "cazymes", "sulf polysaccharides",
    "CAZy-Xylans",                          "cazymes", "xylans",
    "CAZy-Xyloglucan",                      "cazymes", "xyloglucan",
    # Carbon pathways
    "xylose_degradation-isomerase_pathway",          "C degradation", "xylose isomerase",
    "xylose_degradation-oxidoreductase_pathway",     "C degradation", "xylose oxidoreductase",
    "xylose_degradation-xylonate_hydratase_pathway", "C degradation", "xylose xylonate hydratase",
    "xylose_degradation-weimburg_dahms",             "C degradation", "xylose weimburg-dahms",
    "lactose_degradation-all",                       "C degradation", "lactose",
    "galactose_degradation-all",                     "C degradation", "galactose",
    "sucrose_degradation-all",                       "C degradation", "sucrose",
    "fructose_degradation-all",                      "C degradation", "fructose",
    "mannose_degradation-all",                       "C degradation", "mannose",
    "fucose_degradation-all",                        "C degradation", "fucose",
    "galacturonic_acid_degradation-all",             "C degradation", "galacturonic acid",
    "glycerol_degradation-glycerol",                 "C degradation", "glycerol",
    "glycerol_degradation-glycerone",                "C degradation", "glycerone",
    "acetogenesis-simple",                           "C fermentation", "acetogenesis",
    "ethanol_fermentation-all",                      "C fermentation", "ethanol",
    "lactate_fermentation-all",                      "C fermentation", "lactate",
    "propionate_fermentation-succinate_pathway",     "C fermentation", "propionate succinate",
    "propionate_fermentation-acrylate_pathway",      "C fermentation", "propionate acrylate",
    "butanoate_fermentation-succinate_pathway",      "C fermentation", "butanoate succinate",
    "butanoate_fermentation-acetyl_coa_pathway",     "C fermentation", "butanoate acetyl-CoA",
    "Calvin_Benson-all",                             "C fixation", "Calvin-Benson",
    "rTCA-all",                                      "C fixation", "rTCA",
    "Wood_Ljungdahl-acetogen",                       "C fixation", "Wood-Ljungdahl ace",
    "Wood_Ljungdahl-methanogen",                     "C fixation", "Wood-Ljungdahl met",
    "Fuchs_Holo-all",                                "C fixation", "Fuchs-Holo",
    "HPHB_DCHB-DCHB",                                "C fixation", "DCHB",
    "HPHB_DCHB-HPHB",                                "C fixation", "HPHB",
    "reductive_glycine-serine_pathway",              "C fixation", "rGlycine serine",
    "reductive_glycine-acetyl_pathway",              "C fixation", "rGlycine acetyl",
    # Methane pathways
    "hydrogenotrophic_methanogenesis-all",      "methanogenesis", "hydrogenotrophic",
    "acetoclastic_methanogenesis-all",          "methanogenesis", "acetoclastic",
    "methylotrophic_methanogenesis-h2_dep",     "methanogenesis", "H2 dep methylotrophic",
    "methylotrophic_methanogenesis-h2_indep",   "methanogenesis", "H2 indep methylotrophic",
    "carbon_redox-assimilatory_methanotrophy",  "methanotrophy", "assim methanotrophy",
    "carbon_redox-dissimilatory_methanotrophy", "methanotrophy", "methanotrophy",
    # Other
    "chemotaxis-all",        "motility", "chemotaxis",
    "flagella_assembly-all", "motility", "flagella",
    # Nitrogen redox pathways
    "nitrogen_redox-ammonia_oxidation",              "nitrogen", "ammonia ox",
    "nitrogen_redox-assimilatory_nitrate_reduction", "nitrogen", "assim nitrate red",
    "nitrogen_redox-denitrification",                "nitrogen", "denitrification",
    "nitrogen_redox-nitrate_reduction",              "nitrogen", "nitrate red",
    "nitrogen_redox-nitric_oxide_oxidation",         "nitrogen", "nitric oxide ox",
    "nitrogen_redox-nitrite_oxidation",              "nitrogen", "nitrite ox",
    "nitrogen_redox-nitrogen_fixation",              "nitrogen", "nitrogen fixation",
    "nitrogen_redox-anammox",                        "nitrogen", "anammox",
    "nitrogen_redox-comammox",                       "nitrogen", "comammox",
    "urea_degradation-all",                          "nitrogen", "urea degradation",
    # Sulfur redox pathways
    "sulfur_redox-assimilatory_sulfate_reduction",  "sulfur", "assim sulfate red",
    "sulfur_redox-dissimilatory_sulfate_reduction", "sulfur", "sulfate red",
    "sulfur_redox-dissimilatory_sulfur_oxidation",  "sulfur", "sulfur ox",
    # Phosphorous redox pathways
    "phosphorous_redox-phosphite_oxidation",                "phosphorous", "phosphite ox",
    "phosphorous_redox-assimilatory_phosphonate_oxidation", "phosphorous", "assim phosphonate ox",
    "phosphorous_redox-anabolic_phosphate_reduction",       "phosphorous", "anab phosphate red"
  )

  return(d)
}

get_reaction_abundance <- function(product, rel_abund, samples, memory_option = memory_saving) {
  if(memory_option) {
    # Prepare data.tables 
    product_dtbl <- as.data.table(product)
    product_dtbl[, complex_info := NULL]
    product_dtbl <- product_dtbl[, subpathway := paste(pathway, reaction, sep = "-")]
    product_dtbl <- product_dtbl[call == T,]
    
    rel_abund_dtbl <- as.data.table(rel_abund)
    selected_cols <- c("genome", samples)
    rel_abund_dtbl <- rel_abund_dtbl[, ..selected_cols]
    rel_abund_dtbl <- melt(rel_abund_dtbl, id.vars = c("genome"), measure.vars = 2:length(selected_cols), 
                           variable.name = "sample", value.name = "rel_abund") 
    setkey(rel_abund_dtbl, genome, sample)
    setkey(product_dtbl, genome)
    # Join data.tables
    dr1 <- merge(rel_abund_dtbl, product_dtbl, by = c("genome"), all.x = T, allow.cartesian = T)
    d <- as_tibble(dr1)
    
    # Cleanup
    rm(product_dtbl)
    rm(rel_abund_dtbl)
    rm(dr1)
    gc()

  } else {
    #Original code
    d <- rel_abund %>%
      select(genome, all_of(samples)) %>%
      pivot_longer(-genome, names_to = "sample", values_to = "rel_abund") %>%
      left_join(
        product %>%
          select(-complex_info) %>%
          mutate(subpathway = map2_chr(pathway, reaction, str_c, sep = "-")),
        by = "genome", multiple = "all") %>%
      filter(call)
  }
  return(d)
}

get_pathway_abundance <- function(product_refined, rel_abund, samples) {
  d <- rel_abund %>%
    select(genome, all_of(samples)) %>%
    pivot_longer(-genome, names_to = "sample", values_to = "rel_abund") %>%
    left_join(product_refined, by = "genome", multiple = "all") %>%
    filter(call)

  return(d)
}

get_cumulative_abundance <- function(abund) {
  d <- abund %>%
    group_by(subpathway, sample) %>%
    summarise(abund = sum(rel_abund)) %>%
    pivot_wider(names_from = "subpathway", values_from = "abund")

  return(d)
}

get_metapathway_abundance <- function(pathway_abundance, metapathway_groups) {
  d <- pathway_abundance %>%
    left_join(metapathway_groups) %>%
    group_by(metapath, sample, genome) %>%
    summarise(rel_abund = first(rel_abund)) %>%
    group_by(metapath, sample) %>%
    summarise(rel_abund = sum(rel_abund))

  return(d)
}

# Read in network stats
get_network_stats <- function() {
  nw_stats_fp <- list.files(here("data", "network_stats"), pattern = "*.csv", 
                            full.names = T)
  
  tbl_colnames <- read_csv(nw_stats_fp[1], n_max = 0) %>%
    names()
  
  d <- lapply(nw_stats_fp, function(x) {
    read_csv(x, col_names = tbl_colnames, skip = 1, col_types = "ccddddi") %>% # setting column types to ensure correct parsing
      mutate(Habitat = gsub("^network_stats_(.{1,}).csv$","\\1", basename(x))) %>%
      mutate(Habitat = str_to_title(Habitat)) %>%
      mutate(Habitat = factor(Habitat, levels = c("Palsa", "Bog", "Fen"))) %>%
      mutate(Year = 2010 + layer) %>%
      filter(node_id != "node_id") %>%
      select(-node_id)
  })
  
  d <- d %>%
    reduce(bind_rows)
  
  return(d)
}

# Create an outlier and depth filter function (to remove samples deeper than 39.9 cm,
# and the temporal outlier samples)
filter_depth_and_outliers <- function(input_data = input_ra, 
                                      input_data_type = NULL,
                                      filter_outliers = TRUE,
                                      filter_deep = TRUE,
                                      how_deep = 40,
                                      verbose = FALSE) {
  # Function arguments
  # input_data ........ either input_ra, input_counts, or sample_metadata generated by
  #                     setup.R; this function will decide on the basis of length
  #                     whether or not this is a dataframe or a list of dataframes (
  #                     one of which is sample_metadata)
  # input_data_type ... default is unspecified; if set, either "dataframe" or 
  #                     "list of dataframes"; if "dataframe" a format of 
  #                     sample_metadata is assumed. if "list of dataframes"
  #                     a format of input_counts/input_ra is assumed
  # filter_outliers ... default is TRUE; if TRUE temporal.outliers in Outlier_status column
  #                     will be filtered out (otherwise they will be retained).
  # filter_deep ....... default is TRUE; if TRUE "deep" samples (DepthAvg__>39.999)
  #                     will be filtered out (otherwise they will be retained).
  # how_deep .......... default is 40 (DepthAvg >= 40 will be filtered out); optional
  #                     variant to change the depth at which filtering happens. (Included
  #                     for future convenience, but should not be changed for temporal
  #                     paper analysis)
  # verbose ........... default = FALSE; should the function report which samples are filtered?
  # Returns ........... filtered version of input_data. 
  
  # Decide if the input is input_ra, input_counts, or sample_metadata
  # Alternatively you can specify it manually with input_data_type
  if(length(input_data) > 4 | (!is.null(input_data_type) && input_data_type == "dataframe")) {
    sample_metadata_tmp <- input_data
  }
  
  if(length(input_data) <= 4 | (!is.null(input_data_type) && input_data_type == "list of dataframes")) {
    sample_metadata_tmp <- input_data$sample_metadata
  }
  
  # Filter outliers out
  if(filter_outliers) {
    # Check input data for relevant columns
    if(!("Outlier_status" %in% names(sample_metadata_tmp))) {
      stop("Outlier_status column not found in input_data")
    }
    
    if(verbose) {
      samples_to_remove <- sample_metadata_tmp %>%
        filter(Outlier_status %in% c("temporal.outlier"))
      
    }
    
    # Filter temporal.outliers out
    sample_metadata_tmp <- sample_metadata_tmp %>%
      filter(!(Outlier_status %in% c("temporal.outlier")))
    
    if(verbose) {
      writeLines(paste0("Filtered out outliers: ", paste(samples_to_remove$temporal_sample_id, collapse = ", ")))
    }
    
  }
  
  # Filter deep samples out
  if(filter_deep) {
    # Check input data for relevant columns
    if(!("DepthAvg__" %in% names(sample_metadata_tmp))) {
      stop("DepthAvg__ column not found in input_data")
    }
    if(verbose) {
      samples_to_remove <- sample_metadata_tmp %>%
        filter(DepthAvg__ >= how_deep)
      
    }
    
    # Filter samples with DepthAvg__ >39.999 out
    # This follows the same way of assessing depths
    # outlined by Virginia Rich in the basecamp post
    # https://3.basecamp.com/4869263/buckets/18955755/messages/4572257670
    sample_metadata_tmp <- sample_metadata_tmp %>%
      filter(DepthAvg__ < how_deep)
    
    if(verbose) {
      writeLines(paste0("Filtered out deep samples: ", paste(samples_to_remove$temporal_sample_id, collapse = ", ")))
    }
  }
  
  # Return a filtered version of the same input 
  if(length(input_data) > 4 | (!is.null(input_data_type) && input_data_type == "dataframe")) {
    return(sample_metadata_tmp)
  }
  
  if(length(input_data) <= 4 | (!is.null(input_data_type) && input_data_type == "list of dataframes")) {
    # If format is input_ra or input_counts invoke check_sample_metadata to 
    # match filtering
    
    if("taxonomy" %in% names(input_data)) {
      taxa <- input_data$taxonomy
    } else {
      taxa <- NULL
    }
    
    if("tree" %in% names(input_data)) {
      tree <- input_data$tree
    } else {
      tree <- NULL
    }
    
    if(verbose) {
      writeLines("Now generating filtered input object from filtered metadata...")
      input_data_filt <- check_sample_metadata(otu_table = input_data$otu_table, 
                                               sample_metadata = sample_metadata_tmp, 
                                               sample_id_column = "temporal_sample_id",
                                               taxonomy = taxa,
                                               tree = tree)
    } else {
      invisible(capture.output(input_data_filt <- check_sample_metadata(input_data$otu_table, 
                                                 sample_metadata_tmp, 
                                                 sample_id_column = "temporal_sample_id",
                                                 taxonomy = taxa,
                                                 tree = tree)))
      }
    
    return(input_data_filt)
  }
  
}

# Read in module membership
read_in_modules <- function() {
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
  
  return(module_list_comb)
}

# Read in vip members
read_vip_members <- function(vip_fp) {
  # Read in Module VIPs
  if(file.exists(here("data","vip_members_pathways-07-07-22.txt"))) {
    vip_members <- read_tsv(here("data","vip_members_pathways-07-07-22.txt"))
  } else {
    vip_members <- NULL
    writeLines("VIP member file not found, setting vip members to null")
  }
  
  return(vip_members)
}

transform_perc <- function(vec) {
  # Remove 0s for CLR transformation
  # See Smithson & Verkuilen 2006 (https://doi.org/10.1037/1082-989X.11.1.54)
  (vec * (length(vec) - 1) + 0.5) / length(vec)
}

read_aa_frequency <- function() {
  aa_columns <- c("genome", "aa", "freq")

  d <- bind_rows(
      read_tsv(here(v1_aa_frequencies_directory, "Emerge_MAGs_v1_aa_freq.tsv"), col_names = aa_columns) %>% mutate(genome_set = "emerge"),
      read_tsv(here(v1_aa_frequencies_directory, "gtdb_r207_aa_freq.tsv"), col_names = aa_columns) %>% mutate(genome_set = "gtdb")
    ) %>%
    pivot_wider(names_from = aa, values_from = freq, values_fill = 0) %>%
    pivot_longer(cols = -c(genome, genome_set), names_to = "aa", values_to = "freq")

  return(d)
}

read_persistent_modules <- function() {
  d <- read_tsv(here("data", "persistent_module_metabolism.txt")) %>%
    select(genome, Habitat__ = habitat)

  return(d)
}

read_mag_genome_sizes <- function(info, clusters) {
    d <- read_tsv(
          here("data", "MAG_genome_sizes.tsv"),
          col_names = c("genome", "genome_size")
          ) %>%
        full_join(info, by = "genome") %>%
        left_join(clusters, by = "genome") %>%
        mutate(
          completeness = completeness / 100,
          contamination = contamination / 100,
          adj_size = pmap_dbl(
            list(genome_size, completeness, contamination),
            ~ ..1 / ..2 / (1 + ..3)
            )
          ) %>%
        group_by(representative) %>%
        summarise(mean_adj_size = mean(adj_size))

    return(d)
}


calculate_alpha_diversity <- function(otu_tables) {
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

  d <- otu_tables %>%
    filter(gene %in% both_domain_genes) %>%
    select(temporal_sample_id, gene, sequence, coverage) %>%
    group_by(temporal_sample_id, gene) %>%
    nest() %>%
    mutate(nseq = map_dbl(data, nrow)) %>%
    filter(nseq >= 100) %>%
    mutate(
      sample = purrr::map(data,
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

  return(d)
}

alpha_vs_year_correlations <- function(singlem_alpha, sample_metadata) {
  # Habitat correlations
  d1 <- singlem_alpha %>%
    pivot_longer(-c(temporal_sample_id), names_to = "alpha_stat", values_to = "alpha_val") %>%
    left_join(sample_metadata %>% select(temporal_sample_id, Habitat__, Year__, DepthLumping)) %>%
    filter(!is.na(Habitat__)) %>%
    mutate(cat_Year = as.factor(Year__)) %>%
    group_by(Habitat__, alpha_stat) %>% 
    nest() %>%
    mutate(
      cont_year = purrr::map(data, ~ lm(alpha_val ~ Year__, data = .)),
      cat_year = purrr::map(data, ~ lm(alpha_val ~ cat_Year, data = .)),
    ) %>%
    pivot_longer(c(cont_year, cat_year), names_to = "reg_name", values_to = "reg") %>%
    mutate(
      r2 = map_dbl(reg, ~ summary(.)$r.squared %>% round(3)),
      p = map_dbl(reg, ~ anova(.)$`Pr(>F)`[1] %>% round(3)),
      label = str_c(r2, " (", p, ")")
    ) %>%
    mutate(cat = str_c(alpha_stat, reg_name, sep = "-")) %>%
    ungroup() %>%
    select(Habitat__, cat, label) %>%
    pivot_wider(names_from = cat, values_from = label) %>%
    arrange(Habitat__)

  # Habitat-depth correlations
  d2 <- singlem_alpha %>%
    pivot_longer(-c(temporal_sample_id), names_to = "alpha_stat", values_to = "alpha_val") %>%
    left_join(sample_metadata %>% select(temporal_sample_id, Habitat__, Year__, DepthLumping)) %>%
    filter(!is.na(Habitat__)) %>%
    mutate(cat_Year = as.factor(Year__)) %>%
    group_by(Habitat__, DepthLumping, alpha_stat) %>% 
    nest() %>%
    mutate(n_row = map_int(data, nrow)) %>%
    filter(n_row > 2) %>%
    mutate(
      cont_year = purrr::map(data, ~ lm(alpha_val ~ Year__, data = .)),
      cat_year = purrr::map(data, ~ lm(alpha_val ~ cat_Year, data = .)),
    ) %>%
    pivot_longer(c(cont_year, cat_year), names_to = "reg_name", values_to = "reg") %>%
    mutate(
      r2 = map_dbl(reg, ~ summary(.)$r.squared %>% round(3)),
      p = map_dbl(reg, ~ anova(.)$`Pr(>F)`[1] %>% round(3)),
      label = str_c(r2, " (", p, ")")
    ) %>%
    mutate(cat = str_c(alpha_stat, reg_name, sep = "-")) %>%
    ungroup() %>%
    select(Habitat__, DepthLumping, cat, label) %>%
    pivot_wider(names_from = cat, values_from = label) %>%
    arrange(Habitat__, DepthLumping)

  return(list(alpha_habitat = d1, alpha_habitat_depth = d2))
}




################################################################################################
################################################################################################
################################################################################################

genome_info = read_genome_info()
genome_info = read_gtdb2_taxonomy(genome_info = genome_info)
genome_clusters = read_genome_clusters()
genome_sizes <- read_mag_genome_sizes(info = genome_info, clusters = genome_clusters)
trimmed_mean = read_trimmed_mean()
raw_sample_metadata = read_sample_metadata()
read_counts = read_read_counts()
vip_members = read_vip_members()

sample_metadata <- clean_sample_metadata(sample_metadata = raw_sample_metadata)
trees <- read_genome_tree(otu_table = trimmed_mean$rel_abund,
                         visualize = FALSE)
taxonomy <- get_taxonomy(otu_table = trimmed_mean$rel_abund, 
                          genome_info = genome_info)

# Check sample metadata and combine into one output
input_ra <- check_sample_metadata(otu_table = trimmed_mean$rel_abund,
                      sample_metadata = sample_metadata,
                      taxonomy = taxonomy,
                      tree = trees$bac_arch_tree)

input_counts <- check_sample_metadata(otu_table = trimmed_mean$coverage,
                                  sample_metadata = sample_metadata,
                                  taxonomy = taxonomy,
                                  tree = trees$bac_arch_tree)

# Creates a *_all version that contains deep samples but not the outliers
sample_metadata_all <- filter_depth_and_outliers(input_data = input_ra$sample_metadata,
                                                 filter_outliers = TRUE,
                                                 filter_deep = FALSE)
input_ra_all <- filter_depth_and_outliers(input_data = input_ra,
                                          filter_outliers = TRUE,
                                          filter_deep = FALSE,
                                          verbose = TRUE)
input_counts_all <- filter_depth_and_outliers(input_data = input_counts,
                                              filter_outliers = TRUE,
                                              filter_deep = FALSE)

if(memory_saving){rm(sample_metadata_all, input_ra_all, input_counts_all)}

if (dir.exists(v4_singlem_directory)) {
  singlem_otu_tables <- read_singlem_otu_tables()
  singlem_alpha <- calculate_alpha_diversity(singlem_otu_tables)
  alpha_vs_year <- alpha_vs_year_correlations(singlem_alpha, sample_metadata)
}
singlem_appraise_summary <- read_singlem_appraise_summary()

# Contains neither deep samples nor outliers (default option)
sample_metadata_filt <- filter_depth_and_outliers(input_data = input_ra$sample_metadata)
input_ra_filt <- filter_depth_and_outliers(input_data = input_ra, verbose = TRUE)
input_counts_filt <- filter_depth_and_outliers(input_data = input_counts)

# Overwrite sample_metadata, input_counts, and input_ra with *_filt versions
sample_metadata <- sample_metadata_filt
input_counts <- input_counts_filt
input_ra <- input_ra_filt
if(memory_saving){rm(sample_metadata_filt, input_ra_filt, input_counts_filt)}

# Read in module membership
module_list_comb <- read_in_modules()
# Load metaT data
reaction_definitions <- read_reaction_definitions()
pathway_definitions <- read_pathway_definitions(reaction_definitions)
cazy_definitions <- read_cazy_definitions()

if (dir.exists(v5_emerge_metaT_directory) && dir.exists(v4_dram_annotations_directory)) {
  print("Run Metabolic-analysis/05_metaT_analysis/metaT_processing.R to reprocess new data")
}

if (dir.exists(v7_metaT_processed_directory)) {
  metaT_genomes <- read_metaT_genomes()
  metaT_pathways <- read_metaT_pathways()
  metaT_exp_all <- read_metaT_exp()
  metaT_react_all <- read_metaT_react()

  metaT_exp <- sample_metadata %>% select(SampleID__, temporal_sample_id) %>% inner_join(metaT_exp_all) %>% select(-SampleID__)
  metaT_react <- sample_metadata %>% select(SampleID__, temporal_sample_id) %>% inner_join(metaT_react_all) %>% select(-SampleID__)
  if(memory_saving){rm(metaT_exp_all, metaT_react_all)}
}

kegg_brite <- read_kegg_brite()
aa_frequency <- read_aa_frequency()
persistent_modules <- read_persistent_modules()

# Load metabolism data
if (dir.exists(v9_emerge_distillate_directory) && dir.exists(v2_dram_distillate_directory) && dir.exists(v1_manual_methanogen_directory)) {
  emerge_distillate <- read_emerge_distillate()
  emerge_product <- read_emerge_product()
  emerge_product_refined <- read_emerge_product_refined()
  dram_distillate <- read_dram_distillate()
  dram_product <- read_dram_product(definitions = cazy_definitions, distillate = dram_distillate)
  manual_methanogen_calls <- read_manual_methanogen_calls()

  # Distillate: individual gene calls
  distillate <- combine_distillates(emerge_distillate, dram_distillate)
  # Product: reaction calls (based on distillate gene calls)
  product <- combine_products(emerge_product, dram_distillate)
  if(memory_saving){rm(dram_distillate)}
  if(memory_saving){rm(emerge_product)}
  # Product refined: pathway calls (based on reaction calls and individual gene calls)
  product_refined <- combine_products_refined(emerge_product_refined, dram_product, manual_methanogen_calls)
  specialisations <- calculate_specialisations(product_refined, genome_clusters)

  metapathway_groups <- get_metapathway_groups()
  if(memory_saving){mem_use <- gc()}

  reaction_abundance <- get_reaction_abundance(product, trimmed_mean$rel_abund, sample_metadata$temporal_sample_id)
  if(memory_saving){rm(product)}
  if(memory_saving){mem_use <- gc()}
  reaction_abundance_cumu <- get_cumulative_abundance(reaction_abundance)
  pathway_abundance <- get_pathway_abundance(product_refined, trimmed_mean$rel_abund, sample_metadata$temporal_sample_id)
  pathway_abundance_cumu <- get_cumulative_abundance(pathway_abundance)

  metapath_abundance <- get_metapathway_abundance(pathway_abundance, metapathway_groups)
}

if(dir.exists(here("data", "network_stats"))) {
  network_stats <- get_network_stats()
}
if(memory_saving){mem_use <- gc()}
