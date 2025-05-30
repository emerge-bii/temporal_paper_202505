# This script consolidates and plots geochem data for microbiome "outlier" samples vs. everything else.

library(ggplot2)

# Sample Metadata Sheet version in this GitHub repository. 
# TODO: Update to newest version when/if the file below is replaced with the newest version.
samplemetadata <- read.csv("../data/coring_geochem_sequenced_samples_1.0.0.csv", stringsAsFactors=FALSE)

# Subset samples used in this analysis
samplemetadata <- samplemetadata[(samplemetadata$CoreGroup__=="MainAutochamber" & samplemetadata$Month__==7 & 
                                    samplemetadata$Year__>=2011 & samplemetadata$Year__<=2017), ]

# Factor Habitat__ and assign colors
possible_habitats <- c("Palsa", "Bog", "Fen")
samplemetadata$Habitat__ <- ordered(samplemetadata$Habitat__, levels=possible_habitats)
habitat_colors <- c('#703C1B', # palsa
                    '#058000', # sphag
                    '#0001FF')  # erio
names(habitat_colors) <- possible_habitats

# Combined pH column for comparing habitats
samplemetadata$pH__ <- ifelse(is.na(samplemetadata$pH_porewater), samplemetadata$pH_peat, samplemetadata$pH_porewater)

# Import "weird samples" list (downloaded from https://docs.google.com/spreadsheets/d/1f0xB5rd1QY1qxPPKVSXt0itwFomuIpFs6T4ZhLvZCUg/edit#gid=0)
microbiome_outlier_summary <- read.csv("misc_analyses/Sample Outlier Collection.csv", stringsAsFactors=FALSE, na.strings=c("NA", ""))

weird_samples_JGI <- gsub("_", " ", microbiome_outlier_summary$Sample.ID)
weird_samples_ACE <- as.character(na.omit(microbiome_outlier_summary$Sample.ID.2))

is_outlier <- (samplemetadata$metaG_JGI_NovaSeq__ %in% weird_samples_JGI) | (samplemetadata$metaG_ACE__ %in% weird_samples_ACE)
is_outlier[is.na(is_outlier)] <- FALSE

samplemetadata$is_outlier <- is_outlier

plot_compare_outliers <- function(data, xvar, yvar) {
  plot <- ggplot(data, aes_string(x=xvar, y=yvar, color="Habitat__", shape="is_outlier")) + theme_bw() + geom_point(size=3, alpha=0.6) + scale_color_manual(name="Habitat", values=habitat_colors) + scale_shape_manual(values=c(16, 17))
  
  if(grepl("^Depth", yvar)) {
    plot <- plot + scale_y_reverse()
  }
  
  # Show text for outlier samples-- the original weird samples list currently identifies them by metaG_JGI_NovaSeq__, which all of them have, so use this name
  plot <- plot + geom_text(data=subset(data, is_outlier), aes(label=metaG_JGI_NovaSeq__), vjust=-1, fontface='bold') #hjust=-0.1
  
  return(plot)
}

plot_compare_outliers(samplemetadata, "pH__", "DepthAvg__") + xlab('pH (combined porewater and peat)')
plot_compare_outliers(samplemetadata, "O2.percent", "DepthAvg__")
plot_compare_outliers(samplemetadata, "d13C_CH4__", "DepthAvg__")
plot_compare_outliers(samplemetadata, "(d13C_CO2__+1000)/(d13C_CH4__+1000)", "DepthAvg__") + xlab("alphaC")
plot_compare_outliers(samplemetadata, "CO2.mM__/CH4.mM__", "DepthAvg__") + xlab("CO2/CH4 concentration ratio")

# No outlier sample data for: BulkDensity, DOC.mM__, TN.mM__, Sulfate.uM__, Nitrate.uM__, Ammonia.uM__, Phosphate.uM__, C.percent__, N.percent__, d13C_peat__, d15N_peat__
# Too variable: Acetate.uM__ (and other VFAs)
