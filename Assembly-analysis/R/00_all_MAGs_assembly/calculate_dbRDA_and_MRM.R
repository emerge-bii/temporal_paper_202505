#' ## Interpret Assembly Analysis results
#' This step takes the output from 01_calc_betaNTI.R and 02_calc_rcbc.R and runs
#' initial figures and analyses on them.

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
library(cowplot)
library(here)
library(GGally)
library(ecodist)

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

#### Read in data 
#### ====================================================================== ####
# Read in OTU table and env data
input <- input_ra
input$otu_table <- input$otu_table[-1]

# Turn Habitat into a factor
input$sample_metadata <- input$sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, 
                            levels = c("Palsa", "Collapsed Palsa", "Bog", "Fen")))

# Read in BetaNTI and RCBC wide format matrices
# Read in betaNTI results
betaNTI <- read_csv(file = paste0(outputs.fp, "/weighted_bNTI.csv")) %>%
  rename(SampleID = 1) %>%
  column_to_rownames(var = "SampleID")

# Read in RCBC results
RCBC.raw <- read_csv(file = paste0(outputs.fp, "/rcbc_matrix.csv")) %>%
  rename(SampleNumber = 1) %>%
  column_to_rownames(var = "SampleNumber")

# Read in long format
betanull.lf <- read_csv(file = paste0(outputs.fp, "/betanull.lf.csv"))

# Do minimal reformatting to column names of RCBC
RCBC <- RCBC.raw[c(as.character(1:ncol(input$otu_table))),
                  c(as.character(1:ncol(input$otu_table)))]
RCBC <- RCBC.raw # Temporary until sample ids fixed
RCBC <- as.matrix(RCBC)
colnames(RCBC) <- names(input$otu_table)
rownames(RCBC) <- names(input$otu_table)

RCBC_rank <- RCBC
RCBC_rank[RCBC_rank >= 0.98] <- 3
RCBC_rank[RCBC_rank >= -0.98 & RCBC_rank < 0.98] <- 2
RCBC_rank[RCBC_rank < -0.98] <- 1

par(mfrow = c(1, 2))
image(as.matrix(RCBC_rank))
image(as.matrix(RCBC))

betaNTI_rank <- betaNTI
betaNTI_rank[betaNTI_rank >= 2] <- 3
betaNTI_rank[betaNTI_rank >= -2 & betaNTI_rank < 2] <- 2
betaNTI_rank[RCBC_rank < -2] <- 1

par(mfrow = c(1, 2))
image(as.matrix(betaNTI_rank))
image(as.matrix(betaNTI))



# Calculate and plot PCOAs of the data to get a sense
env_data_filt <- input_ra$sample_metadata %>%
  filter(Min_PF_dist_from_sample >= 0) %>%
  mutate(Year__ = as.numeric(as.character(Year__))) %>%
  #mutate(Habitat__ = as.numeric(as.factor(Habitat__))) %>%
  filter_at(.vars = c("DepthLumping", "DepthAvg__"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id, Habitat__, 
         matches(c("DepthLumping", "DepthAvg__")))

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))
RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                    env_data_filt$temporal_sample_id]))

# Figure out the transformations since neitehr betaNTI nor RCBC meet the 3 
# criteria for a distance metric. 
# 1. Symmetry; d(p,q) == d(q,p) Order of the comparison can't matter. <-- this is true of both measures
# 2. The distance two samples is 0 ONLY if they are identical <--- this is true of both measures
# 3. The distance can't be negative. The distance between all samples has to be greater than 0. <--- This one doesn't work; 
# To solve this problem: multiply negative values by -1, and just to make the distances match other distance metrics, normalize
# between 0 and 1. BUT the negative values are important because they tell us something about
# the assembly processes. SO, we will also keep track of which values used to be negative (0),
# and which have always been positive, (1)
# 
# BetaNTI first

betaNTI.filt.trans <- ifelse(betaNTI.filt<0, betaNTI.filt*-1, betaNTI.filt) 
betaNTI.filt.norm <- (betaNTI.filt.trans - min(betaNTI.filt.trans))/(max(betaNTI.filt.trans) - min(betaNTI.filt.trans))
betaNTI.filt.dum <- ifelse(betaNTI.filt<0, 0, 1) 

# Plot to test
betaNTI.filt.trans.lf <- reshape2::melt(betaNTI.filt.trans)
betaNTI.filt.norm.lf <- reshape2::melt(betaNTI.filt.norm)
betaNTI.filt.lf <- reshape2::melt(betaNTI.filt)
betaNTI.filt.dum.lf <- reshape2::melt(betaNTI.filt.dum)
# get x intercepts for significance cutoffs
norm_int <- (2 - min(betaNTI.filt.trans))/(max(betaNTI.filt.trans) - min(betaNTI.filt.trans))
orig_int <- c(2, -2)
trans_int <- 2
bnti_int_data <- data.frame(bnti_type = c("bnti.norm", "bnti.orig", "bnti.orig", "bnti.trans"),
                       intercept = c(norm_int, orig_int, trans_int))

betanti_trans_plot <- left_join(betaNTI.filt.norm.lf, betaNTI.filt.trans.lf, by = c("Var1", "Var2")) %>%
  rename(bnti.norm = value.x, bnti.trans = value.y) %>%
  left_join(betaNTI.filt.lf, by = c("Var1", "Var2")) %>%
  rename(bnti.orig = value) %>%
  left_join(betaNTI.filt.dum.lf, by = c("Var1", "Var2")) %>%
  mutate(VariableHomogenous = ifelse(value == 1, "Variable", "Homogenous"),
         VariableHomogenous = ifelse(Var1 == Var2, "Neither", VariableHomogenous),
         VariableHomogenous = ifelse(bnti.orig > -2 & bnti.orig < 2, "Neither", VariableHomogenous)) %>%
  pivot_longer(contains("bnti"), names_to = "bnti_type", values_to = "bnti_measure")


ggplot(data = betanti_trans_plot, aes(x = bnti_measure)) +
  geom_histogram(aes(fill = VariableHomogenous), bins = 100, alpha = 0.5) +
  facet_wrap(~bnti_type, scale = "free_x") +
  geom_vline(data = bnti_int_data, aes(xintercept = intercept),
             color = "red", linetype = "dashed")
  

# RCBC
RCBC.filt.trans <- ifelse(RCBC.filt<0, RCBC.filt*-1, RCBC.filt) 
RCBC.filt.norm <- (RCBC.filt.trans - min(RCBC.filt.trans))/(max(RCBC.filt.trans) - min(RCBC.filt.trans))
RCBC.filt.dum <- ifelse(RCBC.filt<0, 0, 1) 

# Plot to test
# get x intercepts for significance cutoffs
norm_int <- (0.95 - min(RCBC.filt.trans))/(max(RCBC.filt.trans) - min(RCBC.filt.trans))
orig_int <- c(0.95, -0.95)
trans_int <- 0.95
int_data <- data.frame(rcbc_type = c("rcbc.norm", "rcbc.orig", "rcbc.orig", "rcbc.trans"),
                       intercept = c(norm_int, orig_int, trans_int))

RCBC_trans_plot <- left_join(RCBC.filt.norm.lf, RCBC.filt.trans.lf, by = c("Var1", "Var2")) %>%
  rename(rcbc.norm = value.x, rcbc.trans = value.y) %>%
  left_join(RCBC.filt.lf, by = c("Var1", "Var2")) %>%
  rename(rcbc.orig = value) %>%
  left_join(RCBC.filt.dum.lf, by = c("Var1", "Var2")) %>%
  mutate(VariableHomogenous = ifelse(value == 1, "Variable", "Homogenous"),
         VariableHomogenous = ifelse(Var1 == Var2, "Neither", VariableHomogenous),
         VariableHomogenous = ifelse(rcbc.orig > -0.95 & rcbc.orig < 0.95, "Neither", VariableHomogenous)) %>%
  pivot_longer(contains("rcbc"), names_to = "rcbc_type", values_to = "rcbc_measure")


ggplot(data = RCBC_trans_plot, aes(x = rcbc_measure)) +
  geom_histogram(aes(fill = VariableHomogenous), bins = 100, alpha = 0.5) +
  facet_wrap(~rcbc_type, scale = "free_x") +
  geom_vline(data = int_data, aes(xintercept = intercept),
             color = "red", linetype = "dashed")

# Run a pcoa on the new distance metric
betaNTI.pcoa <- cmdscale(as.dist(betaNTI.filt.norm), k = 2, eig = TRUE, add = TRUE)
betaNTI.pct_ex <- round((betaNTI.pcoa$eig/sum(betaNTI.pcoa$eig)) * 100, 1)
betaNTI.pcoa.mapdata <- betaNTI.pcoa$points %>% data.frame() %>%
  rownames_to_column(var = "temporal_sample_id") %>%
  dplyr::rename(PCOA1 = X1, PCOA2 = X2) %>%
  left_join(input_ra$sample_metadata, by = "temporal_sample_id") %>%
  dplyr::group_by(Habitat__) %>%
  dplyr::mutate(PCOA1HabAvg = mean(PCOA1),
                PCOA2HabAvg = mean(PCOA2)) %>%
  ungroup() %>%
  dplyr::select(temporal_sample_id, starts_with("PCOA1"), starts_with("PCOA2"), PCOA1, PCOA2, everything())
hab_bnti.pcoa.chulls <- plyr::ddply(betaNTI.pcoa.mapdata, ~ Habitat__, 
                                    function(df) df[chull(df$PCOA1, df$PCOA2),])

bnti.pcoa.plot <- ggplot(betaNTI.pcoa.mapdata, aes(x = PCOA1, y = PCOA2)) +
  #geom_point(alpha = 0.5, aes(color = as.factor(DepthLumping))) +
  geom_polygon(data = hab_bnti.pcoa.chulls, 
               aes(x = PCOA1, y = PCOA2, fill = Habitat__), alpha = 0.1) +
  geom_point(alpha = 0.5, aes(color = DepthLumping)) +
  geom_text(data = betaNTI.pcoa.mapdata %>%
              select(PCOA1HabAvg, PCOA2HabAvg, Habitat__) %>% distinct(), 
            aes(x = PCOA1HabAvg, y = PCOA2HabAvg, label = Habitat__),
            color = "grey20") +
  # geom_text(data = betaNTI.pcoa.mapdata %>% 
  #             filter(grepl("717", temporal_sample_id)), 
  #           aes(label = temporal_sample_id)) +
  theme_bw() + 
  scale_color_viridis(name = "DepthLumping", discrete = T, direction = -1) +
  scale_fill_discrete(name = "Habitat") +
  xlab(paste0("PCOA1 (", betaNTI.pct_ex[1], " %)")) + 
  ylab(paste0("PCOA2 (", betaNTI.pct_ex[2], " %)")) +
  ggtitle("Deterministic BetaNTI PCoA Plot")
bnti.pcoa.plot


# Run a pcoa on the new distance metric
# Mask all Heterogenous and "neither" selection
betaNTI.filt.norm.hom <- betaNTI.filt.norm
betaNTI.filt.norm.var <- betaNTI.filt.norm
betaNTI.filt.norm.hom[betaNTI.filt > -2] <- 0
betaNTI.filt.norm.var[betaNTI.filt < 2] <- 0
# Image heatmap
betaNTI.filt.norm.var %>% 
  as.data.frame() %>%
  rownames_to_column("f_id") %>%
  pivot_longer(-c(f_id), names_to = "samples", values_to = "value") %>%
  #filter(samples %in% c("714_S3_5-9", "715_E3_1-5")) %>%
  ggplot(aes(x=samples, y=f_id, fill=value)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

betaNTI.pcoa <- cmdscale(as.dist(betaNTI.filt.norm.var), k = 2, eig = TRUE, add = TRUE)
betaNTI.pct_ex <- round((betaNTI.pcoa$eig/sum(betaNTI.pcoa$eig)) * 100, 1)
betaNTI.pcoa.mapdata <- betaNTI.pcoa$points %>% data.frame() %>%
  rownames_to_column(var = "temporal_sample_id") %>%
  dplyr::rename(PCOA1 = X1, PCOA2 = X2) %>%
  left_join(input_ra$sample_metadata, by = "temporal_sample_id") %>%
  dplyr::group_by(Habitat__) %>%
  dplyr::mutate(PCOA1HabAvg = mean(PCOA1),
                PCOA2HabAvg = mean(PCOA2)) %>%
  ungroup() %>%
  dplyr::select(temporal_sample_id, starts_with("PCOA1"), starts_with("PCOA2"), PCOA1, PCOA2, everything())
hab_bnti.pcoa.chulls <- plyr::ddply(betaNTI.pcoa.mapdata, ~ Habitat__, 
                                    function(df) df[chull(df$PCOA1, df$PCOA2),])

bnti.pcoa.plot <- ggplot(betaNTI.pcoa.mapdata, aes(x = PCOA1, y = PCOA2)) +
  #geom_point(alpha = 0.5, aes(color = as.factor(DepthLumping))) +
  geom_polygon(data = hab_bnti.pcoa.chulls, 
               aes(x = PCOA1, y = PCOA2, fill = Habitat__), alpha = 0.1) +
  geom_point(alpha = 0.5, aes(color = DepthLumping)) +
  geom_text(data = betaNTI.pcoa.mapdata %>%
              select(PCOA1HabAvg, PCOA2HabAvg, Habitat__) %>% distinct(), 
            aes(x = PCOA1HabAvg, y = PCOA2HabAvg, label = Habitat__),
            color = "grey20") +
  geom_text(data = betaNTI.pcoa.mapdata,
            aes(label = temporal_sample_id)) +
  theme_bw() + 
  scale_color_viridis(name = "DepthLumping", discrete = T, direction = -1) +
  scale_fill_discrete(name = "Habitat") +
  xlab(paste0("PCOA1 (", betaNTI.pct_ex[1], " %)")) + 
  ylab(paste0("PCOA2 (", betaNTI.pct_ex[2], " %)")) +
  ggtitle("Deterministic BetaNTI PCoA Plot - Homogenous Selection Only")
bnti.pcoa.plot






RCBC.pcoa <- cmdscale(as.dist(RCBC.filt.norm), k = 2, eig = TRUE)
RCBC.pct_ex <- round((RCBC.pcoa$eig/sum(RCBC.pcoa$eig)) * 100, 1)
RCBC.pcoa.mapdata <- RCBC.pcoa$points %>% data.frame() %>%
  rownames_to_column(var = "temporal_sample_id") %>%
  dplyr::rename(PCOA1 = X1, PCOA2 = X2) %>%
  left_join(input_ra$sample_metadata, by = "temporal_sample_id") %>%
  dplyr::group_by(Habitat__) %>%
  dplyr::mutate(PCOA1HabAvg = mean(PCOA1),
                PCOA2HabAvg = mean(PCOA2)) %>%
  ungroup() %>%
  dplyr::select(temporal_sample_id, starts_with("PCOA1"), starts_with("PCOA2"), PCOA1, PCOA2, everything())
hab_RCBC.pcoa.chulls <- plyr::ddply(RCBC.pcoa.mapdata, ~ Habitat__, 
                                    function(df) df[chull(df$PCOA1, df$PCOA2), ])

RCBC.pcoa.plot <- ggplot(RCBC.pcoa.mapdata, aes(x = PCOA1, y = PCOA2)) +
  geom_point(alpha = 0.5, aes(color = Habitat__)) + 
  geom_polygon(data = hab_RCBC.pcoa.chulls, 
               aes(x = PCOA1, y = PCOA2), 
               alpha = 0, color = "grey80") +
  geom_text(data = RCBC.pcoa.mapdata %>%
              select(PCOA1HabAvg, PCOA2HabAvg, Habitat__) %>% distinct(), 
            aes(x = PCOA1HabAvg, y = PCOA2HabAvg, label = Habitat__),
            alpha = 0.5, color = "grey40") +
#  scale_color_viridis(name = "Sample-WT distance", direction = -1) +
  scale_fill_discrete(name = "Habitat") +
  xlab(paste0("PCOA1 (", RCBC.pct_ex[1], " %)")) + 
  ylab(paste0("PCOA2 (", RCBC.pct_ex[2], " %)"))
RCBC.pcoa.plot

# Mask all Heterogenous and "neither" selection
betaNTI.filt.norm.hom <- betaNTI.filt.norm
betaNTI.filt.norm.var <- betaNTI.filt.norm
betaNTI.filt.norm.hom[betaNTI.filt > -2] <- 0
betaNTI.filt.norm.var[betaNTI.filt < 2] <- 0
# Image heatmap
betaNTI.filt.norm.var %>% 
  as.data.frame() %>%
  rownames_to_column("f_id") %>%
  pivot_longer(-c(f_id), names_to = "samples", values_to = "value") %>%
  ggplot(aes(x=samples, y=f_id, fill=value)) + 
  geom_raster() +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#### ====================================================================== ####
# Useful functions
#### ====================================================================== ####
# data(pyrifos)
# ditch <- gl(12, 1, length=132)
# week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
# dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
# res_rda_pyrifos_inv_BC<-dbrda(pyrifos~interaction(week,dose) + Condition(week), distance = "bray")
# ctrl_pyrifos <- how(plots = Plots(strata = ditch,type = "free"), within = Within(type = "series"), nperm = 99)
# permutest(res_rda_pyrifos_inv_BC,  permutations=ctrl_pyrifos, first=TRUE)

#### ====================================================================== ####

# Set up metadata
#### ====================================================================== ####
env_data <- input_ra$sample_metadata %>%
  select(temporal_sample_id, Habitat__, Year__, DepthAvg__,
         T_air.deg_C, pH_peat, d13C_CH4__, d13C_CO2__, alphaC) %>%
  mutate(Year__ = as.factor(Year__)) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

lapply(seq_along(env_data[,4:ncol(env_data)]), function(x) {
  env_data.df <- env_data %>%
    select(1, x+3)
  env_hist.plot <- env_data.df %>%
    ggplot(aes(x = .data[[names(env_data.df)[2]]])) + 
    geom_histogram()
  
  # Create directory for plots if it doesn't exist
  if(!dir.exists("~/Downloads/histograms")) {
    dir.create("~/Downloads/histograms")}
  
  # Save plots
  plot_filename <- paste0("~/Downloads/histograms/", 
                               names(env_data.df[,2]), "_hist.png")
  # fix special characters
  plot_filename <- gsub("%", "Perc", plot_filename)
  
  ggsave(env_hist.plot, file = plot_filename)
})

# Calculate BC distance matrix for MAGs, for comparison
# Bray-curtis dissimilarities
motu_transformed <- t(sqrt(input$otu_table))
dm_bc <- vegdist(motu_transformed, method = "bray")
dm_bc <- as.matrix(dm_bc)

#### ====================================================================== ####
# Run cross-habitat dbRDA for Deterministic and Stochastic matrices
# Question: To what extent do Habitat and depth explain stochastic and deterministic
# assembly processes?
#### ====================================================================== ####
# Filter data
env_data_filt <- input_ra$sample_metadata %>%
  filter(Min_PF_dist_from_sample >= 0) %>%
  mutate(Year__ = as.numeric(as.character(Year__))) %>%
  #mutate(Habitat__ = as.numeric(as.factor(Habitat__))) %>%
  filter_at(.vars = c("DepthLumping", "DepthAvg__"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id, Habitat__, 
         matches(c("DepthLumping", "DepthAvg__"))) %>%
  mutate(Habitat__ = factor(Habitat__, levels = c("Palsa", "Bog", "Fen")))

ggpairs(env_data_filt, columns = 2:ncol(env_data_filt))


dim(env_data_filt)

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))
RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))

# Run dbrda
rcbc.rda <- dbrda(as.dist(RCBC_rank.filt) ~ Habitat__*DepthAvg__, 
                     data = env_data_filt[,2:ncol(env_data_filt)])
rcbc.rda.sum <- summary(rcbc.rda)
rcbc.rda.sum
anova.cca(rcbc.rda) # is model significant?
anova.cca(rcbc.rda, by = "terms") # what terms are signifcant


# MRM
MRM(as.dist(betaNTI.filt.norm) ~ dist(DepthAvg__) + as.dist(betaNTI.filt.dum),
    data =env_data_filt)

bnti.rda <- dbrda(as.dist(betaNTI.filt.norm) ~ Habitat__*DepthAvg__, 
                  data = env_data_filt, na.action = na.omit)
bnti.rda <- dbrda(as.dist(betaNTI.filt.norm.hom) ~ Habitat__*DepthAvg__, 
                  data = env_data_filt, na.action = na.omit)
plot(bnti.rda)

bnti.rda.sum <- summary(bnti.rda)
bnti.rda.sum
anova.cca(bnti.rda)
anova.cca(bnti.rda, by = "terms")

# compare to permanova
adonis(as.dist(betaNTI.filt.norm) ~ Habitat__*DepthAvg__, 
       data = env_data_filt, na.action = na.omit)


# Prepare ordination for plotting:

test <- summary(rcbc.rda)
test$biplot


ggord(rcbc.rda) +
  scale_x_continuous(limits = c(-1, 2)) +
  scale_y_continuous(limits = c(-2, 0.5))
anova(rcbc.rda)
anova(rcbc.rda, by = "terms", permu = 1000)

ggord(bnti.rda)
plot(bnti.rda)
anova(bnti.rda)
anova(bnti.rda, by = "terms", permu = 1000)

#### ====================================================================== ####
# 
# input$sample_metadata %>%
#   filter(Habitat__ == "Palsa") %>%
#   mutate(WTD = -1*WTD) %>%
#   mutate(ALD_detectable_num = ifelse(ALD_detectable == "No", Min_ALD, NA)) %>%
#   ggplot(aes(y = DepthAvg__, x = Year__)) +
#   geom_point(aes(y = ALD_detectable_num), color = "red", size = 4, shape = 21) +
#   geom_point(aes(y = DepthAvg__), size = 3) +
#   geom_point(aes(y = WTD), color = "steelblue", size = 3) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "darkolivegreen") +
#   geom_point(aes(y = Min_ALD), color = "burlywood4", size = 3) +
#   facet_grid(Habitat__~ Core__) +
#   scale_y_reverse()
# 
# input$sample_metadata %>%
#   mutate(ALD_detectable_num = ifelse(ALD_detectable == "No", Min_ALD, NA)) %>%
#   ggplot(aes(y = DepthAvg__, x = Year__)) + 
#   geom_point(aes(y = Min_ALD), color = "brown", size = 3) +
#   geom_point(aes(y = ALD_detectable_num), color = "red", size = 4, shape = 21) +
#   geom_point(aes(y = DepthAvg__), size = 3) +
#   facet_grid(Habitat__~ Core__) + 
#   scale_y_reverse()

# Question: For each of the three habitats, what env factors are most
# important in explaining assembly distance matrices; are they the same as
# factors that are most important in bc distance matrix?
#### ====================================================================== ####
# Filter data for Palsa
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Palsa") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("DepthAvg__", "Year__", 
                      "ALD", 
                      "T_air.deg_C", "T_soil.deg_C",
                      "Min_PF_dist_from_sample"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "Year__", 
                   "^ALD$", 
                   "T_air.deg_C", "T_soil.deg_C",
                   "Min_PF_dist_from_sample", "DepthLumping"))) %>%
  filter(Min_PF_dist_from_sample >= 0) # filter samples below ALD because they're unevenly sampled in Palsa

# Main missing information comes from Soil T and a couple from Air Temp

palsa.pairs.plot <- ggpairs(env_data_filt, columns = 2:ncol(env_data_filt))
palsa.pairs.plot

dim(env_data_filt)

palsa.samp.size <- nrow(env_data_filt)

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))
RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                    env_data_filt$temporal_sample_id]))
bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt$temporal_sample_id,
                                  env_data_filt$temporal_sample_id]))

# Run dbRDAs:
env_data_filt.rda <- input$sample_metadata %>% 
  filter(Habitat__ == "Palsa") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("DepthAvg__", "Year__", 
                      "ALD", 
                      "T_air.deg_C", "T_soil.deg_C",
                      "Min_PF_dist_from_sample"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "Year__", 
                   "^ALD$", 
                   "T_air.deg_C", "T_soil.deg_C",
                   "Min_PF_dist_from_sample", "DepthLumping"))) %>%
  filter(Min_PF_dist_from_sample >= 0) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

RCBC_rank.filt <- as.matrix(as.dist(RCBC_rank[env_data_filt.rda$temporal_sample_id,
                                    env_data_filt.rda$temporal_sample_id]))
dbrda.palsa.rcbc <- dbrda(as.dist(RCBC_rank.filt) ~ Year__ + Min_PF_dist_from_sample + T_air.deg_C  + T_soil.deg_C + Condition(DepthLumping),
                             data = env_data_filt.rda)
dbrda.palsa.rcbc0 <- dbrda(as.dist(RCBC_rank.filt) ~ 1,
                             data = env_data_filt.rda)

dbrda.palsa.rcbc.ordi <- ordistep(dbrda.palsa.rcbc0, 
                                  scope = formula(dbrda.palsa.rcbc), 
                                  direction = "both", permutations = 1000)
dbrda.palsa.rcbc.ordi$anova

env_data_filt.rdarc <- input$sample_metadata %>% 
  filter(Habitat__ == "Palsa") %>%
  filter(Min_PF_dist_from_sample >= 0) %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("T_soil.deg_C", "T_air.deg_C", "Year__"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("T_soil.deg_C", "T_air.deg_C", "Year__", "DepthLumping"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt.rdarc$temporal_sample_id,
                                          env_data_filt.rdarc$temporal_sample_id]))

dbrda.palsa.rcbc <- dbrda(as.dist(RCBC_rank.filt) ~ T_soil.deg_C + T_air.deg_C + Condition(DepthLumping),
                          data = env_data_filt.rda)
dbrda.palsa.rcbc.aov <- anova.cca(dbrda.palsa.rcbc)
anova.cca(dbrda.palsa.rcbc, by = "terms")

# Bnti
dbrda.palsa.bnti <- capscale(as.dist(betaNTI.filt) ~ Year__ + Min_PF_dist_from_sample + T_air.deg_C  + T_soil.deg_C + Condition(DepthLumping),
                             data = env_data_filt.rda, add = TRUE)

dbrda.palsa.bnti0 <- capscale(as.dist(betaNTI.filt) ~ 1,
                              data = env_data_filt.rda, add = TRUE)

dbrda.palsa.bnti.ordi <- ordistep(dbrda.palsa.bnti0, scope = formula(dbrda.palsa.bnti), 
                                  direction = "both", permutations = 1000)
dbrda.palsa.bnti.ordi$anova

env_data_filt.rdabnti <- input$sample_metadata %>% 
  filter(Habitat__ == "Palsa") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("Min_PF_dist_from_sample"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("Min_PF_dist_from_sample", "DepthLumping"))) %>%
  filter(Min_PF_dist_from_sample >= 0) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt.rdabnti$temporal_sample_id,
                                          env_data_filt.rdabnti$temporal_sample_id]))

dbrda.palsa.bnti <- capscale(as.dist(betaNTI.filt) ~ Min_PF_dist_from_sample + Condition(DepthLumping),
                             data = env_data_filt.rdabnti, add = TRUE)
dbrda.palsa.bnti.aov <- anova.cca(dbrda.palsa.bnti)

# BC
dbrda.palsa.bc <- capscale(as.dist(bc.filt) ~ Year__ + Min_PF_dist_from_sample + T_air.deg_C  + T_soil.deg_C + Condition(DepthLumping),
                             data = env_data_filt.rda, add = TRUE)

dbrda.palsa.bc0 <- capscale(as.dist(bc.filt) ~ 1,
                              data = env_data_filt.rda, add = TRUE)

dbrda.palsa.bc.ordi <- ordistep(dbrda.palsa.bc0, scope = formula(dbrda.palsa.bc), 
                                  direction = "both", permutations = 1000)
dbrda.palsa.bc.ordi$anova

dbrda.palsa.bc <- capscale(as.dist(bc.filt) ~ T_soil.deg_C + Year__ + Condition(DepthLumping),
                             data = env_data_filt.rda, add = TRUE)
dbrda.palsa.bc.aov <- anova.cca(dbrda.palsa.bc)

# Conclusions from dbRDA - distance from permafrost is a significant predictor of 
# both stochastic and deterministic turnover in community structure; by contrast
# community species turnover by b-c distance are best predicted by Soil temperature
# and year. 

# Run BioEnv to select variables
rcbc.palsa.bioenv <- bioenv(as.dist(RCBC.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

betanti.palsa.bioenv <- bioenv(as.dist(betaNTI.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

bc.palsa.bioenv <- bioenv(as.dist(bc.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

rcbc.palsa.bioenv
betanti.palsa.bioenv
bc.palsa.bioenv

# Run MRM to test strength of associations in RCBC
env_data_filt <- input$sample_metadata %>% 
  filter(Min_PF_dist_from_sample >= 0) %>%
  filter(Habitat__ == "Palsa") %>%
  filter_at(.vars = c("DepthAvg__", "T_air.deg_C"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "T_air.deg_C"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                    env_data_filt$temporal_sample_id]))

mrm.palsa.rcbc <- MRM(as.dist(RCBC.filt) ~ dist(DepthAvg__) + dist(T_air.deg_C), 
                data = env_data_filt, 
               nperm = 1000, 
               mrank = TRUE)

mrm.palsa.rcbc$samp.size <- nrow(env_data_filt)

# Run MRM to test strength of associations in BetaNTI
env_data_filt <- input$sample_metadata %>% 
  filter(Min_PF_dist_from_sample >= 0) %>%
  filter(Habitat__ == "Palsa") %>%
  filter_at(.vars = c("ALD"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("^ALD$")))  %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))

mrm.palsa.bnti <- MRM(as.dist(betaNTI.filt) ~ dist(ALD), 
                data = env_data_filt, 
                nperm = 1000, 
                mrank = TRUE) # not sig

mrm.palsa.bnti$samp.size <- nrow(env_data_filt)

# Run MRM to test strength of associations in bray-curtis
env_data_filt <- input$sample_metadata %>%
  filter(Min_PF_dist_from_sample >= 0) %>%
  filter(Habitat__ == "Palsa") %>%
  filter_at(.vars = c("DepthAvg__", "T_air.deg_C"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "T_air.deg_C"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))

mrm.palsa.bc <- MRM(as.dist(bc.filt) ~ dist(DepthAvg__) + dist(T_air.deg_C), 
                      data = env_data_filt, 
                      nperm = 1000, 
                      mrank = TRUE) # not sig

mrm.palsa.bc$samp.size <- nrow(env_data_filt)

mrm.palsa.rcbc
mrm.palsa.bnti
mrm.palsa.bc

# Conclusions: In the palsa, Depth and air temperature are sig. correlated to
# differences in stochastic forces, while none of the environmental variables
# are significantly correlated with phylogenetic dispersion matrix; Regular
# bc matrix is sig. driven by DepthAvg__

#### ====================================================================== ####

#### ====================================================================== ####
# Filter data for Bog
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Bog") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("DepthAvg__", "Year__", "T_air.deg_C", "T_soil.deg_C",
                      "WTD",
                      "alphaC", "pH_porewater"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "Year__", "T_air.deg_C", "T_soil.deg_C",
                   "^WTD$", "alphaC", "pH_porewater")))


ggpairs(env_data_filt, columns = 2:ncol(env_data_filt))

# Lots of strong co-correlations, we will make some executive decisions
# Air temp and year strongly correlated (0.847***), choose air temp b/c more representative of environment
# alphaC and depth correlated (0.793***); Choose alphaC b/c more representative of processes
# soil temp and year strongly correlated (0.630***); choose soil temp b/c more representative
# Soil temp and air temp mildly correlated (0.472*); leave both in
# Air temp and WTD strongly negatively correlated (-0.949***); holding off
# Upon further testing, BIO-ENV seems robust to inclusion of highly correlated variables
# Leaving in.

dim(env_data_filt)

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))
RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                    env_data_filt$temporal_sample_id]))
bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt$temporal_sample_id,
                                   env_data_filt$temporal_sample_id]))


# Run BioEnv to select variables
rcbc.bog.bioenv <- bioenv(as.dist(RCBC.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

betanti.bog.bioenv <- bioenv(as.dist(betaNTI.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

bc.bog.bioenv <- bioenv(as.dist(bc.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

rcbc.bog.bioenv
betanti.bog.bioenv
bc.bog.bioenv

# rcbc 
# Filter data for Bog
env_data_filt.rda <- input$sample_metadata %>% 
  filter(Habitat__ == "Bog") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("DepthAvg__", "Year__", "T_air.deg_C", "T_soil.deg_C",
                      "WTD",
                      "alphaC", "pH_porewater", "DepthLumping"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "Year__", "T_air.deg_C", "T_soil.deg_C",
                   "^WTD$", "alphaC", "pH_porewater", "DepthLumping"))) %>% 
  mutate(across(where(is.numeric), ~ scale(.)[,1]))


betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt.rda$temporal_sample_id,
                                          env_data_filt.rda$temporal_sample_id]))
RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt.rda$temporal_sample_id,
                                    env_data_filt.rda$temporal_sample_id]))
bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt.rda$temporal_sample_id,
                                   env_data_filt.rda$temporal_sample_id]))
  
dbrda.bog.rcbc <- capscale(as.dist(RCBC.filt) ~ Year__ +  T_air.deg_C + T_soil.deg_C + WTD + alphaC + pH_porewater + Condition(DepthLumping),
                           data = env_data_filt.rda, add = TRUE)

dbrda.bog.rcbc0 <- capscale(as.dist(RCBC.filt) ~ 1,
                            data = env_data_filt.rda, add = TRUE)

dbrda.bog.rcbc.ordi <- ordistep(dbrda.bog.rcbc0, scope = formula(dbrda.bog.rcbc), 
                                direction = "both", permutations = 1000) # NOT SIGNIFICANT

dbrda.bog.bnti <- capscale(as.dist(betaNTI.filt) ~ Year__ +  T_air.deg_C + T_soil.deg_C + WTD + alphaC + pH_porewater + Condition(DepthLumping),
                           data = env_data_filt.rda, add = TRUE)

dbrda.bog.bnti0 <- capscale(as.dist(betaNTI.filt) ~ 1,
                            data = env_data_filt.rda, add = TRUE)

dbrda.bog.bnti.ordi <- ordistep(dbrda.bog.bnti0, scope = formula(dbrda.bog.bnti), 
                                direction = "both", permutations = 1000) # NOT SIGNIFICANT
dbrda.bog.bnti.ordi$anova

env_data_filt.rdabnti <- input$sample_metadata %>% 
  filter(Habitat__ == "Bog") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("alphaC"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("alphaC", "DepthLumping"))) %>% 
  mutate(across(where(is.numeric), ~ scale(.)[,1]))


betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt.rdabnti$temporal_sample_id,
                                          env_data_filt.rdabnti$temporal_sample_id]))
dbrda.bog.bnti <- capscale(as.dist(betaNTI.filt) ~ alphaC + Condition(DepthLumping),
                           data = env_data_filt.rdabnti, add = TRUE)

anova.cca(dbrda.bog.bnti)
anova.cca(dbrda.bog.bnti, by = "terms")
# Only alphaC is significant

dbrda.bog.bc <- capscale(as.dist(bc.filt) ~ Year__ +  T_air.deg_C + T_soil.deg_C + WTD + alphaC + pH_porewater + Condition(DepthLumping),
                           data = env_data_filt.rda, add = TRUE)

dbrda.bog.bc0 <- capscale(as.dist(bc.filt) ~ 1,
                            data = env_data_filt.rda, add = TRUE)

dbrda.bog.bc.ordi <- ordistep(dbrda.bog.bc0, scope = formula(dbrda.bog.bc), 
                                direction = "both", permutations = 1000) # NOT SIGNIFICANT
dbrda.bog.bc.ordi$anova

env_data_filt.rdabc <- input$sample_metadata %>% 
  filter(Habitat__ == "Bog") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>% 
  filter_at(.vars = c("alphaC" , "Year__", "T_air.deg_C"),
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("alphaC", "Year__", "T_air.deg_C", "DepthLumping"))) %>% 
  mutate(across(where(is.numeric), ~ scale(.)[,1]))


bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt.rdabc$temporal_sample_id,
                                          env_data_filt.rdabc$temporal_sample_id]))
dbrda.bog.bc <- capscale(as.dist(bc.filt) ~ alphaC + Year__ + T_air.deg_C + Condition(DepthLumping),
                           data = env_data_filt.rdabc, add = TRUE)

anova.cca(dbrda.bog.bc)
anova.cca(dbrda.bog.bc, by = "terms")


# Run MRM to test strength of associations in RCBC
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Bog") %>%
  filter_at(.vars = c("DepthAvg__", "WTD"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "WTD"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

dim(env_data_filt)

RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                    env_data_filt$temporal_sample_id]))

mrm.bog.rcbc <- MRM(as.dist(RCBC.filt) ~ dist(DepthAvg__) + dist(WTD), 
                data = env_data_filt, 
                nperm = 1000, 
                mrank = TRUE)

mrm.bog.rcbc$samp.size <- nrow(env_data_filt)

# Run MRM to test strength of associations in BetaNTI
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Bog") %>%
  filter_at(.vars = c("T_soil.deg_C", "WTD"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("T_soil.deg_C", "^WTD$"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

dim(env_data_filt)

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))

mrm.bog.bnti <- MRM(as.dist(betaNTI.filt) ~ dist(T_soil.deg_C) + dist(WTD), 
                data = env_data_filt, 
                nperm = 1000, 
                mrank = TRUE) # not sig

mrm.bog.bnti$samp.size <- nrow(env_data_filt)

# Run MRM to test strength of associations in BC
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Bog") %>%
  filter_at(.vars = c("DepthAvg__"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

dim(env_data_filt)

bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))

mrm.bog.bc <- MRM(as.dist(bc.filt) ~ dist(DepthAvg__), 
                    data = env_data_filt, 
                    nperm = 1000, 
                    mrank = TRUE) # not sig

mrm.bog.bc$samp.size <- nrow(env_data_filt)


mrm.bog.rcbc
mrm.bog.bnti
mrm.bog.bc
# Conclusions: In the bog, Depth is sig. correlated (.23) to rcbc while none of
# the variables tested are strongly correlated to differences associated with betaNTI,
# DepthAvg__ is strongly correlated to bray-curtis distances
#### ====================================================================== ####

#### ====================================================================== ####
# Filter data for fen
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Fen") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>%
  filter_at(.vars = c("DepthAvg__", 
                      "Year__", 
                      "T_air.deg_C", 
                      "T_soil.deg_C",
                      "WTD", 
                      "alphaC", 
                      "pH_porewater"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", 
                   "Year__", 
                   "T_air.deg_C", 
                   "T_soil.deg_C",
                   "^WTD$", 
                   "alphaC", 
                   "pH_porewater"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))



ggpairs(env_data_filt, columns = 2:ncol(env_data_filt))

# Lots of strong co-correlations, we will make some executive decisions
# Air temp and year strongly correlated (0.860***);
# soil temp and depth mildly correlated (0.461*);
# soil temp and year strongly correlated (0.590***)
# Soil temp and air temp mildly correlated (0.440*); leave both in
# WTD and pH strongly correlated (0.599***)
# Upon further testing, BIO-ENV seems robust to inclusion of highly correlated variables
# Leaving in.

dim(env_data_filt)

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))
RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                    env_data_filt$temporal_sample_id]))
bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt$temporal_sample_id,
                                   env_data_filt$temporal_sample_id]))

# Run BioEnv to select variables
rcbc.fen.bioenv <- bioenv(as.dist(RCBC.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

betanti.fen.bioenv <- bioenv(as.dist(betaNTI.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

bc.fen.bioenv <- bioenv(as.dist(bc.filt) ~ ., data = env_data_filt[,2:ncol(env_data_filt)])

rcbc.fen.bioenv
betanti.fen.bioenv
bc.fen.bioenv

# dbrda
env_data_filt.rda <- input$sample_metadata %>% 
  filter(Habitat__ == "Fen") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>%
  filter_at(.vars = c("DepthAvg__", 
                      "Year__", 
                      "T_air.deg_C", 
                      "T_soil.deg_C",
                      "WTD", 
                      "alphaC", 
                      "pH_porewater"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", 
                   "Year__", 
                   "T_air.deg_C", 
                   "T_soil.deg_C",
                   "^WTD$", 
                   "alphaC", 
                   "pH_porewater", "DepthLumping"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt.rda$temporal_sample_id,
                                    env_data_filt.rda$temporal_sample_id]))
betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt.rda$temporal_sample_id,
                                          env_data_filt.rda$temporal_sample_id]))
bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt.rda$temporal_sample_id,
                                   env_data_filt.rda$temporal_sample_id]))

dbrda.fen.rcbc <- capscale(as.dist(RCBC.filt) ~ Year__ +  T_air.deg_C + T_soil.deg_C + WTD + alphaC + pH_porewater + Condition(DepthLumping),
                           data = env_data_filt.rda, add = TRUE)

dbrda.fen.rcbc0 <- capscale(as.dist(RCBC.filt) ~ 1,
                            data = env_data_filt.rda, add = TRUE)

dbrda.fen.rcbc.ordi <- ordistep(dbrda.fen.rcbc0, scope = formula(dbrda.fen.rcbc), 
                                direction = "both", permutations = 1000) 
dbrda.fen.rcbc.ordi$anova

env_data_filt.rdarc <- input$sample_metadata %>% 
  filter(Habitat__ == "Fen") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>%
  filter_at(.vars = c("alphaC"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("alphaC", "DepthLumping")))  %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt.rdarc$temporal_sample_id,
                                    env_data_filt.rdarc$temporal_sample_id]))

dbrda.fen.rcbc <- capscale(as.dist(RCBC.filt) ~ alphaC + Condition(DepthLumping),
                           data = env_data_filt.rdarc, add = TRUE)
dbrda.fen.rcbc.aov <- anova.cca(dbrda.fen.rcbc)


dbrda.fen.bnti <- capscale(as.dist(betaNTI.filt) ~ Year__ +  T_air.deg_C + T_soil.deg_C + WTD + alphaC + pH_porewater + Condition(DepthLumping),
                           data = env_data_filt.rda, add = TRUE)

dbrda.fen.bnti0 <- capscale(as.dist(betaNTI.filt) ~ 1,
                            data = env_data_filt.rda, add = TRUE)

dbrda.fen.bnti.ordi <- ordistep(dbrda.fen.bnti0, scope = formula(dbrda.fen.bnti), 
                                direction = "both", permutations = 1000) # not significant


# BC
dbrda.fen.bc <- capscale(as.dist(bc.filt) ~ Year__ +  T_air.deg_C + T_soil.deg_C + WTD + alphaC + pH_porewater + Condition(DepthLumping),
                           data = env_data_filt.rda, add = TRUE)

dbrda.fen.bc0 <- capscale(as.dist(bc.filt) ~ 1,
                            data = env_data_filt.rda, add = TRUE)

dbrda.fen.bc.ordi <- ordistep(dbrda.fen.bc0, scope = formula(dbrda.fen.bc), 
                                direction = "both", permutations = 1000) # not significant

env_data_filt.rdabc <- input$sample_metadata %>% 
  filter(Habitat__ == "Fen") %>%
  select_if(~sum(is.na(.)) < 0.60*length(.)) %>%
  filter_at(.vars = c("alphaC"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("alphaC", "DepthLumping")))  %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt.rdabc$temporal_sample_id,
                                    env_data_filt.rdabc$temporal_sample_id]))

dbrda.fen.bc <- capscale(as.dist(RCBC.filt) ~ alphaC + Condition(DepthLumping),
                           data = env_data_filt.rdabc, add = TRUE)
dbrda.fen.bc.aov <- anova.cca(dbrda.fen.bc)


# Conclusions: None of the environmental factors tested are signficant predictors of rcbc or betanti in the 
# fen


# Run MRM to test strength of associations in RCBC
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Fen") %>%
  filter_at(.vars = c("alphaC"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("alphaC"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

dim(env_data_filt)

RCBC.filt <- as.matrix(as.dist(RCBC[env_data_filt$temporal_sample_id,
                                    env_data_filt$temporal_sample_id]))

mrm.fen.rcbc <- MRM(as.dist(RCBC.filt) ~ dist(alphaC), 
                    data = env_data_filt, 
                    nperm = 1000, 
                    mrank = TRUE)

mrm.fen.rcbc$samp.size <- nrow(env_data_filt)

# Run MRM to test strength of associations in BetaNTI
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Fen") %>%
  filter_at(.vars = c("alphaC"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("alphaC"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

dim(env_data_filt)

betaNTI.filt <- as.matrix(as.dist(betaNTI[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))

mrm.fen.bnti <- MRM(as.dist(betaNTI.filt) ~ dist(alphaC), 
                    data = env_data_filt, 
                    nperm = 1000, 
                    mrank = TRUE) # not sig

mrm.fen.bnti$samp.size <- nrow(env_data_filt)

# Run MRM to test strength of associations in BC
env_data_filt <- input$sample_metadata %>% 
  filter(Habitat__ == "Fen") %>%
  filter_at(.vars = c("DepthAvg__","alphaC"), 
            ~!is.na(.)) %>%
  select(temporal_sample_id,
         matches(c("DepthAvg__", "alphaC"))) %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))

dim(env_data_filt)

bc.filt <- as.matrix(as.dist(dm_bc[env_data_filt$temporal_sample_id,
                                          env_data_filt$temporal_sample_id]))

mrm.fen.bc <- MRM(as.dist(bc.filt) ~ dist(DepthAvg__) + dist(alphaC), 
                    data = env_data_filt, 
                    nperm = 1000, 
                    mrank = TRUE) # not sig

mrm.fen.bc$samp.size <- nrow(env_data_filt)

mrm.fen.rcbc
mrm.fen.bnti
mrm.fen.bc

# Conclusions: In the fen, alphaC is strongly correlated to both rcbc and betaNTI,
# while both depth and alpha C are strongly correlated to bc
rcbc.filt.lf <- pivot_longer(RCBC.filt %>% as.data.frame() %>% 
               rownames_to_column(var = "SampleID.site1"),
             cols = !matches("SampleID.site1"), names_to = "SampleID.site2",
             values_to = "rcbc.filt")
bnti.filt.lf <- pivot_longer(betaNTI.filt %>% as.data.frame() %>% 
                             rownames_to_column(var = "SampleID.site1"),
                           cols = !matches("SampleID.site1"), names_to = "SampleID.site2",
                           values_to = "betanti.filt")
bc.filt.lf <- pivot_longer(bc.filt %>% as.data.frame() %>% 
                             rownames_to_column(var = "SampleID.site1"),
                           cols = !matches("SampleID.site1"), names_to = "SampleID.site2",
                           values_to = "bc.filt")

alphaC.dm.lf <- env_data_filt %>% 
  column_to_rownames("temporal_sample_id") %>%
  select(alphaC) %>%
  dist() %>% as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID.site1") %>%
  pivot_longer(.,
               cols = !matches("SampleID.site1"), 
               names_to = "SampleID.site2",
               values_to = "alphaC.dist")

DepthAvg.dm.lf <- env_data_filt %>% 
  column_to_rownames("temporal_sample_id") %>%
  select(DepthAvg__) %>%
  dist() %>% as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID.site1") %>%
  pivot_longer(.,
               cols = !matches("SampleID.site1"), 
               names_to = "SampleID.site2",
               values_to = "DepthAvg__.dist")

MRM.plot.fen.df <- left_join(rcbc.filt.lf, bnti.filt.lf, by = c("SampleID.site1", "SampleID.site2")) %>%
  left_join(bc.filt.lf, by = c("SampleID.site1", "SampleID.site2")) %>%
  left_join(alphaC.dm.lf, by = c("SampleID.site1", "SampleID.site2")) %>%
  left_join(DepthAvg.dm.lf, by = c("SampleID.site1", "SampleID.site2")) %>%
  filter(SampleID.site1 != SampleID.site2) %>%
  pivot_longer(cols = contains("filt"), names_to = "DistanceMeasure", 
               values_to = "Turnover")

ggplot(data = MRM.plot.fen.df, 
       aes(x = Turnover, y = alphaC.dist)) +
  geom_point() +
  facet_wrap(~DistanceMeasure, scales = "free_x")

#### ====================================================================== ####

# Combine and compare MRM results
#### ====================================================================== ####
extract_bioenv_mrm_results <- function(bioenv_result = rcbc.palsa.bioenv,
                                       mrm_result = mrm.palsa.rcbc,
                                       habitat = "palsa",
                                       DM_type = "rcbc") {
  bioenv_var_num <- bioenv_result$models[bioenv_result$whichbest][[1]]$best
  all_bioenv_vars <- bioenv_result$names
  bioenv_var_names_tested <- all_bioenv_vars[bioenv_var_num]
  

  
  mrm_sig <- mrm_result$r.squared
  # rename rownames to remove "dist()"
  rownames(mrm_result$coef) <- gsub("(dist\\()(.{1,})(\\))", "\\2", rownames(mrm_result$coef))
  mrm_coef <- as.data.frame(mrm_result$coef)
  
  # Create final data frame
  mrm_bioenv.df <- data.frame(Vars_tested_with_BioEnv = all_bioenv_vars, 
             Habitat = habitat,
             DM_type = DM_type,
             SampleSize = mrm_result$samp.size,
             MRM_sig_pval = mrm_sig[["pval"]],
             MRM_sig_r2 = mrm_sig[["R2"]]) %>%
    mutate(vars_in_MRM = ifelse(Vars_tested_with_BioEnv %in% bioenv_var_names_tested, 
                                Vars_tested_with_BioEnv,
                                NA),
           sig_in_MRM = ifelse(!is.na(vars_in_MRM), mrm_coef[vars_in_MRM, 2], NA),
           cor_in_MRM = ifelse(!is.na(vars_in_MRM), mrm_coef[vars_in_MRM, 1], NA))
  
  return(mrm_bioenv.df)
}

palsa.rcbc <- extract_bioenv_mrm_results(bioenv_result = rcbc.palsa.bioenv,
                                         mrm_result = mrm.palsa.rcbc,
                                         habitat = "Palsa", DM_type = "rcbc")
palsa.bnti <- extract_bioenv_mrm_results(bioenv_result = betanti.palsa.bioenv,
                                         mrm_result = mrm.palsa.bnti,
                                         habitat = "Palsa", DM_type = "bnti")
palsa.bc <- extract_bioenv_mrm_results(bioenv_result = bc.palsa.bioenv,
                                         mrm_result = mrm.palsa.bc,
                                         habitat = "Palsa", DM_type = "bc")

bog.rcbc <- extract_bioenv_mrm_results(bioenv_result = rcbc.bog.bioenv,
                                         mrm_result = mrm.bog.rcbc,
                                         habitat = "Bog", DM_type = "rcbc")
bog.bnti <- extract_bioenv_mrm_results(bioenv_result = betanti.bog.bioenv,
                                         mrm_result = mrm.bog.bnti,
                                         habitat = "Bog", DM_type = "bnti")
bog.bc <- extract_bioenv_mrm_results(bioenv_result = bc.bog.bioenv,
                                       mrm_result = mrm.bog.bc,
                                       habitat = "Bog", DM_type = "bc")

fen.rcbc <- extract_bioenv_mrm_results(bioenv_result = rcbc.fen.bioenv,
                                         mrm_result = mrm.fen.rcbc,
                                         habitat = "Fen", DM_type = "rcbc")
fen.bnti <- extract_bioenv_mrm_results(bioenv_result = betanti.fen.bioenv,
                                         mrm_result = mrm.fen.bnti,
                                         habitat = "Fen", DM_type = "bnti")
fen.bc <- extract_bioenv_mrm_results(bioenv_result = bc.fen.bioenv,
                                       mrm_result = mrm.fen.bc,
                                       habitat = "Fen", DM_type = "bc")

MRM_ecodist_results.df  <- bind_rows(palsa.rcbc, palsa.bnti, palsa.bc, 
          bog.rcbc, bog.bnti, bog.bc,
          fen.rcbc, fen.bnti, fen.bc) %>%
  mutate(fill_yesno = ifelse(sig_in_MRM < 0.05, TRUE, FALSE),
         included_in_model = ifelse(!is.na(vars_in_MRM), TRUE, FALSE)) %>%
  mutate(Habitat = factor(Habitat, levels = c("Palsa", "Bog", "Fen")),
         DM_type = factor(DM_type, levels = c("bnti", "rcbc", "bc")),
         DM_type_pretty = ifelse(DM_type == "bnti", "Deterministic",
                                 ifelse(DM_type == "rcbc", "Stochastic",
                                        "Bray-Curtis"))) %>%
  mutate(y_names = paste0(DM_type_pretty, "\nspearman rho = ", 
                          round(as.numeric(MRM_sig_r2), digits = 3),
                   "\n(n = ", SampleSize, ")"),
         fill_names = ifelse(MRM_sig_pval <= 0.05 & !is.na(vars_in_MRM), 
                             paste0(round(cor_in_MRM, digits = 3),"\n", 
                                    sig_in_MRM),
                             ""),
         fill_names = ifelse(MRM_sig_pval > 0.05 & !is.na(vars_in_MRM), 
                             paste0("NS"), fill_names),
         y_names = ifelse(MRM_sig_pval > 0.05, 
                             paste0("MRM not significant"), y_names))


MRM_bioenv_plot <- ggplot(MRM_ecodist_results.df, 
       aes(x = fct_reorder(Vars_tested_with_BioEnv, included_in_model, 
                           .desc = TRUE, .fun = sum), 
           y = y_names)) +
  geom_tile(aes(fill = fill_yesno), color = "black", size = 0.5) +
  geom_text(aes(label = fill_names, color = fill_yesno), size = rel(2.5)) +
  scale_color_manual(values = c("grey30", "black")) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = c(0.1, 0.4)) +
  scale_fill_brewer(name = "",
                    limits = c(TRUE, FALSE),
                    breaks = c(1),
                    labels = c("Best predictors"),
                    na.value = "transparent",
                    palette = "Greens", direction = -1) +
  facet_wrap(~Habitat, ncol = 3, scales = "free", 
             strip.position = "bottom") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0,  size = rel(1)),
        axis.text.y = element_text(angle = 0, hjust = 0, size = rel(1)),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.03,0.06,0.03,0.03), "npc"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.switch.pad.wrap = unit(c(0.03,0.03,0.03,0.03), "npc"),
        #strip.placement = "outside",
        legend.position = "none")

# Create BioEnv Plot legend
dummy_legend_df <- data.frame(x = c(1,2), y = c(1,1), 
                              fill = "coefficient\np-value", 
                              sig = factor(c("significant", "non-significant")))

bioenv_legend <- ggplot(dummy_legend_df, aes(x = fct_reorder(sig, x), y = y)) + 
  geom_tile(aes(fill = sig), color = "black", size = 0.5) +
  scale_fill_brewer(name = "",
                    limits = c("significant", "non-significant"),
                    breaks = c(1),
                    na.value = "transparent",
                    palette = "Greens", direction = -1) +
  geom_text(aes(label = fill, color = sig), size = rel(2)) +
  scale_color_manual(values = c("grey30", "black")) +
  scale_x_discrete(position = "bottom", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("") +
  facet_wrap(~ fct_reorder(sig, x), scales = "free", drop = TRUE) +
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(c(0.18), "npc"),
        legend.position = "none")


ggsave(MRM_bioenv_plot, filename = paste0(figures.fp, "/MRM_results.pdf"),
       width = 14, height = 4)
#### ====================================================================== ####


# Get a sense for how correlated betaNTI and rcbc values are
#### ====================================================================== ####
# Scale both matrices
min_dist <- min(betaNTI, na.rm = T)
max_dist <- max(betaNTI, na.rm = T)
betaNTI.scl <- apply(betaNTI, MARGIN = c(1,2), 
      FUN = function(x) {((x-min_dist)/(max_dist-min_dist))})
min_dist <- min(as.matrix(RCBC), na.rm = T)
max_dist <- max(as.matrix(RCBC), na.rm = T)
RCBC.scl <- apply(as.matrix(RCBC), MARGIN = c(1,2), 
                  FUN = function(x) {((x-min_dist)/(max_dist-min_dist))})

mantel(as.dist(betaNTI.scl), as.dist(RCBC.scl)) # mantel stat r = 0.479; sig = 0.001

mantel(as.dist(betaNTI.scl), as.dist(RCBC.scl)) # mantel stat r = 0.479; sig = 0.001


# plot BetaNTI against RCBC

ggplot(betanull.lf %>%
         filter(DepthLumping.Site1 == DepthLumping.Site2), 
       aes(x = BetaNTI, y = RCBC.nona)) +
  geom_point(aes(color = Assembly_Process))# +
  geom_smooth(method = "lm")
  
  ggplot(betanull.lf, 
         aes(x = BetaNTI, y = RCBC.nona)) +
    geom_point(aes(color = Assembly_Process))

#### ====================================================================== ####


