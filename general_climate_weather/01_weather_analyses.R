#' ## Weather analysis

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
library(viridis)
library(cowplot)
library(lubridate)
library(broom)
library(broom.mixed)
library(MuMIn)
library(here)


# Load required data
source(here("setup.R"))

## Set up input and output directories
#+ directory creation, eval = FALSE
# setting up names of file paths to stay organized
outputs.fp <- here("general_climate_weather", "outputs")
figures.fp <- here("general_climate_weather", "figures")

if (!dir.exists(outputs.fp)) {dir.create(outputs.fp)}
if (!dir.exists(figures.fp)) {dir.create(figures.fp)}

#### ====================================================================== ####

#### Read in data 
#### ====================================================================== ####
# Read in ANS weather data
ANS_weather.raw <- readxl::read_excel(here("data", "WeatherData", "Abisko_1913-2016.xls"), sheet = 2) %>%
  rename(AirTempC = `Avg_AirTemperature (°C)`,
         PrecipMM = `Precipitation (mm)`) %>%
  mutate(year = year(Date),
         month = month(Date)) %>%
  mutate(temporalstudy = ifelse(year %in% c(2011:2017), TRUE, FALSE)) %>%
  mutate(PrecipMM = ifelse(PrecipMM == -9999, NA, PrecipMM))

# Read in Swedish Polar Secretariate data (Station 188800 is the manual observations since 1913; there's another station that's digital since 1985)
ANS_weather_temp.raw <- read_delim(here("data", "smhi-opendata_2_188800_20230613_170724.csv"), 
                               skip = 9, delim = ";") %>% # warns about rows 15 forward because they don't have "Tidsutsnitt:" column
  rename(Date = `Representativt dygn`,
         AirTempC = Lufttemperatur)

ANS_weather_temp.raw <- read_delim(here("data", "smhi-opendata_1_188800_20230613_173232.csv"), 
                                   skip = 9, delim = ";") %>% # warns about rows 15 forward because they don't have "Tidsutsnitt:" column
  rename(Date = `Datum`,
         Time = `Tid (UTC)`,
         AirTempC = Lufttemperatur) %>%
  #mutate(year = year(Date)) %>%
  group_by(Date) %>%
  summarize(AirTempC = mean(AirTempC, na.rm = T)) %>%
  filter(Date <= "2018-12-31")
# Rainfall
ANS_weather_rain.raw <- read_delim(here("data", "smhi-opendata_5_188800_20230613_165112.csv"), 
                                   skip = 9, delim = ";") %>% # warns about rows 15 forward because they don't have "Tidsutsnitt:" column
  rename(Date = `Representativt dygn`,
         PrecipMM = Nederbördsmängd)

# Check quality:
# ANS_weather_rain.raw %>%
#   ggplot(aes(y = PrecipMM, x = Date)) +
#   geom_point(aes(color = Kvalitet))
#   
# ANS_weather_rain.raw %>% 
#   ggplot(aes(x = PrecipMM)) + 
#   geom_histogram()
# 
# ANS_weather_temp.raw %>%
#   ggplot(aes(y = AirTempC, x = Date)) +
#   geom_point() +
#   geom_point(data = ANS_weather_temp.raw %>% 
#                mutate(year = year(Date)) %>%
#                filter(year <= 2018) %>% # 2020 and 2019 are incomplete
#                group_by(year) %>%
#                summarize(AirTempC = mean(AirTempC)) %>%
#                mutate(Date = as_date(paste0(year, "-01-01"))),
#                color = "red") +
#   # geom_smooth(data = ANS_weather_temp.raw %>% 
#   #              mutate(year = year(Date)) %>%
#   #              filter(year <= 2018) %>% # 2020 and 2019 are incomplete
#   #              group_by(year) %>%
#   #              summarize(AirTempC = mean(AirTempC)) %>%
#   #              mutate(Date = as_date(paste0(year, "-01-01"))),
#              # color = "blue", method = "loess") +
#   geom_line(data = ANS_weather_temp.raw %>% 
#                 mutate(year = year(Date)) %>%
#                 filter(year <= 2018) %>% # 2020 and 2019 are incomplete
#                 group_by(year) %>%
#                 summarize(AirTempC = mean(AirTempC)) %>%
#                 mutate(Date = as_date(paste0(year, "-01-01"))),
#               color = "white")

# ANS_weather_temp.raw %>% 
#   mutate(data_source = "ANSWebsite") %>%
#   bind_rows(ANS_weather.raw %>% mutate(data_source = "EMERGE")) %>%
#   ggplot(aes(y = AirTempC, x = Date)) + 
#   geom_point(aes(color = Kvalitet)) +
#   facet_wrap(~data_source)
# 
# ANS_weather_temp.raw %>%
#   rename(ANSWebsiteAvg = AirTempC) %>%
#   full_join(ANS_weather.raw %>% mutate(EMERGEAirTempC = AirTempC), by = "Date") %>%
#   mutate(DiffANS_EMERGEAirTemp = ANSWebsiteAvg - EMERGEAirTempC) %>%
#   ggplot(aes(y = DiffANS_EMERGEAirTemp, x = Date)) +
#   geom_point() +
#   geom_point(data = ANS_weather_temp.raw %>% mutate(Date = lubridate::as_datetime(Date)),
#              aes(y = AirTempC), color = "red", shape = 21) +
#   geom_point(data = ANS_weather.raw,
#              aes(y = AirTempC), color = "blue", shape = 21)
# 
# 
# ANS_weather_rain.raw %>% 
#   ggplot(aes(x = PrecipMM)) + 
#   geom_histogram()
# There's a discrepancy in the EMERGE-db data and the ANS data in temperature. Have asked Suzanne 
# if she can add the next couple of years; for now temporarily using the ANS data directly
ANS_weather.raw <- full_join(ANS_weather_temp.raw, ANS_weather_rain.raw, by = "Date") %>%
  select(Date, AirTempC, PrecipMM) %>%
  mutate(year = year(Date),
         month = month(Date)) %>%
  mutate(temporalstudy = ifelse(year %in% c(2011:2017), TRUE, FALSE))

# Setup plotting parameters
#### ====================================================================== ####
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
shape_habitat  <- c(15, 16, 17)
fill_habitat   <- colour_habitat
#### ====================================================================== ####

# Our 7 years in context of temperature
#### ====================================================================== ####
# Years layered plot
years_layered_df <- ANS_weather.raw %>%
  filter(year %in% c(1913:2019)) %>%
  mutate(JulianDate = yday(Date),
         Date = as.Date(Date),
         WkAirTempC = zoo::rollmean(AirTempC, 7, na.pad = T, align = "right")) %>% # weekly air temperature
  mutate(JulianDate = as.Date("2021-12-31") + JulianDate)

years_layered_plot <- ggplot(years_layered_df, 
                             aes(x = JulianDate, y = AirTempC)) +
  geom_rect(xmin = as.Date("2022-07-01"), xmax = as.Date("2022-07-31"), 
            ymax = Inf, ymin = -Inf, alpha = 0.2, fill = "grey80") + # July
  xlab("Day of Year") +
  ylab("Temperature (Celsius, smoothed)") +
  geom_smooth(aes(color = year, group = year), se = F, alpha = 0.4) +
  geom_smooth(data = years_layered_df %>% 
              filter(year %in% c(2011:2017)),
            aes(x = JulianDate, y = AirTempC,
                group = year),se = F,
                color = "firebrick", linetype = "dashed",
            linewidth = 0.5, alpha = 0.4) +
  scale_color_viridis(name = "Year", direction = -1) +
  scale_x_date(date_breaks = "1 month", date_labels="%B") +
  theme_bw()
years_layered_plot
ggsave(filename = here(figures.fp, "layered_years_plot.png"), 
       plot = years_layered_plot, height = 5, width = 12)

extremes_ANS <- ANS_weather.raw %>%
  group_by(year) %>%
  mutate(yearMax = max(AirTempC),
         yearMin = min(AirTempC)) %>%
  ungroup() %>%
  mutate(avgyearMax = mean(yearMax, na.rm = T),
         sdyearMax = sd(yearMax, na.rm = T),
         avgyearMin = mean(yearMin, na.rm = T),
         sdyearMin = sd(yearMin, na.rm = T)) %>%
  select(contains("year")) %>%
  distinct()


#MAT
MAT_plot <- ANS_weather.raw %>%
  group_by(year) %>%
  summarize(MAT = mean(AirTempC, na.rm = T)) %>%
  filter(year <= 2017) %>%
  ggplot(aes(x = year, y = MAT)) +
  # Based on Gin's suggestion - this is the Varner 2021 Paper range:
  geom_rect(xmin = 1970, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80") +
  # This is the climate normal
  # geom_rect(xmin = 1992, xmax = 2017,
  #           ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80") +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60") +
  annotate(geom = "text", label = "Range covered in \nVarner et al. 2022",
           x = 1971, y = -Inf, hjust = 0, vjust = -0.10) +
  annotate(geom = "text", label = "This study",
           x = 2007, y = 2, hjust = 1, vjust = 1) +
  annotate(geom = "curve", arrow = arrow(length = unit(0.4, "lines")),
           x = 2007, xend = 2012, y = 2, yend = 2.11, curvature = -0.2) +
  geom_point() +
  geom_smooth(method = "loess", se = F) + # Note - gam vs. loess smoothing changes shape of fit line - imply different things...take with salt
  geom_line() + 
  ylab("Mean Annual Temperature") +
  xlab("Year (1913-2017)") +
  theme_bw()
MAT_plot
ggsave(filename = here(figures.fp, "MAT_plot.png"), 
       plot = MAT_plot, height = 5, width = 12)

# Seasonal MAT
seasonal_MT_plot <- ANS_weather.raw %>%
  filter(year <=2017) %>%
  mutate(JulianDate = yday(Date),
         Date = as.Date(Date), 
         Season = ifelse(month %in% c(12,1,2), "DJF",
                         ifelse(month %in% c(3:5), "MAM",
                                ifelse(month %in% c(6:8), "JJA",
                                       "SON"))),
         Season = factor(Season, levels = c("DJF", "MAM", "JJA", "SON")),
         SeasonLabels = fct_recode(Season, `Winter (DJF)` = "DJF",
                                   `Spring (MAM)` = "MAM",
                                   `Summer (JJA)` = "JJA",
                                   `Fall (SON)` = "SON")) %>%
  group_by(year, Season, SeasonLabels) %>%
  summarize(SeasonMean = mean(AirTempC)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = SeasonMean, group = Season)) +
  facet_grid(~SeasonLabels, margin = TRUE) +
  # Based on Gin's suggestion - this is the Varner 2021 Paper range:
  geom_rect(xmin = 1970, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80") +
  # This is the climate normal
  # geom_rect(xmin = 1992, xmax = 2017,
  #           ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80") +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60") +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  geom_line(aes(group = Season)) + 
  scale_y_continuous(breaks = seq(from = -20, to = 20, by = 1)) +
  ylab("Seasonal Mean Temperature") +
  xlab("Year (1913-2017)") +
  theme_bw()
seasonal_MT_plot

ggsave(filename = here(figures.fp, "seasonal_MT_plot.png"), 
       plot = seasonal_MT_plot, height = 5, width = 12)

# Days above freezing
DAF_plot <- ANS_weather.raw %>%
  filter(year <=2017) %>%
  group_by(year) %>%
  mutate(AboveFreezing = ifelse(AirTempC > 0, 1, 0),
         DaysAboveFreezing = max(cumsum(AboveFreezing))) %>%
  select(year, DaysAboveFreezing) %>%
  ggplot(aes(x = year, y = DaysAboveFreezing)) +
  # Based on Gin's suggestion - this is the Varner 2021 Paper range:
  geom_rect(xmin = 1970, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80") +
  # This is the climate normal
  # geom_rect(xmin = 1992, xmax = 2017,
  #           ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80") +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60") +
  geom_point() + 
  geom_line() +
  geom_smooth() + 
  ylab("Days above freezing") +
  xlab("Year (1913-2017)") +
  theme_bw()
DAF_plot

ggsave(filename = here(figures.fp, "DAF_plot.png"), 
       plot = DAF_plot, height = 5, width = 12)



# Days from last spring freeze to first fall freeze
growingSeason_df <- ANS_weather.raw %>%
  group_by(year) %>%
  mutate(DOY = yday(Date),
         #WkAirTempC = zoo::rollmean(AirTempC, 7, na.pad = TRUE, align = "right"),
         BelowFreezing = ifelse(AirTempC < 0, 1, 0),
         SpringFreeze = ifelse(BelowFreezing & DOY < 182, DOY, NA),
         FallFreeze = ifelse(BelowFreezing & DOY > 182, DOY, NA),
         LastSpFreeze = max(SpringFreeze, na.rm = T),
         FirstFlFreeze = min(FallFreeze, na.rm = T)) %>%
  mutate(across(all_of(c("LastSpFreeze", "FirstFlFreeze")), 
                ~as.Date("2021-12-31") + .x)) %>%
  select(year,LastSpFreeze, FirstFlFreeze) %>%
  distinct() %>% ungroup() %>%
  mutate(GrowingSeasonLength = FirstFlFreeze - LastSpFreeze,
         GrowingSeasonLength = as.numeric(GrowingSeasonLength))

growingSeason_plot <- ggplot(growingSeason_df) +
  geom_rect(xmin = 1992, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.8, fill = "grey80") +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60") +
  geom_segment(aes(x = year, xend = year, y = LastSpFreeze, yend = FirstFlFreeze),
               color = "grey30") +
  geom_point(aes(x = year, y = LastSpFreeze), color = rgb(0.2,0.7,0.1,0.5)) +
  geom_point(aes(x = year, y = FirstFlFreeze), color = rgb(0.7,0.2,0.1,0.5)) +
  coord_flip(xlim = c(1912,2017), expand = F) +
  ylab("Growing Season Length\n(last spring freeze to first fall freeze)") +
  xlab("Year (1913-2016)") +
  theme_bw()
growingSeason_plot
ggsave(filename = here(figures.fp, "GrowingSeason_plot.png"), 
       plot = growingSeason_plot, height = 12, width = 5)

growingSeason_separate_plot <- ggplot(growingSeason_df) +
  geom_rect(xmin = 1992, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.8, fill = "grey80") +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60") +
#  geom_segment(aes(x = year, xend = year, y = LastSpFreeze, yend = FirstFlFreeze),
#               color = "grey30") +
  geom_line(aes(x = year, y = LastSpFreeze), color = rgb(0.2,0.7,0.1,0.5)) +
  geom_line(aes(x = year, y = FirstFlFreeze), color = rgb(0.7,0.2,0.1,0.5)) +
#  coord_flip() +
  ylab("Growing Season\nLast Spring Freeze (green); First Fall Freeze (red)") +
  xlab("Year (1913-2016)") +
  theme_bw()
growingSeason_separate_plot

growingSeason_length_plot <- ggplot(growingSeason_df, 
                                    aes(x = year, y = GrowingSeasonLength)) +
  geom_rect(xmin = 1992, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.8, fill = "grey80") +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60") +
  geom_point() +
  geom_line() + 
  geom_smooth(se = F) +
  ylab("Growing Season Length (days)") +
  xlab("Year (1913-2016)") +
  theme_bw()
growingSeason_length_plot
ggsave(filename = here(figures.fp, "GrowingSeasonLength_plot.png"), 
       plot = growingSeason_length_plot, height = 5, width = 12)

#### ====================================================================== ####

# What is a hot/cold year
#### ====================================================================== ####
hotcold_plot <- input$sample_metadata %>%
  select(Year__, Habitat__, T_air.deg_C, contains("mean_AirTemperature")) %>%
  pivot_longer(cols = matches("air"), names_to = "TemperatureMeasure", values_to = "Temperature") %>% 
  mutate(TemperatureMeasure = gsub("mean_AirTemperature", "", TemperatureMeasure),
         TemperatureMeasure = gsub("_", "", TemperatureMeasure),
         TemperatureMeasure = factor(TemperatureMeasure, levels = c("Tair.degC", "samplingdate", "7d", "14d",
                                                                    "21d", "28d", "growing", "allgrowing"))) %>%
  group_by(Habitat__, TemperatureMeasure) %>%
  mutate(MaxMin = ifelse(Temperature == max(Temperature, na.rm = T), "max", 
                         ifelse(Temperature == min(Temperature, na.rm = T), "min",
                                NA))) %>%
  ungroup() %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
  ggplot(aes(x = Year__, y = Temperature)) + 
  geom_point(aes(color = MaxMin), size = rel(2)) + 
  scale_color_manual(name = "Hottest and Coldest Year", values = c("blue", "red"), breaks = c("min", "max")) +
  facet_grid(TemperatureMeasure~Habitat__, switch = "y") +
  xlab("Year") + ylab("Temperature") +
  theme_bw() + 
  theme(axis.title = element_text(size = rel(2)), 
        strip.text = element_text(size = rel(1)))
hotcold_plot
ggsave(paste0(figures.fp, "hotcold_stordalen.png"), plot = hotcold_plot, device = "png", dpi = 300,
       height = 11, width = 10)

#### ====================================================================== ####

# Ground temperature related
#### ====================================================================== ####
soil_temp_depth.plot <- input$sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels),
         Year__ = factor(Year__)) %>%
  ggplot( aes(x = DepthAvg__, y = T_soil.deg_C, group = Year__, fill = Year__, color = Year__)) + 
  geom_point() +
  facet_wrap(~Habitat__) +
  geom_smooth(method = "lm") +
  ylab("Ground Temperature (°C)") + xlab("Average Sample Depth (cm)") +
  scale_color_manual(values = RColorBrewer::brewer.pal(7, "YlOrBr"), name = "Year") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(7, "YlOrBr"), name = "Year") +
  theme_bw()
ggsave(paste0(figures.fp, "soil_depth_plot.png"), plot = soil_temp_depth.plot, device = "png", dpi = 300,
       height = 7, width = 14)


# Sideways
soil_temp_depth_sideways.plot <- input$sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels),
         Year__ = factor(Year__)) %>%
  ggplot( aes(x = DepthAvg__, y = T_soil.deg_C, group = Year__, fill = Year__, color = Year__)) + 
  geom_point() +
  facet_wrap(Habitat__~., ncol = 1, strip.position = "right") +
  geom_smooth(method = "lm") +
  ylab("Ground Temperature (°C)") + xlab("Average Sample Depth (cm)") +
  scale_x_reverse() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_color_manual(values = RColorBrewer::brewer.pal(7, "YlOrBr"), name = "Year") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(7, "YlOrBr"), name = "Year") +
  theme_bw()
ggsave("~/Downloads/soil_depth_plot_sideways.png", plot = soil_temp_depth_sideways.plot, device = "png", dpi = 300,
       height = 10, width = 8)





Intercept_plot.df <- sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
  dplyr::select(temporal_sample_id, Habitat__, Year__, DepthAvg__, T_soil.deg_C) %>%
  filter(!is.na(T_soil.deg_C)) %>%
  dplyr::select(-temporal_sample_id) %>% # removes misleading extra "observations" of CH4 flux that actually are just the same across samples
  group_by(Habitat__) %>%
  nest() %>% #pluck(2,3)
  mutate(lmer.model = purrr::map(data, ~lmer(T_soil.deg_C~DepthAvg__ + (1|Year__), data = .))) %>% # random intercept (nearly equivalent to random slope + rand intercept in all habitats)
  mutate(lmer.tidied = purrr::map(lmer.model, tidy), # %>%  coef()
         lmer.model.intest = purrr::map(lmer.model, ~coef(.)[[1]] %>% 
                                          data.frame() %>% 
                                          rownames_to_column(var = "Year") %>% 
                                          rename(Intercept = X.Intercept.))) %>% 
  unnest(lmer.model.intest) %>%
  # Prepare to calculate relationship between fitted intercepts and year
  mutate(Year = as.numeric(Year)) %>%
  rename(LMERInt = Intercept) %>%
  group_by(Habitat__) %>%
  nest() %>%
  mutate(lm.model = purrr::map(data, ~lm(LMERInt ~ Year, data = .)),
         lm.tided = purrr::map(lm.model, tidy),
         lm.model.rsqr = map_dbl(lm.model, ~glance(.)$adj.r.squared),
         lm.model.pval = map_dbl(lm.model, ~glance(.)$p.value)) %>%
  unnest(data) %>%
  mutate(Sig = ifelse(lm.model.pval < 0.05, T, F))


Intercept_plot.df %>% 
  dplyr::select(Habitat__, data, lmer.model, lmer.tidied) %>%
  # undo the unnest so there's only one copy per habitat
  nest() %>% #pluck(2,1)
  mutate(origdata = purrr::map(data, ~pluck(.,1,1)),
         lmer.model = purrr::map(data, ~pluck(.,2,1)),
         lmer.tidied = purrr::map(data, ~pluck(.,3,1))) %>%
  dplyr::select(-data) %>%
  unnest(lmer.tidied) %>%
  mutate(lwr = estimate - std.error, upr = estimate + std.error)


Intercept_plot <- Intercept_plot.df %>%
  ggplot(aes(x = Year, y = LMERInt)) + 
  geom_point(aes(fill = Habitat__), stat = "identity", shape = 21, size = 2) +
  facet_wrap(~Habitat__, scales = "free_y") +
  geom_smooth(data = Intercept_plot.df %>% filter(Sig),
              method = "lm", color = "black", aes(group = Habitat__)) +
  geom_text(data = Intercept_plot.df %>%
              dplyr::select(Habitat__, lm.model.rsqr, lm.model.pval, Sig) %>%
              filter(Sig) %>% distinct() %>%
              mutate(label = paste0("\n   p = ", round(lm.model.pval,2), "\n",
                                    "   r2 = ", round(lm.model.rsqr, 2))),
            aes(label = label), color = "black",
            x = -Inf, y = Inf, vjust = 1, hjust = 0) +
  scale_fill_manual(breaks = habitat_levels, values = colour_habitat) +
  theme_bw() +
  ylab("Intercept") +
  ggtitle("Fitted Y-intercept of Ground Temperature Values and Depth over time")

Intercept_plot
ggsave("~/Downloads/PooledSoilInterceptTime.png", plot = Intercept_plot, device = "png", dpi = 300,
       height = 7, width = 14)


AvgTemp_Time.plot.df <- sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
  dplyr::select(Habitat__, Year__, DepthAvg__, T_soil.deg_C) %>%
  filter(!is.na(T_soil.deg_C)) %>% # removes misleading extra "observations" of CH4 flux that actually are just the same across samples
  group_by(Habitat__, Year__) %>%
  mutate(AverageSoilTempAcrossDepth = mean(T_soil.deg_C)) %>%
  dplyr::select(Habitat__, Year__, AverageSoilTempAcrossDepth) %>%
  distinct() %>%
  group_by(Habitat__) %>%
  nest() %>%
  mutate(lm.model = purrr::map(data, ~lm(AverageSoilTempAcrossDepth ~ Year__, data = .)),
         lm.tided = purrr::map(lm.model, tidy),
         lm.model.rsqr = map_dbl(lm.model, ~glance(.)$adj.r.squared),
         lm.model.pval = map_dbl(lm.model, ~glance(.)$p.value)) %>%
  unnest(data) %>%
  mutate(Sig = ifelse(lm.model.pval < 0.05, T, F))

AvgTemp_Time.plot <- AvgTemp_Time.plot.df  %>%
  ggplot(aes(x = Year__, y = AverageSoilTempAcrossDepth)) +
  geom_point(aes(fill = Habitat__), shape = 21, size = 2) +
  facet_wrap(~Habitat__, scales = "free_y") +
  geom_smooth(data = AvgTemp_Time.plot.df %>% filter(Sig),
              method = "lm", color = "black", aes(group = Habitat__)) +
  geom_text(data = AvgTemp_Time.plot.df %>%
              dplyr::select(Habitat__, lm.model.rsqr, lm.model.pval, Sig) %>%
              filter(Sig) %>% distinct() %>%
              mutate(label = paste0("\n   p = ", round(lm.model.pval,2), "\n",
                                    "   r2 = ", round(lm.model.rsqr, 2))),
            aes(label = label), color = "black",
            x = -Inf, y = Inf, vjust = 1, hjust = 0) +
  scale_fill_manual(breaks = habitat_levels, values = colour_habitat) +
  xlab("Year") + ylab("Average Soil Temperature (across 0-39 cm)") +
  theme_bw() +
  ggtitle("Ground Temperature (Averaged Across Depth) Over Time")

AvgTemp_Time.plot 

ggsave("~/Downloads/AverageSoilTempTime.png", plot = AvgTemp_Time.plot, device = "png", dpi = 300,
       height = 7, width = 14)


Temp_Time_depth.plot.df <- sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
  dplyr::select(Habitat__, Year__, DepthAvg__, DepthLumping, T_soil.deg_C) %>%
  filter(!is.na(T_soil.deg_C)) %>% # removes misleading extra "observations" of CH4 flux that actually are just the same across samples
  dplyr::select(Habitat__, Year__, DepthLumping, T_soil.deg_C) %>%
  group_by(Habitat__, DepthLumping) %>%
  nest() %>%
  mutate(lm.model = purrr::map(data, ~lm(T_soil.deg_C ~ Year__, data = .)),
         lm.tidied = purrr::map(lm.model, tidy),
         lm.model.rsqr = map_dbl(lm.model, ~glance(.)$adj.r.squared),
         lm.model.pval = map_dbl(lm.model, ~glance(.)$p.value),
         lm.model.yr_est = map_dbl(lm.tidied, ~.$estimate[2]),
         lm.model.yr_int = map_dbl(lm.tidied, ~.$estimate[1])) %>%
  mutate(Sig = ifelse(lm.model.pval < 0.05, T, F),
         plotlab = paste0("\n r2 = ", round(lm.model.rsqr, 2),
                          ", p = ",round(lm.model.pval, 2),
                          "\n slope = ", round(lm.model.yr_est, 2))) %>%
  unnest(data)

# Horizontal Plot;
Temp_Time_depth.plot <- Temp_Time_depth.plot.df %>%
  ggplot(aes(x = Year__, y = T_soil.deg_C)) +
  geom_point(aes(shape = DepthLumping), size = rel(5)) +
  geom_segment(data = Temp_Time_depth.plot.df %>% 
                 mutate(MinYear = min(Year__), MaxYear = max(Year__)) %>%
                 select(-Year__, -T_soil.deg_C) %>% distinct() %>%
                mutate(lm.model.yr_est = ifelse(Sig == T, lm.model.yr_est, NA),
                       lm.model.yr_int = ifelse(Sig == T, lm.model.yr_int, NA)),
              aes(x = MinYear, xend = MaxYear,
                  y = lm.model.yr_int + lm.model.yr_est*MinYear,
                  yend = lm.model.yr_int + lm.model.yr_est*MaxYear,
                  group = DepthLumping, linetype = DepthLumping,
                  color = Habitat__), linewidth = rel(2)) +
  geom_smooth(data = Temp_Time_depth.plot.df %>% filter(Sig == T),
              method = "lm", aes(group = DepthLumping, fill = Habitat__),
              alpha = 0.2, color = NA, linetype = "blank") +
  # remove model output text
  # geom_text(data = Temp_Time_depth.plot.df %>% filter(Sig == T) %>% 
  #             dplyr::select(DepthLumping, Habitat__, plotlab) %>% 
  #             distinct(),
  #           aes(label = plotlab), size = rel(2.5),
  #           x = -Inf, y = Inf, vjust = 0.8, hjust = 0) +
  facet_grid(~Habitat__) +
  scale_x_continuous(limits = c(2011, 2017), breaks = seq(from = 2012, to = 2017, by = 2)) +
  coord_fixed(ratio = 0.29) +
  scale_shape_manual(name = "Depth", breaks = c("0-9", "10-19", "20-29", "30-39"), 
                       values = c(1,2,3,4)) +
  scale_linetype_manual(name = "Depth", breaks = c("0-9", "10-19", "20-29", "30-39"),
                        values = c(1,2,3,4)) +
  scale_color_manual(name = "Habitat", breaks = habitat_levels, values = colour_habitat) +
  scale_fill_manual(name = "Habitat", breaks = habitat_levels, values = colour_habitat) +
  guides(linetype = guide_legend(keywidth = unit(5, "cm"), override.aes = list(fill = NA, color = "black", linewidth = rel(1)))) +
  theme_bw() +
  theme(axis.title = element_text(size = rel(3)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(2.5)),
        axis.text.y = element_text(size = rel(2.5)),
        strip.text = element_blank(),
        legend.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(2)),
        panel.spacing = unit(1.3, "lines"),
        strip.background = element_rect(fill = "white", color = NA),
        plot.margin = unit(c(1,1,1,1), "lines")) +
  xlab("Year") + ylab("Soil temperature (°C)")
Temp_Time_depth.plot

Temp_Time_depth.plot_legend <- get_legend(Temp_Time_depth.plot)

Temp_Time_depth.plot_no_legend <- Temp_Time_depth.plot + theme(legend.position = "none")

ggsave(paste0(figures.fp, "SoilTempDepthTime_horizontal.png"), plot = Temp_Time_depth.plot_no_legend, device = "png", dpi = 300,
       height = 7, width = 10)

# Vertical
Temp_Time_depth.plot_v <- Temp_Time_depth.plot.df %>%
  ggplot(aes(x = Year__, y = T_soil.deg_C)) +
  geom_point(aes(shape = DepthLumping), size = rel(5)) +
  geom_segment(data = Temp_Time_depth.plot.df %>% 
                 mutate(MinYear = min(Year__), MaxYear = max(Year__)) %>%
                 select(-Year__, -T_soil.deg_C) %>% distinct() %>%
                 mutate(lm.model.yr_est = ifelse(Sig == T, lm.model.yr_est, NA),
                        lm.model.yr_int = ifelse(Sig == T, lm.model.yr_int, NA)),
               aes(x = MinYear, xend = MaxYear,
                   y = lm.model.yr_int + lm.model.yr_est*MinYear,
                   yend = lm.model.yr_int + lm.model.yr_est*MaxYear,
                   group = DepthLumping, linetype = DepthLumping,
                   color = Habitat__), linewidth = rel(2)) +
  geom_smooth(data = Temp_Time_depth.plot.df %>% filter(Sig == T),
              method = "lm", aes(group = DepthLumping, fill = Habitat__),
              alpha = 0.2, color = NA, linetype = "blank") +
  # remove model output text
  # geom_text(data = Temp_Time_depth.plot.df %>% filter(Sig == T) %>% 
  #             dplyr::select(DepthLumping, Habitat__, plotlab) %>% 
  #             distinct(),
  #           aes(label = plotlab), size = rel(2.5),
  #           x = -Inf, y = Inf, vjust = 0.8, hjust = 0) +
  facet_grid(Habitat__~.) +
  scale_shape_manual(name = "Depth", breaks = c("0-9", "10-19", "20-29", "30-39"), 
                     values = c(1,2,3,4)) +
  scale_linetype_manual(name = "Depth", breaks = c("0-9", "10-19", "20-29", "30-39"),
                        values = c(1,2,3,4)) +
  scale_color_manual(name = "Habitat", breaks = habitat_levels, values = colour_habitat) +
  scale_fill_manual(name = "Habitat", breaks = habitat_levels, values = colour_habitat) +
  guides(linetype = guide_legend(keywidth = unit(5, "cm"), override.aes = list(fill = NA, color = "black", linewidth = rel(1)))) +
  theme_bw() +
  theme(axis.title = element_text(size = rel(2)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        strip.text = element_blank(),
        legend.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(2)),
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = "white", color = NA)) +
  xlab("Year") + ylab("Soil temperature (°C)")
Temp_Time_depth.plot_v
ggsave(paste0(figures.fp, "SoilTempDepthTime_vertical.png"), plot = Temp_Time_depth.plot_v, device = "png", dpi = 300,
       height = 15, width = 7)




#### ====================================================================== ####


# Ground temperature model comparison
#### ====================================================================== ####

# lay out competing models
m1 <- "T_soil.deg_C ~ DepthAvg__"
m2 <- "T_soil.deg_C ~ DepthAvg__ + (1|Year__)"
m3 <- "T_soil.deg_C ~ DepthAvg__ + (DepthAvg__|Year__)"
# Swapping out Year for air temperature metrics
m4 <- "T_soil.deg_C ~ DepthAvg__ + T_air.deg_C"
m5 <- "T_soil.deg_C ~ DepthAvg__ + samplingdate_mean_AirTemperature"
m6 <- "T_soil.deg_C ~ DepthAvg__ + mean_AirTemperature_7d"
m7 <- "T_soil.deg_C ~ DepthAvg__ + mean_AirTemperature_14d"
m8 <- "T_soil.deg_C ~ DepthAvg__ + mean_AirTemperature_21d"
m9 <- "T_soil.deg_C ~ DepthAvg__ + mean_AirTemperature_28d"
m10 <- "T_soil.deg_C ~ DepthAvg__ + mean_AirTemperature_growing"
m11 <- "T_soil.deg_C ~ DepthAvg__ + mean_AirTemperature_all_growing"

modelnames = c("DepthOnly", "DepthYrRdInt", "DepthYrRdIntSl", 
              "DepthTairSampTime", "DepthTairSampDate", "DepthTair7d",
              "DepthTair14d", "DepthTair21d", "DepthTair28d", 
              "DepthTairGrowing", "DepthTairAllGrowing")
modelnamedlist <- c(m1, m2, m3, m4, m5, m6,
                    m7, m8, m9, m10, m11)
names(modelnamedlist) <- modelnames

# Create tibble with models
model_frame <- tibble(modelname = modelnames,
                          model = as.list(modelnamedlist)) %>%
  mutate(model = purrr::map(model, as.formula),
         modelType = ifelse(modelname %in% c("DepthYrRdInt", "DepthYrRdIntSl"), "Mixed", "LM"))

  
# Create a model frame for each habitat's analysis
model_frame <- bind_rows(model_frame %>% mutate(Habitat__ = "Palsa"),
                         model_frame %>% mutate(Habitat__ = "Bog"),
                         model_frame %>% mutate(Habitat__ = "Fen")) %>%
  mutate(Habitat__ = factor(Habitat__, levels =  habitat_levels))

# Prepare metadata
Tsoil_data <- sample_metadata %>%
  mutate(Habitat__ = factor(Habitat__, levels = habitat_levels)) %>%
  dplyr::select(temporal_sample_id, Habitat__, Year__, DepthAvg__, T_soil.deg_C,
                T_air.deg_C, samplingdate_mean_AirTemperature, mean_AirTemperature_7d,
                mean_AirTemperature_14d, mean_AirTemperature_21d, 
                mean_AirTemperature_28d, mean_AirTemperature_growing,
                mean_AirTemperature_all_growing) %>%
  # Remove any entries where the data is na
  filter(!is.na(T_soil.deg_C)) %>%
#  filter(if_all(T_air.deg_C:mean_AirTemperature_all_growing, ~!is.na(.))) %>% 
  dplyr::select(-temporal_sample_id) %>% 
  group_by(Habitat__) %>%
  nest()

# Combine sample data and model frame
TsoilAirModComp.df <- left_join(Tsoil_data, model_frame, by = "Habitat__") %>%
  # Since T_air.deg_C is the only other column with NAs, filter out NAs only in 
  # those specific models
  mutate(data = purrr::map(data, ~filter(.x, !is.na(T_air.deg_C)))) %>%
  # mutate(data = ifelse(modelname == "DepthTairSampTime", 
  #                      purrr::map(data, ~filter(.x, !is.na(T_air.deg_C))),
  #                      data)) %>% 
  # Run an LMER model for those with fixed effects, or LM model for others
  # note - DepthYrRdIntSl model is singular fit in some cases
  mutate(modelFit = purrr::map2(data, model, ~{if(grepl("Year__", format(.y))) { # use format to convert formula to character; not as.character
                                                      lme4::lmer(.y, data = .x)} else {
                                                        lm(.y, data = .x)
                                                      }
                                               })) %>%
  mutate(mod.tidied = purrr::map(modelFit, broom::tidy),
         mod.AIC = map_dbl(modelFit, ~broom::glance(.)$AIC),
         mod.BIC = map_dbl(modelFit, ~broom::glance(.)$BIC))

TsoilAirModComp.df %>% arrange(Habitat__, mod.AIC) %>% View()

TsoilAirModComp.df %>% filter(Habitat__ == "Palsa") %>% pluck("modelFit") %>% model.sel(.)
TsoilAirModComp.df %>% filter(Habitat__ == "Bog") %>% pluck("modelFit") %>% model.sel(.)
TsoilAirModComp.df %>% filter(Habitat__ == "Fen") %>% pluck("modelFit") %>% model.sel(.)

# Palsa models
palsa_mod_comp <- TsoilAirModComp.df %>% filter(Habitat__ == "Palsa") %>% 
  pluck("modelFit")
names(palsa_mod_comp) <-TsoilAirModComp.df %>% filter(Habitat__ == "Palsa") %>% 
  pluck("modelname") 
palsa_mod_comp <- model.sel(palsa_mod_comp)

palsa_mod_comp %>% subset(., delta<5) 
palsa_mod_comp %>%
  # uses Royall's 1/8th rule for strength of evidence
  subset(., 1/8 < weight/max(.$weight)) 
palsa_mod_comp %>%
  subset(., cumsum(.$weight) <= 0.95)

palsa.sel.table<-as.data.frame(palsa_mod_comp)[,c(11:16)]


# Bog Models:
bog_mod_comp <- TsoilAirModComp.df %>% filter(Habitat__ == "Bog") %>% 
  pluck("modelFit")
names(bog_mod_comp) <-TsoilAirModComp.df %>% filter(Habitat__ == "Bog") %>% 
  pluck("modelname") 
bog_mod_comp <- model.sel(bog_mod_comp)

bog_mod_comp %>% subset(., delta<5) 
bog_mod_comp %>%
  # uses Royall's 1/8th rule for strength of evidence
  subset(., 1/8 < weight/max(.$weight)) 
bog_mod_comp %>%
  subset(., cumsum(.$weight) <= 0.95)

bog.sel.table<-as.data.frame(bog_mod_comp)[,c(11:16)]


# Fen Models:
fen_mod_comp <- TsoilAirModComp.df %>% filter(Habitat__ == "Fen") %>% 
  pluck("modelFit")
names(fen_mod_comp) <-TsoilAirModComp.df %>% filter(Habitat__ == "Fen") %>% 
  pluck("modelname") 
fen_mod_comp <- model.sel(fen_mod_comp)

fen_mod_comp %>% subset(., delta<5) 
fen_mod_comp %>%
  # uses Royall's 1/8th rule for strength of evidence
  subset(., 1/8 < weight/max(.$weight)) 
fen_mod_comp %>%
  subset(., cumsum(.$weight) <= 0.95)

fen.sel.table<-as.data.frame(fen_mod_comp)[,c(11:16)]



palsa.sel.table
bog.sel.table
fen.sel.table

#         model.rsqr = map_dbl(modelFit, ~broom::glance(.)$adj.r.squared),
#         model.pval = map_dbl(modelFit, ~broom::glance(.)$p.value),
#         model.yr_est = map_dbl(mod.tidied, ~.$estimate[2])) %>%
  mutate(Sig = ifelse(lm.model.pval < 0.05, T, F),
         plotlab = paste0("\n r2 = ", round(lm.model.rsqr, 2),
                          ", p = ",round(lm.model.pval, 2),
                          "\n slope = ", round(lm.model.yr_est, 2)))



#### ====================================================================== ####

# Our 7 years in context of precipitation
#### ====================================================================== ####
#MAP
MAP_df <- ANS_weather.raw %>%
  mutate(DOY = yday(Date),
         Season = ifelse(month %in% c(12,1,2), "DJF",
                         ifelse(month %in% c(3:5), "MAM",
                                ifelse(month %in% c(6:8), "JJA",
                                       "SON"))),
         Season = factor(Season, levels = c("DJF", "MAM", "JJA", "SON")),
         SeasonLabels = fct_recode(Season, `Winter (DJF)` = "DJF",
                                   `Spring (MAM)` = "MAM",
                                   `Summer (JJA)` = "JJA",
                                   `Fall (SON)` = "SON")) %>%
  mutate(JulianDate = as.Date("2021-12-31") + DOY) %>%
  group_by(year) %>%
  mutate(MAP = mean(PrecipMM, na.rm = T),
         NACount = ifelse(is.na(PrecipMM), 1, 0),
         NASum = cumsum(NACount),
         halfNAs = ifelse((max(NASum)/365) > 0.5, T, F),
         CumPrecip = cumsum(ifelse(is.na(PrecipMM), 0, 1))) 
MAP_plot <- ggplot(MAP_df, aes(x = year, y = MAP)) +
  geom_rect(xmin = 1992, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80") +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60") +
  geom_point(data = MAP_df %>%
               select(year, MAP, halfNAs) %>%
               distinct(), 
             aes(x = year, y = MAP, 
                 shape = halfNAs)) +
  scale_shape_manual(values = c(19,21), guide = "none") +
  annotate("curve", xend = 1956, x = 1960, yend = 3.8, y = 3, curvature = 0.2,
           arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", x = 1960, y = 3, hjust = 0.2, vjust = 1,
           label = "open circle indicates \nyear with > 50% NAs") +
  geom_smooth() +
  geom_line() + 
  ylab("Mean Annual Precipitation (mm)") +
  xlab("Year (1913-2016)") +
  theme_bw()
MAP_plot
ggsave(filename = here(figures.fp, "MAP_plot.png"), 
       plot = MAP_plot, height = 5, width = 12)


# Seasonal MAP
seasonal_MP_plot <- MAP_df %>%
  ungroup() %>%
  group_by(year, Season) %>%
  mutate(SeasonMean = mean(PrecipMM, na.rm = T)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = SeasonMean, group = Season)) +
  facet_wrap(~SeasonLabels, scales = "free_y") +
  geom_rect(xmin = 1992, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 1, fill = "grey80", color = NA) +
  geom_rect(xmin = 2011, xmax = 2017,
            ymax = Inf, ymin = -Inf, alpha = 0.02, fill = "grey60", color = NA) +
  geom_point() +
  geom_smooth() +
  geom_line(aes(group = Season)) + 
  ylab("Seasonal Mean Precip") +
  xlab("Year (1913-2016)") +
  theme_bw()
seasonal_MP_plot

ggsave(filename = here(figures.fp, "seasonal_MP_plot.png"), 
       plot = seasonal_MP_plot, height = 5, width = 12)

# Cumulative precipitation
CumPrecip_plot <- ggplot(MAP_df, 
                         aes(x = JulianDate, y = CumPrecip)) +
#  geom_point(aes(color = year)) +
  geom_line(aes(color = year, group = year)) +
  scale_color_viridis(name = "Year", direction = -1) +
  scale_x_date(date_breaks = "1 month", date_labels="%B") +
  ylab("Cumulative Precipitation (mm)") +
  xlab("Day of Year") +
  theme_bw()
CumPrecip_plot


# Cumulative precipitation
WkPrecip_plot <- MAP_df %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(wkAvgPrecip = zoo::rollmean(PrecipMM, 7, na.pad = T, align = "right",)) %>%
  ungroup() %>%
  ggplot(aes(x = JulianDate, y = wkAvgPrecip)) +
  geom_line(aes(color = year, group = year)) +
  scale_color_viridis(name = "Year", direction = -1) +
  scale_x_date(date_breaks = "1 month", date_labels="%B") +
  ylab("Precip ove Time (rolling week average, mm)") +
  xlab("Day of Year") +
  theme_bw()
WkPrecip_plot


PrecipTime <- ggplot(MAP_df, 
                         aes(x = JulianDate, y = PrecipMM)) +
  #  geom_point(aes(color = year)) +
  geom_point(aes(color = year)) +
  geom_bar(aes(color = year, group = year), stat = "identity", position = "dodge") +
  scale_color_viridis(name = "Year", direction = -1) +
  scale_x_date(date_breaks = "1 month", date_labels="%B") +
  ylab("Precipitation (mm)") +
  xlab("Date") +
  theme_bw()
PrecipTime
ggsave(filename = here(figures.fp, "PrecipOverYear_plot.png"), 
       plot = PrecipTime, height = 5, width = 12)
#### ====================================================================== ####