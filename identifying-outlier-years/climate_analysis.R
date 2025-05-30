library(ggplot2)
library(tsibble)
library(lubridate)
library(dplyr)

source('functions_climate_analysis.R')


# Import Sample Metadata Sheet and extract the relevant sampling dates ----
# use the same version as used elsewhere in this repository (remember to update when this updates)
samplemetadata <- read.csv("../data/coring_geochem_sequenced_samples_1.0.0.csv")

# additional dates from Depth-Info cached query (for weather timeseries plot)
coringinfo <- read.csv("data/cached-queries/cached_query_Depth-Info_20220120.csv", stringsAsFactors = FALSE, na.strings=c("NaN", "", "No core"))

sampling_dates_all <- sort(unique(as.Date(coringinfo$Date__)))  # ALL sampling dates (simple vector) - currently this is only used for the weather timeseries plot

# Get sampling dates for autochamber sites only, July 2011-2017, labeled for each site in a small dataframe
sampling_AC_July_2011to2017 <- unique(samplemetadata[(samplemetadata$CoreGroup__=="MainAutochamber" & 
                                                    samplemetadata$Month__==7 & 
                                                    samplemetadata$Year__>=2011 & samplemetadata$Year__<=2017), 
                                                 c("Site__", "Date__", "Year__")])
sampling_AC_July_2011to2017 <- na.omit(sampling_AC_July_2011to2017)  # remove one row with NAs
sampling_AC_July_2011to2017$Date <- as.Date(sampling_AC_July_2011to2017$Date__)
sampling_AC_July_2011to2017 <- sampling_AC_July_2011to2017[order(sampling_AC_July_2011to2017$Date), c("Date", "Date__", "Year__", "Site__")]  # sort by Date and reorder columns
row.names(sampling_AC_July_2011to2017) <- 1:nrow(sampling_AC_July_2011to2017)  # reindex

sampling_dates_July_2011to2017 <- unique(sampling_AC_July_2011to2017$Date)  # July sampling dates only (same format as sampling_dates_all)

# Import and clean WTD, ALD data ----

# TODO: Also import this data into the graph DB, using some of the logic below (except removing <, >, other "flags") to standardize.
# (Note: May also use a different "formatted" version for the DB with fewer changes to column names, and then put edited names (with units added) in "standardized" version.)

wtd_ald <- read.csv("data/WTD_ALD/Active_Layer_Water_Table_03-17_formatted.csv", stringsAsFactors = FALSE, na.strings=c(NA, ""))

# remove false 0 values (should be NA) in palsa WTD
for (wtd_palsa_col in grep("WTD_palsa", names(wtd_ald), value=TRUE)) {
  corr_wtds <- wtd_ald[[wtd_palsa_col]]
  corr_wtds[corr_wtds==0] <- NA  # for some reason this works even if character "0"
  wtd_ald[[wtd_palsa_col]] <- corr_wtds
}

# Remove other special characters from data columns (all cols after DOY).
# Specific examples were identified based on visual inspection.
# NOTE THAT MANY OF THESE REMOVE INFO ABOUT UNCERTAINTY, SO THE RESULTING "CLEANED" DATA SHOULD BE USED WITH CAUTION

wtd_ald_orig <- wtd_ald  # make copy with original characters included
first_col_index <- which(names(wtd_ald)=="DOY")+1

for (colname in names(wtd_ald)[first_col_index:ncol(wtd_ald)]) {
  newcol <- wtd_ald[[colname]]
  
  # remove or replace specific strings
  newcol <- gsub("ice", "", newcol, fixed=TRUE, useBytes=TRUE)
  newcol <- gsub("Fr.+ga niklas, ", "", newcol, useBytes=TRUE)  # use regex to ensure match of Swedish "å"
  newcol <- gsub("?", "", newcol, fixed=TRUE, useBytes=TRUE)
  newcol <- gsub("--", "-", newcol, fixed=TRUE, useBytes=TRUE)  # one erroneous double "-"
  newcol <- gsub(" ", "", newcol, fixed=TRUE, useBytes=TRUE)
  
  if(grepl("^ALD", colname)) {
    newcol[grepl("inf", newcol, fixed=TRUE)] <- -130  # use deepest "detection limit" for ALD in this sheet
    newcol[grepl("<", newcol, fixed=TRUE)] <- -130
    newcol[grepl(">", newcol, fixed=TRUE)] <- -130  # some "< -X" values erroneously labeled as "> X"
    
    # TODO - Another (possible?) correction: Many of the bog ALDs are "missing", but many of these should probably be changed to -130. Patrick's file just says "later- missing depths either not measured or not measureable," which unfortunately is ambiguous, but the average bog ALDs calculated without these missing values seem too shallow.
    
    # finally, convert to numeric, and then make then all negative (as is the convention for this sheet; but should it be positive to match coring sheet data?)
    newcol <- -1*abs(as.numeric(newcol))
  }
  
  if(grepl("^WTD", colname)) {  # for WTD, just remove <, > and use the value
    newcol <- gsub("<", "", newcol, fixed=TRUE, useBytes=TRUE)
    newcol <- gsub(">", "", newcol, fixed=TRUE, useBytes=TRUE)
    newcol <- as.numeric(newcol)
  }
  
  # finally, update the column in the dataframe
  wtd_ald[[colname]] <- newcol
}

# Calculate average columns for bog
# (i.e., the average across different locations within the bog; see IsoGenieSite_AL_WTD_MapsVisualNotes_200310.pdf for a visual. Other habitats don't have multiple locations, except for Palsa which has "Palsa center" in 2017 only-- but we will ignore that for inter-year consistency)
wtd_ald$WTD_bog_avg <- apply(wtd_ald[, c("WTD_bog_246", "WTD_bog_E", "WTD_bog_F")], 1, mean, na.rm=TRUE)
wtd_ald$ALD_bog_avg <- apply(wtd_ald[, c("ALD_bog_246", "ALD_bog_E", "ALD_bog_F")], 1, mean, na.rm=TRUE)

# Re-assign Date column to "Date" type based on formatted Date__
wtd_ald <- within(wtd_ald, Date <- as.Date(Date__))

# Assign Month, Day, Year 
wtd_ald$Month <- month(wtd_ald$Date)
wtd_ald$Day <- day(wtd_ald$Date)
wtd_ald$Year <- year(wtd_ald$Date)

# Convert 2011-2017 data to tsibble and fill gaps
# (FYI: data from 2006 and prior has duplicate dates)
wtd_ald_tsbl <- as_tsibble(wtd_ald[wtd_ald$Year %in% 2011:2017, ], index=Date)
wtd_ald_tsbl <- fill_gaps(wtd_ald_tsbl, .full = TRUE)
# Re-assign Month, Day, Year (to fill these for NA rows)
wtd_ald_tsbl$Month <- month(wtd_ald_tsbl$Date)
wtd_ald_tsbl$Day <- day(wtd_ald_tsbl$Date)
wtd_ald_tsbl$Year <- year(wtd_ald_tsbl$Date)


# Import ANS weather data ----

weather_ans_15 <- read.csv("data/ANS/ANS_Daily_Wx_Jul84_Mar15rev150918.csv", stringsAsFactors = FALSE)
weather_ans_17 <- read.csv("data/ANS/ANS_Daily_Wx_Jul84_Dec17.txt", stringsAsFactors = FALSE)

weather_ans <- weather_ans_17  # use version that goes to 2017

# correct columns formatted as text, converting non-numeric text strings (structured like "9 9.564") to NA
weather_ans$SoilTemperature_20cm <- as.numeric(weather_ans$SoilTemperature_20cm)
weather_ans$WindDirection <- as.numeric(weather_ans$WindDirection)

# convert Month, Day, Year columns to Date format
weather_ans <- within(weather_ans, Date <- as.Date(paste(Year, Month, Day, sep="-")))


# Clean up duplicate and missing dates in ANS data, and convert to tsibble ----
# Make a copy of the raw weather_ans DF prior to "cleaning"
weather_ans_raw <- weather_ans

# Make sure successive dates are always separated by 1 day

unique(diff(weather_ans$Date))  # Should all be 1, but at least one duplicate date (diff=0) was identified
duplicate_dates <- weather_ans$Date[which(diff(weather_ans$Date)==0)]  # Currently just one date: 2015-06-09

for (date in as.list(duplicate_dates)) {  # as.list ensures that the Date format is preserved; not sure why this should be necessary, but it is
  date_df <- weather_ans[weather_ans$Date==date, ]
  new_row <- sapply(date_df, consolidate_duplicates, simplify=FALSE)
  
  # remove the original duplicate rows and append the new_row
  weather_ans <- weather_ans[-as.numeric(row.names(date_df)), ]
  weather_ans <- rbind(weather_ans, new_row)
}

# get dates back in order again, then re-index
weather_ans <- weather_ans[order(weather_ans$Date), ]
rownames(weather_ans) <- 1:nrow(weather_ans)

# See which date diffs are > 1 (as an FYI; converting to tsibble will correct these).
diffs <- diff(weather_ans$Date)
date_gaps <- diffs[diffs != 1]
names(date_gaps) <- weather_ans[diffs != 1, 'Date']  # these date names refer to the last date before the gap
print(date_gaps)

# Convert to tsibble and fill implicit gaps (gaps in dates, denoted by date_gaps) with NAs
weather_ans_tsbl <- as_tsibble(weather_ans, index=Date)
weather_ans_tsbl <- fill_gaps(weather_ans_tsbl, .full = TRUE)
# Re-assign Month, Day, Year (to fill these for NA rows)
weather_ans_tsbl$Month <- month(weather_ans_tsbl$Date)
weather_ans_tsbl$Day <- day(weather_ans_tsbl$Date)
weather_ans_tsbl$Year <- year(weather_ans_tsbl$Date)
# TODO: Interpolate gaps for calculating averages (and flag them as such)

# Import WeatherHawk weather data ----
weatherhawk_10min <- read.csv("data/WeatherHawk/WeatherHawk_Stordalen_MetData_2013-2019_10_min.csv", stringsAsFactors = FALSE)
weatherhawk_10min[grep("^X", names(weatherhawk_10min))] <- NULL
weatherhawk_10min$SeqDay.1 <- NULL

weatherhawk_10min$DateTime_s <- weatherhawk_10min$Date  # copied datetime strings to new column, so that Date (without time) can be its own column
weatherhawk_10min$Time <- substr(weatherhawk_10min$DateTime_s, start=12, stop=16)
weatherhawk_10min$Date <- as.Date(substr(weatherhawk_10min$DateTime_s, start=1, stop=10))  # OVERWRITES ORIGINAL "Date" COLUMN
weatherhawk_10min$Year <- weatherhawk_10min$year  # capitalized version for merging with same column name in ANS dataset
weatherhawk_10min$DateTime <- ymd_hm(weatherhawk_10min$DateTime_s) # convert datetime strings to POSIXct

# Clean up duplicate rows in weatherhawk_10min and convert to tsibble

duplicated_wh <- duplicates(weatherhawk_10min, index=DateTime)

if(nrow(duplicated_wh)==0) {
  # do nothing
  
} else if (unique(duplicated_wh$DateTime_s)=="2014-10-22 07:30") {  # just one known duplicated time identified in initial version of WeatherHawk data
  new_row <- sapply(duplicated_wh, consolidate_duplicates, simplify=FALSE)
  
  # remove the original duplicate rows and append the new_row
  weatherhawk_10min_cleaned <- 
    rbind(weatherhawk_10min[weatherhawk_10min$DateTime_s != "2014-10-22 07:30", ], 
          new_row)
  
  # re-order and re-index
  weatherhawk_10min_cleaned <- weatherhawk_10min_cleaned[order(weatherhawk_10min_cleaned$DateTime), ]
  rownames(weatherhawk_10min_cleaned) <- 1:nrow(weatherhawk_10min_cleaned)
  
} else {
  stop("New duplicated time(s) identified in weatherhawk_10min. Check duplicated_wh to see which rows are duplicated before proceeding further with corrections.")
}

# convert to tsibble and fill implicit gaps
weatherhawk_10min_tsbl <- as_tsibble(weatherhawk_10min_cleaned, index=DateTime)
weatherhawk_10min_tsbl <- fill_gaps(weatherhawk_10min_tsbl, .full = TRUE)

# Weatherhawk daily weather summaries (with same column names as the ANS weather dataframe) ----
# TODO: THIS IS SLOPPY AS IT DOESN'T ACCOUNT FOR DATES WITH UNEVEN TIME COVERAGE; figure out a better way of calculating dailies!
weatherhawk_daily <- weatherhawk_10min %>% group_by(Date, Year) %>% summarize(AirTemperature = mean(Tair_C, na.rm=TRUE), Maximum_AirTemperature = max(Tair_C, na.rm=TRUE), Precipitation.mm. = sum(Pptn_mm, na.rm=TRUE))

# Merged dataframe for ***DATES >= 2010*** with both ANS and WeatherHawk (WH) data, for plotting the whole timeseries ----
weather <- merge(x=weather_ans[weather_ans$Year>2009,], y=weatherhawk_daily, by="Date", all=TRUE, suffixes=c("_ANS", "_WH"))

# Plot long-term climate timeseries ----

# Overall climate timeseries plot
ts_plot <- plot_climate_data()  # show the plot by simply running "ts_plot" in the console

# Boxplots of years 2010+, omitting last year of ANS data and first & last years of WH data (due to incomplete year)
last_yr_ans <- year(max(weather_ans$Date))
first_yr_wh <- year(min(weatherhawk_daily$Date))
last_yr_wh <- year(max(weatherhawk_daily$Date))

#boxplot(AirTemperature ~ Year, data=weather_ans[weather_ans$Year >= 2010 & weather_ans$Year <= 2014, ], ylim=c(-36, 26), main="ANS")
#boxplot(AirTemperature ~ Year, data=weatherhawk_daily[weatherhawk_daily$Year >=2015 & weatherhawk_daily$Year <= 2018, ], ylim=c(-36, 26), main="Stordalen")


# For year comparison: Create columns for day of growing season (Day.Grow) and days to July sampling for each site (Day.P, Day.S, Day.E) ----
start_grow_md <- "06-01"  # define "growing season" (somewhat arbitrarily) as June 1 thru Sept 30
end_grow_md <- "09-30"
autochamber_sites <- unique(sampling_AC_July_2011to2017$Site__)

weather_ans_tsbl$Day.Grow <- NA
weather_ans_tsbl$Day.P <- NA
weather_ans_tsbl$Day.S <- NA
weather_ans_tsbl$Day.E <- NA

wtd_ald_tsbl$Day.Grow <- NA
wtd_ald_tsbl$Day.P <- NA
wtd_ald_tsbl$Day.S <- NA
wtd_ald_tsbl$Day.E <- NA

for(year in 2011:2017) {
  start_grow <- as.Date(paste0(year, "-", start_grow_md))
  end_grow <- as.Date(paste0(year, "-", end_grow_md))
  growing <- seq(start_grow, end_grow, by="days")
  
  # Note: The code below assumes that any implicit gaps have been filled with NAs (as is now the case for both datasets)
  weather_ans_tsbl[weather_ans_tsbl$Date %in% growing, "Day.Grow"] <- 1:length(growing)
  wtd_ald_tsbl[wtd_ald_tsbl$Date %in% growing, "Day.Grow"] <- 1:length(growing)
  
  for(site in autochamber_sites) {
    sampling_date <- sampling_AC_July_2011to2017[sampling_AC_July_2011to2017$Year__==year & 
                                                   sampling_AC_July_2011to2017$Site__==site, 
                                                 "Date"]
    if(length(sampling_date) != 1) {
      stop(paste0("No unique sampling date found for ", site, " in ", year, "."))
    }
    
    # Fill in days from sampling for this site
    days_from_sampling <- growing - sampling_date
    weather_ans_tsbl[weather_ans_tsbl$Date %in% growing, paste0("Day.", substr(site,1,1))] <- days_from_sampling
    wtd_ald_tsbl[wtd_ald_tsbl$Date %in% growing, paste0("Day.", substr(site,1,1))] <- days_from_sampling
    
    # While here, fill in sampling_AC_July_2011to2017 with the T avg, min, & max for this sampling date 
    sampling_AC_July_2011to2017[sampling_AC_July_2011to2017$Date==sampling_date, "samplingdate_mean_AirTemperature"] <- 
      weather_ans_tsbl[weather_ans_tsbl$Date == sampling_date, "AirTemperature"]
    sampling_AC_July_2011to2017[sampling_AC_July_2011to2017$Date==sampling_date, "samplingdate_min_AirTemperature"] <- 
      weather_ans_tsbl[weather_ans_tsbl$Date == sampling_date, "Minimim_AirTemperature"]  # misspelling from original dataset
    sampling_AC_July_2011to2017[sampling_AC_July_2011to2017$Date==sampling_date, "samplingdate_max_AirTemperature"] <- 
      weather_ans_tsbl[weather_ans_tsbl$Date == sampling_date, "Maximum_AirTemperature"]
  }
}

# Exploratory plots of the 28-day and since-June-1 data for the Sphagnum Autochamber Site
t_comparison_28d_S <- year_comparison(data=weather_ans_tsbl, variable="AirTemperature", site_code="S", start_date=-28)
t_comparison_June1_S <- year_comparison(data=weather_ans_tsbl, variable="AirTemperature", site_code="S", start_date="06-01")
t_comparison_28d_S[['Plot']]
t_comparison_June1_S[['Plot']]

wtd_comparison_28d_S <- year_comparison(data=wtd_ald_tsbl, variable="WTD_bog_avg", site_code="S", start_date=-28)
wtd_comparison_June1_S <- year_comparison(data=wtd_ald_tsbl, variable="WTD_bog_avg", site_code="S", start_date="06-01")
wtd_comparison_28d_S[['Plot']]
wtd_comparison_June1_S[['Plot']]

ald_comparison_july_P <- year_comparison(data=wtd_ald_tsbl, variable="ALD_palsa_135", site_code="P", start_date="07-01", end_date="07-31")
ald_comparison_july_P[['Plot']]
ald_comparison_july_S <- year_comparison(data=wtd_ald_tsbl, variable="ALD_bog_avg", site_code="S", start_date="07-01", end_date="07-31")
ald_comparison_july_S[['Plot']]
ald_comparison_july_E <- year_comparison(data=wtd_ald_tsbl, variable="ALD_fen_78", site_code="E", start_date="07-01", end_date="07-31")
ald_comparison_july_E[['Plot']]

# saw some weird-looking #s for fen growing season ALD; so look at this too to be sure (it appears OK, just showing abrupt thaw in June, after comparing to the original data file)
ald_comparison_allgrowing_E <- year_comparison(data=wtd_ald_tsbl, variable="ALD_fen_78", site_code="E", start_date="06-01", end_date="09-30")
ald_comparison_allgrowing_E[['Plot']]

# Generate temperature & WTD summary stats for *all* sites and time intervals (7, 14, 21, and 28 days; and since June 1) ----

intervals <- c("7d", "14d", "21d", "28d", "july", "growing", "all_growing")
wtd_intervals <- intervals[!grepl('^7d', intervals) & !grepl('^14d', intervals)]  # same intervals as for T, but removing short intervals

temperature_variable <- "AirTemperature"
summary_statistic_types <- c("n", "median", "mean", "sd")
t_summary_colnames <- paste0(summary_statistic_types, "_", temperature_variable)
wtd_summary_colnames <- paste0(summary_statistic_types, "_", "WTD")  # manually define 1 column name suffix for WTD 
ald_summary_colnames <- paste0(summary_statistic_types, "_", "ALD") 
t_wtd_ald_summaries <- NULL

for (interval in intervals) {

  interval_df <- sampling_AC_July_2011to2017
  interval_df$Interval <- interval
  names(interval_df)[names(interval_df)=="Date"] <- "Sampling_Date"
  
  if(grepl("[0-9]d$", interval)) {
    start_date_arg <- -1*as.numeric(gsub("d", "", interval))
    interval_df$Start_Date <- interval_df$Sampling_Date + start_date_arg  # added (not subtracted) because it's negative
    end_date_arg <- FALSE
    interval_df$End_Date <- interval_df$Sampling_Date
    
  } else if(interval=="growing") {
    start_date_arg <- "06-01"
    interval_df$Start_Date <- paste0(interval_df$Year__, "-", start_date_arg)
    end_date_arg <- FALSE
    interval_df$End_Date <- interval_df$Sampling_Date
    
  } else if(interval=="all_growing") {
    start_date_arg <- "06-01"
    interval_df$Start_Date <- paste0(interval_df$Year__, "-", start_date_arg)
    end_date_arg <- "09-30"
    interval_df$End_Date <- paste0(interval_df$Year__, "-", end_date_arg)
    
  } else if(interval=="july") {
    start_date_arg <- "07-01"
    interval_df$Start_Date <- paste0(interval_df$Year__, "-", start_date_arg)
    end_date_arg <- "07-31"
    interval_df$End_Date <- paste0(interval_df$Year__, "-", end_date_arg)
    
  } else {
    stop(paste0("Invalid interval string: ", interval))
  }
  
  # allocate empty columns for summary stats
  for (colname in t_summary_colnames) {
    interval_df[[colname]] <- NA
  }
  for (colname in wtd_summary_colnames) {
    interval_df[[colname]] <- NA
  }
  for (colname in ald_summary_colnames) {
    interval_df[[colname]] <- NA
  }
  
  # re-order rows and columns
  interval_df <- interval_df[order(interval_df$Site__, interval_df$Year__), c("Site__", "Year__", "Interval", "Start_Date", "End_Date", "Sampling_Date", "samplingdate_mean_AirTemperature", "samplingdate_min_AirTemperature", "samplingdate_max_AirTemperature", t_summary_colnames, wtd_summary_colnames, ald_summary_colnames)]
  
  # loop through sites, then calculate and fill in summary stats
  for(site in autochamber_sites) {
    # start with T
    site_interval_t_results <- year_comparison(data=weather_ans_tsbl, variable=temperature_variable, 
                                             site_code=substr(site, 1, 1), start_date=start_date_arg,
                                             end_date=end_date_arg)
    site_interval_t_summary <- site_interval_t_results$Summary
    # fill values (NOTE: ROWS AND SUMMARY STATISTIC COLUMNS MUST BE IN THE SAME ORDER, which they should be here, due to ordering by year, and due to how t_summary_colnames and summary_statistic_types are defined (static).)
    interval_df[interval_df$Site__==site, t_summary_colnames] <- site_interval_t_summary[, summary_statistic_types]
    
    
    # Do the same for WTD (WTD variable needs to be defined)
    # (note: there are actually two WTD columns for palsa, WTD_palsa_135 and WTD_palsa_center, but both have all NAs)
    wtd_variable <- switch(substr(site, 1, 1), "P"="WTD_palsa_135", "S"="WTD_bog_avg", "E"="WTD_fen_78")
    site_interval_wtd_results <- year_comparison(data=wtd_ald_tsbl, variable=wtd_variable, 
                                               site_code=substr(site, 1, 1), start_date=start_date_arg,
                                               end_date=end_date_arg)
    site_interval_wtd_summary <- site_interval_wtd_results$Summary
    # fill values (NOTE: ROWS AND SUMMARY STATISTIC COLUMNS MUST BE IN THE SAME ORDER, which they should be here, due to ordering by year, and due to how t_summary_colnames and summary_statistic_types are defined (static).)
    interval_df[interval_df$Site__==site, wtd_summary_colnames] <- site_interval_wtd_summary[, summary_statistic_types]
    
    
    # Finally do the same for ALD (same process as for WTD)
    # (note: there are actually two ALD columns for palsa, ALD_palsa_135 and ALD_palsa_center, but ALD_palsa_center only exists for 2017 so can't be compared across years)
    ald_variable <- switch(substr(site, 1, 1), "P"="ALD_palsa_135", "S"="ALD_bog_avg", "E"="ALD_fen_78")
    site_interval_ald_results <- year_comparison(data=wtd_ald_tsbl, variable=ald_variable, 
                                                 site_code=substr(site, 1, 1), start_date=start_date_arg,
                                                 end_date=end_date_arg)
    site_interval_ald_summary <- site_interval_ald_results$Summary
    # fill values (NOTE: ROWS AND SUMMARY STATISTIC COLUMNS MUST BE IN THE SAME ORDER, which they should be here, due to ordering by year, and due to how t_summary_colnames and summary_statistic_types are defined (static).)
    interval_df[interval_df$Site__==site, ald_summary_colnames] <- site_interval_ald_summary[, summary_statistic_types]
  }
  
  # combine into main summary table
  t_wtd_ald_summaries <- rbind(t_wtd_ald_summaries, interval_df)
}

# Generate long plots of entire growing season T and WTD (bog site) for all years.

# The below assumes that interval=="all_growing" and site=="Sphagnum Autochamber Site" on the last iteration of the above loop.
all_T <- site_interval_t_results$Data
all_WTD <- site_interval_wtd_results$Data

# subset of July sampling dates for Sphagnum site only
sampling_dates_July_2011to2017_S <- unique(sampling_AC_July_2011to2017[sampling_AC_July_2011to2017$Site__=="Sphagnum Autochamber Site", "Date"])

# make long plots
plot_long_T <- long_plot(all_T, variable="AirTemperature", sampling_dates=sampling_dates_July_2011to2017)
plot_long_WTD_S <- long_plot(all_WTD, variable="WTD_bog_avg", sampling_dates=sampling_dates_July_2011to2017_S)

# export results
write.csv(t_wtd_ald_summaries, file='results/t_wtd_ald_summaries_July2011-2017samplings_version-20240506.csv', row.names=FALSE)

# Generate % of time below WTD for each interval for each sample. ----

# Start with metadata sheet, selecting the subsets of relevant data.
wtd_summaries <- samplemetadata[(samplemetadata$CoreGroup__=="MainAutochamber" & 
                                   samplemetadata$Month__==7 & 
                                   samplemetadata$Year__>=2011 & samplemetadata$Year__<=2017), 
                                c('DepthID', 'SampleID__', 'SampleID_old__', 'Month__', 'Year__', 'Date__', 
                                  'Site__', 'Core__', 'DepthAvg__', 'ALD.cm__', 'WTD.cm_neg_is_below_sfc__', 
                                  'T_air.deg_C', 'T_soil.deg_C')]

# Generate nested list of WTD measurements for each site, interval, and year, to be used when filling wtd_summaries.

wtd_measurement_list <- vector(mode="list", length=length(autochamber_sites))
names(wtd_measurement_list) <- autochamber_sites

for(site in autochamber_sites) {
  wtd_measurement_list[[site]] <- vector(mode="list", length=length(wtd_intervals))
  names(wtd_measurement_list[[site]]) <- wtd_intervals
  
  # Choose the appropriate WTD column name
  # (note: there are actually two WTD columns for palsa, WTD_palsa_135 and WTD_palsa_center, but both have all NAs)
  wtd_variable <- switch(substr(site, 1, 1), "P"="WTD_palsa_135", "S"="WTD_bog_avg", "E"="WTD_fen_78")
  
  for(interval in wtd_intervals) {
    
    if(grepl("[0-9]d$", interval)) {
      start_date_arg <- -1*as.numeric(gsub("d", "", interval))
      end_date_arg <- FALSE
      
    } else if (interval=="growing") {
      start_date_arg <- "06-01"
      end_date_arg <- FALSE
      
    } else if(interval=="all_growing") {
      start_date_arg <- "06-01"
      end_date_arg <- "09-30"
      
    } else if(interval=="july") {
      start_date_arg <- "07-01"
      end_date_arg <- "07-31"
      
    } else {
      stop(paste0("Invalid interval string: ", interval))
    }
    
    # Get relevant yearly data
    site_interval_data <- year_comparison(data=wtd_ald_tsbl, variable=wtd_variable, 
                                          site_code=substr(site, 1, 1), start_date=start_date_arg,
                                          end_date=end_date_arg)$Data
    
    yearly_list <- lapply(setNames(2011:2017, as.character(2011:2017)), 
                          function(x) site_interval_data[site_interval_data$Year==x, wtd_variable])
    yearly_list <- lapply(yearly_list, function(x) as.vector(x[!is.na(x)]))  # drop extra dimensions and NAs
    
    wtd_measurement_list[[site]][[interval]] <- yearly_list
  }
}

# fill in wtd_summaries
for(interval in wtd_intervals) {
  wtd_summaries[[paste0('n_WTD_', interval)]] <- 
    apply(wtd_summaries, 1, pct_time_below_wtd, interval=interval, wtd_measurement_list=wtd_measurement_list, return_value='n')
  
  wtd_summaries[[paste0('pct_time_below_WTD_', interval)]] <- 
    apply(wtd_summaries, 1, pct_time_below_wtd, interval=interval, wtd_measurement_list=wtd_measurement_list, return_value='pct')
}

# export results
write.csv(wtd_summaries, file='results/wtd_summaries_July2011-2017samples_version-20240506.csv', row.names=FALSE)