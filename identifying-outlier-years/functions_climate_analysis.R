# Function for aggregating timeseries data (tsibble) across coarser timescales (e.g., from 10-min data to daily summaries), in a way that accounts for missing data by requiring a minimum number of observations to produce a non-NA result.
aggreggate_timeseries <- function(data, na_max, na_max_consec, fun) {
  # data            input data (tsibble)
  # na_max          maximum allowed number of total NAs
  # na_max_consec   maximum allowed number of consecutive NAs
  # fun             aggregation function (mean or sum)
  
  # Started this function and kept as placeholder, but may not end up using it as-is.
  
}

# Function for consolidating duplicate rows (to be used with sapply)
consolidate_duplicates <- function(dup_data) {
  # dup_data      vector containing values from a single column of a "duplicates" dataframe (i.e., containing duplicated data for single date-time index)
  
  # first get unique non-NA values
  if(all(is.na(dup_data))) {
    unique <- NA
  } else {
    # If there's data there, remove NAs to get only the unique data-containing values -- this eliminates the need to use na.rm=TRUE below.
    dup_data <- na.omit(dup_data)
    unique <- unique(dup_data)
  }
  
  # if only 1 unique value, return the value; otherwise calculate the mean if it's numeric
  if(length(unique)==1) {
    return(unique)
  } else if(is.numeric(unique)) {
    return(mean(unique))
  } else {
    print(unique)
    stop("Could not consolidate these values.")
  }
}


# Function for calculating % of time below WTD for a given sample, along with WTD measurement counts used for calculating this.
# (to be used with apply over rows of dataframe)
pct_time_below_wtd <- function(row, interval, wtd_measurement_list, return_value) {
  # row                     Single row of a dataframe (subset of Sample Metadata Sheet). Must include:
  #                         Site__, Year__, DepthAvg__.
  # interval                Time interval (string) for subsetting wtd_measurement_list.
  # wtd_measurement_list    Nested list of WTD measurements from different sites, time intervals, and years.
  # return_value            Which value to return. Allowed options are "n" (for WTD measuremnet counts over 
  #                         the specified interval), or "pct" (for % of time below WTD over the interval)
  
  wtd_values <- wtd_measurement_list[[row['Site__']]][[interval]][[as.character(row['Year__'])]]
  meas_count <- length(wtd_values)
  
  if(return_value=="n") {
    return(meas_count)
    
  } else if(return_value=="pct") {
    if(grepl("Palsa Autochamber Site", row['Site__']) & meas_count==0) {
      percent_below_wtd <- 0
    } else {
      percent_below_wtd <- sum(-1*as.numeric(row['DepthAvg__']) <= wtd_values, na.rm=TRUE)/meas_count  # inclusive of depth=WTD, for consistency with Carrie's method
    }
    return(percent_below_wtd)
    
  } else {
    stop(paste0("Invalid return_value type: ", return_value))
  }
}

# "Wrapper" function for plotting and generating summary statistics for climate data over a given date range 
# ahead of each sampling date, compared between years
year_comparison <- function(data, variable, site_code, start_date, end_date=FALSE) {
  # data              Data timeseries (tsibble).
  # variable          Which variable to plot (string).
  # site_code         Letter code of the site (string: P, S, or E) for which to get date ranges.
  # start_date        Start for the date range to plot/summarize. Can be either a specific day of year  
  #                   (string, MM-DD), or # of days (negative number) prior to the sampling date. All dates  
  #                   within the resulting range must be within the defined growing season (June 1 - Sept 30).
  #                   The type of start date specified (exact date string, or # of days) will also affect how 
  #                   the plot's x-axis is defined.
  # end_date          End of the date range to plot/summarize (string, MM-DD). Used to optionally override 
  #                   the default (end_date=FALSE) use of the sampling date as the end date. This option 
  #                   currently only works if start_date is provided in this same format (string, MM-DD), 
  #                   and end_date must be within the defined growing season (June 1 - Sept 30).
  
  day_site_colname <- paste0("Day.", site_code)  # Day.P, Day.S, or Day.E
  days_from_sampling <- data[, day_site_colname]
  
  # Subset the relevant dates.
  # The checks below assume that start_date is already in a valid format.
  
  if(is.numeric(start_date)) {  # number of days from sampling
    data_to_compare <- 
      data[days_from_sampling >= start_date & days_from_sampling <= 0 & !is.na(days_from_sampling), ]
    
    x_axis_col <- day_site_colname
    
  } else if(is.character(start_date)) {  # specific start date
    start_grow_day <- as.Date(paste0("2010-", start_date)) - as.Date("2010-05-31")  # convert start_date to day of growing season (using arbitrary year)
    
    if(end_date==FALSE) {
      data_to_compare <- 
        data[data$Day.Grow >= start_grow_day & days_from_sampling <= 0 & !is.na(days_from_sampling), ]
      
    } else if(is.character(end_date)) {
      end_grow_day <- as.Date(paste0("2010-", end_date)) - as.Date("2010-05-31")
      data_to_compare <- 
        data[data$Day.Grow >= start_grow_day & data$Day.Grow <= end_grow_day & !is.na(data$Day.Grow), ]
      
    } else {
      stop("Invalid end_date format.")
    }
    
    x_axis_col <- "Day.Grow"
    
  } else {
    return("Invalid start_date format.")
  }
  
  # Copy of Year column as factor (for plotting)
  data_to_compare$Year_f <- factor(data_to_compare$Year)
  
  # Colorblind-friendly palette (from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/, reordered to look more continuous, and removing the hard-to-see yellow)
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cbbPalette <- c("#000000", "#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7")
  
  # For comparison plot, omit data with NAs so that points connect.
  comparison_plot <- ggplot(data_to_compare[!is.na(data_to_compare[, variable]), ], aes_string(x=x_axis_col, y=variable, group="Year_f", color="Year_f")) +
    theme_bw() + geom_line() + geom_point() + scale_color_manual(values=cbbPalette) 
  #+ scale_color_manual(values=c("#AA0000", "#FF3300", "#FFAA00", "#CCCC00", "#00FF00", "#0088AA", "#0000FF")) #alternate custom rainbow scale
  
  # Generate summary statistics for each year.
  # This only works with as.data.frame() to remove the tsibble attribute. Integer grouping variables are OK.
  summary_stats <- as.data.frame(data_to_compare) %>% group_by(Year) %>% 
    summarise(site_code=site_code,
              start_date=start_date,
              variable=variable,
              n=sum(!is.na(.data[[variable]])),
              median=median(.data[[variable]], na.rm=TRUE),
              mean=mean(.data[[variable]], na.rm=TRUE),
              sd=sd(.data[[variable]], na.rm=TRUE))
  
  return(list(Data=data_to_compare, Summary=summary_stats, Plot=comparison_plot))
}


# Function for making long plots of growing season data (separate panels for each year, arranged horizontally), annotated with sampling dates.
# This function is designed to work with dataframes exported via the "Data" item in the year_comparison function above.
long_plot <- function(data, variable, sampling_dates=FALSE) {
  
  plot <- ggplot(data, aes_string(x="Day.Grow", y=variable)) + theme_bw() + geom_line() + geom_point() + facet_wrap(vars(Year_f), nrow=1)
  
  if(is.Date(sampling_dates)) {
    
    sampling_dates_df <- data[data$Date %in% sampling_dates, c("Year_f", "Date", "Day.Grow")]
    
    plot <- plot + geom_vline(aes_string(xintercept="Day.Grow"), data=sampling_dates_df, linetype = "dashed", color="blue")
  }
  return(plot)
  
}


# "Wrapper" function for plotting the climate data defined in climate_analysis.R, used for specifying which 
# subsets of data to plot. Default is to plot everything.
plot_climate_data <- function(plot_ANS=TRUE, plot_WH=TRUE, 
                              plot_Tmax=TRUE, plot_Tmean=TRUE,
                              plot_precip=TRUE, plot_samplings=TRUE) {
  
  colors_ts_plot <- c("orange", "red", "gray", "black", "turquoise", "navyblue")
  names(colors_ts_plot) <- c("Max Daily T (ANS)",
                             "Max Daily T (Stordalen)",
                             "Mean Daily T (ANS)",
                             "Mean Daily T (Stordalen)",
                             "Precipitation (ANS)",
                             "Precipitation (Stordalen)")
  
  ts_plot <- ggplot(weather, aes(x=Date)) + theme_bw() 
  
  if(plot_samplings) {
    # sampling dates (porewater or solid phase)
    # filter to only include years within the weather date range
    end_weather <- as.Date(paste0(year(max(weather$Date)), "-12-31"))
    sampling_dates_in_range <- sampling_dates_all[sampling_dates_all <= end_weather]
    
    ts_plot <- ts_plot + annotate("point", x=sampling_dates_in_range, y=30, color="blue")
  }
  
  if(plot_ANS) {
    
    if(plot_Tmax) {
      ts_plot <- ts_plot + geom_line(aes(y=Maximum_AirTemperature_ANS, color="Max Daily T (ANS)"))
    }
    if(plot_Tmean) {
      ts_plot <- ts_plot + geom_line(aes(y=AirTemperature_ANS, color="Mean Daily T (ANS)"))
    }
    if(plot_precip) {
      ts_plot <- ts_plot + geom_line(aes(y=Precipitation.mm._ANS - 50, color="Precipitation (ANS)"))
    }
    
  }
  
  if(plot_WH) {
    if(plot_Tmax) {
      ts_plot <- ts_plot + geom_line(aes(y=Maximum_AirTemperature_WH, color="Max Daily T (Stordalen)"))
    }
    if(plot_Tmean) {
      ts_plot <- ts_plot + geom_line(aes(y=AirTemperature_WH, color="Mean Daily T (Stordalen)"))
    }
    if(plot_precip) {
      ts_plot <- ts_plot + geom_line(aes(y=Precipitation.mm._WH - 50, color="Precipitation (Stordalen)"))
    }
  }
  
  ts_plot <- ts_plot + scale_x_date(date_breaks = "1 year", date_labels = "%Y", name="Year")
  
  # TODO: Remove unused axes and legend entries based on what's plotted.
  ts_plot <- ts_plot + scale_y_continuous(name = "Air Temperature (\u00b0C)", 
                                          sec.axis = sec_axis(trans = ~.+50, name = "Precipitation (mm)"))
  ts_plot <- ts_plot + scale_color_manual(values = colors_ts_plot)
    
  ts_plot <- ts_plot + labs(color = "Climate Data")
  
  return(ts_plot)
}