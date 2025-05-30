---
Author:  Suzanne Hodgkins, EMERGE Institute
Date of this version:  April 25, 2022
---

# Temperature and water-table-depth summaries

(filepaths below are relative to `temporal_paper/identifying-outlier-years` unless otherwise noted)

## OVERVIEW

This directory (`results/`) contains the
results of running the script: `climate_analysis.R`
and exporting the temperature (T) and water table depth (WTD) summaries to CSV.
Descriptions of each of these summary files are given below.

## t_wtd_summaries_July2011-2017samplings.csv

This file gives summary statistics for the following:

- mean daily air temperature ("AirTemperature", from the file `data/ANS/ANS_Daily_Wx_Jul84_Dec17.txt`) measured at the Abisko Scientific Research Station (ANS)
- water table depths ("WTD", from the file `data/WTD_ALD/Active_Layer_Water_Table_03-17_formatted.csv`) measured manually by Patrick Crill's group at the autochamber sites at Stordalen

...over various time intervals defined relative to the sampling date for each site (defined separately for each site, due to the sites having slightly different sampling dates in some years).

The specific time intervals used are defined as follows:

- 7d (7 days prior to the sampling date, plus the sampling date itself)
- 14d (14 days prior to the sampling date, plus the sampling date itself)
- 21d (21 days prior to the sampling date, plus the sampling date itself)
- 28d (28 days prior to the sampling date, plus the sampling date itself)
- growing (Time from beginning of growing season (defined as June 1) until (and including) the sampling date)
- all_growing (Entire growing season (June 1 – Sept. 30))

For clarity, the start and end dates for each time interval (inclusive) are also given under the columns Start_Date and End_Date (End_Date=Sampling_Date for all intervals except all_growing).

Summary statistics for each interval include: measurement count (n), median (median), mean (mean), and standard deviation (sd), and are given under the column names beginning with these statistic labels.

**IMPORTANT NOTE**: For temperature, these statistics are calculated based on the average temperature measured on each day, meaning that the standard deviations do NOT account for within-day temperature variation. To get values that account for full diurnal T variation, we would need the hourly data, which we currently don't have for the years being analyzed (2011-2017). For now, to provide short-term (1 day) T variation context for each sampling date, the within-day mean, minimum, and maximum air temperatures for the sampling date only (taken directly from the corresponding row & columns in the source ANS data file above) are provided in the columns samplingdate_mean_AirTemperature, samplingdate_min_AirTemperature, and samplingdate_max_AirTemperature.

## wtd_summaries_July2011-2017samples.csv

This file gives the percent of time that the DepthAvg__ for each sample was at or below the WTD, and the number of WTD measurements used for calculating this (from Patrick Crill's manual WTD measurements in the file `data/WTD_ALD/Active_Layer_Water_Table_03-17_formatted.csv`), for the longer subset of the same time intervals (≥21 days) used for the temperature & WTD summaries above. (Note that unlike T, WTD is usually measured only every few days; therefore the shorter time intervals were omitted due to low n.)

The first few columns are taken directly from the Sample Metadata Sheet (using the version of the file at `temporal_paper/data/coring_geochem_sequenced_samples_1.0.0.csv`), for the samples collected in July of 2011-2017 from the main autochamber sites. The last set of columns include the following, with the time interval labels (definitions are the same as for temperature summaries) appended at the end of each column name:

- n_WTD_* (Number of WTD measurements used in calculation)
- pct_time_below_WTD_* (Fraction of measured WTDs over the given time interval that were at or above the DepthAvg__ for each sample, which equates to the fraction of measurement timepoints during which the given sample was at or below the WTD. This is the same method used for calculating "% Time below water table" in Figure 6 of Singleton et al. (2018). For palsa sites, this is automatically set to 0 based on the lack of a water table at all timepoints in the analysis.)

## Newer versions of the above files (ending in "_version-20240506.csv")

On May 6 2024, the script was amended to also include the "july" time interval (including all July dates for each year), and to also include "active layer depth" (ALD) in the summary statistics file (now renamed to "t_wtd_ald_summaries_..."). To distinguish from the original results, these newer files have the version date appended to the filename, although the same summary statistics should be identical between versions.

**IMPORTANT NOTE**: The summary statistics for bog ALD treat "missing" values as truly missing, even though many of them were likely deeper than detection, so the calculated average bog ALDs are likely too shallow. Patrick's original file just has a note saying "later- missing depths either not measured or not measureable", which is ambiguous. This needs to be corrected in a future version.