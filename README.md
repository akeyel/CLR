# CLR
Supporting R code for the manuscript "SARS-CoV-2 vaccine breakthrough by Omicron and Delta variants: comparative assessments with New York State genomic surveillance data"

## Index to files:
00_Delta_Prelim_Settings.R:            Settings for Delta analysis
00_Omicron_Prelim_Settings.R:          Settings for Omicron analysis
00_Descriptive_Omicron.R:              Basic descriptive analysis of the data set
01_ConditionalLogistic_Omicron.R:      Main Omicron analysis
02_ConditionalLogistic_AGE_Omicron.R:  Omicron analysis without age matching
03_ConditionalLogistic_VACCINATED_Omicron.R: Omicron analysis, restricted to vaccinated only
04_Overview_Figure_CondLogistic.R      Creates the stacked bar chart that shows regional and temporal distribution
05_Descriptive_Delta.R                 Basic descriptive analysis of the Delta data set
06_ConditionalLogistic_Delta.R         Main Delta analysis
07_ConditionalLogistic_AGE_Delta.R     Delta analysis without age matching
08_PowerAnalysis.R                     Power analysis for Delta analysis
DOH_COVID_hlpr.R                       File with helper functions