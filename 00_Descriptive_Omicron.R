# Goal: Analyze Omicron data using the case-control approach used for Delta and VOC's

# Created 2021-12-15
# Updated 2022-03-22

# Author: Alexander "Sasha" Keyel <alexander.keyel@health.ny.gov

# Load required standard libraries
library(lubridate)

# Set working directory
setwd("C:/hri/DOH_COVID") # Adjust as necessary
# Load helper functions
source("R/DOH_COVID_hlpr.R") # Requires dfmip package, devtools::install_github('akeyel/dfmip')

# Load consistent set of labels for all Omicron analyses
source("R/00_Omicron_Prelim_Settings.R") 

# Data Pre-processing Instructions
# Data were in .xlsx format.
# A find-and-replace was used to replace ',' with ';' 101770 replacements
# 6495 (formerly 8187) blank vaccine type entries were changed to "None" for 'Not Vaccinated'
# Days since vaccine were left as blank - preserves numeric nature of this column, and NA is appropriate here
# 10309 (formerly 10356) blank booster type entries were changed to "None" for "Not Boosted"
# A .csv file was created by saving the Excel file as a .csv
DOH.data.file = "Data/omicron_delta_2022_03_02.csv"

# NOTES:
# In three cases, the same patient was included in both the case and control data sets
# (unvaccinated, then vaccinated infection)
# Note that they would not be matched to themselves, must be within 6 days, and the vaccination series takes longer than that.
# Identifier 1745 + 1749
# Identifier 1558 + 4720
# Identifier 1294 + 4281

# Read in the input data file
new.data = read.csv(DOH.data.file)

# Adjust data set based on any errors from the data.setup.v2 function below
# head(new.data) # Make sure position did not change
colnames(new.data)[7] = "sequence_result"
colnames(new.data)[11] = "booster_type" # Change to be consistent with expected formatting
colnames(new.data)[17] = "Vaccine.date.5" # Change to be consistent with others

table(new.data$Age)
# Convert ages in months to years
month.index = grep(" month", new.data$Age)
new.data$Age_Months = NA
new.data$Age_Months[month.index] = sapply(new.data$Age[month.index], dfmip::splitter, ' ', 1, 0)
new.data$Age[month.index] = floor(new.data$Age_Months[month.index] / 12)
new.data$Age_Months = NULL # Drop added intermediate field

# Convert ages in days to years (assume all ages in day are less than 1 year)
day.index = grep(' day', new.data$Age) # grep will match to day and days with this syntax. In this case, that is desirable.
new.data$Age[day.index] = 0
new.data$Age = as.numeric(new.data$Age)

# This block is left over from a previous version of the data. Not sure why the data pattern changed.
# # Ages in 200's: Drop the leading 2 to get the age in months
# new.data$Age[new.data$Age >= 200 & new.data$Age < 300] = floor((new.data$Age[new.data$Age >= 200 & new.data$Age < 300] - 200)/ 12)

# # Ages in 300's: Drop the leading 3 to get the age in days
# new.data$Age[new.data$Age >= 300 & new.data$Age < 400] = 0 # Age in days would need to be > 365 to be > 1 year

# Correct sex entries
new.data$Sex[new.data$Sex == "m"] = "M"
new.data$Sex[new.data$Sex == "f"] = "F"

#### FILTER DATA BASED ON NUMBER OF DOSES.
# Convert date fields to date in R
new.data$collection_date_dateformat = sapply(new.data$collection_date, as.Date, format = c("%m/%d/%Y"))
new.data$vd1 = sapply(new.data$Vaccine.date.1, as.Date, format = c("%m/%d/%Y"))
new.data$vd2 = sapply(new.data$Vaccine.date.2, as.Date, format = c("%m/%d/%Y"))
new.data$vd3 = sapply(new.data$Vaccine.date.3, as.Date, format = c("%m/%d/%Y"))
new.data$vd4 = sapply(new.data$Vaccine.date.4, as.Date, format = c("%m/%d/%Y"))
new.data$vd5 = sapply(new.data$Vaccine.date.5, as.Date, format = c("%m/%d/%Y"))


# Check that all individuals with Vaccine type 'None' have a vaccination date ON or AFTER their infection date (if applicable)
# Individuals with a positive test on date of vaccination would have been infected prior to vaccination.
test = new.data$Identifier[!is.na(new.data$vd1) & new.data$vaccine_type == "None" & new.data$vd1 < new.data$collection_date_dateformat]
length(test)
# Two were vaccinated on the date they were infected (now no longer pulled out)
# 9 were vaccinated once with vaccine type listed as 'None', these should be dropped

# Do any of these have a second dose before infection?
test2 = new.data$Identifier[!is.na(new.data$vd2) & new.data$vaccine_type == "None" & new.data$vd2 <= new.data$collection_date_dateformat]
# No, none of these individuals remained in the data set

# Remove or flag individuals with one dose of Pfizer or Moderna at time of infection
test3 = new.data$Identifier[is.na(new.data$vd2) & new.data$vaccine_type == "Pfizer"] # No individuals with only a single Pfizer shot
test4 = new.data$Identifier[is.na(new.data$vd2) & new.data$vaccine_type == "Moderna"] # No individuals with only a single Moderna shot

# Changed criterion - if individual was infected prior to vaccination, there is no problem here as long as they are treated as unvaccinated.
test5 = new.data$Identifier[new.data$vaccine_type == "Pfizer" & new.data$collection_date_dateformat > new.data$vd1 & new.data$collection_date_dateformat <= (new.data$vd2 + 14)]
test6 = new.data$Identifier[new.data$vaccine_type == "Moderna" & new.data$collection_date_dateformat > new.data$vd1 & new.data$collection_date_dateformat <= (new.data$vd2 + 14)] # individual 1500 needs to be removed - only 13 days, perhaps why sequence.result was 'None'
test6b = new.data$Identifier[new.data$vaccine_type == "Janssen" & new.data$collection_date_dateformat > new.data$vd1 & new.data$collection_date_dateformat <= (new.data$vd1 + 14)]

# Remove or flag individuals with 4+ doses at time of infection
test7 = new.data$Identifier[!is.na(new.data$vd4) & new.data$collection_date_dateformat >= (new.data$vd4)]
# 16 need to be removed from the data set 

# Remove or flag individuals with 3rd shot < 135 days after 2nd shot (cross-check against individuals with 3rd shot prior to August)
test8 = new.data$Identifier[!is.na(new.data$vd3) & new.data$vd3 <= (new.data$vd2 + 135)] # 58 individuals
shot3.file = sprintf("%s/Individuals_3rd_shot_less_than_135_days%s.csv", out.folder, figure.label)
cat("IDS:,", file = shot3.file, append = FALSE, sep = ',')
cat(test8, file = shot3.file, sep = ',', append = TRUE)
cat(sprintf("\nN:,%s", length(test8)), sep = ',', append = TRUE, file = shot3.file)

# Remove individuals with a 2nd shot <135 days after a 1st shot of Janssen vaccine
test8b = new.data$Identifier[new.data$vaccine_type == "Janssen" & !is.na(new.data$vd2) & new.data$vd2 <= (new.data$vd1 + 135)]
# 4 individuals

# Get date number for Aug 1, 2021
boost.start = as.Date("8/1/2021", format = "%m/%d/%Y")
test9 = new.data$Identifier[!is.na(new.data$vd3) & new.data$vd3 < as.numeric(boost.start)]
length(test9) # Only 52 by the Aug. 1 criterion. Use 135 day based criterion instead.

# Drop individuals with a mixed primary vaccine series
mixed = new.data$Identifier[tolower(new.data$vaccine_type) == "mixed"]

# Rounding to remove the fractional days. 
# Days.between.collection.date.and.booster.date is missing in the updated data, but was calculated below. Change field name in calculation below.
#new.data$Days.between.collection.date.and.booster.date = floor(new.data$Days.between.collection.date.and.booster.date)
# Subtracting in R shows that floor is the appropriate way to round. Matched with R's calculations, but I didn't need an extra
# step to process the negatives, so I left it with the Excel calculation.
new.data$Days.between.collection.date.and.booster.date = new.data$collection_date_dateformat - new.data$vd3
new.data$Days.between.collection.date.and.booster.date[new.data$vaccine_type == "Janssen"] = (new.data$collection_date_dateformat[new.data$vaccine_type == "Janssen"]
                                                       - new.data$vd2[new.data$vaccine_type == "Janssen"])

new.data$Days.between.collection.date.and.booster.date = new.data$collection_date_dateformat - new.data$vd3 
new.data$Days.between.collection.date.and.booster.date[new.data$vaccine_type == "Janssen"] = (new.data$collection_date_dateformat[new.data$vaccine_type == "Janssen"]
                                                        - new.data$vd2[new.data$vaccine_type == "Janssen"])

# Do not consider an individual boosted if the collection date is before the booster shot
new.data$Days.between.collection.date.and.booster.date[new.data$Days.between.collection.date.and.booster.date < 1] = NA
new.data$time_between_vaccination_and_booster = new.data$Days.between.collection.date.and.vaccination.complete.date -
  new.data$Days.between.collection.date.and.booster.date

# Drop individuals with a lineage assigned of 'None'
nones = new.data$Identifier[new.data$sequence_result == "None"] #5 entries had a 'None' value here

# Remove pending individuals from the analysis
#pending.index = grep("pending", tolower(in.data$sequence_result)) # None were lowercase, but figured I'd do it anyhow # 146
new.data$sequence_result[new.data$sequence_result == "Pending "] = "Pending" # Remove trailing space
pendings = new.data$Identifier[new.data$sequence_result == "Pending"]

# EXCLUDE INDIVIDUALS IDENTIFIED IN TESTS ABOVE
# 8299 added - identified as a duplicate by visual inspection of the ID's highlighted in the possible duplicates table generated below.
duplicates = c("8299")
drops = unique(c(test, test6, test6b, test7, test8, test8b, mixed, nones, pendings, duplicates)) # Dropping 232 records

drop.df = new.data[new.data$Identifier %in% drops, ]
dropped.rec.file = sprintf("%s/DroppedRecords%s.csv", out.folder, figure.label)
write.table(drop.df, file = dropped.rec.file, sep = ',', col.names = TRUE, row.names = FALSE)
new.data = new.data[!new.data$Identifier %in% drops, ] #10633 records

# correct inconsistent labeling of controls
new.data$case_control[new.data$case_control == "Control "] = "Control"
new.data$case_control[new.data$case_control == "Case "] = "Case"

# correct missing regions with information
# Kings/Brooklyn, so did not match a region, assign to NYC
new.data$Economic_region[new.data$Economic_region == "#N/A"] = "New York City"

# Add VACCINE_DOSES field
# NOTE THAT SINGLE-SHOT PFIZER & MODERNA INDIVIDUALS WERE EXCLUDED
new.data$VACCINE_DOSES = 0
new.data$VACCINE_DOSES[new.data$case_control == "Case" & new.data$vaccine_type != "Janssen"] = 2
new.data$VACCINE_DOSES[new.data$case_control == "Case" & new.data$vaccine_type == "Janssen"] = 1
new.data$VACCINE_DOSES[!is.na(new.data$Days.between.collection.date.and.booster.date) & new.data$vaccine_type != "Janssen"] = 3
new.data$VACCINE_DOSES[!is.na(new.data$Days.between.collection.date.and.booster.date) & new.data$vaccine_type == "Janssen"] = 2

# Look Days Post Vaccination and Days Post Boost. I suspect controls who were subsequently vaccinated were causing the problem.
new.data$Days.between.collection.date.and.vaccination.complete.date = floor(new.data$Days.between.collection.date.and.vaccination.complete.date)

# Correct input data based on any errors from the below function
in.data = data.setup.v2(new.data, start.year)
nrow(in.data)
View(in.data) 

table(in.data$sequence_result)

#Classify strains as omicron or not
# #original.delta = grep("B.1.617.2", in.data$seqeunce_result) 
# #derivative.delta = grep("AY", in.data$seqeunce_result)
#all.focal = sort(c(original.delta, derivative.delta))
original.omicron = grep("B.1.1.529", in.data$sequence_result)
derivative.omicron = grep("BA", in.data$sequence_result)
all.focal = c(original.omicron, derivative.omicron)

#### ANALYSIS SETUP ####
in.seed = 20220106 # Ensure repeatability of results

##### Examine descriptive information, basic QC ######
# Create an overview table about the entire data set
# Check age
summary(in.data$Age)
age.out = matrix(summary(in.data$Age))
row.names(age.out) = names(summary(in.data$Age))
colnames(age.out) = c("Age")
age.file = sprintf("%s/AgeTable%s.csv", out.folder, figure.label)
cat("Statistic,", file = age.file) # moves header over to account for row names
suppressWarnings(write.table(age.out, sep = ',',
            file = age.file,
            col.names = TRUE, append = TRUE)) # Gives a warning about appending to an existing table

# Make age histogram
hist(in.data$Age)
tiff(filename = sprintf("%s/AgeHistogram%s.tif", out.folder, figure.label), res = c(300),
     compression = c('lzw'), height = 1200, width = 1200)
hist(in.data$Age)
dev.off()

# Check sex
table(in.data$Sex)
mf.file = sprintf("%s/SexTable%s.csv", out.folder, figure.label)
write.table(table(in.data$Sex), sep = ',',
                             file = mf.file,
                             col.names = TRUE, row.names = FALSE, append = FALSE)

# Check region
table(in.data$Economic_region)
region.file = sprintf("%s/RegionTable%s.csv", out.folder, figure.label)
write.table(table(in.data$Economic_region), sep = ',', file = region.file,
                  col.names = TRUE, row.names = FALSE, append = FALSE)

# Check case/control status and consistency of labels
table(in.data$case_control)

# Check sequence information
seq.info = table(in.data$sequence_result)
seq.info = sort(seq.info, decreasing = TRUE) # put strains with most representation at the top
seq.info
seq.file = sprintf("%s/SequenceTable%s.csv", out.folder, figure.label)
write.table(seq.info, sep = ',', file = seq.file, col.names = TRUE, row.names = FALSE,
            append = FALSE)

# Classify strains
in.data$IS_FOCAL = 0
in.data$IS_FOCAL[all.focal] = 1
table(in.data$IS_FOCAL)
om.file = sprintf("%s/FocalTable%s.csv", out.folder, figure.label)
write.table(table(in.data$IS_FOCAL), sep = ',', file = om.file, col.names = TRUE,
            row.names = FALSE, append = FALSE)

# Check number of doses
table(in.data$VACCINE_DOSES)
dose.file = sprintf("%s/VaccineDosesTable%s.csv", out.folder, figure.label)
write.table(table(in.data$VACCINE_DOSES), sep = ',', file = dose.file, col.names = TRUE, row.names = FALSE,
            append = FALSE)

# Check vaccine type information
table(in.data$vaccine_type)
vax.type.file = sprintf("%s/VaccineTypeTable%s.csv", out.folder, figure.label)
write.table(table(in.data$vaccine_type), sep = ',', file = vax.type.file, col.names = TRUE, row.names = FALSE,
            append = FALSE)

# Check booster type information
table(in.data$booster_type)
boost.type.file = sprintf("%s/BoosterTypeTable%s.csv", out.folder, figure.label)
write.table(table(in.data$booster_type), sep = ',', file = boost.type.file, col.names = TRUE, row.names = FALSE,
            append = FALSE)

# Output the cleaned in.data object for use in other analyses
out.data.file = sprintf("Data/cleaned_data%s.csv", figure.label)
write.table(in.data, out.data.file, sep = ',', col.names = TRUE, row.names = FALSE,
            append = FALSE)

# Check for potential duplicated entries
dup.check.df = in.data[ , c(2,3,4,5,6,7,8,9,10,11)]
n.possible.duplicates = sum(duplicated(dup.check.df))
possible.dups = in.data[duplicated(dup.check.df) | duplicated(dup.check.df, fromLast = TRUE), ]
View(possible.dups)
duplicate.table = sprintf("%s/possible_duplicates%s.csv", out.folder, figure.label)
write.table(possible.dups, file = duplicate.table, sep = ',', row.names = FALSE,
            col.names = TRUE, append = FALSE)

