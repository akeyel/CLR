# Goal: Generate descriptive information for the Delta vs. other strains analysis

# Created 2021-12-15

# Author: Alexander "Sasha" Keyel <alexander.keyel@health.ny.gov

stop("This script may no longer be needed now that Delta and Omicron data sets are merged into one file")

#### THIS SECTION NEEDS UPDATING ####

# NOTES:
# Set working directory
setwd("C:/hri/DOH_COVID")
# Load helper functions
source("R/DOH_COVID_hlpr.R") # Requires dfmip package, devtools::install_github('akeyel/dfmip')
in.seed = 20211130 #20220112 # Ensure repeatability of results

# Create a label for file outputs
focal.label = "Delta" # Indicator for what strain is being analyzed. Used in text labels.
out.folder = "Results/Delta"
figure.label = "_delta"
sim.data = 0 # Set to 0 when running for real
start.year = 2021 # Earliest year of a collection date in the data set
DOH.data.file = "Data/Final_dataset_2021_12_07.csv"

# Read in the input data file
new.data = read.csv(DOH.data.file)

# Adjust data set based on any errors from the data.setup.v2 function below
# Correct data set irregularities
new.data$sequence_result = new.data$seqeunce_result
new.data$seqeunce_result = NULL # Get rid of typo
# Add columns for consistency purposes. These are NA for the delta data set.
new.data$booster_type = NA
new.data$Days.between.collection.date.and.booster.date = NA

if (sim.data == 1){
  new.data = fix.sim.data.1(new.data)
}

# Correct input data based on any errors from the below function
in.data = data.setup.v2(new.data, start.year)
nrow(in.data)
View(in.data) 

# Check for duplicated entries
dup.check.df = in.data[ , c(2,3,4,5,6,7,8,9,10,11)]
n.possible.duplicates = sum(duplicated(dup.check.df))
possible.dups = in.data[duplicated(dup.check.df) | duplicated(dup.check.df, fromLast = TRUE), ]
# Note switch of in.data object for dup.check.df. They should be in the same order
# But this way the identifiers are in the final file.
View(possible.dups)
duplicate.table = sprintf("%s/possible_duplicates%s.csv", out.folder, figure.label)
write.table(possible.dups, file = duplicate.table, sep = ',', row.names = FALSE,
            col.names = TRUE, append = FALSE)

warning("Please check sequence_result column manually to ensure that sequences are correct and in the correct format")
table(in.data$sequence_result)

# Classify into focal strain or other strains
# #For Delta Select 1.617.2 and everything that has an AY (NOTE: if AY is anywhere, it will grab it!)
original.delta = grep("B.1.617.2", in.data$sequence_result) 
derivative.delta = grep("AY", in.data$sequence_result)
all.focal = sort(c(original.delta, derivative.delta))


#### END OF SECTION ####

#### ANALYSIS SETUP ####
# Load packages and source helper code
library(exact2x2)
library(survival)

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

# Check sequence information
seq.info = table(in.data$sequence_result)
seq.info = sort(seq.info, decreasing = TRUE) # put strains with most representation at the top
seq.info
seq.file = sprintf("%s/SequenceTable%s.csv", out.folder, figure.label)
write.table(seq.info, sep = ',', file = seq.file, col.names = TRUE, row.names = FALSE,
            append = FALSE)

# Classify strains
#Classify strains as omicron or not (old code for delta included as a template and if a delta classification is still of interest)
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


