# Goal: Analyze Delta data using a case-control approach

# Created 2021-12-15
# Updated 2022-03-28

# Author: Alexander "Sasha" Keyel <alexander.keyel@health.ny.gov

# Load required packages
library(survival)

setwd("C:/hri/DOH_COVID")
source("R/DOH_COVID_hlpr.R") # Requires dfmip package, devtools::install_github('akeyel/dfmip')
in.seed = 20211130 #20220112 # Ensure repeatability of results

source("R/00_Delta_Prelim_Settings.R")

# Read in the input data file
DOH.data.file = sprintf("Data/cleaned_data%s.csv", figure.label)
if (!file.exists(DOH.data.file)){
  stop("THIS ANALYSIS REQUIRES THE CLEANED DATA SET CREATED IN
       05_Descriptive_Delta.Rmd or 05_Descriptive_Delta.R. If you have created
       that file and are still seeing this message, please check the file paths.")
}
in.data = read.csv(DOH.data.file)

# Check data object to ensure the correct object loaded and there are no irregularities in the data set
View(in.data)

##### Match cases in flipped design & perform conditional logistic regression #####
date.offset = 6
age.bins = c(4, 11,17,29,  49, 69, 89, 110) # Assumes no one over 110

case.options = case.control.options.v4(in.data, in.seed, date.offset, age.bins)

nrow(case.options[case.options$N_EXACT > 0, ]) # Number with at least one exact match
# Sort in.data so that those with the fewest matches are matched first.
# Sort on Male/Female mismatch first
case.options = case.options[order(case.options$N_NO_MF), ]
# Then sort on exact matches. This makes this the first sort level, secondarily sorted by number of M/F matches
case.options = case.options[order(case.options$N_EXACT), ]
case.options$SORT = seq(1,nrow(case.options))

in.data.2 = merge(in.data, case.options[ , c("Identifier", "SORT")],
                  by = "Identifier", all.x = TRUE)
in.data.2 = in.data.2[order(in.data.2$SORT), ]

strata.df = match.focal.v2(in.data.2, in.seed, date.offset, age.bins,
                        removal = TRUE)

# Add a time-since-vaccination binary variable
strata.df$TimePostVax = 0
strata.df$TimePostVax[strata.df$Days.between.collection.date.and.vaccination.complete.date < 90] = 1 # Does not affect NA values, so those remain zero
strata.df$TimePostVax[strata.df$Days.between.collection.date.and.vaccination.complete.date >= 90] = 2

strata.df$TimePostBoost = 0
strata.df$TimePostBoost[strata.df$Days.between.collection.date.and.booster.date < 90] = 1 # Does not affect NA values, so those remain zero
strata.df$TimePostBoost[strata.df$Days.between.collection.date.and.booster.date >= 90] = 2

# Make demographic summary table for the flipped design
table.label = "_ConditionalLogistic"
# UPDATED 2022-04-29 - changed make.flipped.demographic.table, and selectively ran code to just update that piece.
#CL.file = sprintf("%s/01%s_stratadf%s.csv", out.folder, table.label, figure.label)
#strata.df = read.csv(CL.file)
make.flipped.demographic.table(strata.df, out.folder, figure.label, focal.label, table.label)
# Note: There will be a "'" before 10-19; this is to keep Excel from putting it into Date format, and that will need to be manually removed.

# Test for significant differences in age (paired) and sex (unpaired)
# Sort data to ensure a paired test
focal.cases = strata.df[strata.df$IS_FOCAL == 1, ]
other.cases = strata.df[strata.df$IS_FOCAL == 0, ]

focal.cases = focal.cases[order(focal.cases$stratum), ]
other.cases = other.cases[order(other.cases$stratum), ]
if (nrow(focal.cases) != nrow(other.cases)){stop("this assumes a 1:1 matched design. Need to change the t-test if not doing a 1:1 match")}
age.t.test = t.test(focal.cases$Age, other.cases$Age, paired = TRUE, alternative = 'two.sided')
age.t.test # Display results to screen
age.t.test.file = sprintf("%s/age_t_test_results%s%s.csv", out.folder, table.label, figure.label)
age.t.test.results = data.frame(t_statistic = age.t.test$statistic,
                                df = age.t.test$parameter, p_value = age.t.test$p.value,
                                mean_difference = age.t.test$estimate)

if (age.t.test$p.value < 0.05){
  warning(sprintf("There is a significant difference in ages of %s years between
                  the case and controls (p = %s)", age.t.test$estimate, age.t.test$p.value))
}

write.table(age.t.test.results, file = age.t.test.file, sep = ',',
            row.names = FALSE, col.names = TRUE,
            append = FALSE)

# Note: assumption of normality is violated, but as long as cases are above 5 in each group, this should be fine
mf.t.test = t.test(as.numeric(as.factor(focal.cases$Sex)),as.numeric(as.factor(other.cases$Sex)),
                   alternative = 'two.sided', paired = FALSE)
mf.t.test

mf.t.test.file = sprintf("%s/sex_t_test_results%s%s.csv", out.folder, table.label, figure.label)
mf.t.test.results = data.frame(t_statistic = mf.t.test$statistic,
                                df = mf.t.test$parameter, p_value = mf.t.test$p.value,
                                group_mean_1 = mf.t.test$estimate[1], group_mean_2 = mf.t.test$estimate[2])

if (mf.t.test$p.value < 0.05){
  warning(sprintf("There is a significant difference in distribution of sexes between
                  the case and controls (p = %s)", mf.t.test$p.value))
}

write.table(mf.t.test.results, file = mf.t.test.file, sep = ',',
            row.names = FALSE, col.names = TRUE,
            append = FALSE)

##### PERFORM CONDITIONAL LOGISTIC REGRESSION ON THE FLIPPED CASE/CONTROL DESIGN #####
# Add an indicator variable, to make 0, 'no vaccine' the reference state
strata.df$vaccine_type_code = 0
strata.df$vaccine_type_code[strata.df$vaccine_type == "Pfizer"] = 1
strata.df$vaccine_type_code[strata.df$vaccine_type == "Moderna"] = 2
strata.df$vaccine_type_code[strata.df$vaccine_type == "Janssen"] = 3

# Add a time since last dose category (coded on the indicator variable not on the days)
strata.df$TimePostLastDose = strata.df$TimePostVax
# If booster is not 0, then use the time from the booster shot.
strata.df$TimePostLastDose[strata.df$TimePostBoost > 0] = strata.df$TimePostBoost[strata.df$TimePostBoost > 0]

# Add a simplified time post vaccination to reflect significance testing
strata.df$TimePostVaxSimple = 0
strata.df$TimePostVaxSimple[strata.df$Days.between.collection.date.and.vaccination.complete.date >= 90] = 1

# Add a vaccination indicator variable
strata.df$IS.VAX = as.factor(strata.df$case_control)
strata.df$IS.VAX = relevel(strata.df$IS.VAX, ref = "Control")

# Patch on 2022-03-28 - just read in existing strata.df without re-matching everything
# need to run CL.file = ... line first
#strata.df = read.csv(CL.file)

# A variable similar to vaccine type, but with the MRNA vaccines merged
strata.df$MRNA_MERGED = 0
strata.df$MRNA_MERGED[strata.df$IS.VAX == "Case"] = 1
strata.df$MRNA_MERGED[strata.df$vaccine_type == "Janssen"] = 2

# Write matched data set for record-keeping purposes
CL.file = sprintf("%s/01%s_stratadf%s.csv", out.folder, table.label, figure.label)
write.table(strata.df, file = CL.file, append = FALSE, row.names = FALSE, col.names = TRUE,
            sep = ',')

# Examine individual models
vax.type.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + strata(stratum))
vax.type = survival::clogit(vax.type.formula, data=strata.df)

tpv.formula = as.formula(IS_FOCAL ~ as.factor(TimePostVax) + strata(stratum))
time.post.vax = survival::clogit(tpv.formula, data=strata.df)

tpv2.formula = as.formula(IS_FOCAL ~ as.factor(TimePostVaxSimple) + strata(stratum))
time.post.vax2 = survival::clogit(tpv2.formula, data=strata.df)

# Examine basic 2-way interactions 

# This interaction works, only one parameter is unestimated due to duplication
vaxtype.x.timepostvax.formula = IS_FOCAL ~ as.factor(vaccine_type_code)*as.factor(TimePostVax) + strata(stratum)
vaxtype.x.timepostvax = survival::clogit(vaxtype.x.timepostvax.formula, data = strata.df)

# Look at full model with no interactions
full.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + as.factor(TimePostVax) +
                            strata(stratum))
full.model = survival::clogit(full.formula, data=strata.df)

has.vax.formula = as.formula(IS_FOCAL ~ IS.VAX + strata(stratum)) # Here, case_control is not related to the case and control in IS_FOCAL, but is an indicator for vaccination status.
has.vax = survival::clogit(has.vax.formula, data = strata.df)

# Add number of doses - this was a top model from the age analysis
# Removed - not cleanly interpretable. Replaced with a merged variable for the mRNA vaccines
#ndoses.formula = as.formula(IS_FOCAL ~ VACCINE_DOSES + strata(stratum))
#ndoses = survival::clogit(ndoses.formula, data=strata.df)

MRNA.formula = as.formula(IS_FOCAL ~ as.factor(MRNA_MERGED) + strata(stratum))
MRNA = survival::clogit(MRNA.formula, data=strata.df)

# Calculate AIC scores and identify the model with the lowest AIC
models = c("VaccineType", "TimePostVaccination", "FullModel",
           "VaccineTypeXTimePostVaccination", "Greater90daysPostVax",
           "HasVaccine", "VaccineType_mRNA_merged")
model.objects = list(vax.type, time.post.vax, full.model, vaxtype.x.timepostvax,
                     time.post.vax2, has.vax, MRNA) #ndoses
formula.objects = list(vax.type.formula, tpv.formula, full.formula, vaxtype.x.timepostvax.formula,
                       tpv2.formula, has.vax.formula, MRNA.formula) #ndoses.formula

aic.results = data.frame(model = models, AIC = NA, K = NA, LogLikelihood = NA,
                         coefficients = NA,
                         coefficient.Odds.Ratios = NA,
                         coefficient.Odds.Ratios.lower = NA,
                         coefficient.Odds.Ratios.upper = NA,
                         p.value = NA, coefficient.p.values = NA)
for (i in 1:length(models)){
  model = models[i]
  model.object = model.objects[[i]]
  model.summary = summary(model.object)
  model.AIC = extractAIC(model.object)[2]
  model.K = extractAIC(model.object)[1]
  model.loglik = model.object$loglik[2]
  aic.results$AIC[i] = model.AIC
  aic.results$K[i] = model.K
  aic.results$LogLikelihood[i] = model.loglik
  aic.results$p.value[i] = model.summary$logtest[3] # 3rd element is the p-value, logtest corresponds to the Likelihood Ratio test
  aic.results$coefficients[i] = paste(rownames(model.summary$coefficients), collapse = '; ')
  aic.results$coefficient.p.values[i] = paste(round(model.summary$coefficients[ ,5], 4), collapse = '; ')
  aic.results$coefficient.Odds.Ratios[i] = paste(round(model.summary$coefficients[ ,2],6), collapse = '; ')
  aic.results$coefficient.Odds.Ratios.lower[i] = paste(round(model.summary$conf.int[ ,3],3), collapse = '; ')
  aic.results$coefficient.Odds.Ratios.upper[i] = paste(round(model.summary$conf.int[ ,4],3), collapse = '; ')
  
}
aic.results$SORT.ORDER = seq(1:nrow(aic.results)) # Get an indicator for the R object associated with each model
aic.results$DeltaAIC = aic.results$AIC - min(aic.results$AIC)
aic.results = aic.results[order(aic.results$DeltaAIC), ]

aic.file = sprintf("%s/AIC_RESULTS%s%s.csv", out.folder, table.label, figure.label)
write.table(aic.results, aic.file, sep = ',', row.names = FALSE, 
            col.names = TRUE)

stop("Please examine the model results, and add additional models as desired based on the
     results")

# Takes the model with the lowest delta AIC.
# This ignores ties - so it may be desirable to change this manually
obj.index = aic.results$SORT.ORDER[1] 
model.object = model.objects[[obj.index]] 
formula.object = formula.objects[[obj.index]]
make.leverage.plot(model.object, formula.object, out.folder, figure.label, table.label)

