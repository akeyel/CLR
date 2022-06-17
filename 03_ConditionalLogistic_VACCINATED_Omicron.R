# Goal: Analyze Omicron data using a matched case-control approach
# Exclude unvaccinated individuals, in order to look at time post vaccination in a continuous manner

# Created 2021-12-15
# Updated 2022-03-23

# Author: Alexander "Sasha" Keyel <alexander.keyel@health.ny.gov

# Load required packages
library(survival)

setwd("C:/hri/DOH_COVID")
source("R/DOH_COVID_hlpr.R") # Requires dfmip package, devtools::install_github('akeyel/dfmip')
in.seed = 20220106 # Ensure repeatability of results

# Create a label for file outputs
source("R/00_Omicron_Prelim_Settings.R")

# Read in the input data file
DOH.data.file = sprintf("Data/cleaned_data%s.csv", figure.label)
if (!file.exists(DOH.data.file)){
  stop("THIS ANALYSIS REQUIRES THE CLEANED DATA SET CREATED IN
       00_Descriptive_Omicron.Rmd or 00_Descriptive_Omicron.R. If you have created
       that file and are still seeing this message, please check the file paths.")
}
in.data = read.csv(DOH.data.file)

# Check data object to ensure the correct object loaded and there are no irregularities in the data set
View(in.data)

##### DROP UNVACCINATED INDIVIDUALS PRIOR TO MATCHING #####
vax.data = in.data[in.data$VACCINE_DOSES > 0, ]
View(vax.data)

##### Match cases in flipped design & perform conditional logistic regression #####
date.offset = 6
age.bins = c(4, 11,17,29,  49, 69, 89, 110) # Assumes no one over 110

case.options = case.control.options.v4(vax.data, in.seed, date.offset, age.bins)

nrow(case.options[case.options$N_EXACT > 0, ]) # Number with at least one exact match
# Sort in.data so that those with the fewest matches are matched first.
# Sort on Male/Female mismatch first
case.options = case.options[order(case.options$N_NO_MF), ]
# Then sort on exact matches. This makes this the first sort level, secondarily sorted by number of M/F matches
case.options = case.options[order(case.options$N_EXACT), ]
case.options$SORT = seq(1,nrow(case.options))

in.data.2 = merge(vax.data, case.options[ , c("Identifier", "SORT")],
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

strata.df$IS.BOOST = 0
strata.df$IS.BOOST[strata.df$VACCINE_DOSES == 3] = 1 # If 3 doses, it's boosted (unless an upstream subsetting criterion changes)
strata.df$IS.BOOST[strata.df$VACCINE_DOSES == 2 & strata.df$vaccine_type == "Janssen"] = 1
# Could have done booster_type != 'none' as an alternate way of subsetting. Also, it does not look like there are any true Janssen boosters.

# Make demographic summary table for the flipped design
# UPDATED 2022-04-29
#CL.vax.file = sprintf("%s/02%s_stratadf%s.csv", out.folder, table.label, figure.label)
#strata.df = read.csv(CL.vax.file)
table.label = "_ConditionalLogistic_VAXONLY"
make.flipped.demographic.table(strata.df, out.folder, figure.label, focal.label,
                               table.label, incl.booster = 1)
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


# Note: assumption of normality is violated, but as long as cases are above 5 in each group, I think this is OK.
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

##### ANALYSIS 2: PERFORM CONDITIONAL LOGISTIC REGRESSION ON THE FLIPPED CASE/CONTROL DESIGN #####
# Add an indicator variable, to make 0, 'no vaccine' the reference state
strata.df$vaccine_type_code = 0
strata.df$vaccine_type_code[strata.df$vaccine_type == "Pfizer"] = 1
strata.df$vaccine_type_code[strata.df$vaccine_type == "Moderna"] = 2
strata.df$vaccine_type_code[strata.df$vaccine_type == "Janssen"] = 3

# Add a time since last dose category (coded on the indicator variable not on the days)
strata.df$TimePostLastDose = strata.df$TimePostVax
# If booster is not 0, then use the time from the booster shot.
strata.df$TimePostLastDose[strata.df$TimePostBoost > 0] = strata.df$TimePostBoost[strata.df$TimePostBoost > 0]

# Calculate days since last dose
strata.df$DaysPostLastDose = mapply(min, strata.df$Days.between.collection.date.and.vaccination.complete.date,
                                    strata.df$Days.between.collection.date.and.booster.date, na.rm = TRUE)

# Add an indicator variable for the Janssen vaccine
strata.df$JANSSEN = 0
strata.df$JANSSEN[strata.df$vaccine_type == "Janssen"] = 1

# Write matched data set for record-keeping purposes
CL.vax.file = sprintf("%s/02%s_stratadf%s.csv", out.folder, table.label, figure.label)
write.table(strata.df, file = CL.vax.file, append = FALSE, row.names = FALSE, col.names = TRUE,
            sep = ',')

# Examine individual models
vax.type.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + strata(stratum))
vax.type = survival::clogit(vax.type.formula, data=strata.df)

tpv.formula = as.formula(IS_FOCAL ~ Days.between.collection.date.and.vaccination.complete.date
                         + strata(stratum))
time.post.vax = survival::clogit(tpv.formula, data=strata.df)

has.boost.formula = as.formula(IS_FOCAL ~ IS.BOOST + strata(stratum))
has.boost = survival::clogit(has.boost.formula, data = strata.df)

# Plot time since vaccination vs. probability
warning("Plots are not currently being output. Code is mainly here as a template to adapt when conducting the analysis for real")
# https://stat.ethz.ch/R-manual/R-patched/library/survival/html/predict.coxph.html
expected.prob = predict(time.post.vax, data = strata.df, type = 'expected')
plot(strata.df$Days.between.collection.date.and.vaccination.complete.date,
     expected.prob) # Absolute measure of risk.

est.risk = predict(time.post.vax, data = strata.df, type = 'risk')
plot(strata.df$Days.between.collection.date.and.vaccination.complete.date,
     est.risk)

# Remove number of vaccine doses. This variable is difficult to interpret in its current formulation
#nd.formula = as.formula(IS_FOCAL ~ VACCINE_DOSES + strata(stratum))
#n.doses = survival::clogit(nd.formula, data=strata.df)

# THIS MODEL IS NOT COMPARABLE, BECAUSE THE SAMPLE SIZE CHANGES.
#tpb.formula = as.formula(IS_FOCAL ~ Days.between.collection.date.and.booster.date + strata(stratum))
#time.post.boost = survival::clogit(tpb.formula, data=strata.df)

tpd.formula = as.formula(IS_FOCAL ~ DaysPostLastDose + strata(stratum))
time.post.dose = survival::clogit(tpd.formula, data=strata.df)
expected.prob = predict(time.post.dose, data = strata.df, type = 'expected')
plot(strata.df$DaysPostLastDose, expected.prob) # Absolute measure of risk.

# Examine basic 2-way interactions for 5 models = 10 combinations (will code out by hand, for more combinations, will need to do programmatically)
vaxtype.x.timepostvax.formula = IS_FOCAL ~ as.factor(vaccine_type_code)*DaysPostLastDose + strata(stratum)
vaxtype.x.timepostvax = survival::clogit(vaxtype.x.timepostvax.formula, data = strata.df)
# Interactions are all coming out as NA in estimation. I think this is because both code for unvaccinated?

# Dropping the VACCINE_DOSES variable from the analysis
# This interaction works, only one parameter is unestimated due to duplication
#vaxtype.x.ndoses.formula = IS_FOCAL ~ as.factor(vaccine_type_code)*VACCINE_DOSES + strata(stratum)
#vaxtype.x.ndoses = survival::clogit(vaxtype.x.ndoses.formula, data = strata.df)

# Look at interaction between time post dose and booster
timepostdose.x.boost.formula = IS_FOCAL ~ as.factor(IS.BOOST)*DaysPostLastDose + strata(stratum)
timepostdose.x.boost = survival::clogit(timepostdose.x.boost.formula, data = strata.df)

# Try additive model rather than interaction
timepostdose.plus.boost.formula = IS_FOCAL ~ as.factor(IS.BOOST) + DaysPostLastDose + strata(stratum)
timepostdose.plus.boost = survival::clogit(timepostdose.plus.boost.formula, data = strata.df)

# Look at full model with no interactions (time post boost not included - NA values for unboosted individuals) 
full.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + 
                            Days.between.collection.date.and.vaccination.complete.date +
                            DaysPostLastDose +
                            strata(stratum)) # Dropped: VACCINE_DOSES + 
full.model = survival::clogit(full.formula, data=strata.df)

# Custom model - examine vaccine type + time since last dose, without interaction
vax.type.plus.time.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + DaysPostLastDose + 
                                          strata(stratum))
vax.type.plus.time = survival::clogit(vax.type.plus.time.formula, data = strata.df)

# Custom model - time since last dose plus number of doses, without an interaction
#tpd.and.ndoses.formula = as.formula(IS_FOCAL ~ DaysPostLastDose + VACCINE_DOSES + strata(stratum))
#time.post.dose.and.ndoses = survival::clogit(tpd.and.ndoses.formula, data=strata.df)

# Add a model to simulate the number of doses model that had been the top model, but with clearer interpretability
vax.status.plus.janssen.formula = as.formula(IS_FOCAL ~ IS.BOOST + JANSSEN + strata(stratum))
vax.status.plus.janssen = survival::clogit(vax.status.plus.janssen.formula, data = strata.df)

# Add a model to see if we can improve the top model
vax.status.plus.vaccine.type.formula = as.formula(IS_FOCAL ~ IS.BOOST + as.factor(vaccine_type_code) + strata(stratum))
vax.status.plus.vaccine.type = survival::clogit(vax.status.plus.vaccine.type.formula, data = strata.df)

# Add a model to see if we can improve the top model
timepostdose.plus.Janssen.formula = IS_FOCAL ~ JANSSEN + DaysPostLastDose + strata(stratum)
timepostdose.plus.Janssen = survival::clogit(timepostdose.plus.Janssen.formula, data = strata.df)

# Add a Janssen-only model for contrast
Janssen.formula = IS_FOCAL ~ JANSSEN + strata(stratum)
Janssen = survival::clogit(Janssen.formula, data = strata.df)

# Calculate AIC scores and identify the model with the lowest AIC
models = c("VaccineType", "TimePostVaccination", 
           "TimePostDose", "FullModel",
           "VaccineTypeXTimePostDose",
           "VaccineType+TimePostDose",
           "IsBoosted", "TimePostDoseXIsBoosted",
           "TimePostDosePlusBoosted", "Vax Status + Janssen",
           "Vax Status + VaccineType",
           "TimePostDose + Janssen", "Janssen") # "TimePostBooster",
# "NumberVaccineDoses","VaccineTypeXNumberVaccineDoses",           "TimePostDose+NumberVaccineDoses",
model.objects = list(vax.type, time.post.vax,
                     time.post.dose, full.model, 
                     vaxtype.x.timepostvax, vax.type.plus.time,
                     has.boost, timepostdose.x.boost,
                     timepostdose.plus.boost, vax.status.plus.janssen,
                     vax.status.plus.vaccine.type,
                     timepostdose.plus.Janssen, Janssen) # time.post.boost,
#  n.doses, vaxtype.x.ndoses, time.post.dose.and.ndoses,
formula.objects = list(vax.type.formula, tpv.formula, 
                       tpd.formula, full.formula,
                       vaxtype.x.timepostvax.formula,
                       vax.type.plus.time.formula,
                       has.boost.formula, timepostdose.x.boost.formula,
                       timepostdose.plus.boost.formula,
                       vax.status.plus.janssen.formula,
                       vax.status.plus.vaccine.type.formula,
                       timepostdose.plus.Janssen.formula,
                       Janssen.formula) # tpb.formula,
# nd.formula, vaxtype.x.ndoses.formula, tpd.and.ndoses.formula,
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

##### THIS SECTION REQUIRES UPDATING #####
stop("Please examine the model results, and add additional models as desired based on the
     results")

# Takes the model with the lowest delta AIC.
# This ignores ties - so it may be desirable to change this manually
obj.index = aic.results$SORT.ORDER[1] 
model.object = model.objects[[obj.index]] 
formula.object = formula.objects[[obj.index]]
make.leverage.plot(model.object, formula.object, out.folder, figure.label, table.label)


## SUPPLEMENTAL: REMOVE BOOSTED INDIVIDUALS FROM ANALYSIS, AND CHECK TIME POST LAST DOSE EFFECT
no.boost = strata.df[strata.df$IS.BOOST == 0, ]
strata.info = table(no.boost$stratum)
partial.strata = strata.info[strata.info < 2]
no.boost.2 = no.boost[!no.boost$stratum %in% names(partial.strata), ] # 152 records

tpd.formula.no.boost = as.formula(IS_FOCAL ~ DaysPostLastDose + strata(stratum))
time.post.dose.no.boost = survival::clogit(tpd.formula.no.boost, data=no.boost.2)
tpd.AIC = extractAIC(time.post.dose.no.boost)[2]

tpd.plus.Janssen.formula.no.boost = as.formula(IS_FOCAL ~ DaysPostLastDose + JANSSEN + strata(stratum))
time.post.dose.plus.Janssen.no.boost = survival::clogit(tpd.plus.Janssen.formula.no.boost, data=no.boost.2)

tpdJ.AIC = extractAIC(time.post.dose.plus.Janssen.no.boost)[2]

Janssen.formula.no.boost = as.formula(IS_FOCAL ~ JANSSEN + strata(stratum))
Janssen.no.boost = survival::clogit(Janssen.formula.no.boost, data=no.boost.2)

J.AIC = extractAIC(Janssen.no.boost)[2]

plot(strata.df$DaysPostLastDose, strata.df$IS.BOOST, col = (strata.df$IS_FOCAL + 1))

plot(strata.df$DaysPostLastDose, strata.df$IS_FOCAL, col = strata.df$IS.BOOST + 1)

## LOOK AT ONLY BOOSTED INDIVIDUALS, AND CHECK TIME POST LAST DOSE EFFECT
boost = strata.df[strata.df$IS.BOOST == 1, ]
strata.info.2 = table(boost$stratum)
# Drop partial strata where half of the match is boosted and half is not.
partial.strata.2 = strata.info.2[strata.info.2 < 2] 
boost.2 = boost[!boost$stratum %in% names(partial.strata.2), ] # 12 samples meet the criteria (6 pairs)

tpd.formula.boost = as.formula(IS_FOCAL ~ DaysPostLastDose + strata(stratum))
time.post.dose.boost = survival::clogit(tpd.formula.no.boost, data=boost.2)
time.post.dose.boost
# Not significant