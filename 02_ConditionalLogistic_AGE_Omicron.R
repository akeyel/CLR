# Goal: Examine age effects on Omicron infection using a matched case-control approach.

# Created 2021-12-15
# Updated 2022-03-28

# Author: Alexander "Sasha" Keyel <alexander.keyel@health.ny.gov

# Load required packages
library(survival)

setwd("C:/hri/DOH_COVID")
source("R/DOH_COVID_hlpr.R") # Requires dfmip package, devtools::install_github('akeyel/dfmip')
in.seed = 20220106 # Ensure repeatability of results

source("R/00_Omicron_Prelim_Settings.R")
DOH.data.file = sprintf("Data/cleaned_data%s.csv", figure.label)
if (!file.exists(DOH.data.file)){
  stop("THIS ANALYSIS REQUIRES THE CLEANED DATA SET CREATED IN
       00_Descriptive_Omicron.Rmd or 00_Descriptive_Omicron.R. If you have created
       that file and are still seeing this message, please check the file paths.")
}
in.data = read.csv(DOH.data.file)
# Check data object to ensure the correct object loaded and there are no irregularities in the data set
View(in.data)


##### Match cases in flipped design & perform conditional logistic regression #####
# REMOVE AGE AS A MATCHING CRITERION
date.offset = 6
table.label = "_ConditionalLogistic_Age"

# No pre-sorting, because it includes age as a matching criterion.
warning("Data are unsorted. Pre-sorting the data may increase the number of successful matches")
strata.df = match.focal(in.data, in.seed, date.offset, age.offset,
                        removal = TRUE, age.offset.override = NA,
                        exclude.age = TRUE)

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
# Updated on 2022-04-29
#CL.vax.file = sprintf("%s/03%s_stratadf%s.csv", out.folder, table.label, figure.label)
#strata.df = read.csv(CL.vax.file)

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


# Note: assumption of normality is violated, but as long as cases are above 5 in each group, this should be OK.
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

##### ANALYSIS 3: PERFORM CONDITIONAL LOGISTIC REGRESSION ON THE FLIPPED CASE/CONTROL DESIGN #####
# Add an indicator variable, to make 0, 'no vaccine' the reference state
strata.df$vaccine_type_code = 0
strata.df$vaccine_type_code[strata.df$vaccine_type == "Pfizer"] = 1
strata.df$vaccine_type_code[strata.df$vaccine_type == "Moderna"] = 2
strata.df$vaccine_type_code[strata.df$vaccine_type == "Janssen"] = 3

# Add a time since last dose category (coded on the indicator variable not on the days)
strata.df$TimePostLastDose = strata.df$TimePostVax
# If booster is not 0, then use the time from the booster shot.
strata.df$TimePostLastDose[strata.df$TimePostBoost > 0] = strata.df$TimePostBoost[strata.df$TimePostBoost > 0]

strata.df$IS.VAX = as.factor(strata.df$case_control)
strata.df$IS.VAX = relevel(strata.df$IS.VAX, ref = "Control")

# Add a combined variable that includes three categories: 'Unvaccinated', 'Vaccinated', and 'Vaccinated and Boosted'
strata.df$VAX.STATUS = 0
strata.df$VAX.STATUS[strata.df$case_control == "Case"] = 1
strata.df$VAX.STATUS[strata.df$VACCINE_DOSES == 3] = 2 # If 3 doses, it's boosted (unless an upstream subsetting criterion changes)
strata.df$VAX.STATUS[strata.df$VACCINE_DOSES == 2 & strata.df$vaccine_type == "Janssen"] = 2
strata.df$VAX.STATUS = as.factor(strata.df$VAX)

# Add an indicator variable for the Janssen vaccine
strata.df$JANSSEN = 0
strata.df$JANSSEN[strata.df$vaccine_type == "Janssen"] = 1

# Add an indicator variable for an old booster shot (>90 days) (added based on preliminary analyses showing an effect here)
strata.df$OldBoost = 0
strata.df$OldBoost[strata.df$TimePostBoost == 2] = 1

# Load strata.df to patch in an age.bin model & re-run the lower part of the code
# strata.df = read.csv(sprintf("%s/03%s_stratadf%s.csv", out.folder, table.label, figure.label))
# Flipping the order, so that 90+ can be reference bin, and bins increase in decreasing age order
strata.df$AGE.BIN = 0
age.bins = sort(c(4, 11,17,29,  49, 69, 89, 110), decreasing = TRUE) # Assumes no one over 110
count = 0
bin.upper = 110
for (i  in 1:length(age.bins)){
  bin.upper = age.bins[i]
  if (i != length(age.bins)){
    bin.lower = age.bins[i + 1] + 1
  }else{
    bin.lower = 0
  }
  strata.df$AGE.BIN[strata.df$Age >= bin.lower & strata.df$Age <= bin.upper] = count
  
  # Update lower age limit
  #bin.lower = age + 1
  # Update count
  count = count + 1
}

# Write matched data set for record-keeping purposes
CL.vax.file = sprintf("%s/03%s_stratadf%s.csv", out.folder, table.label, figure.label)
write.table(strata.df, file = CL.vax.file, append = FALSE, row.names = FALSE, col.names = TRUE,
            sep = ',')


# Examine individual model with age
age.formula = as.formula(IS_FOCAL ~ Age + strata(stratum))
age.result = survival::clogit(age.formula, data=strata.df)

# Examine basic 2-way interactions for 5 models
vaxtype.x.age.formula = as.formula(IS_FOCAL ~ Age*as.factor(vaccine_type_code) + strata(stratum))
vaxtype.x.age = survival::clogit(vaxtype.x.age.formula, data=strata.df)

tpv.x.age.formula = as.formula(IS_FOCAL ~ Age*as.factor(TimePostVax) + strata(stratum))
timepostvax.x.age = survival::clogit(tpv.x.age.formula, data=strata.df)

#nd.x.age.formula = as.formula(IS_FOCAL ~ Age*VACCINE_DOSES + strata(stratum))
#ndoses.x.age = survival::clogit(nd.x.age.formula, data=strata.df)

tpb.x.age.formula = as.formula(IS_FOCAL ~ Age*as.factor(TimePostBoost) + strata(stratum))
timepostboost.x.age = survival::clogit(tpb.x.age.formula, data=strata.df)

tpd.x.age.formula = as.formula(IS_FOCAL ~ Age*as.factor(TimePostLastDose) + strata(stratum))
timepostdose.x.age = survival::clogit(tpd.x.age.formula, data=strata.df)

# Look at full model with no interactions
full.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + as.factor(TimePostVax) +
                            VACCINE_DOSES + as.factor(TimePostBoost) + Age +
                            strata(stratum))
full.model = survival::clogit(full.formula, data=strata.df)

# Custom
# Add age + vaccine doses WITHOUT interaction (which was non-sig.)
#nd.and.age.formula = as.formula(IS_FOCAL ~ Age + VACCINE_DOSES + strata(stratum))
#ndoses.and.age = survival::clogit(nd.and.age.formula, data=strata.df)

boost.plus.vax.formula = as.formula(IS_FOCAL ~ as.factor(VAX.STATUS) + Age + strata(stratum))
boost.plus.vax = survival::clogit(boost.plus.vax.formula, data = strata.df)

# Add best model without age, to test for an age effect
boost.plus.vax.no.age.formula = as.formula(IS_FOCAL ~ as.factor(VAX.STATUS) + strata(stratum))
boost.plus.vax.no.age = survival::clogit(boost.plus.vax.no.age.formula, data = strata.df)

# Add a model to simulate the number of doses model that had been the top model, but with clearer interpretability
vax.status.plus.janssen.formula = as.formula(IS_FOCAL ~ as.factor(VAX.STATUS) + JANSSEN + strata(stratum))
vax.status.plus.janssen = survival::clogit(vax.status.plus.janssen.formula, data = strata.df)

vax.status.plus.janssen.plus.age.formula = as.formula(IS_FOCAL ~ as.factor(VAX.STATUS) + JANSSEN + Age + strata(stratum))
vax.status.plus.janssen.plus.age = survival::clogit(vax.status.plus.janssen.plus.age.formula, data = strata.df)

# Add a model with an indicator for old booster shots - data inspection suggests this is an important variable.
#vax.status.plus.oldboost.plus.age.formula = as.formula(IS_FOCAL ~ VAX.STATUS + OldBoost + Age + strata(stratum))
#vax.status.plus.oldboost.plus.age = survival::clogit(vax.status.plus.oldboost.plus.age.formula, data = strata.df)
# Excluding this model from final presentation - not appreciably better than the next best model, and it is more difficult to explain.
# A bit surprised - from the demographic table, I expected this to be quite significant in the final results.

# Add a model to examine age by bin, rather than as a linear effect
age.bin.formula = as.formula(IS_FOCAL ~ as.factor(AGE.BIN) + strata(stratum))
bin.model = survival::clogit(age.bin.formula, data = strata.df)

age.bin.status.formula = as.formula(IS_FOCAL ~ as.factor(AGE.BIN) + as.factor(VAX.STATUS) + strata(stratum))
bin.status.model = survival::clogit(age.bin.status.formula, data = strata.df)

# Looking at plot, would a log-transform of age perform better? With an indicator for Age Bin 3, which does not fit the line?
# need to do the small offset to keep the age 0's from crashing the model
log.age.formula = as.formula(IS_FOCAL ~ log(Age + 0.001) + as.factor(VAX.STATUS) + strata(stratum))
log.age.model = survival::clogit(log.age.formula, data = strata.df)
log.age.AIC = extractAIC(log.age.model)[2]
bva.AIC = extractAIC(boost.plus.vax)[2]

# Try a 1/age model
strata.df$INVERSE_AGE = 1/(strata.df$Age + 0.001) 
inverse.age.formula = as.formula(IS_FOCAL ~ INVERSE_AGE + as.factor(VAX.STATUS) + strata(stratum))
inverse.age.model = survival::clogit(inverse.age.formula, data = strata.df)
inv.age.AIC = extractAIC(inverse.age.model)[2]

# Try indicator for 0-4 and 18 - 29, with age as an otherwise continuous variable
strata.df$AGE0_4 = 0
strata.df$AGE0_4[strata.df$Age < 5] = 1
strata.df$AGE18_29 = 0
strata.df$AGE18_29[strata.df$Age >= 18 & strata.df$Age < 30] = 1
custom.age.formula = as.formula(IS_FOCAL ~ Age + as.factor(AGE0_4) + as.factor(AGE18_29) + as.factor(VAX.STATUS) + strata(stratum))
custom.age.model = survival::clogit(custom.age.formula, data = strata.df)
custom.age.AIC = extractAIC(custom.age.model)[2]


# Calculate AIC scores and identify the model with the lowest AIC
models = c("Age", "VaccineType x Age", "TimePostVaccination x Age",
           "TimePostBooster x Age", "TimePostDose x Age", "FullModel",
           "Age + Vax + Boost", "Vax + Boost",
           "Vax + Boost + Janssen", "Vax + Boost + Janssen + Age",
           "Custom Age")
# "Age + Vaccine Doses",  "NumberVaccineDoses x Age","Vax + Boost + Age + Boost Age"
model.objects = list(age.result, vaxtype.x.age, timepostvax.x.age, timepostboost.x.age,
                     timepostdose.x.age, full.model, boost.plus.vax,
                     boost.plus.vax.no.age, vax.status.plus.janssen,
                     vax.status.plus.janssen.plus.age, custom.age.model)
#ndoses.x.age, ndoses.and.age, vax.status.plus.oldboost.plus.age
formula.objects = list(age.formula, vaxtype.x.age.formula, tpv.x.age.formula, 
                       tpb.x.age.formula, tpd.x.age.formula, full.formula,
                       boost.plus.vax.formula,
                       boost.plus.vax.no.age.formula,
                       vax.status.plus.janssen.formula, vax.status.plus.janssen.plus.age.formula,
                       custom.age.formula)
# nd.x.age.formula, nd.and.age.formula, vax.status.plus.oldboost.plus.age.formula
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


####   Add an examination of what happens if under 11's are dropped from analysis ####
strata.df.12plus = strata.df[strata.df$Age > 11, ]
#589
# No age matching, so there will be partial drops - need to drop the unmatched pair, or redo the matching process.
to.keep = names(table(strata.df.12plus$stratum)[table(strata.df.12plus$stratum) > 1])
strata.df.12plus.fixed = strata.df.12plus[strata.df.12plus$stratum %in% to.keep, ] # 560 after removing partial drops.
#boost.plus.vax.formula = as.formula(IS_FOCAL ~ IS.VAX + IS.BOOST + Age + strata(stratum)) Same as line defining it above.
boost.plus.vax.12plus = survival::clogit(boost.plus.vax.formula, data = strata.df.12plus.fixed)
boost.plus.vax.12plus

model.summary = summary(boost.plus.vax.12plus)
model.summary
# Age remains HIGHLY sig. after the <12's are removed.
# In fact, the OR's are stronger for all variables 

#### Look at pattern of decrease in OR by age

# Make a plot of the binned age results from above to assess linearity of decrease in OR
bin.OR = summary(bin.status.model)$coefficients[ ,2]
bin.summary = summary(bin.status.model)
bin.OR.upper = bin.summary$conf.int[ ,4]
bin.OR.lower = bin.summary$conf.int[ ,3]

bin.OR.df = data.frame(BIN.OR = bin.OR[1:7], BIN.OR.lower = bin.OR.lower[1:7], BIN.OR.upper = bin.OR.upper[1:7], bin = seq(1,7))
bin.OR.df = bin.OR.df[order(bin.OR.df$bin, decreasing = TRUE), ]
bin.plot = sprintf("%s/BinPlot%s%s.tif", out.folder, table.label, figure.label)
text.size = 2
tiff(bin.plot, width = 2400, height = 2400, res = c(300), compression = c('lzw'))
par(mar = c(7,6,1,1))
plot(seq(1,7), bin.OR.df$BIN.OR, xaxt = 'n',
     xlab = "", ylab = "Odds Ratio relative to 90+ age group",
     cex = text.size, pch = 16, cex.axis = text.size, cex.lab = text.size,
     ylim = c(0,100))
axis(side = 1, at = seq(1,7), labels = c("0-4","5-11","12-17","18-29","30-49","50-69","70-89"),
     cex.axis = text.size, cex.lab = text.size, las = 2)
mtext("Age Bin", side = 1, line = 5.5, cex = text.size)

# Add confidence intervals to show that they primarily overlap 0
arrows(seq(1,7), bin.OR.df$BIN.OR, seq(1,7), bin.OR.df$BIN.OR.upper,
       angle = 90)
arrows(seq(1,7), bin.OR.df$BIN.OR, seq(1,7), bin.OR.df$BIN.OR.lower,
       angle = 90)

segments(0,1,8,1, lwd = 6, lty = 2)

dev.off()


##### Predict specific effects for the terms without stratum-specific effects
terms.prob.2 = predict(custom.age.model, data = strata.df, type = 'terms')
plot(strata.df$Age, terms.prob.2[ ,1] + terms.prob.2[ ,2] + terms.prob.2[ ,3])

# Not getting error information out of 'predict'. Add +/- 1 SE
cam.summary = summary(custom.age.model)

# Doing +/- 1 SE instead of 95% CI
# For 0/1 factors, set 0.5 as the reference level. This way half the error is attributed to each series instead of all to one or the other.

continuous_variable_approach = 'standard'
continuous_variable_approach = 'per_one_unit'

# Check that this factor is consistent with another variable
#factor.2 = (log(cam.summary$conf.int[3,4]) - cam.summary$coefficients[3,1]) / cam.summary$coefficients[3,3]
factor = 1 # For SE; saves re-coding the variables
ref = 0.5 # Change reference to 0.5 instead of 0 or 1. Then we can add the error to both factor levels.
# Only for factors)
# The survival package does not have an option to set the reference level to 0.5.
# This does not appear to be standard statistical practice
# However, assigning all error to only one series makes no conceptual sense to me.

# Formulation based on error being associated with a single-unit change
Age_LogOdds = strata.df$Age *cam.summary$coefficients[1,1]

if (continuous_variable_approach == "per_one_unit"){
  Age_Lower = strata.df$Age *cam.summary$coefficients[1,1] - cam.summary$coefficients[1,3]*factor
  Age_Upper = strata.df$Age *cam.summary$coefficients[1,1] + cam.summary$coefficients[1,3]*factor
}

if (continuous_variable_approach == "standard"){
  # Formulation based on statistical package output
  terms.prob.3 = predict(custom.age.model, data = strata.df, type = 'terms', se.fit = TRUE)
  age.se = terms.prob.3[[2]][ ,1] # Second list item is the standard error, first column corresponds to age
  Age_Lower = Age_LogOdds - age.se*factor
  Age_Upper = Age_LogOdds + age.se*factor
}

A0_4_LogOdds = strata.df$AGE0_4 *cam.summary$coefficients[2,1]
A0_4_Lower = A0_4_LogOdds - cam.summary$coefficients[2,3]*factor*ref
A0_4_Upper = A0_4_LogOdds + cam.summary$coefficients[2,3]*factor*ref

A_18_LogOdds = strata.df$AGE18_29 *cam.summary$coefficients[3,1]
A_18_Lower = A_18_LogOdds - cam.summary$coefficients[3,3]*factor*ref
A_18_Upper = A_18_LogOdds + cam.summary$coefficients[3,3]*factor*ref

# Need to add joint SE, since all will be either 0 or 1 for both variables - so error applies to both.
Vax_Boost_LogOdds = rep(0, nrow(strata.df))
Vax_Boost_LogOdds[strata.df$IS.VAX == "Case" & strata.df$IS.BOOST == 0] = cam.summary$coefficients[4,1]
Vax_Boost_LogOdds[strata.df$IS.BOOST == 1] = cam.summary$coefficients[5,1]

Vax_Boost_LogOdds_Lower = Vax_Boost_LogOdds - cam.summary$coefficients[4,3]*factor*ref - cam.summary$coefficients[5,3]*factor*ref
Vax_Boost_LogOdds_Upper = Vax_Boost_LogOdds + cam.summary$coefficients[4,3]*factor*ref + cam.summary$coefficients[5,3]*factor*ref

Offset = terms.prob.2[1,1] - Age_LogOdds[1]

# Check that individual pieces were correctly calculated
strata.df$Predict_Age_LogOdds = terms.prob.2[ , 1]
strata.df$Age_LogOdds = Age_LogOdds + Offset

mini.df = strata.df[ ,c("Age", "Predict_Age_LogOdds", "Age_LogOdds")]
mini.df$Diff = round(mini.df$Predict_Age_LogOdds - mini.df$Age_LogOdds, 4)
# This shows them being the same

strata.df$LogOdds_custom_age = Age_LogOdds + A0_4_LogOdds + A_18_LogOdds + Vax_Boost_LogOdds + Offset
strata.df$LogOdds.Lower_custom_age = Age_Lower + A0_4_Lower + A_18_Lower + Vax_Boost_LogOdds_Lower + Offset
strata.df$LogOdds.Upper_custom_age = Age_Upper + A0_4_Upper + A_18_Upper + Vax_Boost_LogOdds_Upper + Offset

# Check LogOdds against sum of terms - these should be the same
#colnames(terms.prob.2)
#[1] "Age"                   "as.factor(AGE0_4)"     "as.factor(AGE18_29)"   "as.factor(VAX.STATUS)"
plot(strata.df$LogOdds_custom_age, terms.prob.2[ ,1] + terms.prob.2[ ,2] + terms.prob.2[ ,3] + terms.prob.2[,4])

# Check age
plot(Age_LogOdds + Offset, terms.prob.2[ ,1]) # It's 1:1, but it's offset by 1.5.
plot(A0_4_LogOdds, terms.prob.2[,2]) # This checks out
plot(A_18_LogOdds, terms.prob.2[,3]) # This checks out
plot(Vax_Boost_LogOdds, terms.prob.2[ ,4])

# Export visualization of odds of getting Omicron for the age model plus vaccination status
age.vax.file = sprintf("%s/OddsPlots%s%s_%s.tif", out.folder, table.label,
                       figure.label, continuous_variable_approach)
tiff(age.vax.file, res = 300, compression = c('lzw'), height = 2400, width = 3600)
par(mfrow = c(1,2))
par(mar = c(5,5,0,0))
# Full model without stratum effects
plot.text = 2
x.limits = c(0,95)
y.limits = c(-4,6.8)

plot(strata.df$Age, strata.df$LogOdds_custom_age,
     xlab = "Age", ylab = "ln(Odds of Omicron)",
     col = strata.df$VAX.STATUS + 1, pch = strata.df$VAX.STATUS,
     cex.lab = plot.text, cex.axis = plot.text, cex = plot.text,
     xlim = x.limits, ylim = y.limits)
legend(par('usr')[1], par('usr')[3] + 2.4, legend = c("Unvaccinated", "Vaccinated", "Boosted"),
       col = seq(1,3), pch = seq(15,17), cex = plot.text*0.95)
segments(0, 0, 100, 0, lty = 2, lwd = 4)

# Add for unvaccinated
unvax.df = strata.df[strata.df$IS.VAX != "Case", ]
unvax.df = unvax.df[order(unvax.df$Age), ]
lines(unvax.df$Age, unvax.df$LogOdds.Lower_custom_age,
      xlim = x.limits, ylim = y.limits, lty = 2, lwd = 2)
lines(unvax.df$Age, unvax.df$LogOdds.Upper_custom_age,
      xlim = x.limits, ylim = y.limits, lty = 2, lwd = 2)

# Add for Vaccinated
vax.df = strata.df[strata.df$IS.VAX == "Case" & strata.df$IS.BOOST == 0, ]
vax.df = vax.df[order(vax.df$Age), ]
lines(vax.df$Age, vax.df$LogOdds.Lower_custom_age,
      xlim = x.limits, ylim = y.limits, col = 2, lty = 2, lwd = 2)
lines(vax.df$Age, vax.df$LogOdds.Upper_custom_age,
      xlim = x.limits, ylim = y.limits, col = 2, lty = 2, lwd = 2)

boo.df = strata.df[strata.df$IS.BOOST == 1, ]
boo.df = boo.df[order(boo.df$Age), ]
lines(boo.df$Age, boo.df$LogOdds.Lower_custom_age,
      xlim = x.limits, ylim = y.limits, col = 3, lty = 2, lwd = 2)
lines(boo.df$Age, boo.df$LogOdds.Upper_custom_age,
      xlim = x.limits, ylim = y.limits, col = 3, lty = 2, lwd = 2)

# Try replacing lines with polygons to see if it's better


polygon(c(unvax.df$Age, rev(unvax.df$Age)), c(unvax.df$LogOdds.Lower_custom_age, rev(unvax.df$LogOdds.Upper_custom_age)),
        col = 'black', density = 8, angle = 90)

polygon(c(vax.df$Age, rev(vax.df$Age)), c(vax.df$LogOdds.Lower_custom_age, rev(vax.df$LogOdds.Upper_custom_age)),
        col = 'red', density = 8, angle = 0)


polygon(c(boo.df$Age, rev(boo.df$Age)), c(boo.df$LogOdds.Lower_custom_age, rev(boo.df$LogOdds.Upper_custom_age)),
        col = 'green', density = 8, angle = 45)


# Re-add points, so that points are in front of lines
par(new = TRUE)
plot(strata.df$Age, strata.df$LogOdds_custom_age,
     xlab = "Age", ylab = "ln(Odds of Omicron)",
     col = strata.df$VAX.STATUS + 1, pch = strata.df$VAX.STATUS + 15,
     cex.lab = plot.text, cex.axis = plot.text, cex = plot.text,
     xlim = x.limits, ylim = y.limits)




plot(strata.df$Age, exp(terms.prob.2[ ,1] + terms.prob.2[ ,2] + terms.prob.2[ ,3] + terms.prob.2[,4]),
     xlab = "Age", ylab = "Odds of Having Omicron",
     col = strata.df$VAX.STATUS + 1, pch = strata.df$VAX.STATUS,
     cex.lab = plot.text, cex.axis = plot.text, cex = plot.text)
segments(0, 1, 100, 1, lty = 2, lwd = 4)
polygon(c(-1,95,95,-1), c(0,0,25,25), cex = 2) # Go 1 outside data range to be inclusive

dev.off()

# Make an inset that can be pasted in using PowerPoint (Or illustrator)
age.vax.file = sprintf("%s/OddsPlotInset%s%s.tif", out.folder, table.label, figure.label)
tiff(age.vax.file, res = 300, compression = c('lzw'), height = 2400, width = 1300)
par(mar = c(5,5,0,0))
plot(strata.df$Age, exp(terms.prob.2[ ,1] + terms.prob.2[ ,2] + terms.prob.2[ ,3] + terms.prob.2[,4]),
     xlab = "Age", ylab = "Odds of Having Omicron", ylim = c(0,25),
     col = strata.df$VAX.STATUS + 1, pch = strata.df$VAX.STATUS,
     cex.lab = 2, cex.axis = 2, cex = 2)
segments(0, 1, 100, 1, lty = 2, lwd = 4)

dev.off()

