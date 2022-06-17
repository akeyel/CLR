# Goal: Analyze Omicron data using the case-control approach used for Delta and VOC's

# Created 2021-12-15
# Updated 2022-03-28

# Author: Alexander "Sasha" Keyel <alexander.keyel@health.ny.gov

# Load required packages
library(survival)

setwd("C:/hri/DOH_COVID")
source("R/DOH_COVID_hlpr.R") # Requires dfmip package, devtools::install_github('akeyel/dfmip')
in.seed = 20220106 # Ensure repeatability of results

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
# REMOVE AGE AS A MATCHING CRITERION
date.offset = 6
table.label = "_ConditionalLogistic_Age"

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

                        
# Make demographic summary table for the flipped design
# Updated 2022-04-29
#CL.file = sprintf("%s/08%s_stratadf%s.csv", out.folder, table.label, figure.label)
#strata.df = read.csv(CL.file)
make.flipped.demographic.table(strata.df, out.folder, figure.label, focal.label, table.label)
# Note: There will be a "'" before 10-19; this is to keep Excel from putting it into Date format, and that will need to be manually removed.
# I couldn't figure out how to get Excel to make it silently disappear the way it does when you enter it in Excel directly.

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

##### ANALYSIS 3: PERFORM CONDITIONAL LOGISTIC REGRESSION ON THE FLIPPED CASE/CONTROL DESIGN #####
# Add an indicator variable, to make 0, 'no vaccine' the reference state
strata.df$vaccine_type_code = 0
strata.df$vaccine_type_code[strata.df$vaccine_type == "Pfizer"] = 1
strata.df$vaccine_type_code[strata.df$vaccine_type == "Moderna"] = 2
strata.df$vaccine_type_code[strata.df$vaccine_type == "Janssen"] = 3

# Add a time since last dose category (coded on the indicator variable not on the days)
strata.df$TimePostLastDose = strata.df$TimePostVax # No boosters in this data set

# Add a vaccination indicator variable
strata.df$IS.VAX = as.factor(strata.df$case_control)
strata.df$IS.VAX = relevel(strata.df$IS.VAX, ref = "Control")

# Add indicator for Janssen vaccine, to try to make number of doses effect more interpretable
#strata.df$JANSSEN = 0
#strata.df$JANSSEN[strata.df$vaccine_type == "Janssen"] = 1
# This was not as interpretable - reframed below

# A variable similar to vaccine type, but with the MRNA vaccines merged
strata.df$MRNA_MERGED = 0
strata.df$MRNA_MERGED[strata.df$IS.VAX == "Case"] = 1
strata.df$MRNA_MERGED[strata.df$vaccine_type == "Janssen"] = 2

# Read in strata.df file to go back and edit analyses after the code was initially run
# strata.df = read.csv(CL.file)

# Add an age bin category. Start with oldest individuals, so the 90+ will be the reference category
strata.df$AGE.BIN = 0
#age.bins = sort(c(4, 11,17,29,  49, 69, 89, 110), decreasing = TRUE) # Assumes no one over 110
age.bins = sort(c(17,29,  49, 69, 89, 110), decreasing = TRUE) # Assumes no one over 110
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
  count = count + 1
}
nrow(strata.df[strata.df$Age < 12, ]) # Only 8 under 12, but getting model convergence
nrow(strata.df[strata.df$Age < 5, ]) # 3 under 5, also leading to model convergence
# The low sample size in the lower bins is causing problems. Going to do an 'under 18' bin.


# Write matched data set for record-keeping purposes
CL.file = sprintf("%s/08%s_stratadf%s.csv", out.folder, table.label, figure.label)
write.table(strata.df, file = CL.file, append = FALSE, row.names = FALSE, col.names = TRUE,
            sep = ',')

# Examine individual model with age
age.formula = as.formula(IS_FOCAL ~ Age + strata(stratum))
age.result = survival::clogit(age.formula, data=strata.df)

# Examine basic 2-way interactions for 5 models
vaxtype.x.age.formula = as.formula(IS_FOCAL ~ Age*as.factor(vaccine_type_code) + strata(stratum))
vaxtype.x.age = survival::clogit(vaxtype.x.age.formula, data=strata.df)

#tpv.x.age.formula = as.formula(IS_FOCAL ~ Age*as.factor(TimePostVax) + strata(stratum))
#timepostvax.x.age = survival::clogit(tpv.x.age.formula, data=strata.df)

#nd.x.age.formula = as.formula(IS_FOCAL ~ Age*VACCINE_DOSES + strata(stratum))
#ndoses.x.age = survival::clogit(nd.x.age.formula, data=strata.df)

tpd.x.age.formula = as.formula(IS_FOCAL ~ Age*as.factor(TimePostLastDose) + strata(stratum))
timepostdose.x.age = survival::clogit(tpd.x.age.formula, data=strata.df)

# Look at full model with no interactions
full.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + as.factor(TimePostLastDose) +
                            VACCINE_DOSES + Age +
                            strata(stratum))
full.model = survival::clogit(full.formula, data=strata.df)

# Removed number of doses. Not very interpretable, as it is a mix of vaccine type and vaccination status.
#nd.formula = as.formula(IS_FOCAL ~ VACCINE_DOSES + strata(stratum))
#ndoses = survival::clogit(nd.formula, data=strata.df)

#nd.plus.age.formula = as.formula(IS_FOCAL ~ Age + VACCINE_DOSES + strata(stratum))
#ndoses.plus.age = survival::clogit(nd.plus.age.formula, data=strata.df)

# Look at vaccination status or not, as a cleaner model
has.vax.formula = as.formula(IS_FOCAL ~ IS.VAX + strata(stratum)) # Here, case_control is not related to the case and control in IS_FOCAL, but is an indicator for vaccination status.
has.vax = survival::clogit(has.vax.formula, data = strata.df)

# Look at vaccination status or not, as a cleaner model
age.x.vax.formula = as.formula(IS_FOCAL ~ Age*IS.VAX + strata(stratum)) # Here, case_control is not related to the case and control in IS_FOCAL, but is an indicator for vaccination status.
age.x.vax = survival::clogit(age.x.vax.formula, data = strata.df)

# Vaccine Type alone
vaxtype.formula = as.formula(IS_FOCAL ~ as.factor(vaccine_type_code) + strata(stratum))
vaxtype = survival::clogit(vaxtype.formula, data=strata.df)

# Vaccine Type + Age
vaxtype.plus.age.formula = as.formula(IS_FOCAL ~ Age + as.factor(vaccine_type_code) + strata(stratum))
vaxtype.plus.age = survival::clogit(vaxtype.plus.age.formula, data=strata.df)

# Examine age bins plus vaccination status
age.bin.status.formula = as.formula(IS_FOCAL ~ as.factor(AGE.BIN) + IS.VAX + strata(stratum))
bin.status.model = survival::clogit(age.bin.status.formula, data = strata.df)
bin.status.AIC = extractAIC(bin.status.model)[2]

# BELOW TWO MODELS WERE REPLACED WITH A MERGED MRNA - WILL MAKE THE STATISTICS EASIER TO INTERPRET
# Neither these nor those were included in the final analysis.
# There was not much of an improvement, and vaccine type was cleaner to interpret.

# Add a Janssen vaccine indicator plus age to replace number of vaccine doses
#j.formula = as.formula(IS_FOCAL ~ JANSSEN + IS.VAX + strata(stratum))
#janssen = survival::clogit(j.formula, data=strata.df)
#janssen.AIC = extractAIC(janssen)[2]

# Add a Janssen vaccine indicator plus age to replace number of vaccine doses
#j.age.formula = as.formula(IS_FOCAL ~ JANSSEN + IS.VAX + Age + strata(stratum))
#janssen.age = survival::clogit(j.age.formula, data=strata.df)
#janssen.age.AIC = extractAIC(janssen.age)[2]

MRNA.formula = as.formula(IS_FOCAL ~ as.factor(MRNA_MERGED) + strata(stratum))
MRNA = survival::clogit(MRNA.formula, data=strata.df)
MRNA.AIC = extractAIC(MRNA)[2]

# Add a Janssen vaccine indicator plus age to replace number of vaccine doses
MRNA.age.formula = as.formula(IS_FOCAL ~ as.factor(MRNA_MERGED) + Age + strata(stratum))
MRNA.age = survival::clogit(MRNA.age.formula, data=strata.df)
MRNA.age.AIC = extractAIC(MRNA.age)[2]


# Also, no bin 0 individuals in this analysis, so it's using 70 - 89 as the reference group
#NOTE: Bin approach is not currently being included in the AIC comparison.
# It is not at all competitive based on AIC score

# Calculate AIC scores and identify the model with the lowest AIC
models = c("Age", "VaccineType x Age",
           "TimePostDose x Age", "FullModel",
           "HasVaccine", "HasVaccine x Age",
           "VaccineType", "VaccineType + Age",
           "VaccineType_mRNA_merged", "VaccineType_mRNA_merged + Age")
# "NumberVaccineDoses x Age", "NumberVaccineDoses", "NumberVaccineDoses + Age",
model.objects = list(age.result, vaxtype.x.age, 
                     timepostdose.x.age, full.model, 
                     has.vax, age.x.vax, vaxtype, vaxtype.plus.age,
                     MRNA, MRNA.age)
# ndoses.x.age, ndoses, ndoses.plus.age,
formula.objects = list(age.formula, vaxtype.x.age.formula, 
                       tpd.x.age.formula, full.formula, has.vax.formula, age.x.vax.formula,
                       vaxtype.formula, vaxtype.plus.age.formula,
                       MRNA.formula, MRNA.age.formula)
# nd.x.age.formula, nd.formula, nd.plus.age.formula, 
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
obj.index = aic.results$SORT.ORDER[2] # Changed to take vaccine type
# Fewer parameters than the top model with only a marginal reduction in AIC
# Top model is likely overfit
model.object = model.objects[[obj.index]] 
formula.object = formula.objects[[obj.index]]
make.leverage.plot(model.object, formula.object, out.folder, figure.label, table.label)

## Make a plot to examine the bin results (not included in final manuscript)
bin.OR = summary(bin.status.model)$coefficients[ ,2]
bin.OR.df = data.frame(BIN.OR = bin.OR[1:4], bin = seq(1,4))
bin.OR.df = bin.OR.df[order(bin.OR.df$bin, decreasing = TRUE), ]

plot(seq(1,4), bin.OR.df$BIN.OR, xaxt = 'n',
     xlab = "Bin", ylab = "Odds Ratio relative to 70-89 age group",
     ylim = c(0,4))
axis(side = 1, at = seq(1,4), labels = c("0-17","18-29","30-49","50-69"))


