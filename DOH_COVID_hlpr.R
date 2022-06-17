# Helper functions for DOH_COVID.R

#' Data Setup
#' 
#' Function to hide the mechanics of setting up the data set to streamline code and enable reuse
#' 
#' 
data.setup = function(DOH.data.file, drop.2020 = FALSE){
  
  in.data = read.csv(DOH.data.file)
  
  # Add year, month, day columns
  in.data$YEAR = sapply(in.data$collection_date, dfmip::splitter, "/", 3)
  in.data$MONTH = sapply(in.data$collection_date, dfmip::splitter, "/", 1)
  in.data$DAY = sapply(in.data$collection_date, dfmip::splitter, "/", 2)
  
  # OK. We have records in 2020, so can't just use day of year. Need an index.
  # Where do I have code for that? Mesonet folder.
  
  if (max(in.data$YEAR) > 2021){stop("2022 and beyond are not yet supported by this script. Please fix the date index to properly order the dates")}
  
  # Convert to day of year
  in.data$DOY = mapply(dfmip:::get.DOY, in.data$YEAR, in.data$MONTH, in.data$DAY)
  in.data$DATE_INDEX = in.data$DOY
  # Code change on 2021-12-07. This function is used in the initial analysis, and the change should be backwards compatible.
  if (drop.2020 == FALSE){
    in.data$DATE_INDEX[in.data$YEAR == 2021] = in.data$DATE_INDEX[in.data$YEAR == 2021] + 366
  }
  
  # Correct data set irregularities (these may no longer apply - developed for the initial data set)
  in.data$Sex = gsub("f", "F", in.data$Sex)
  in.data$Age = gsub("2 months", "0", in.data$Age)
  in.data$Age = gsub("5 days", "0", in.data$Age)
  in.data$Age = as.numeric(in.data$Age)
  
  return(in.data)
}

#' Data Setup
#' 
#' Function to hide the mechanics of setting up the data set to streamline code and enable reuse
#' v2 adds support for 2022, and includes data checks for standardized entries
#' # Reading the file is now external to the file, to allow custom data corrections with gsub and colnames
#' 
#' 
data.setup.v2 = function(new.data, start.year){

  # Check for required column names
  check.formatting(new.data)
  
  # Check required columns for required values
  check.column.values(new.data)
  
  # Add year, month, day columns
  new.data$YEAR = sapply(new.data$collection_date, dfmip::splitter, "/", 3)
  new.data$MONTH = sapply(new.data$collection_date, dfmip::splitter, "/", 1)
  new.data$DAY = sapply(new.data$collection_date, dfmip::splitter, "/", 2)
  
  new.data$INDEX_ADJUST = sapply(new.data$YEAR, calculate.adjustment, start.year)
  
  # Convert to day of year
  new.data$DOY = mapply(dfmip:::get.DOY, new.data$YEAR, new.data$MONTH, new.data$DAY)
  new.data$DATE_INDEX = new.data$DOY + new.data$INDEX_ADJUST
  new.data$INDEX_ADJUST = NULL # Remove column that is no longer needed

  return(new.data)
}

#' Adjust for a leap year
#' 
#' @year The current year (excluded from adjustment)
#' @start.year the starting year
#' 
#' @return Returns an adjustment for the number of leap years from start.year to present year minus 1.
#' 
calculate.adjustment = function(year, start.year){
  # Get a date index adjustment
  year.adjust = year - start.year
  index.adjust = year.adjust * 365
  
  # Correct for leap years
  day.adjust = 0
  if (year.adjust > 0){
    test = seq(start.year, year)
    test = test[1:(length(test) - 1)] # Exclude the present year from consideration - we're adding past years
    
    remainder = test %% 4
    day.adjust = length(remainder[remainder == 0])
  }
  
  # Adjust index.adjust for leap years
  index.adjust = index.adjust + day.adjust
  
  return(index.adjust)
}

#' Check that formatting matches what was used to develop the code
#' 
#' @new.data The input data object
#' 
check.formatting = function(new.data){
  required.fields = c("Identifier", "collection_date", "Age", "Sex", 
                      "Economic_region", "sequence_result", "vaccine_type",
                      "Days.between.collection.date.and.vaccination.complete.date",
                      "booster_type",
                      "Days.between.collection.date.and.booster.date")

  #Loop through required fields and check that they are in the new.data names
  field.names = names(new.data)
  missing.fields = c()
  for (field in required.fields){
    if (!field %in% field.names){
      missing.fields = c(missing.fields, field)
    }
  }
  
  if (length(missing.fields) > 0){
    stop(sprintf("The following fields are missing or contain a typo: %s", paste(missing.fields, collapse = ", ")))
  }
  
  optional.fields = c("Mutations", "case_control")

}

#' Check format values
#' 
check.column.values = function(new.data){
  
  problem.entries = c()
  # Check Age Column
  if (min(new.data$Age, na.rm = TRUE) < 0){
    problem.entries = c(problem.entries, "Age:ValuesBelow0")
  }
  if (max(new.data$Age, na.rm = TRUE) > 110){
    problem.entries = c(problem.entries, "Age:ValueAbove110")
  }
  if (typeof(new.data$Age) == "character"){
    problem.entries = c(problem.entries, "Age:CharacterFormat")
  }
  
  # Check Sex Column
  # Should be M, F, U based on previous data sets
  s.values = c("M","F","U")
  sexes = unique(new.data$Sex)
  for (s in sexes){
    if (!s %in% s.values){
      problem.entries = c(problem.entries, sprintf("Sex:%s", s))
    }
  }
  
  # Check Economic_region Column
  r.values = c("Long Island", "Western New York", "Mid-Hudson", "Capital Region",
                "Central New York","Finger Lakes", "Southern Tier",
               "Mohawk Valley", "North Country", "New York City")
  regions = unique(new.data$Economic_region)
  for (r in regions){
    if (!r %in% r.values){
      problem.entries = c(problem.entries, sprintf("Region:%s", r))
    }
  }
  
  # Check vaccine_type column
  v.values = c("Pfizer", "Moderna", "Janssen", "None")
  v.types = unique(new.data$vaccine_type)
  for (v in v.types){
    if (!v %in% v.values){
      problem.entries = c(problem.entries, sprintf("VaccineType:%s", v))
    }
  }

  
  # If booster data is missing, skip the booster columns
  # Check if booster_type has only a single entry
  skip_booster = 0
  if (length(unique(new.data$booster_type)) == 1){
    # Confirm the single entry is NA
    if (is.na(unique(new.data$booster_type))){
      skip_booster = 1
    }
  }
  
  if (skip_booster == 0){
    # Check booster_type column
    b.types = unique(new.data$booster_type)
    for (b in b.types){
      if (!b %in% v.values){
        problem.entries = c(problem.entries, sprintf("BoosterType:%s", b))
      }
    }
  }
  
  # Check sequence_result column
  # NOT IMPLEMENTED, but probably want to ensure there is omicron in your sample data
  
  # Check Days.between.collection.date.and.vaccination.complete.date
  if (min(new.data$Days.between.collection.date.and.vaccination.complete.date, na.rm = TRUE) < 1){
    problem.entries = c(problem.entries, "DaysPostVax:ValuesBelow1")
  }
  if (max(new.data$Days.between.collection.date.and.vaccination.complete.date, na.rm = TRUE) > 450){
    problem.entries = c(problem.entries, "DaysPostVax:ValueAbove450") # Seems implausible - Vaccine in Dec. 2020, so shouldn't be much above 426 by end of January
  }
  if (typeof(new.data$Days.between.collection.date.and.vaccination.complete.date) == "character"){
    problem.entries = c(problem.entries, "DaysPostVax:CharacterFormat")
  }
  
  if (skip_booster == 0){
    # Check Days.between.collection.date.and.booster.date
    if (min(new.data$Days.between.collection.date.and.booster.date, na.rm = TRUE) < 1){
      problem.entries = c(problem.entries, "DaysPostBoost:ValuesBelow1")
    }
    if (max(new.data$Days.between.collection.date.and.booster.date, na.rm = TRUE) > 200){
      warning("Entries with >200 days post-booster. These should be examined carefully")
      #problem.entries = c(problem.entries, "DaysPostBoost:ValueAbove200") # Seems implausible - Boosters starting in August? 2021, so shouldn't be much above 5 months (30 * 5 = 150) by end of January
    }
    if (typeof(new.data$Days.between.collection.date.and.booster.date) == "character"){
      problem.entries = c(problem.entries, "DaysPostBoost:CharacterFormat")
    }
  }
  
  # Note: there is no check for the time_between_vaccination_and_booster column!
  if (length(new.data$time_between_vaccination_and_booster[!is.na(new.data$time_between_vaccination_and_booster) & new.data$time_between_vaccination_and_booster <= 30] > 0)){
    warning("Less than 30 days between vaccination and booster. Should these individuals be excluded from the analysis?")
  }
  
  # Throw an error if something is wrong
  if (length(problem.entries) > 0){
    stop(sprintf("Data entries did not match pre-planned formatting. Please fix the following entries:
         %s", paste(problem.entries, collapse = ', ')))
  }
}


fix.sim.data.1 = function(new.data){
  new.data$Identifier = seq(1,nrow(new.data)) # Arbitrary, since this is simulation data and there is no linked information
  new.data$collection_date = new.data$COLDATE
  new.data$COLDATE = NULL
  new.data$Age = new.data$AGE
  new.data$AGE = NULL
  new.data$Sex = new.data$SEX
  new.data$SEX = NULL
  new.data$Economic_region = new.data$REGION
  new.data$REGION = NULL
  new.data$sequence_result = new.data$STRAIN
  new.data$STRAIN = NULL
  new.data$vaccine_type = new.data$VACCINE_TYPE
  new.data$VACCINE_TYPE = NULL
  new.data$booster_type = new.data$BOOST_TYPE
  new.data$BOOST_TYPE = NULL
  
  new.data$Days.between.collection.date.and.vaccination.complete.date =
    date.difference(new.data$collection_date, new.data$VACCINE_DATE)
  new.data$Days.between.collection.date.and.booster.date =
    date.difference(new.data$collection_date, new.data$BOOST_DATE)

  return(new.data)
}

date.difference = function(dates1, dates2){

  differences = mapply(calc.date.difference, dates1, dates2)
  
  return(differences)
}

# Test code for calc.date.difference
#calc.date.difference("09/01/2021", "08/31/2021")
#calc.date.difference("12/31/2021", "01/01/2021")
#calc.date.difference("12/31/2020", "01/01/2020")
#calc.date.difference("01/01/2021", "12/31/2020")
#calc.date.difference("01/01/2022", "12/31/2021")
#calc.date.difference("01/01/2022", "12/1/2021")
#calc.date.difference("01/01/2022", "11/30/2021")
#calc.date.difference("01/02/2021", "01/02/2020")
#calc.date.difference("01/02/2022", "01/02/2021")
#calc.date.difference("03/01/2021", "03/01/2020")
#calc.date.difference("03/01/2022", "03/01/2021")
#calc.date.difference("2021-03-01", "2021-02-28")
#calc.date.difference("03/01/2021", "2021-02-28")
#calc.date.difference(NA, "2021-02-28")
#calc.date.difference("03/01/2022", NA)

#' Calculate the difference in days between two dates in MM/DD/YYYY format. Should also work for M/D/YY format
#' WILL NOT WORK FOR 1900 and 2100
calc.date.difference = function(x, y){
  
  # If one value is NA, return NA
  if (is.na(x) | is.na(y)){
    days = NA
  }else{
    # Split date apart based on /. If NOT 3 parts, throw an error.
    x.parts = strsplit(x, '/')[[1]]
    if (length(x.parts) != 3){ stop("Problem with input date 1. Dates should be in MM/DD/YYYY")}
    x.month = as.numeric(x.parts[1])
    x.day = as.numeric(x.parts[2])
    x.year = as.numeric(x.parts[3])
    
    y.parts = strsplit(y, "/")[[1]]
    if (length(y.parts) != 3){ stop("Problem with input date 2. Dates should be in MM/DD/YYYY")}
    y.month = as.numeric(y.parts[1])
    y.day = as.numeric(y.parts[2])
    y.year = as.numeric(y.parts[3])
    
    # If in same year, just take the difference
    if (x.year == y.year){
      days = dfmip:::get.DOY(x.year, x.month, x.day) - dfmip:::get.DOY(y.year, y.month, y.day)
    }else{
      x.doy = dfmip:::get.DOY(x.year, x.month, x.day)
      y.doy = dfmip:::get.DOY(y.year, y.month, y.day)
      
      # Get difference in years
      year.diff = x.year - y.year - 1 # Minus 1 drops the y.year, because that will be handled differently
      
      # Check for intermediate leap years and add a day for those
      years = seq(y.year, x.year)
      if (length(years > 2)){
        years = years[2:(length(years) - 1)] # remove endpoint years
        leap.adjust = length(years[years %% 4 == 0])
      }
      
      days.in.year = 365
      if (y.year %% 4 == 0){ days.in.year = 365 }
      y.adjust = days.in.year - y.doy # count the days FOLLOWING the y.doy between y.doy and x.doy.
      days = x.doy + year.diff * 365 + leap.adjust + y.adjust
    }
  }
  
  return(days)
}

#' Match Cases
#' 
#' ORDER OF REMOVAL COULD MATTER SUBSTANTIALLY - create a second version that can look at potential matches
#' 
case.control.match = function(in.data, in.seed, date.offset, age.offset){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  case.df = in.data[in.data$case_control == "case", ] # 608
  case.df$CONTROL_ID = NA # Initialize a column to store control identity
  case.df$CONTROL_SEQ = NA # Initialize a column to store control strain information
  control.df = in.data[in.data$case_control == "control ", ] # 5292
  
  no.match.vec = c()
  
  # Loop through cases and find matches
  for (i in 1:nrow(case.df)){
    case.id = case.df$Identifier[i]
    case.age = case.df$Age[i]
    case.sex = case.df$Sex[i]
    case.region = case.df$Economic_region[i]
    case.date = case.df$DATE_INDEX[i]
    
    # Set a no-match indicator
    no.match = 0
    
    # Match based on Date (+/- 3 days)
    #date.offset = 3
    date.match.index = control.df$DATE_INDEX >= (case.date - date.offset) & control.df$DATE_INDEX <= (case.date + date.offset)
    date.matches = control.df[date.match.index, ]
    # For i = 1, there are 131 potential matches.
    
    if (nrow(date.matches) == 0){ no.match = 1 }
    
    # Match based on Region (# Kustin say municipality, geographic region, and sector, whatever that means. For us it means region.)
    if (no.match == 0){
      region.matches = date.matches[date.matches$Economic_region == case.region, ] # 55 that match on date and region
      if (nrow(region.matches) == 0){ no.match = 1 }
    }
    
    # Match based on Age (+/- 10 years)
    if (no.match == 0){
      #age.offset = 10
      age.match.index = region.matches$Age >= (case.age - age.offset) & region.matches$Age <= (case.age + age.offset)
      age.matches = region.matches[age.match.index, ] # 21 that match on age, date, and region for i = 1
      if (nrow(age.matches) == 0){ no.match = 1 }    
    }
    
    # Match based on sex (M vs. F)
    if (no.match == 0){
      sex.matches = age.matches[age.matches$Sex == case.sex, ] # 10 that match on all criteria
      if (nrow(sex.matches) == 0){ no.match = 1 }
    }
    
    if (no.match == 0){
      # Select an identifier at random from remaining cases (if more than one case)
      if (nrow(sex.matches) == 1){
        control.id = sex.matches$Identifier[1]
        control.seq = sex.matches$seqeunce_result[1]
      }else{
        control.id = sample(sex.matches$Identifier, 1) # This will not work if there is only one row.
        control.seq = sex.matches$seqeunce_result[sex.matches$Identifier == control.id] 
      }
      case.df$CONTROL_ID[i] = control.id
      case.df$CONTROL_SEQ[i] = control.seq
      
      # Remove the control from the control data base to ensure each case is matched to only one control.
      control.df = control.df[control.df$Identifier != control.id, ]
    }else{
      no.match.vec = c(no.match.vec, case.id)
    }
  }
  
  return(list(case.df, control.df, no.match.vec))
}


#' Match Cases (Development version)
#' 
#' ORDER OF REMOVAL COULD MATTER SUBSTANTIALLY - create a second version that can look at potential matches
#' 
#' # age.offset.override added 12/20/21. Should be a backwards compatible change
#' # Changed to not be backwards compatible on 2022-01-12 - seqeunce_result changed to sequence_result
#' 
case.control.match.dev = function(in.data, in.seed, date.offset, age.offset,
                                  removal = TRUE, age.offset.override = NA){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  case.df = in.data[tolower(in.data$case_control) == "case", ] # 608
  case.df$CONTROL_ID = NA # Initialize a column to store control identity
  case.df$CONTROL_SEQ = NA # Initialize a column to store control strain information
  case.df$MF_FLAG = NA
  control.df = in.data[tolower(in.data$case_control) == "control" | in.data$case_control == "control ", ] # 5292
  
  no.match.vec = c()
  
  # Loop through cases and find matches
  for (i in 1:nrow(case.df)){
    case.id = case.df$Identifier[i]
    case.age = case.df$Age[i]
    case.sex = case.df$Sex[i]
    case.region = case.df$Economic_region[i]
    case.date = case.df$DATE_INDEX[i]
    
    # If case sex is not known, DO NOT match on sex
    unknown.sex = 0
    if (case.sex == "0"){ unknown.sex = 1 }
    
    # Set a no-match indicator
    no.match = 0
    
    # Match based on Date (+/- 3 days)
    #date.offset = 3
    date.match.index = control.df$DATE_INDEX >= (case.date - date.offset) & control.df$DATE_INDEX <= (case.date + date.offset)
    date.matches = control.df[date.match.index, ]
    # For i = 1, there are 131 potential matches.
    
    if (nrow(date.matches) == 0){ no.match = 1 }
    
    # Match based on Region (# Kustin say municipality, geographic region, and sector, whatever that means. For us it means region.)
    if (no.match == 0){
      region.matches = date.matches[date.matches$Economic_region == case.region, ] # 55 that match on date and region
      if (nrow(region.matches) == 0){ no.match = 1 }
    }
    
    # Match based on Age (+/- 10 years)
    if (no.match == 0){
      #age.offset = 10
      
      # Change to a fixed case age structure if override set to 1.
      # This allows us to match +/- 5 years for younger people and +/- 10 years for older
      # without a discontinuity. This will OVERRIDE the input age.offset, NOT adjust it.
      if (!is.na(age.offset.override)){
        age.offset = get.age.offset(case.age)
      }
      
      age.match.index = region.matches$Age >= (case.age - age.offset) & region.matches$Age <= (case.age + age.offset)
      age.matches = region.matches[age.match.index, ] # 21 that match on age, date, and region for i = 1
      if (nrow(age.matches) == 0){ no.match = 1 }    
    }
    
    mf.flag = 1
    # Match based on sex (M vs. F)
    if (no.match == 0 & unknown.sex == 0){
      sex.matches = age.matches[age.matches$Sex == case.sex, ] # 10 that match on all criteria
      if (nrow(sex.matches) > 0){ mf.flag = 0 }
    }
    
    # Remove the matching sex requirement if not a match. Still require a match on the above criteria
    if (mf.flag == 1 & no.match == 0){
      # Just select from the age match data set - renamed because the sex.matches name was already in use in the code below.
      sex.matches = age.matches
    }
    
    
    if (no.match == 0){
      # Select an identifier at random from remaining cases (if more than one case)
      if (nrow(sex.matches) == 1){
        control.id = sex.matches$Identifier[1]
        control.seq = sex.matches$sequence_result[1]
      }else{
        control.id = sample(sex.matches$Identifier, 1) # This will not work if there is only one row.
        control.seq = sex.matches$sequence_result[sex.matches$Identifier == control.id] 
      }
      case.df$CONTROL_ID[i] = control.id
      case.df$CONTROL_SEQ[i] = control.seq
      case.df$MF_FLAG[i] = mf.flag # Flag records where an exact sex match was not used
      
      if (removal == TRUE){
        # Remove the control from the control data base to ensure each case is matched to only one control.
        control.df = control.df[control.df$Identifier != control.id, ]
      }
    }else{
      no.match.vec = c(no.match.vec, case.id)
    }
  }
  
  return(list(case.df, control.df, no.match.vec))
}


# unit tests for get.age.offset
# get.age.offset(1)  # 5
# get.age.offset(10) # 5
# get.age.offset(20) # 5
# get.age.offset(21) # 5
# get.age.offset(22) # 6
# get.age.offset(23) # 6
# get.age.offset(24) # 7
# get.age.offset(25) # 7
# get.age.offset(26) # 8
# get.age.offset(27) # 8
# get.age.offset(28) # 9
# get.age.offset(29) # 9
# get.age.offset(30) # 10
# get.age.offset(31) # 10
# get.age.offset(42) # 10
# get.age.offset(63) # 10
# get.age.offset(112) # 10
# 
# 
#' Create a gradient age offset based on case age
#' 
#' @param case.age age of the case individual
#' 
#' @return age.offset Cases will be matched to controls that are +/- the age.offset.
#' So an age.offset of 5 and a case.age of 20, will match with 15 - 25 year olds.
#' 
get.age.offset = function(case.age){
  age.offset = NA
  if (case.age <= 21){  age.offset = 5  }
  if (case.age > 21 & case.age <= 23){ age.offset = 6 }
  if (case.age > 23 & case.age <= 25){ age.offset = 7 }
  if (case.age > 25 & case.age <= 27){ age.offset = 8 }
  if (case.age > 27 & case.age <= 29){ age.offset = 9 }
  if (case.age > 29){ age.offset = 10 }
  
  return(age.offset)
}

#' Match Cases (Development version)
#' 
#' Identify all records that match a given case - look at options for order of selection
#' 
case.control.options = function(in.data, in.seed, date.offset, age.offset){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  case.df = in.data[in.data$case_control == "case", ] # 608
  case.df$CONTROL_OPTIONS = NA # Initialize a column to store control identity
  case.df$MATCH_FAIL = NA
  #case.df$CONTROL_SEQ = NA # Initialize a column to store control strain information
  control.df = in.data[in.data$case_control == "control ", ] # 5292
  
  no.match.vec = c()
  
  # Loop through cases and find matches
  for (i in 1:nrow(case.df)){
    case.id = case.df$Identifier[i]
    case.age = case.df$Age[i]
    case.sex = case.df$Sex[i]
    case.region = case.df$Economic_region[i]
    case.date = case.df$DATE_INDEX[i]
    
    # Set a no-match indicator
    no.match = ""
    
    # Match based on Date (+/- 3 days)
    #date.offset = 3
    date.match.index = control.df$DATE_INDEX >= (case.date - date.offset) & control.df$DATE_INDEX <= (case.date + date.offset)
    date.matches = control.df[date.match.index, ]
    # For i = 1, there are 131 potential matches.
    
    if (nrow(date.matches) == 0){ no.match = "D" }
    
    # Match based on Region (# Kustin say municipality, geographic region, and sector, whatever that means. For us it means region.)
    if (no.match == ""){
      region.matches = date.matches[date.matches$Economic_region == case.region, ] # 55 that match on date and region
      if (nrow(region.matches) == 0){ no.match = "R" }
    }
    
    # Match based on Age (+/- 10 years)
    if (no.match == ""){
      #age.offset = 10
      age.match.index = region.matches$Age >= (case.age - age.offset) & region.matches$Age <= (case.age + age.offset)
      age.matches = region.matches[age.match.index, ] # 21 that match on age, date, and region for i = 1
      if (nrow(age.matches) == 0){ no.match = "A" }    
    }
    
    # Match based on sex (M vs. F)
    if (no.match == ""){
      sex.matches = age.matches[age.matches$Sex == case.sex, ] # 10 that match on all criteria
      if (nrow(sex.matches) == 0){ no.match = "S" }
    }
    
    if (no.match == ""){
      # Select an identifier at random from remaining cases (if more than one case)
      case.df$CONTROL_OPTIONS[i] = paste(sex.matches$Identifier, collapse = ",")
      
    }
    case.df$MATCH_FAIL[i] = no.match
      
  }
  
  return(case.df)
}


#' Match Cases (Development version)
#' 
#' Identify all records that match a given case - look at options for order of selection
#' This version looks at relaxing each criterion to give a full picture of options.
#' Not really sure how we'll SELECT from all those options, but... it at least gives us
#' the options on the table.
#' 
#' age.override added on 12/20/21 to provide a means of matching with more restricted values for younger ages.
#' 
case.control.options.v2 = function(in.data, in.seed, date.offset, age.offset,
                                   age.offset.override = NA){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  case.df = in.data[tolower(in.data$case_control) == "case", ] # 608
  case.df$EXACT_OPTIONS = NA # Initialize a column to store control identity
  case.df$NO_MF_OPTIONS = NA
  case.df$NO_AGE_OPTIONS = NA
  case.df$NO_REG_OPTIONS = NA
  case.df$NO_DATE_OPTION = NA
  case.df$N_EXACT = NA
  case.df$N_NO_MF = NA
  case.df$N_NO_AGE = NA
  case.df$N_NO_REG = NA
  case.df$N_NO_DATE = NA
  control.df = in.data[tolower(in.data$case_control) == "control" | tolower(in.data$case_control) == "control ", ] # 5292
  
  case.df = update.options(case.df, control.df, "exact", date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, "mf", date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, 'age', date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, 'region', date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, 'date', date.offset, age.offset, age.offset.override)
  
  return(case.df)
}


#' Create a function to allow relaxation of a selected criterion
update.options = function(case.df, control.df, option.type,
                          date.offset, age.offset, age.offset.override){

  # Loop through cases and find matches
  for (i in 1:nrow(case.df)){
    case.id = case.df$Identifier[i]
    case.age = case.df$Age[i]
    case.sex = case.df$Sex[i]
    case.region = case.df$Economic_region[i]
    case.date = case.df$DATE_INDEX[i]
    
    # Copy control.df into an object that can be selectively modified
    option.df = control.df

    if (option.type != 'date'){
      # Match based on Date (+/- date.offset)
      date.match.index = option.df$DATE_INDEX >= (case.date - date.offset) & option.df$DATE_INDEX <= (case.date + date.offset)
      option.df = option.df[date.match.index, ]
    }

    # Match based on Region
    if (option.type != 'region'){
      option.df = option.df[option.df$Economic_region == case.region, ]
    }

    # Match based on Age (+/- date.offset years)
    if (option.type != 'age'){
      # If we want a tiered system for matching, we can override the age matching flexibility and replace it with
      # a fixed tiered system.
      if (!is.na(age.offset.override)){
        age.offset = get.age.offset(case.age)
      }
      
      age.match.index = option.df$Age >= (case.age - age.offset) & option.df$Age <= (case.age + age.offset)
      option.df = option.df[age.match.index, ]
    }

    # Match based on sex (M vs. F)
    if (option.type != 'mf'){
      option.df = option.df[option.df$Sex == case.sex, ]
    }
    
    if (option.type == 'exact') {
      update.column = "EXACT_OPTIONS"
      update.column.2 = "N_EXACT"
    }
    if (option.type == "mf" )   {
      update.column = "NO_MF_OPTIONS"
      update.column.2 = "N_NO_MF"
    }
    if (option.type == "date")  {
      update.column = "NO_DATE_OPTION"
      update.column.2 = "N_NO_DATE"
    }
    if (option.type == "region"){
      update.column = "NO_REG_OPTIONS"
      update.column.2 = "N_NO_REG"
    }
    if (option.type == "age")   {
      update.column = "NO_AGE_OPTIONS"
      update.column.2 = "N_NO_AGE"
    }
    
    # Update the appropriate column
    case.df[[update.column]][i] = paste(option.df$Identifier, collapse = ",")
    case.df[[update.column.2]][i] = length(option.df$Identifier)
      
  }
  return(case.df)
}

#' Create a function to allow relaxation of a selected criterion
#' 
#' v4 corresponds to case.control.options.v4
#' 
update.options.v4 = function(case.df, control.df, option.type,
                          date.offset, age.bins){
  
  # Loop through cases and find matches
  for (i in 1:nrow(case.df)){
    case.id = case.df$Identifier[i]
    case.age = case.df$Age[i]
    case.sex = case.df$Sex[i]
    case.region = case.df$Economic_region[i]
    case.date = case.df$DATE_INDEX[i]
    
    # Copy control.df into an object that can be selectively modified
    option.df = control.df
    
    if (option.type != 'date'){
      # Match based on Date (+/- date.offset)
      date.match.index = option.df$DATE_INDEX >= (case.date - date.offset) & option.df$DATE_INDEX <= (case.date + date.offset)
      option.df = option.df[date.match.index, ]
    }
    
    # Match based on Region
    if (option.type != 'region'){
      option.df = option.df[option.df$Economic_region == case.region, ]
    }
    
    # Match based on Age (+/- date.offset years)
    if (option.type != 'age'){

      # Test case ages
      #case.age = 0  # 0 4
      #case.age = 4  # 0 4
      #case.age = 5  # 5 11
      #case.age = 11 # 5 11
      #case.age = 17 # 12 17
      #case.age = 28 # 18 29
      #case.age = 90 # 90 110
      #case.age = 111 # 111 110 # This will cause problems, but should be prevented upstream.
      
      # identify bin of focal individual
      bin.min = 0
      for (j in 1:length(age.bins)){
        bin.max = age.bins[j]
        
        # Stop the loop when max and min bound the case age.
        if (case.age >= bin.min & case.age <= bin.max){
          break
        }
        
        # otherwise increment the bins
        bin.min = bin.max + 1
      }
      
      age.match.index = option.df$Age >= bin.min & option.df$Age <= bin.max
      option.df = option.df[age.match.index, ]
    }
    
    # Match based on sex (M vs. F)
    if (option.type != 'mf'){
      option.df = option.df[option.df$Sex == case.sex, ]
    }
    
    if (option.type == 'exact') {
      update.column = "EXACT_OPTIONS"
      update.column.2 = "N_EXACT"
    }
    if (option.type == "mf" )   {
      update.column = "NO_MF_OPTIONS"
      update.column.2 = "N_NO_MF"
    }
    if (option.type == "date")  {
      update.column = "NO_DATE_OPTION"
      update.column.2 = "N_NO_DATE"
    }
    if (option.type == "region"){
      update.column = "NO_REG_OPTIONS"
      update.column.2 = "N_NO_REG"
    }
    if (option.type == "age")   {
      update.column = "NO_AGE_OPTIONS"
      update.column.2 = "N_NO_AGE"
    }
    
    # Update the appropriate column
    case.df[[update.column]][i] = paste(option.df$Identifier, collapse = ",")
    case.df[[update.column.2]][i] = length(option.df$Identifier)
    
  }
  return(case.df)
}



#' Convert to a power calculator for a given scale factor, group1 efficacy and group2 efficacy
#' 
#' For a McNemar's Exact test, one-sided.
#' 
#' @param wt The vaccine efficacy for the wild-type of the virus
#' @param vox The vaccine efficacy for the variant of concern virus strain
#' @param p.wt The proportion of circulating strains that are 'wild-type'
#' @param p.voc The proportion of circulating strains that are 'variant of concern'
#' @param n.vax.cases The number of vaccinated cases used in the analysis
#' 
estimate.power = function(wt, voc, p.wt, p.voc,
                          n.vax.cases){
  
  # Scale.factor should be irrelevant - I think you get the same answer no matter what,
  # since it is deterministic, and n.vax.cases controls the sample size in the
  # actual power analysis.
  # scale.factor is number of control COVID cases. It cancels out in the math, but it is helpful for conceptualization
  scale.factor = 100
  
  #scale.factor = n.unvaccinated.infections
  wt.cases = (1 - wt)*scale.factor*p.wt # These are break-through cases
  voc.cases = (1 - voc) * scale.factor*p.voc
  total.cases = wt.cases + voc.cases
  # Note that this is the total vaccinated cases expected for 100 unvaccinated cases.
  # n.vax.cases is the total vaccinated cases used to calculate statistical power.
  
  a = wt.cases * p.wt # (wt & wt)
  b = wt.cases * p.voc # wt vax & voc unvax
  c = voc.cases * p.wt # voc vax & wt unvax
  d = voc.cases * p.voc # voc & voc
  
  # so pb in this scenario is:
  pb = b/total.cases #2.5 / 25 # 0.1
  # pc is:
  pc = c / total.cases # 10/25 # 0.4
  
  # Calculate statistical power
  result = powerPaired2x2(pb, pc, n.vax.cases, alternative = c('one.sided'))
  
  return(c(result$power, pb, pc, b, c))
  
}


#' Match Cases (age bin version)
#' 
#' ORDER OF REMOVAL COULD MATTER SUBSTANTIALLY - create a second version that can look at potential matches
#' Inclusive of bin start, exclusive of bin end (i.e. 10 - 20 would be 10 - 19, because it is exclusive of 20)
#' 
case.control.match.bin = function(in.data, in.seed, date.offset, bin.width,
                                  bin.start, removal = TRUE){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  case.df = in.data[tolower(in.data$case_control) == "case", ] # 608
  case.df$CONTROL_ID = NA # Initialize a column to store control identity
  case.df$CONTROL_SEQ = NA # Initialize a column to store control strain information
  case.df$MF_FLAG = NA
  control.df = in.data[tolower(in.data$case_control) == "control" | in.data$case_control == "control ", ] # 5292
  
  no.match.vec = c()
  
  # Create age bins
  max.age = max(case.df$Age)
  bins = seq(bin.start,(max.age + bin.width), bin.width) # Create bins from bin start to the upper limit beyond the maximum age, using bin-width as an interval.
  case.df$Bin = sapply(case.df$Age, assign.bin, bins)
  control.df$Bin = sapply(control.df$Age, assign.bin, bins)
  
  # Loop through cases and find matches
  for (i in 1:nrow(case.df)){
    case.id = case.df$Identifier[i]
    case.age = case.df$Age[i]
    case.sex = case.df$Sex[i]
    case.region = case.df$Economic_region[i]
    case.date = case.df$DATE_INDEX[i]
    case.bin = case.df$Bin[i]
    
    # If case sex is not known, DO NOT match on sex
    unknown.sex = 0
    if (case.sex == "0"){ unknown.sex = 1 }
    
    # Set a no-match indicator
    no.match = 0
    
    # Match based on Date (+/- 3 days)
    #date.offset = 3
    date.match.index = control.df$DATE_INDEX >= (case.date - date.offset) & control.df$DATE_INDEX <= (case.date + date.offset)
    date.matches = control.df[date.match.index, ]
    # For i = 1, there are 131 potential matches.
    
    if (nrow(date.matches) == 0){ no.match = 1 }
    
    # Match based on Region (# Kustin say municipality, geographic region, and sector, whatever that means. For us it means region.)
    if (no.match == 0){
      region.matches = date.matches[date.matches$Economic_region == case.region, ] # 55 that match on date and region
      if (nrow(region.matches) == 0){ no.match = 1 }
    }
    
    # Match based on Age (bin)
    if (no.match == 0){
      #age.offset = 10
      age.matches = region.matches[region.matches$Bin == case.bin, ] # 21 that match on age, date, and region for i = 1
      if (nrow(age.matches) == 0){ no.match = 1 }    
    }
    
    mf.flag = 1
    # Match based on sex (M vs. F)
    if (no.match == 0 & unknown.sex == 0){
      sex.matches = age.matches[age.matches$Sex == case.sex, ] # 10 that match on all criteria
      if (nrow(sex.matches) > 0){ mf.flag = 0 }
    }
    
    # Remove the matching sex requirement if not a match. Still require a match on the above criteria
    if (mf.flag == 1 & no.match == 0){
      # Just select from the age match data set - renamed because the sex.matches name was already in use in the code below.
      sex.matches = age.matches
    }
    
    
    if (no.match == 0){
      # Select an identifier at random from remaining cases (if more than one case)
      if (nrow(sex.matches) == 1){
        control.id = sex.matches$Identifier[1]
        control.seq = sex.matches$seqeunce_result[1]
      }else{
        control.id = sample(sex.matches$Identifier, 1) # This will not work if there is only one row.
        control.seq = sex.matches$seqeunce_result[sex.matches$Identifier == control.id]
      }
      case.df$CONTROL_ID[i] = control.id
      case.df$CONTROL_SEQ[i] = control.seq
      case.df$MF_FLAG[i] = mf.flag # Flag records where an exact sex match was not used
      
      if (removal == TRUE){
        # Remove the control from the control data base to ensure each case is matched to only one control.
        control.df = control.df[control.df$Identifier != control.id, ]
      }
    }else{
      no.match.vec = c(no.match.vec, case.id)
    }
  }
  
  return(list(case.df, control.df, no.match.vec))
}

# Test code for bin.index function
#bins = seq(10,110,10)
#assign.bin(9, bins) #NA
#assign.bin(11, bins) # 1
#assign.bin(15, bins) # 1
#assign.bin(20, bins) # 2
#assign.bin(111, bins) #NA
#assign.bin(109, bins) # 10

#' Assign a given age to a bin
#' 
#' @param Age The input age in years
#' @param bins A list of bins to use for age classification
#' 
#' @return bin.index The index for the corresponding bin. For example 1 would be the first element from bins, etc.
#' 
assign.bin = function(Age, bins){
  
  bin.index = NA
  for (i in 1:(length(bins) - 1)){
    bin.upper = bins[i + 1]
    bin.lower = bins[i]
    
    if (Age >= bin.lower & Age < bin.upper){
      bin.index = i 
    }
  }
  return(bin.index)
}

#' Create a function for looking at prevalence over time
#' 
#' Assumes a DATE_INDEX column is present
#' Uses the date from the case, not the date of the matched control (for ease of programming)
#' So a matched control will be assigned to the same week as its associated case.
#' 
#' @param in.df the data.frame to consider
#' @param bin.start the index date to start using for binning
#' @param bin.size the number of days to include in the bin. 
#' 
#' 
make.bin.prevalence.plot = function(in.df, bin.start, bin.size,
                                    case.column.name, control.column.name,
                                    b, c, p, OR,
                                    this.mutation, issue.warning = TRUE){
  
  if(issue.warning == TRUE){
    warning("Prevalence is calculated based on the case/control data. These
          may not be representative of the population as a whole")
  }
  # Calculate prevalence for each week in the period covered by the data
  min.date = min(in.df$DATE_INDEX)
  max.date = max(in.df$DATE_INDEX)
  bin.breaks = seq(bin.start, max.date + bin.size, bin.size)
  
  sample.size.vec = c()
  prevalence.vec = c()
  delta.prev.vec = c()
  
  for (j in 1:(length(bin.breaks)-1)){
    bin.lower = bin.breaks[j]
    bin.upper = bin.breaks[j + 1]
    
    df.subset = in.df[in.df$DATE_INDEX >= bin.lower & in.df$DATE_INDEX < bin.upper, ]
    
    sample.size.vec = c(sample.size.vec, (nrow(df.subset)*2)) # *2 to account for a case and a control from each record
    
    # If there is no sample size, just assign to NA
    if (nrow(df.subset) == 0){
      prevalence.vec = c(prevalence.vec, NA)
      delta.prev.vec = c(delta.prev.vec, NA)
      # Otherwise, calculate prevalence from case and control
    }else{
      # Data needs to be 0/1, otherwise this won't work.
      occurrences = sum(df.subset[[case.column.name]]) + sum(df.subset[[control.column.name]])
      prevalence = occurrences / (nrow(df.subset) * 2) # Divide total by number of rows times two- because there is a case and a control in each row
      prevalence.vec = c(prevalence.vec, prevalence)
      
      # Calculate mutations prevalence among delta strains
      case.occurrences = sum(df.subset[[case.column.name]][df.subset$CASE_DELTA == 1])
      control.occurrences = sum(df.subset[[control.column.name]][df.subset$CTRL_DELTA == 1])
      delta.occurrences = case.occurrences + control.occurrences
      case.n = length(df.subset[[case.column.name]][df.subset$CASE_DELTA == 1])
      control.n = length(df.subset[[control.column.name]][df.subset$CTRL_DELTA == 1])
      delta.n = case.n + control.n
      
      if (delta.n == 0){
        delta.prev.vec = c(delta.prev.vec, NA)
      }else{
        delta.prevalence = delta.occurrences / delta.n
        delta.prev.vec = c(delta.prev.vec, delta.prevalence)
      }
    }
  }
  
  # Plot calculated prevalences
  par(mar = c(4,4,2,4))
  
  # Plot sample size on right axis to give some gauge in onfidence in prevalence.
  plot(bin.breaks[1:length(bin.breaks) - 1], sample.size.vec, xaxt = 'n', yaxt = 'n',
       xlab = "", ylab = "", col = 'white')
  lines(bin.breaks[1:length(bin.breaks) - 1], sample.size.vec, col = 'blue')
  axis(side = 4)
  mtext(sprintf("Sample Size (per %s days)", bin.size), line = 2, side = 4,
        col = 'blue')
  
  par(new = TRUE)
  
  plot(bin.breaks[1:length(bin.breaks) - 1], prevalence.vec, ylim =c(0,1),
       xlab = "Date Index", ylab = "Mutation Prevalence")
  lines(bin.breaks[1:length(bin.breaks) - 1], prevalence.vec)

  # Add prevelence relative to Delta strains
  par(new = TRUE)
  plot(bin.breaks[1:length(bin.breaks) - 1], delta.prev.vec, ylim = c(0,1),
       xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", cex = 0.3, col = 'red', pch = 16) # Make it really small and invisible
  lines(bin.breaks[1:length(bin.breaks) - 1], delta.prev.vec, col = 'red') # Now make it show up as lines
  mtext("Prevalence in Delta strains", side = 2, col = 'red', line = 2)
  
  mtext(this.mutation, side = 3, line = -1)
  mtext(sprintf("Odds Ratio: %.3f, p = %.3f, b = %s, c = %s", OR, p, b, c))
    
  return(list(prevalence.vec, sample.size.vec))
  
}


#' Match focal group with non-focal group
#' 
#' Order of matching matters, and is NOT controlled for directly in this function
#' Created based on case.control.match.dev function, but adapted to flipped design
#' # age.offset.override added 12/20/21. Should be a backwards compatible change
#' 
match.focal = function(in.data, in.seed, date.offset, age.offset,
                                  removal = TRUE, age.offset.override = NA,
                       exclude.age = FALSE){
  # Set random number to ensure reproducibility
  set.seed(in.seed)

  strata.df = in.data
  strata.df$stratum = NA
  # Make each delta case a stratum
  n.strata = nrow(strata.df[strata.df$IS_FOCAL == 1, ])
  strata.df$stratum[strata.df$IS_FOCAL == 1] = seq(1, n.strata)
  
  # Create a flag for non-matches based on MF
  strata.df$MF_FLAG = 0

  # Set up a data frame of controls to draw from
  control.df = strata.df[strata.df$IS_FOCAL == 0, ]

  no.match.vec = c()
  
  # Loop through strata matching
  for (i in 1:n.strata){
    stratum = strata.df[strata.df$stratum == i, ]
    stratum = stratum[!is.na(stratum$stratum), ]
    # The below works, because there is only one entry in stratum. If there are multiple entries in stratum, this part of the code will fail
    # but if there are multiple entries in stratum, the code has already failed.
    case.id = stratum$Identifier
    case.age = stratum$Age
    case.sex = stratum$Sex
    case.region = stratum$Economic_region
    case.date = stratum$DATE_INDEX
    
    # If case sex is not known, DO NOT match on sex
    unknown.sex = 0
    if (case.sex == "0"){ unknown.sex = 1 }
    
    # Set a no-match indicator
    no.match = 0
    
    # Match based on Date
    date.match.index = control.df$DATE_INDEX >= (case.date - date.offset) & control.df$DATE_INDEX <= (case.date + date.offset)
    date.matches = control.df[date.match.index, ]

    if (nrow(date.matches) == 0){ no.match = 1 }
    
    # Match based on Region (# Kustin say municipality, geographic region, and sector, whatever that means. For us it means region.)
    if (no.match == 0){
      region.matches = date.matches[date.matches$Economic_region == case.region, ]
      if (nrow(region.matches) == 0){ no.match = 1 }
    }
    
    # Skip the age matching, but keep the variable name for down-stream code compatibility
    if (exclude.age == TRUE){
      if (no.match == 0){
        age.matches = region.matches
      }
    # Otherwise, match based on age.
    }else{
      # Match based on Age
      if (no.match == 0){
        
        # Change to a fixed case age structure if override set to 1.
        # This allows us to match +/- 5 years for younger people and +/- 10 years for older
        # without a discontinuity. This will OVERRIDE the input age.offset, NOT adjust it.
        if (!is.na(age.offset.override)){
          age.offset = get.age.offset(case.age)
        }
        
        age.match.index = region.matches$Age >= (case.age - age.offset) & region.matches$Age <= (case.age + age.offset)
        age.matches = region.matches[age.match.index, ] 
        if (nrow(age.matches) == 0){ no.match = 1 }    
      }
      
    }
    
    
    mf.flag = 1
    # Match based on sex (M vs. F)
    if (no.match == 0 & unknown.sex == 0){
      sex.matches = age.matches[age.matches$Sex == case.sex, ]
      if (nrow(sex.matches) > 0){ mf.flag = 0 }
    }
    
    # Remove the matching sex requirement if not a match. Still require a match on the above criteria
    if (mf.flag == 1 & no.match == 0){
      # Just select from the age match data set - renamed because the sex.matches name was already in use in the code below.
      sex.matches = age.matches
    }
    
    
    if (no.match == 0){
      # Select an identifier at random from remaining cases (if more than one case)
      if (nrow(sex.matches) == 1){
        control.id = sex.matches$Identifier[1]
      }else{
        control.id = sample(sex.matches$Identifier, 1) # This will not work if there is only one row.
      }
      # Assign the control entry to the matched stratum
      strata.df$stratum[strata.df$Identifier == control.id] = i
      strata.df$MF_FLAG[strata.df$Identifier == control.id] = mf.flag # Flag records where an exact sex match was not used

      if (removal == TRUE){
        # Remove the control from the control data base to ensure each case is matched to only one control.
        control.df = control.df[control.df$Identifier != control.id, ]
      }
    }else{
      no.match.vec = c(no.match.vec, case.id)
    }
  }
  
  # Discard any strata with no matches
  strat.table = table(strata.df$stratum) 
  keeps = names(strat.table)[strat.table > 1]
  out.strata = strata.df[strata.df$stratum %in% keeps, ]

  return(out.strata)
}


#' Match focal group with non-focal group
#' 
#' Order of matching matters, and is NOT controlled for directly in this function
#' Created based on case.control.match.dev function, but adapted to flipped design
#' v2 switches to a bin-based approach
#' 
match.focal.v2 = function(in.data, in.seed, date.offset, age.bins,
                       removal = TRUE,
                       exclude.age = FALSE){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  strata.df = in.data
  strata.df$stratum = NA
  # Make each delta case a stratum
  n.strata = nrow(strata.df[strata.df$IS_FOCAL == 1, ])
  strata.df$stratum[strata.df$IS_FOCAL == 1] = seq(1, n.strata)
  
  # Create a flag for non-matches based on MF
  strata.df$MF_FLAG = 0
  
  # Set up a data frame of controls to draw from
  control.df = strata.df[strata.df$IS_FOCAL == 0, ]
  
  no.match.vec = c()
  
  # Loop through strata matching
  for (i in 1:n.strata){
    stratum = strata.df[strata.df$stratum == i, ]
    stratum = stratum[!is.na(stratum$stratum), ]
    # The below works, because there is only one entry in stratum. If there are multiple entries in stratum, this part of the code will fail
    # but if there are multiple entries in stratum, the code has already failed.
    case.id = stratum$Identifier
    case.age = stratum$Age
    case.sex = stratum$Sex
    case.region = stratum$Economic_region
    case.date = stratum$DATE_INDEX
    
    # If case sex is not known, DO NOT match on sex
    unknown.sex = 0
    if (case.sex == "0"){ unknown.sex = 1 }
    
    # Set a no-match indicator
    no.match = 0
    
    # Match based on Date
    date.match.index = control.df$DATE_INDEX >= (case.date - date.offset) & control.df$DATE_INDEX <= (case.date + date.offset)
    date.matches = control.df[date.match.index, ]
    
    if (nrow(date.matches) == 0){ no.match = 1 }
    
    # Match based on Region (# Kustin say municipality, geographic region, and sector, whatever that means. For us it means region.)
    if (no.match == 0){
      region.matches = date.matches[date.matches$Economic_region == case.region, ]
      if (nrow(region.matches) == 0){ no.match = 1 }
    }
    
    # Skip the age matching, but keep the variable name for down-stream code compatibility
    if (exclude.age == TRUE){
      if (no.match == 0){
        age.matches = region.matches
      }
      # Otherwise, match based on age.
    }else{
      # Match based on Age
      if (no.match == 0){
        
        # identify bin of focal individual
        bin.min = 0
        for (j in 1:length(age.bins)){
          bin.max = age.bins[j]
          
          # Stop the loop when max and min bound the case age.
          if (case.age >= bin.min & case.age <= bin.max){
            break
          }
          
          # otherwise increment the bins
          bin.min = bin.max + 1
        }
        
        age.match.index = region.matches$Age >= bin.min & region.matches$Age <= bin.max
        age.matches = region.matches[age.match.index, ] 
        
        if (nrow(age.matches) == 0){ no.match = 1 }    
      }
      
    }
    
    
    mf.flag = 1
    # Match based on sex (M vs. F)
    if (no.match == 0 & unknown.sex == 0){
      sex.matches = age.matches[age.matches$Sex == case.sex, ]
      if (nrow(sex.matches) > 0){ mf.flag = 0 }
    }
    
    # Remove the matching sex requirement if not a match. Still require a match on the above criteria
    if (mf.flag == 1 & no.match == 0){
      # Just select from the age match data set - renamed because the sex.matches name was already in use in the code below.
      sex.matches = age.matches
    }
    
    
    if (no.match == 0){
      # Select an identifier at random from remaining cases (if more than one case)
      if (nrow(sex.matches) == 1){
        control.id = sex.matches$Identifier[1]
      }else{
        control.id = sample(sex.matches$Identifier, 1) # This will not work if there is only one row.
      }
      # Assign the control entry to the matched stratum
      strata.df$stratum[strata.df$Identifier == control.id] = i
      strata.df$MF_FLAG[strata.df$Identifier == control.id] = mf.flag # Flag records where an exact sex match was not used
      
      if (removal == TRUE){
        # Remove the control from the control data base to ensure each case is matched to only one control.
        control.df = control.df[control.df$Identifier != control.id, ]
      }
    }else{
      no.match.vec = c(no.match.vec, case.id)
    }
  }
  
  # Discard any strata with no matches
  strat.table = table(strata.df$stratum) 
  keeps = names(strat.table)[strat.table > 1]
  out.strata = strata.df[strata.df$stratum %in% keeps, ]
  
  return(out.strata)
}


#' Match Cases (Development version)
#' 
#' Identify all records that match a given case - look at options for order of selection
#' This version looks at relaxing each criterion to give a full picture of options.
#' Not really sure how we'll SELECT from all those options, but... it at least gives us
#' the options on the table.
#' v3 differs from v2 in that #IS_FOCAL is used to define cases and controls.
#' NOTE: THIS IS NOT SUITABLE FOR THE MCNEMAR'S TEST, where case and control needs 
#' to be defined based on vaccination status
#' 
#' age.override added on 12/20/21 to provide a means of matching with more restricted values for younger ages.
#' 
case.control.options.v3 = function(in.data, in.seed, date.offset, age.offset,
                                   age.offset.override = NA){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  case.df = in.data[in.data$IS_FOCAL == 1, ] #CHANGED THIS LINE FROM .v2 to generalize it!
  case.df$EXACT_OPTIONS = NA # Initialize a column to store control identity
  case.df$NO_MF_OPTIONS = NA
  case.df$NO_AGE_OPTIONS = NA
  case.df$NO_REG_OPTIONS = NA
  case.df$NO_DATE_OPTION = NA
  case.df$N_EXACT = NA
  case.df$N_NO_MF = NA
  case.df$N_NO_AGE = NA
  case.df$N_NO_REG = NA
  case.df$N_NO_DATE = NA
  control.df = in.data[in.data$IS_FOCAL == 0, ] # ALSO CHANGED THIS LINE
  
  case.df = update.options(case.df, control.df, "exact", date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, "mf", date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, 'age', date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, 'region', date.offset, age.offset, age.offset.override)
  case.df = update.options(case.df, control.df, 'date', date.offset, age.offset, age.offset.override)
  
  return(case.df)
}

#' Match Cases (Development version)
#' 
#' Identify all records that match a given case - look at options for order of selection
#' This version looks at relaxing each criterion to give a full picture of options.
#' v3 differs from v2 in that #IS_FOCAL is used to define cases and controls.
#' v4 changes to use pre-set age bins rather than an age offset.
#' 
case.control.options.v4 = function(in.data, in.seed, date.offset, age.bins){
  # Set random number to ensure reproducibility
  set.seed(in.seed)
  
  case.df = in.data[in.data$IS_FOCAL == 1, ] #CHANGED THIS LINE FROM .v2 to generalize it!
  case.df$EXACT_OPTIONS = NA # Initialize a column to store control identity
  case.df$NO_MF_OPTIONS = NA
  case.df$NO_AGE_OPTIONS = NA
  case.df$NO_REG_OPTIONS = NA
  case.df$NO_DATE_OPTION = NA
  case.df$N_EXACT = NA
  case.df$N_NO_MF = NA
  case.df$N_NO_AGE = NA
  case.df$N_NO_REG = NA
  case.df$N_NO_DATE = NA
  control.df = in.data[in.data$IS_FOCAL == 0, ] # ALSO CHANGED THIS LINE
  
  case.df = update.options.v4(case.df, control.df, "exact", date.offset, age.bins)
  case.df = update.options.v4(case.df, control.df, "mf", date.offset, age.bins)
  case.df = update.options.v4(case.df, control.df, 'age', date.offset, age.bins)
  case.df = update.options.v4(case.df, control.df, 'region', date.offset, age.bins)
  case.df = update.options.v4(case.df, control.df, 'date', date.offset, age.bins)
  
  return(case.df)
}


#' Make a demographic summary table for the vaccinated/unvaccinated case/control design
#' 
#'
make.demographic.table = function(in.data.2, match.df, out.folder, figure.label,
                                  table.label){ 
  demographic.table = data.frame(Table.Label = NA, control.count = NA, control.percent = NA,
                                 vax.count = NA, vax.percent = NA)
  
  for.merge = in.data.2[ , c(1,3,4,5)]
  if (paste(colnames(for.merge), collapse = ' ') != "Identifier Age Sex Economic_region"){
    stop("Column names do not match required column order. Columns were indexed by order, so a change in column input order could cause probelms.")
  }
  colnames(for.merge) = c("CONTROL_ID", "CONTROL_AGE", "CONTROL_MF", "CONTROL_REGION")
  
  info.df = merge(match.df, for.merge, by = "CONTROL_ID")
  total.count = nrow(info.df)
  
  # Add age information
  age.limits = c(9,19, 29, 39, 49, 59, 69, 79, 89, 99, 109)
  if (max(match.df$Age) > max(age.limits)){ stop("a sample age exceeded the maximum age bin included")}
  
  lower.age = -1 # Set to be -1 to allow age 0 individuals to be appropriately binned
  for (upper.age in age.limits){
    
    new.label = sprintf("%s-%s", lower.age + 1, upper.age)
    
    # Get control statistics
    control.count = length(info.df$CONTROL_AGE[info.df$CONTROL_AGE<= upper.age & info.df$CONTROL_AGE > lower.age])
    control.percent = round(control.count / total.count,3) * 100
    
    # Get vaccinated statistics
    vax.count = length(info.df$Age[info.df$Age <= upper.age & info.df$Age > lower.age])
    vax.percent = round(vax.count / total.count, 3) * 100
    
    new.row = c(new.label, control.count, control.percent, vax.count, vax.percent)
    demographic.table = rbind(demographic.table, new.row)
    
    # Set the previous upper limit to be the lower limit.
    lower.age = upper.age
    #  count = count + 1
  }
  
  # Sex
  new.label = "Male"
  control.count = length(info.df$CONTROL_MF[info.df$CONTROL_MF == "M"])
  control.percent = round(control.count / total.count,3) * 100
  vax.count = length(info.df$Sex[info.df$Sex == "M"])
  vax.percent = round(vax.count / total.count, 3) * 100
  new.row = c(new.label, control.count, control.percent, vax.count, vax.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  new.label = "Female"
  control.count = length(info.df$CONTROL_MF[info.df$CONTROL_MF == "F"])
  control.percent = round(control.count / total.count,3) * 100
  vax.count = length(info.df$Sex[info.df$Sex == "F"])
  vax.percent = round(vax.count / total.count, 3) * 100
  new.row = c(new.label, control.count, control.percent, vax.count, vax.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  new.label = "Unknown"
  control.count = length(info.df$CONTROL_MF[info.df$CONTROL_MF != "F" & info.df$CONTROL_MF != "M"])
  control.percent = round(control.count / total.count,3) * 100
  vax.count = length(info.df$Sex[info.df$Sex != "F" & info.df$Sex != "M"])
  vax.percent = round(vax.count / total.count, 3) * 100
  new.row = c(new.label, control.count, control.percent, vax.count, vax.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  # Region
  regions = unique(info.df$Economic_region)
  for (region in regions){
    new.label = region
    control.count = length(info.df$CONTROL_REGION[info.df$CONTROL_REGION == region])
    control.percent = round(control.count / total.count,3) * 100
    vax.count = length(info.df$Economic_region[info.df$Economic_region == region])
    vax.percent = round(vax.count / total.count, 3) * 100
    new.row = c(new.label, control.count, control.percent, vax.count, vax.percent)
    demographic.table = rbind(demographic.table, new.row)
  }
  
  demographic.table = demographic.table[2:nrow(demographic.table), ]
  
  demo.out = sprintf("%s/DemographicTable%s%s.csv", out.folder, table.label, figure.label)
  write.table(demographic.table, file = demo.out, sep = ',', row.names = FALSE,
              col.names = TRUE, append = FALSE)
  
  return(info.df)
}


#' Make a demographic summary table for the flipped case/control design
#' 
#' @param strata.df The input data frame with records linked by stratum column
#' @param out.folder The folder in which to put the output data
#' @param figure.label An analysis-specific label
make.flipped.demographic.table = function(strata.df, out.folder, figure.label,
                                          focal.label, table.label, incl.booster = 0){
  demographic.table = data.frame(Table.Label = NA, OtherStrain.count = NA, OtherStrain.percent = NA)
  focal.count.field = sprintf("%s.count", focal.label)
  focal.percent.field = sprintf("%s.percent", focal.label)
  demographic.table[[focal.count.field]] = NA
  demographic.table[[focal.percent.field]] = NA

  focal.cases = strata.df[strata.df$IS_FOCAL == 1, ]
  other.cases = strata.df[strata.df$IS_FOCAL == 0, ]
  total.count = nrow(focal.cases) # should equal number of rows in other.cases if doing 1:1 match. Otherwise, recode
  if (nrow(focal.cases) != nrow(other.cases)){
    stop("number of rows in the case control should match or this function will miscalculate the percentages")
  }
  
  # Add age information
  #age.limits = c(9,19, 29, 39, 49, 59, 69, 79, 89, 99, 109)
  age.limits = c(4,11,17,29,49,69,89,109)
  if (max(strata.df$Age) > max(age.limits)){ stop("a sample age exceeded the maximum age bin included")}
  
  lower.age = -1 # Set to be -1 to allow age 0 individuals to be appropriately binned
  for (upper.age in age.limits){
    
    new.label = sprintf("%s-%s", lower.age + 1, upper.age)
    if (new.label == "10-19"){ new.label = "'10-19" } # Add a leading ' so this field is interpreted correctly in Excel
    
    # Get other strain statistics
    other.count = length(other.cases$Age[other.cases$Age<= upper.age & other.cases$Age > lower.age])
    other.percent = round(other.count / total.count,3) * 100
    
    # Get vaccinated statistics
    focal.count = length(focal.cases$Age[focal.cases$Age <= upper.age & focal.cases$Age > lower.age])
    focal.percent = round(focal.count / total.count, 3) * 100
    
    new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
    demographic.table = rbind(demographic.table, new.row)
    
    # Set the previous upper limit to be the lower limit.
    lower.age = upper.age
    #  count = count + 1
  }
  
  # Sex
  new.label = "Male"
  other.count = length(other.cases$Sex[other.cases$Sex == "M"])
  other.percent = round(other.count / total.count,3) * 100
  focal.count = length(focal.cases$Sex[focal.cases$Sex == "M"])
  focal.percent = round(focal.count / total.count, 3) * 100
  new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  new.label = "Female"
  other.count = length(other.cases$Sex[other.cases$Sex == "F"])
  other.percent = round(other.count / total.count,3) * 100
  focal.count = length(focal.cases$Sex[focal.cases$Sex == "F"])
  focal.percent = round(focal.count / total.count, 3) * 100
  new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  new.label = "Unknown"
  other.count = length(other.cases$Sex[other.cases$Sex != "F" & other.cases$Sex != "M"])
  other.percent = round(other.count / total.count,3) * 100
  focal.count = length(focal.cases$Sex[focal.cases$Sex != "F" & focal.cases$Sex != "M"])
  focal.percent = round(focal.count / total.count, 3) * 100
  new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  # Region
  regions = sort(unique(focal.cases$Economic_region)) # Should be same for other.cases object - matched on region!
  for (region in regions){
    new.label = region
    other.count = length(other.cases$Economic_region[other.cases$Economic_region == region])
    other.percent = round(other.count / total.count,3) * 100
    focal.count = length(focal.cases$Economic_region[focal.cases$Economic_region == region])
    focal.percent = round(focal.count / total.count, 3) * 100
    new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
    demographic.table = rbind(demographic.table, new.row)
  }
  
  # Vaccination status (not matched) (now covered by 'None')
  #new.label = "Unvaccinated"
  #other.count = nrow(other.cases[is.na(other.cases$Days.between.collection.date.and.vaccination.complete.date), ])
  #other.percent = round(other.count / total.count,3) * 100
  #focal.count = nrow(focal.cases[is.na(focal.cases$Days.between.collection.date.and.vaccination.complete.date), ])
  #focal.percent = round(focal.count / total.count, 3) * 100
  #new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  #demographic.table = rbind(demographic.table, new.row)
  
  vaccine.types = unique(strata.df$vaccine_type)
  for (vaccine.type in vaccine.types){
    # Add vaccine type breakdown
    new.label = sprintf("%s Vaccine", vaccine.type)
    other.count = nrow(other.cases[other.cases$vaccine_type == vaccine.type, ])
    other.percent = round(other.count / total.count,3) * 100
    focal.count = nrow(focal.cases[focal.cases$vaccine_type == vaccine.type, ])
    focal.percent = round(focal.count / total.count, 3) * 100
    new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
    demographic.table = rbind(demographic.table, new.row)
    
  }
  
  new.label = "Vaccinated <90 days"
  other.count = nrow(other.cases[!is.na(other.cases$Days.between.collection.date.and.vaccination.complete.date) &
                                   other.cases$Days.between.collection.date.and.vaccination.complete.date < 90, ])
  other.percent = round(other.count / total.count,3) * 100
  focal.count = nrow(focal.cases[!is.na(focal.cases$Days.between.collection.date.and.vaccination.complete.date) &
                                   focal.cases$Days.between.collection.date.and.vaccination.complete.date < 90, ])
  focal.percent = round(focal.count / total.count, 3) * 100
  new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  new.label = "Vaccinated 90+ days"
  other.count = nrow(other.cases[!is.na(other.cases$Days.between.collection.date.and.vaccination.complete.date) &
                                   other.cases$Days.between.collection.date.and.vaccination.complete.date >= 90, ])
  other.percent = round(other.count / total.count,3) * 100
  focal.count = nrow(focal.cases[!is.na(focal.cases$Days.between.collection.date.and.vaccination.complete.date) &
                                   focal.cases$Days.between.collection.date.and.vaccination.complete.date >= 90, ])
  focal.percent = round(focal.count / total.count, 3) * 100
  new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  # Add booster status breakdown

  booster.types = unique(strata.df$booster_type)
  for (booster.type in booster.types){
    # Add Booster type breakdown
    new.label = sprintf("%s Booster", booster.type)
    boost.indicator = 1
    if (booster.type == "None"){
      boost.indicator = 0
    }
    
    other.count = nrow(other.cases[other.cases$booster_type == booster.type & other.cases$IS.BOOST == boost.indicator, ])
    other.percent = round(other.count / total.count,3) * 100
    focal.count = nrow(focal.cases[focal.cases$booster_type == booster.type & focal.cases$IS.BOOST == boost.indicator, ])
    focal.percent = round(focal.count / total.count, 3) * 100
    new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
    demographic.table = rbind(demographic.table, new.row)
  }
  
  new.label = "Boosted <90 days"
  other.count = nrow(other.cases[!is.na(other.cases$Days.between.collection.date.and.booster.date) &
                                   other.cases$Days.between.collection.date.and.booster.date < 90, ])
  other.percent = round(other.count / total.count,3) * 100
  focal.count = nrow(focal.cases[!is.na(focal.cases$Days.between.collection.date.and.booster.date) &
                                   focal.cases$Days.between.collection.date.and.booster.date < 90, ])
  focal.percent = round(focal.count / total.count, 3) * 100
  new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  new.label = "Boosted 90+ days"
  other.count = nrow(other.cases[!is.na(other.cases$Days.between.collection.date.and.booster.date) &
                                   other.cases$Days.between.collection.date.and.booster.date >= 90, ])
  other.percent = round(other.count / total.count,3) * 100
  focal.count = nrow(focal.cases[!is.na(focal.cases$Days.between.collection.date.and.booster.date) &
                                   focal.cases$Days.between.collection.date.and.booster.date >= 90, ])
  focal.percent = round(focal.count / total.count, 3) * 100
  new.row = c(new.label, other.count, other.percent, focal.count, focal.percent)
  demographic.table = rbind(demographic.table, new.row)
  
  
  demographic.table = demographic.table[2:nrow(demographic.table), ]
  
  demo.out = sprintf("%s/DemographicTable%s%s.csv", out.folder, table.label, figure.label)
  write.table(demographic.table, file = demo.out, sep = ',', row.names = FALSE,
              col.names = TRUE, append = FALSE)
  
}

#' Make a plot showing the Odds Ratio for the model with all data, and for each stratum
#' deleted.
#' 
make.leverage.plot = function(model.object, formula.object, out.folder,
                              figure.label, table.label){
  # Define the analysis object
  my.result = summary(model.object)
  n.parameters = nrow(my.result$coefficients)
  
  OR.df = data.frame(Coefficient = NA, OR = NA, OR.lower = NA, OR.upper = NA, stratum.removed = NA)
  
  # Take estimate from full analysis
  OR.matrix = my.result$conf.int
  
  for (j in 1:n.parameters){
    this.record = c(row.names(OR.matrix)[j], OR.matrix[j, 1], OR.matrix[j,3], OR.matrix[j,4], 'NONE')
    OR.df = rbind(OR.df, this.record)
  }
  
  strata = unique(strata.df$stratum)
  # Loop through and repeat analysis deleting one stratum at a time
  for (k in 1:length(strata)){
    this.strata = strata[k]
    delete.strata = strata.df[strata.df$stratum != this.strata, ]
    my.analysis = survival::clogit(formula.object, data=delete.strata)
    my.result = summary(my.analysis)
    
    OR.matrix = my.result$conf.int
    
    for (j in 1:n.parameters){
      this.record = c(row.names(OR.matrix)[j], OR.matrix[j, 1], OR.matrix[j,3], OR.matrix[j,4], this.strata)
      OR.df = rbind(OR.df, this.record)
    }
  }
  
  OR.df = OR.df[2:nrow(OR.df), ]
  
  # Make sure numeric columns are numeric
  OR.df$OR = as.numeric(OR.df$OR)
  OR.df$OR.lower = as.numeric(OR.df$OR.lower)
  OR.df$OR.upper = as.numeric(OR.df$OR.upper)
  
  # Save strata ID's to a table to look up if any strata are of particular interest
  OR.file = sprintf("%s/OddsRatioLeverage%s%s.csv", out.folder, table.label, figure.label)
  write.table(OR.df, file = OR.file, sep = ',', row.names = FALSE, col.names = TRUE,
              append = FALSE)
  
  
  # Plot change based on stratum for each coefficient. Sort based on lowest to highest
  
  for (coef.name in unique(OR.df$Coefficient)){
    
    # Limit data to just a single coefficient
    OR.df.coef = OR.df[OR.df$Coefficient == coef.name, ]
    # Create objects for the row with all data, and for all other rows
    first = OR.df.coef[OR.df.coef$stratum == "NONE", ]
    remaining = OR.df.coef[OR.df.coef$stratum != "NONE", ]
    
    inf.index = grep("Inf", OR.df.coef$OR.upper)
    if (length(inf.index) > 0){
      warning("At least one model did not converge. Infinite parameter estimate was removed for plotting purposes")
      OR.df.coef = OR.df.coef[-inf.index, ]
    }
    
    # Get unique values and count of each unique value to avoid massive plots with repetitive data
    unique.values = table(round(remaining$OR, 4)) # How do we do CI for these? Esp. as it could differ (but doesn't in this example)
    # Added rounding, because with continuous input variables, one gets a truly awful plot with these inputs
    w.adjust = (length(unique.values) + 1) / 4 # adjust width of the plot depending on how many values are being plotted.

    influence.plot = sprintf("%s/InfluencePlot_BestModel_%s%s%s.tif", out.folder,
                             coef.name, table.label, figure.label)
    tiff(filename = influence.plot, width = 1200 * w.adjust, height = 1200, compression = c('lzw'),
         res = 300)
    
    # Set x and y limits
    x.limits = c(1, length(unique.values) + 1)
    y.limits = as.numeric(c(min(OR.df.coef$OR.lower), max(OR.df.coef$OR.upper)))
    bar.length = 0.03
    
    # Plot the 'NONE' row first - with no deletions
    plot(1, first$OR, xlim = x.limits, ylim = y.limits, xaxt = 'n', xlab = "Stratum", ylab = "Odds Ratio and 95% CI")
    arrows(1, first$OR, 1, first$OR.upper, length = bar.length, angle = 90)
    arrows(1, first$OR, 1, first$OR.lower, length = bar.length, angle = 90)
    
    # Sort the remaining data in order of OR estimate (do we need increasing/decreasing depending on if it is positive or negative?)
    #remaining = remaining[order(remaining$OR, decreasing = FALSE), ]
    #remaining$PLOT.ORDER = seq(1,nrow(remaining)) # Add a plot order for plotting. Note that it starts at 1, because 0 was used for the full analysis
    
    # Plot remaining data
    par(new = TRUE)
    #plot(remaining$PLOT.ORDER, remaining$OR, xlim = x.limits, ylim = y.limits,
    #     xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    #arrows(remaining$PLOT.ORDER, remaining$OR, remaining$PLOT.ORDER, remaining$OR.upper, length = bar.length, angle = 90)
    #arrows(remaining$PLOT.ORDER, remaining$OR, remaining$PLOT.ORDER, remaining$OR.lower, length = bar.length, angle = 90)
    #axis(1, at = c(-1,remaining$PLOT.ORDER), labels = c("Full", remaining$stratum.removed), las = 2, cex.axis = 0.8)
    r.index = seq(2, length(unique.values) + 1)
    ORs = as.numeric(names(unique.values))
    OR.uppers = sapply(ORs, find.OR.CI, remaining, 'upper')
    OR.lowers = sapply(ORs, find.OR.CI, remaining, 'lower')
    plot(r.index, ORs,
         xlim = x.limits, ylim = y.limits, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    arrows(r.index, ORs, r.index, OR.uppers, length = bar.length, angle = 90)
    arrows(r.index, ORs, r.index, OR.lowers, length = bar.length, angle = 90)
    
    axis(1, at = c(1,r.index), labels = c("All", sprintf("(%s)", unique.values)), las = 2, cex.axis = 0.8)
    
    # Add line for OR of 1
    segments(0,1,x.limits[2],1, lty = 2)
    
    dev.off()
  }
  
  
}


#' Helper function to pull back an upper odds ratio value
#' 
find.OR.CI = function(x, remaining, direction){
  test = grep(x, round(remaining$OR, 4)) # Add rounding to match rounding applied to unique.values
  
  # If no results, give NA (this should not happen)
  if (length(test) == 0){
    out = NA
  }
  
  # Assume all results are the same
  if (direction == 'upper'){
    out = remaining$OR.upper[test[1]] # pull the first value, assume all upper CI's are equal for the same point estimate
  }
  if (direction == 'lower'){
    out = remaining$OR.lower[test[1]] # Same assumption as above
  }
  return(out)
}

#' Convert Demographic table information into a bar chart
#' 
demographic.table.to.bar.chart = function(demographic.table.file, rows.start, rows.end){
  table.info = read.csv(demographic.table.file)
  
  table.info = table.info[rows.start:rows.end, ]

  # Reformat for side-by-side bar plot
  plot.info = matrix(c(table.info[ ,c(2)], table.info[,c(4)]), ncol = 1)
  plot.info = data.frame(COUNT = plot.info, GROUP = "firebrick3")
  plot.info$GROUP[rows.start:rows.end] = 'coral' # Set the initial rows to be group 0
  diff = rows.end - rows.start
  stop.point = (rows.end - rows.start) * 2 + 1
  plot.info$SORT[rows.start:rows.end] = seq(1,stop.point, 2)
  plot.info$SORT[(rows.end + 1):(rows.end + diff + 1)] = seq(2,stop.point + 1, 2)
  length(plot.info$SORT[(rows.end + 1):(rows.end + diff + 1)])
  plot.info = plot.info[order(plot.info$SORT), ]    
  
  x = barplot(plot.info$COUNT, beside = TRUE, col = plot.info$GROUP, ylab = "Count")
  legend(18, 50, legend = c("Other Strain", "Omicron"), fill = c('coral', 'firebrick3'))
  x.alt = x[seq(1,length(x),2)] # Get every other one, as labels will be the same for the paired bars
  axis.labels = table.info$Table.Label
  axis.labels[2] = '10-19' # Remove leading "'", which was added to improve table behavior in Excel
  axis(1, at = x.alt, labels = axis.labels, las = 2)
  mtext("Age", side = 1, line = 4)
  
}

