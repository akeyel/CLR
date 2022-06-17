# Goal is to make an overview figure showing the strain dominance for the
# sequence data, and indicate the key mixing periods that define the three
# analyses used in the draft manuscript

# Author: A.C. Keyel <alexander.keyel@health.ny.gov>

# Created 2022-01-13 from existing code from DOH_COVID_final.R
# (that script created 2021-11-29.)

# Update: 2022-03-14
# Changed to be a single panel with space for a map
# Color code to match the map
# Change to a stacked bar graph by economic region

# Read in Delta data set
strata.file = "Results/Delta/01_ConditionalLogistic_stratadf_delta.csv"
if (!file.exists(strata.file)){
  stop("THIS ANALYSIS REQUIRES THE MATCHED DATA SET CREATED IN
       06_ConditionalLogistic_Delta.R. If you have created
       that file and are still seeing this message, please check the file paths.")
}

strata.df = read.csv(strata.file)
strata.df$IS_DELTA = strata.df$IS_FOCAL
strata.df$IS_OMICRON = 0

prelim.label = "" # "_Prelim"
if (prelim.label != ""){
  warning("Using preliminary Omicron data")
}

omicron.file = sprintf("Results/Omicron%s/01_ConditionalLogistic_stratadf_omicron%s.csv", prelim.label, tolower(prelim.label))

# Read Omicron data
if (!file.exists(omicron.file)){
  stop("THIS ANALYSIS REQUIRES THE MATCHED DATA SET CREATED IN
       01_ConditionalLogistic_Omicron.R. If you have created
       that file and are still seeing this message, please check the file paths.")
}

omicron.df = read.csv(omicron.file)
omicron.df$IS_OMICRON = omicron.df$IS_FOCAL
omicron.df$IS_DELTA = 0

# Combine delta and omicron data frames for ease of processing
delta.for.merge = strata.df[strata.df$IS_DELTA == 1, c(1,2,5,31,42,43)]
omicron.for.merge = omicron.df[omicron.df$IS_OMICRON == 1, c(1,2,5, 31,43,42)]

# Check that correct fields were extracted
target.names = c("Identifier","collection_date", "Economic_region","DATE_INDEX","IS_DELTA", "IS_OMICRON")
if (paste(names(delta.for.merge), collapse = " ") !=  paste(target.names, collapse = " ")){
  stop("Please fix the fields extracted")
}

cases = rbind(delta.for.merge, omicron.for.merge)
  
# Make the plot

# Make a stacked bar plot by economic region
# No indicator for Omicron/Delta this time around - the mixing periods do not overlap for the studies
stacked.data = data.frame(WEEK = NA, ECONOMIC_REGION = NA, CASE_CONTROL_PAIRS = NA)

regions = unique(cases$Economic_region)
first.day = 1
week.start = first.day
week.number = 1


while(week.start < (max(cases$DATE_INDEX) + 6)){
  week.end = week.start + 7
  
  # Add a record for all case-control pairs for each economic region
  for (region in regions){
    # Subset to this region
    region.cases = cases[cases$Economic_region == region, ]
    
    # Count case-control pairs
    n.pairs.delta = sum(region.cases$IS_DELTA[region.cases$DATE_INDEX >= week.start & region.cases$DATE_INDEX < week.end])
    n.pairs.omicron = sum(region.cases$IS_OMICRON[region.cases$DATE_INDEX >= week.start & region.cases$DATE_INDEX < week.end])
    n.pairs = n.pairs.delta + n.pairs.omicron
    
    new.record = c(week.number, region, n.pairs)
    stacked.data = rbind(stacked.data, new.record)
  }
  
  # Update week.start and week.number for the next week.
  week.start = week.end
  week.number = week.number + 1
}

# Drop leading NA
stacked.data = stacked.data[2:nrow(stacked.data), ]

# Reformat stacked.data data frame to format needed for bar plot command
n.weeks = 57 # #1-57, otherwise the indexing doesn't work, but only want 10 - 57 inclusive
n.regions = length(regions)
zeros = rep(0, n.weeks * n.regions)
new.format = matrix(zeros, ncol = n.weeks, nrow = n.regions)

for (i in 1:nrow(stacked.data)){
  col.index = as.numeric(stacked.data$WEEK[i])
  this.region = stacked.data$ECONOMIC_REGION[i]
  row.index = grep(this.region, regions)
  this.value = as.numeric(stacked.data$CASE_CONTROL_PAIRS[i])
  
  new.format[row.index, col.index] = this.value
}

distribution.figure = sprintf("Results/Temporal_Overview_base%s_v2.tif", prelim.label)
tiff(file = distribution.figure, res = 300, compression = c('lzw'),
     height = 2000, width = 3000)

par(mar = c(6,6,2,0))
y.limits = c(0.60) 
x.xlim = c(10,57) 

# Add date windows for VOC/DELTA/OMICRON
delta.min = min(cases$DATE_INDEX[cases$IS_DELTA == 1])
delta.min.date = cases$collection_date[cases$DATE_INDEX == delta.min][1]
delta.max = max(cases$DATE_INDEX[cases$IS_DELTA == 1])
delta.max.date = cases$collection_date[cases$DATE_INDEX == delta.max][1]
omicron.min =  min(cases$DATE_INDEX[cases$IS_OMICRON == 1])
omicron.max = max(cases$DATE_INDEX[cases$IS_OMICRON == 1])
omicron.min.date = cases$collection_date[cases$DATE_INDEX == omicron.min][1]
omicron.max.date = cases$collection_date[cases$DATE_INDEX == omicron.max][1]

delta.label = sprintf("    Delta: %s - %s", delta.min.date, delta.max.date)
omicron.label = sprintf("Omicron: %s - %s", omicron.min.date, omicron.max.date)

x.label = sprintf("Week (%s January 2021 = Week 1)", first.day)

Lon.Isl = hsv(20/360,0.42,0.92) #20
Mid.Hud = hsv(230/360,0.41,0.89) # Changed to blue
Fin.Lak = hsv(24/360,0.33,0.94)
Cen.NY =  hsv(27/360,0.24,0.97)
Wes.NY =  hsv(0/360,0.95,0.77)
Cap.Reg = hsv(31/360,0.16,1.00)
Sou.Tie = hsv(230/360,0.86,0.79) # Changed to blue
Nor.Cou = hsv(6/360,0.77,0.82)
Moh.Val = hsv(13/360,0.59,0.87)
NYC = hsv(10/360,0.68,0.84)


# Use background shading to indicate Delta vs. Omicron mixing periods
region.colors = c(Lon.Isl, Mid.Hud, Fin.Lak, Cen.NY,
                  Wes.NY, Cap.Reg, Sou.Tie, Nor.Cou,
                  Moh.Val, NYC) # This works as expected, but will require checking that the correct regions are assigned the correct colors

# But they should be in the same order
x = barplot(new.format,
            col = region.colors,
            xlim = c(13.9,67.9),
            xpd = FALSE,
            ylab = "Case-Control Pairs",
            cex.axis = 2,
            cex.lab = 2) #xlab = "Week",

# Add hatching
region.colors[1] = 'black'
region.colors[4] = 'black'
region.colors[9] = 'black'
region.density = c(30,-1,-1,30,-1,-1,-1,-1,30,-1)
region.angle = c(90,NA,NA,0,NA,NA,NA,NA,45,NA)
par(new = TRUE)
x = barplot(new.format,
            col = region.colors,
            density = region.density,
            angle = region.angle,
            xlim = c(13.9,67.9),
            xpd = FALSE,
            ylab = "Case-Control Pairs",
            cex.axis = 2,
            cex.lab = 2) #xlab = "Week",



# xlimits were manually adjusted based on the barplots x coordinates.
# Add shaded polygons for Delta and Omicron
# Or... maybe just put a dividing bar in?
segments(56,0,56,70, lty = 2, lwd = 4)
# Labels with week
#axis(side = 1, at = x, labels = seq(1,ncol(new.format)))

# Labels with Date as an alternative
#as.Date(329, origin = "2020-12-31")# Avoid 0-based day of year
#as.Date(1, origin = "2020-12-31")# Avoid 0-based day of year
# Convert weeks to dates
weeks = seq(1,ncol(new.format))
doys = weeks * 7 - 7 + 1
# dates = sapply(doys, as.Date, origin = "2020-12-31")
# This gives me 18628...
dates = c()
for (doy in doys){
  this.date = as.character(as.Date(doy, origin = '2020-12-31')) # This is character
  # Drop year (should be obvious with a reference in the text)
  date.parts = strsplit(this.date, split = "-")[[1]]
  this.date2 = paste(c(date.parts[2], date.parts[3]), collapse = "-")
  dates = c(dates, this.date2)
}
axis(side = 1, at = x, labels = dates, las = 2,
     cex.axis = 2)
dev.off()


# Combine with economic region map inset (done outside of R)

