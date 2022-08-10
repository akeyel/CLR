# Goal: Quantify rates of vaccination by age for the study period to put
# the study results into the context of vaccination rates.

# Created 2022-07-19
# Updated 2022-08-09

# Author: Alexander "Sasha" Keyel <alexander.keyel@health.ny.gov

# Download CDC Vaccination rate data
# https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc
# Data set downloaded as .csv format on 19 July 2022

# Read in CDC Data
vax.data = read.csv("C:/hri/DOH_COVID/Data/CDC/COVID-19_Vaccinations_in_the_United_States_Jurisdiction.csv")
View(vax.data)

# Subset to NYS
# NOTE: This data set will not allow us to go finer than state-level data
# The data set includes NYC, which was not included in the main analysis
# NYC data is handled separately from the rest of the state.

ny.data = vax.data[vax.data$Location == "NY", ]

# Test code to ensure the make_index function is working properly
#make_index(2020,12,28) # 28
#make_index(2021,1,1)   # 32
#make_index(2021,2,28)  # 90
#make_index(2021,12,31) # 396
#make_index(2022,1,1)   # 397

make_index = function(year, month, day){
  # Set 12/1/2020 as 1
  # 1/1/2021 is 32
  if (year == 2020){ year.adj = 0 }
  if (year == 2021){ year.adj = 31}
  if (year == 2022){ year.adj = 365 + 31}
  month.days =     c(31, 28, 31, 30,   31,  30,  31,  31,  30,  31,  30,  31)
  cum.month.days = c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
  
  month.adj = cum.month.days[as.numeric(month)]
  if (year == 2020){ month.adj = 0   }
  day.adj = day
  index = year.adj + month.adj + day.adj
  
  return(index)
}

# Try an approach with continuous data rather than selected dates
ny.data$month = sapply(ny.data$Date, dfmip::splitter, "/", 1, 0)
ny.data$day = sapply(ny.data$Date, dfmip::splitter, "/", 2, 0)
ny.data$year = sapply(ny.data$Date, dfmip::splitter, "/", 3, 0)
ny.data$plot_index = mapply(make_index, ny.data$year, ny.data$month, ny.data$day)
ny.data2 = ny.data[ny.data$plot_index > 94 & ny.data$plot_index < 429, ] # 94 was where the first complete series began 428 is Feb. 1 and prevents too much extra space being used outside the study period

# Delta emergence: 3/19/21, 8/15/21;  Omicron Emergence: 11/28/21, 1/24/22 
key.years = c(2021, 2022, 2021, 2021)
key.months = c(11, 1, 3, 8)
key.days = c(28, 24, 19, 15)
divider.dates = mapply(make_index, key.years, key.months, key.days)
divider.labels = c("11/28/21", "1/24/22", "3/19/21", "8/15/21")

# Get index locations for months, to have human-readable x-axis labels
index.months = c(seq(3,12), 1, 2)
index.years = c(rep(2021, length(index.months) - 2), 2022, 2022)
index.days = rep(1, length(index.months))
index.dates = mapply(make_index, index.years, index.months, index.days)

# Make the plot of vaccinations over time by vaccine type
tiff("C:/hri/DOH_COVID/Figures/Vaccination_over_time.tif",
     height = 2400, width = 4800, res = 600, compression = c('lzw'))

par(mfrow = c(1,3))
par(oma = c(4,6,2,2))
par(mar = c(0,0,0,1))

y.limits = c(0, max(ny.data2$Series_Complete_Pfizer_18Plus - ny.data2$Series_Complete_Pfizer_65Plus, na.rm = TRUE) / 1000000)
x.limits = c(94, 429)
line.width = 3

scale.factor = 1.5

## Pfizer
# Initialize the plot window
plot(1,1,
     ylim = y.limits, col = 'white', xpd = NA,
     xlim = x.limits, xlab = "Month",
     ylab = "Fully Vaccinated\n(millions)",
     xaxt = 'n',
     cex = 0.01,
     cex.axis = scale.factor, cex.lab = scale.factor)

# Add the background polygons
polygon(x = c(divider.dates[1], divider.dates[2], divider.dates[2], divider.dates[1]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')
polygon(x = c(divider.dates[3], divider.dates[4], divider.dates[4], divider.dates[3]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')

# Add the data series and the x-axis
lines(ny.data2$plot_index, ny.data2$Series_Complete_Pfizer_65Plus / 1000000,
      lty = 1, lwd = line.width)
axis(side = 1, at = index.dates, labels = index.months,
     cex.lab = scale.factor, cex.axis = scale.factor)

remainder = ny.data2$Series_Complete_Pfizer_18Plus - ny.data2$Series_Complete_Pfizer_65Plus
remainder[remainder < 0 | is.na(remainder)] = 0
lines(ny.data2$plot_index, remainder  / 1000000,
      lty = 2, lwd = line.width)

remainder = ny.data2$Series_Complete_Pfizer_12Plus - ny.data2$Series_Complete_Pfizer_18Plus
remainder[remainder < 0 | is.na(remainder)] = 0
lines(ny.data2$plot_index, remainder  / 1000000,
      lty = 3, lwd = line.width)

remainder = ny.data2$Series_Complete_Pfizer_5Plus - ny.data2$Series_Complete_Pfizer_12Plus
remainder[remainder < 0 | is.na(remainder)] = 0
lines(ny.data2$plot_index, remainder  / 1000000,
      lty = 4, lwd = line.width)

# Add a label to the plot
text(180, 5.6, "Pfizer", cex = 1.5, xpd = NA)

## Moderna
# Initialize the plot
plot(480,1,
     ylim = y.limits, col = 'white', pch = 2, ylab = "", yaxt = 'n', xpd = NA,
     xlim = x.limits, xlab = "Month", xaxt = 'n',
     cex.lab = scale.factor, cex.axis = scale.factor)
axis(side = 1, at = index.dates, labels = index.months,
     cex.lab = scale.factor, cex.axis = scale.factor)

# Add the background polygons
polygon(x = c(divider.dates[1], divider.dates[2], divider.dates[2], divider.dates[1]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')
polygon(x = c(divider.dates[3], divider.dates[4], divider.dates[4], divider.dates[3]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')

# Plot the data series (only 2 because <18 were not eligible for this vaccine)
lines(ny.data2$plot_index, ny.data2$Series_Complete_Moderna_65Plus / 1000000,
      lty = 1, lwd = line.width)

remainder = ny.data2$Series_Complete_Moderna_18Plus - ny.data2$Series_Complete_Moderna_65Plus
remainder[remainder < 0 | is.na(remainder)] = 0
lines(ny.data2$plot_index, remainder / 1000000,
      lty = 2, lwd = line.width)

# Add plot label
text(180, 5.6, "Moderna", cex = 1.5, xpd = NA)


## Janssen
# Initialize plot window
plot(480,1,
     ylim = y.limits, col = 'white', pch = 3, ylab = "", yaxt = 'n', xpd = NA,
     xlim = x.limits, xlab = "Month", xaxt = 'n',
     cex.lab = scale.factor, cex.axis = scale.factor)
axis(side = 1, at = index.dates, labels = index.months,
     cex.lab = scale.factor, cex.axis = scale.factor)

# Add background polygons
polygon(x = c(divider.dates[1], divider.dates[2], divider.dates[2], divider.dates[1]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')
polygon(x = c(divider.dates[3], divider.dates[4], divider.dates[4], divider.dates[3]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')

# Plot data
lines(ny.data2$plot_index, ny.data2$Series_Complete_Janssen_65Plus / 1000000,
      lty = 1, lwd = line.width)

remainder = ny.data2$Series_Complete_Janssen_18Plus - ny.data2$Series_Complete_Janssen_65Plus
remainder[remainder < 0 | is.na(remainder)] = 0
lines(ny.data2$plot_index, remainder / 1000000,
      lty = 2, lwd = line.width)

# Add a legend for all three panels
legend(divider.dates[3], par('usr')[4] - 1, legend = c("5 - 11", "12 - 17", "18 - 64", "65+"),
       lty = seq(4,1,-1), lwd = 2)

# Add Janssen plot label
text(180, 5.6, "Janssen", cex = 1.5, xpd = NA)

dev.off()

#### Create plot for change in booster over time ####

tiff("C:/hri/DOH_COVID/Figures/Boosters_over_time.tif", height = 3600, width = 2400,
     res = 600, compression = c('lzw'))

par(mfrow = c(2,1))
par(oma = c(0,0,0,0))
par(mar = c(4,5,1,1))

#y.limits = c(0, max(ny.data2$Additional_Doses_18Plus, na.rm = TRUE)/1000000)
y.limits = c(0,4.1)
x.limits = c(305,428)

## Boosters by age
plot(480,1,
     ylim = y.limits, col = 'white', pch = 3, ylab = "# Boosted (millions)",
     xlim = x.limits, xlab = "Month", xaxt = 'n',
     cex.lab = scale.factor, cex.axis = scale.factor)
axis(side = 1, at = index.dates, labels = index.months,
     cex.lab = scale.factor, cex.axis = scale.factor)

# Add background polygons
polygon(x = c(divider.dates[1], divider.dates[2], divider.dates[2], divider.dates[1]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')
polygon(x = c(divider.dates[3], divider.dates[4], divider.dates[4], divider.dates[3]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')

# 18 - 64
lines(ny.data2$plot_index, (ny.data2$Additional_Doses_18Plus - ny.data2$Additional_Doses_65Plus) / 1000000,
      lty = 2, lwd = 3)

# 65+
lines(ny.data2$plot_index, ny.data2$Additional_Doses_65Plus / 1000000,
      lty = 1, lwd = 3)

legend(par('usr')[1], par('usr')[4], legend = c("18-64", "65+"), lty = c(2,1), lwd = 3)

## Boosters by type

plot(400,1,
     ylim = y.limits, col = 'white', pch = 3, ylab = "# Boosted (millions)", xpd = NA,
     xlim = x.limits, xlab = "Month", xaxt = 'n',
     cex.lab = scale.factor, cex.axis = scale.factor)
axis(side = 1, at = index.dates, labels = index.months,
     cex.lab = scale.factor, cex.axis = scale.factor)

# Add background polygons
polygon(x = c(divider.dates[1], divider.dates[2], divider.dates[2], divider.dates[1]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')
polygon(x = c(divider.dates[3], divider.dates[4], divider.dates[4], divider.dates[3]),
        y = c(par('usr')[3], par('usr')[3], par('usr')[4], par('usr')[4]),
        col = 'gray85')

# Pfizer
points(ny.data2$plot_index, (ny.data2$Additional_Doses_Pfizer) / 1000000,
       pch = 3)

# Moderna
points(ny.data2$plot_index, ny.data2$Additional_Doses_Moderna / 1000000,
       pch = 2)

# Janssen
points(ny.data2$plot_index, ny.data2$Additional_Doses_Janssen / 1000000,
       pch = 1)

legend(par('usr')[1], par('usr')[4], legend = c("Pfizer", "Moderna", "Janssen"), pch = c(3,2,1))

dev.off()

