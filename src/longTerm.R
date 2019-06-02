### IMPORTS ####
library(gdata)
library(xts)
library(dplyr)
library(forecast)

# ------------------------- load the data-------------------------------
# Set path
homePath = "/Users/myfiles/Documents/EPFL/M_II/TSE/SwissEnergyGrid/"
dataPath = paste(homePath, "data/", sep="")
outPath = paste(homePath, "res/Robjects/", sep="")

energyTS = readRDS(paste(outPath, "energyTS.RDS", sep=""))


# ------------------------- remove leap days -------------------------------
# remove leap days
remove_leap_days <- function(ts){
  idx = index(ts)
  select <- (format(idx, "%d")=="29" & format(idx, "%m")=="02") | as.numeric(format(idx, "%Y"))<2009 | as.numeric(format(idx, "%Y"))>=2019
  result = ts[!select]
}
consumption = remove_leap_days(energyTS$TotConsumption)


# ------------------------- obtain peak values -----------------------------
dailymax = apply.daily(consumption, max)
plot(dailymax)

# ------------------------- fit sinus and weekend --------------------------
days = index(dailymax)
t <- 1:length(dailymax)
per = 365

weekend_days = c(0, 6)
workday = as.numeric(!days$wday%in%weekend_days)
weday = as.numeric(days$wday%in%weekend_days)
regr <- lm(dailymax ~ sin(2*pi/per*t)+cos(2*pi/per*t) + 
                sin(2*pi/per*2*t)+cos(2*pi/per*2*t) +
                sin(2*pi/per*3*t)+cos(2*pi/per*3*t) +
                sin(2*pi/per*4*t)+cos(2*pi/per*4*t) +
                workday)

regr_res = dailymax - fitted(regr)
plot(regr_res)
