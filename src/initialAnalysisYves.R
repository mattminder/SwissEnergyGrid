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

# the bars tell us that there is no trend

# ------------------------- remove leap days -------------------------------
# remove leap days
remove_leap_days <- function(ts){
  idx = index(ts)
  select <- (format(idx, "%d")=="29" & format(idx, "%m")=="02") | as.numeric(format(idx, "%Y"))<2009 | as.numeric(format(idx, "%Y"))>=2019
  result = ts[!select]
}
consumption = remove_leap_days(energyTS$TotConsumption)

# plot
plot(consumption)
# plot ~2 weeks of data
plot(consumption[2300:(4000)])

# add frequency attribute
attr(consumption, 'frequency') <- 35040 # 15min, yearly frequency

# trend seasonality residual decomposition
plot(stl(consumption, "per", robust = TRUE))
# strong seasonality

# ------------------------- periodigram and acf -------------------------------
# acf of ts does not make much sense since it is clearly not stationnary
acf(consumption)

# periodigram
per = spectrum(consumption)
plot(per)
# get frequency with highest amplitude
per$freq[which.max(per$spec)]

# ------------------------- fit sarima (000)(000) with yearly regressors -------------------------------
t = 1:length(consumption)
s_yearly = matrix(sin(2*pi*t/(24*4*365)), ncol = 1)
c_yearly = matrix(cos(2*pi*t/(24*4*365)), ncol = 1)

regressor = cbind(c_yearly, s_yearly)
Arima(as.ts(consumption), order = c(0, 0, 0), seasonal = c(0, 0, 0), xreg = regressor)
