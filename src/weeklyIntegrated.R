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

# Remove first value (incomplete week)
weeklyInt = apply.weekly(energyTS$TotConsumption, sum)[-1]

# Remove value after 1. january 2019
weeklyInt = weeklyInt[as.numeric(format(index(weeklyInt), "%Y")) < 2019]

plot(weeklyInt)

# ------------------------- winter holidays ----------------------------
# define as winter holidays weeks that end between december 25th and 
# january 8th
ix = index(weeklyInt)
day = as.numeric(format(ix, "%d"))
month = as.numeric(format(ix, "%m"))

holidays <- (month==12 & day > 24) | (month==1 & day < 9)

# ------------------------ mean temperature ---------------------------
weeklyT = apply.weekly(energyTS$Temperature, mean)[-1]
weeklyT = weeklyT[as.numeric(format(index(weeklyT), "%Y")) < 2019]

par(mfrow=c(2,1))
plot(-weeklyT)
plot(weeklyInt)
par(mfrow=c(1,1))

# -------------------------- regression ---------------------------
reg <- lm(weeklyInt ~ holidays + weeklyT)
res <- weeklyInt - fitted(reg)
plot(res)
plot(diff(res))

spectrum(res)
par(mfrow= c(2,1))
Acf(res)
Pacf(res)
par(mfrow=c(1,1))

# ------------------------- arima fitting --------------------------
res_arima <- arima(res, c(1, 1, 4))
tsdiag(res_arima)
cpgram(res_arima$residuals)
std_res = res_arima$residuals / sd(res_arima$residuals)

qqnorm(std_res)
qqline(std_res)
