### IMPORTS ####
library(gdata)
library(xts)
library(dplyr)
library(forecast)

# ------------------------- load the data-------------------------------
# Set path
homePath = "/Users/yvesrychener/Studium/TimeSeries/SwissEnergyGrid/"
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

# ------------------------- MODEL 1 --------------------------------
# Regress out winter holidays and temperature, subsequent ARIMA
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

# ------------------------- MODEL 2 --------------------------------
# Regress out winter holidays and fit sinusoidal with multiple 
# harmonics of period 
# -------------------------- regression ---------------------------
per = 365.25/7
t = 1:length(weeklyInt)
reg2 <- lm(weeklyInt ~ holidays + 
            sin(2*pi*t/per)+cos(2*pi*t/per) + 
            sin(2*pi*t*2/per)+cos(2*pi*t*3/per) + 
            sin(2*pi*t*3/per)+cos(2*pi*t*4/per) + 
            sin(2*pi*t*4/per)+cos(2*pi*t*5/per))
res2 <- weeklyInt - fitted(reg2)
plot(res2)
plot(diff(res2))

spectrum(res2)
par(mfrow= c(2,1))
Acf(res2)
Pacf(res2)
par(mfrow=c(1,1))

# ------------------------- arima fitting --------------------------
res2_arima <- arima(res2, c(1, 1, 4))
tsdiag(res2_arima)
cpgram(res2_arima$residuals)
std_res2 = res2_arima$residuals / sd(res2_arima$residuals)

qqnorm(std_res2)
qqline(std_res2)

# ------------------------- Model 3 --------------------------
#  Sarima without regressors
# ------------------------- ACF & PACF --------------------------
spectrum(weeklyInt)
par(mfrow= c(2,1))
Acf(diff(weeklyInt, lag = 52))
Pacf(diff(weeklyInt, lag = 52))
par(mfrow=c(1,1))
# ------------------------- Fitting model --------------------------
res_sarima = arima(weeklyInt, order=c(0, 0, 3), seasonal = list(order = c(0, 1, 2), period = 52))

tsdiag(res_sarima)
cpgram(res_sarima$residuals)
std_res_sarima = res_sarima$residuals / sd(res_sarima$residuals)

qqnorm(std_res_sarima)
qqline(std_res_sarima)