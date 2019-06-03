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
mod1blank <- arima(weeklyInt/100, 
              xreg=cbind(holidays, weeklyT))

par(mfrow=c(2,1))
Acf(mod1blank$residuals)
Acf(diff(mod1blank$residuals))

Pacf(mod1blank$residuals)
Pacf(diff(mod1blank$residuals))
par(mfrow=c(1,1))



mod11 <- arima(weeklyInt/100, order=c(1, 0, 4), 
              xreg=cbind(holidays, weeklyT))
tsdiag(mod11)

cpgram(mod11$residuals)
qqnorm(mod11$residuals)
qqline(mod11$residuals)


mod12 <- arima(weeklyInt/100, order=c(4, 1, 1), 
              xreg=cbind(holidays, weeklyT))
tsdiag(mod12)

cpgram(mod12$residuals)
qqnorm(mod12$residuals)
qqline(mod12$residuals)


# ------------------------- MODEL 2 --------------------------------
# Regress out winter holidays and fit sinusoidal with multiple 
# harmonics of period 
# -------------------------- regression ---------------------------
yearPer = 365.25/7
t = 1:length(weeklyInt)

sy1 <- sin(2*pi*t/yearPer)
cy1 <- cos(2*pi*t/yearPer)
sy2 <- sin(2*pi*t*2/yearPer)
cy2 <- cos(2*pi*t*2/yearPer)
sy3 <- sin(2*pi*t*3/yearPer)
cy3 <- cos(2*pi*t*3/yearPer)
sy4 <- sin(2*pi*t*4/yearPer)
cy4 <- cos(2*pi*t*4/yearPer)
sy5 <- sin(2*pi*t*5/yearPer)
cy5 <- cos(2*pi*t*5/yearPer)


mod2blank <- arima(weeklyInt/100,
                   xreg = cbind(sy1, cy1, sy2, cy2, 
                                sy3, cy3, sy4, cy4, sy5, cy5))

par(mfrow=c(2,1))
Acf(mod2blank$residuals)
Acf(diff(mod2blank$residuals))

Pacf(mod2blank$residuals)
Pacf(diff(mod2blank$residuals))
par(mfrow=c(1,1))

mod2 <- arima(weeklyInt/100, order = c(6, 1, 4),
                   xreg = cbind(sy1, cy1, sy2, cy2, 
                                sy3, cy3, sy4, cy4, sy5, cy5))

tsdiag(mod2)

cpgram(mod2$residuals)
qqnorm(mod2$residuals)
qqline(mod2$residuals)


# ------------------------- Model 3 --------------------------
#  Sarima without regressors
# ------------------------- ACF & PACF --------------------------
spectrum(weeklyInt)
par(mfrow= c(2,1))
Acf(weeklyInt, lag.max = 200)
Pacf(weeklyInt, lag.max = 200)

Acf(diff(weeklyInt, lag = 52), lag.max = 200)
Pacf(diff(weeklyInt, lag = 52), lag.max = 200)


Acf(diff(weeklyInt))
Pacf(diff(weeklyInt))
par(mfrow=c(1,1))
# ------------------------- Fitting model --------------------------
mod3 = arima(weeklyInt, order=c(3, 0, 0), seasonal = list(order = c(0, 1, 2), period = 52))

tsdiag(mod3)
cpgram(mod3$residuals)

qqnorm(mod3$residuals)
qqline(mod3$residuals)
