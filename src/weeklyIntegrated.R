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
plotPath = paste(homePath, "res/plots/", sep="")

energyTS = readRDS(paste(outPath, "energyTS.RDS", sep=""))

# Remove first value (incomplete week)
weeklyInt = apply.weekly(energyTS$TotConsumption, sum)[-1]

# Remove value after 1. january 2019
weeklyInt = weeklyInt[as.numeric(format(index(weeklyInt), "%Y")) < 2018]

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
weeklyT = weeklyT[as.numeric(format(index(weeklyT), "%Y")) < 2018]

pdf(paste(plotPath, "weeklyInt+Temperature.pdf", sep=""), height=4, width=8)
par(mfrow=c(2,1))
plot(weeklyInt)
plot(-weeklyT)
par(mfrow=c(1,1))
dev.off()

attr(weeklyInt, 'frequency') <- 52
spectrum(weeklyInt)

# ------------------------- MODEL 1 --------------------------------
# Regress out winter holidays and temperature, subsequent ARIMA
# -------------------------- regression ---------------------------
mod1blank <- arima(weeklyInt/10e6, 
              xreg=cbind(holidays, weeklyT))



pdf(paste(plotPath, "Weekly_TempFit_before.pdf", sep=""), height=8, width=8,
    useDingbats = F)
par(mfrow=c(2,2))
Acf(mod1blank$residuals)
Acf(diff(mod1blank$residuals))

Pacf(mod1blank$residuals)
Pacf(diff(mod1blank$residuals))
par(mfrow=c(1,1))
dev.off()


mod11 <- arima(weeklyInt/10e6, order=c(1, 0, 4), 
              xreg=cbind(holidays, weeklyT))
pdf(paste(plotPath, "Weekly_TempFit_diag.pdf", sep=""), height=8, width=4,
    useDingbats = F)
tsdiag(mod11)
dev.off()


pdf(paste(plotPath, "Weekly_TempFit_diag2.pdf", sep=""), height=8, width=4,
    useDingbats = F)
par(mfrow=c(2,1))
cpgram(mod11$residuals)
qqnorm(mod11$residuals)
qqline(mod11$residuals)
dev.off()

mod12 <- arima(weeklyInt/10e6, order=c(4, 1, 1), 
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


mod2blank <- arima(weeklyInt/10e6,
                   xreg = cbind(sy1, cy1, sy2, cy2, 
                                sy3, cy3, sy4, cy4, sy5, cy5))

pdf(paste(plotPath, "Weekly_SinCosFit_before.pdf", sep=""), height=8, width=8)
par(mfrow=c(2,2))
Acf(mod2blank$residuals)
Acf(diff(mod2blank$residuals))

Pacf(mod2blank$residuals)
Pacf(diff(mod2blank$residuals))
par(mfrow=c(1,1))
dev.off()

mod2 <- arima(weeklyInt/10e6, order = c(4, 1, 4),
                   xreg = cbind(sy1, cy1, sy2, cy2, 
                                sy3, cy3, sy4, cy4, sy5, cy5))

pdf(paste(plotPath, "Weekly_SinCosFit_diag.pdf", sep=""), height=8, width=4,
    useDingbats=F)
tsdiag(mod2)
dev.off()

pdf(paste(plotPath, "Weekly_SinCosFit_diag2.pdf", sep=""), height=8, width=4,
    useDingbats=F)
par(mfrow=c(2,1))
cpgram(mod2$residuals)
qqnorm(mod2$residuals)
qqline(mod2$residuals)
dev.off()
par(mfrow=c(1,1))


# Forecasting
t_forecast = (length(weeklyInt)+1):(length(weeklyInt)+52)
sy1_f <- sin(2*pi*t_forecast/yearPer)
cy1_f <- cos(2*pi*t_forecast/yearPer)
sy2_f <- sin(2*pi*t_forecast*2/yearPer)
cy2_f <- cos(2*pi*t_forecast*2/yearPer)
sy3_f <- sin(2*pi*t_forecast*3/yearPer)
cy3_f <- cos(2*pi*t_forecast*3/yearPer)
sy4_f <- sin(2*pi*t_forecast*4/yearPer)
cy4_f <- cos(2*pi*t_forecast*4/yearPer)
sy5_f <- sin(2*pi*t_forecast*5/yearPer)
cy5_f <- cos(2*pi*t_forecast*5/yearPer)
predictions_2 = forecast(mod2, h = 52, xreg = cbind(sy1_f, cy1_f, sy2_f, cy2_f, 
                                                   sy3_f, cy3_f, sy4_f, cy4_f, sy5_f, cy5_f))
yl = 'Forecast'
plot(predictions_2, include = 104, ylab = yl)
pdf(paste(plotPath, "Weekly_arima_fourier_forecast.pdf", sep=""), height=8, width=4,
    useDingbats=F)
plot(predictions_2, include = 104, ylab = yl)
dev.off()

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
mod3 = arima(weeklyInt/10e6, order=c(1, 0, 2), seasonal = list(order = c(0, 1, 2), period = 52))

tsdiag(mod3)
cpgram(mod3$residuals)

par(mfrow=c(1,1))
qqnorm(mod3$residuals)
qqline(mod3$residuals)

pdf(paste(plotPath, "Weekly_sarima.pdf", sep=""), height=8, width=4,
    useDingbats=F)
tsdiag(mod3)
dev.off()

par(mfrow=c(1,1))
predictions_3 = forecast(mod3, h = 52)
plot(predictions_3, include = 104, ylab = yl)
yl = 'Forecast'
pdf(paste(plotPath, "Weekly_sarima_forecast.pdf", sep=""), height=8, width=4,
    useDingbats=F)
plot(predictions_3, include = 104, ylab = yl)
dev.off()

# ------------------------- Model 4 --------------------------


mod4blank <- arima(weeklyInt/10e6,
                   xreg=cbind(holidays, weeklyT))

Acf(mod4blank$residuals)
Acf(diff(mod4blank$residuals, 52))
Pacf(diff(mod4blank$residuals, 52))
Acf(diff(mod4blank$residuals))
Pacf(diff(mod4blank$residuals))

mod4 <- arima(weeklyInt/10e9,
              xreg=cbind(holidays, weeklyT),
              seasonal=list(order=c(1,1,1), period=52),
              order=c(4,1,1))
              #order=c(1,0,4))
tsdiag(mod4)
par(mfrow=c(1,1))
cpgram(mod4$residuals)
qqnorm(mod4$residuals)
qqline(mod4$residuals)
