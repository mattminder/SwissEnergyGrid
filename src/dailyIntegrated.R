### IMPORTS ####
library(gdata)
library(xts)
library(dplyr)
library(forecast)
library(astsa)
library(tis)
# ------------------------- load the data-------------------------------
# Set path
homePath = "/Users/myfiles/Documents/EPFL/M_II/TSE/SwissEnergyGrid/"
dataPath = paste(homePath, "data/", sep="")
outPath = paste(homePath, "res/Robjects/", sep="")

energyTS = readRDS(paste(outPath, "energyTS.RDS", sep=""))

dailyInt = apply.daily(energyTS$TotConsumption, sum)

# Remove leap days
remove_leap_days <- function(ts){
  idx = index(ts)
  select <- (format(idx, "%d")=="29" & format(idx, "%m")=="02") | 
    as.numeric(format(idx, "%Y"))<2009 | as.numeric(format(idx, "%Y"))>=2018
  result = ts[!select]
}

dailyInt = remove_leap_days(dailyInt)
attr(dailyInt, 'frequency') <- 365


# ----------------- define holidays -----------------------
# Holiday season between 24.12 and 2.1
# + August first, easter monday, good friday, pentecoste
# ascension
ix = index(dailyInt)
day = as.numeric(format(ix, "%d"))
month = as.numeric(format(ix, "%m"))
year = as.numeric(format(ix, "%y"))

ndays = length(ix)

weekendIx = ix$wday %in% c(0,6)

christmasIx = (day %in% c(24, 25, 26, 27, 28, 29, 30, 31) & month == 12) |
  (day%in% c(1,2) & month == 1)

easterIx = isEaster(ix)
whereEaster = which(easterIx)

goodfridayIx = rep(F, ndays)
goodfridayIx[whereEaster - 2] = T

easterMondayIx = rep(F, ndays)
easterMondayIx[whereEaster + 1] = T

pentecostIx = rep(F, ndays)
pentecostIx[whereEaster + 50] = T

ascensionIx = rep(F, ndays)
ascensionIx[whereEaster + 39] = T

ascensionFridayIx = rep(F, ndays)
ascensionFridayIx[whereEaster + 40] = T

firstOfAugustIx = (day == 1 & month == 8)

holidayIx = christmasIx | goodfridayIx | pentecostIx  | 
  ascensionIx | firstOfAugustIx | easterMondayIx | ascensionFridayIx

freedayIx = weekendIx | holidayIx

saturdayIx = ix$wday == 6
sundayIx = ix$wday == 0 | holidayIx

# -------------------- temperature ------------------------------
dailyT = apply.daily(energyTS$Temperature, mean)[-1]
dailyT = remove_leap_days(dailyT)


# -------------------- sine regressors ---------------------------
weekPer = 7
yearPer = 365
t = 1:length(dailyInt)

sw1 <- sin(2*pi*t/weekPer)
cw1 <- cos(2*pi*t/weekPer)
sw2 <- sin(2*pi*t*2/weekPer)
cw2 <- cos(2*pi*t*2/weekPer)
sw3 <- sin(2*pi*t*3/weekPer)
cw3 <- cos(2*pi*t*3/weekPer)
sw4 <- sin(2*pi*t*4/weekPer)
cw4 <- cos(2*pi*t*4/weekPer)
sw5 <- sin(2*pi*t*5/weekPer)
cw5 <- cos(2*pi*t*5/weekPer)

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


# --------------------- MODEL 1 ---------------------------------
# ---------- regress out weekly, yearly and holidays ------------

Acf(diff(dailyInt, 7), lag.max = 100)
Pacf(diff(dailyInt, 7), lag.max = 100)

mod1blank <- arima(dailyInt/10e9, #order=c(1, 0, 3), 
              xreg=cbind(saturdayIx, sundayIx, dailyT,
                         sw1,cw1,sw2,cw2,sw3,cw3,#sw4,cw4,sw5,cw5,
                         sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))

Acf(diff(mod1blank$residuals), lag.max = 100)
Pacf(diff(mod1blank$residuals), lag.max = 100)


mod1 <- arima(dailyInt/100, order=c(3, 1, 1), 
              xreg=cbind(saturdayIx, sundayIx, dailyT, sw1, cw1, sw2, cw2,
                         sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))
tsdiag(mod1)

# Still seasonality within residuals

par(mfrow=c(1,1))
qqnorm(mod1$residuals)
qqline(mod1$residuals)


# --------------------- MODEL 2 ---------------------------------
# ---------- regress out yearly and holidays, sarima ------------

mod2 <- arima(dailyInt/100, order=c(1, 0, 3), 
              seasonal = list(order=c(2, 1, 2), period=7),
              xreg=cbind(saturdayIx, sundayIx, dailyT, 
                         sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))
tsdiag(mod2)


par(mfrow=c(1,1))
qqnorm(mod2$residuals)
qqline(mod2$residuals)



# --------------------- MODEL 3 ---------------------------------
# ------------- set weekends and holidays to NA  ----------------
dailyInt_noWE <- dailyInt
dailyInt_noWE[saturdayIx] <- NA
dailyInt_noWE[sundayIx] <- NA

mod3blank <- arima(dailyInt_noWE/100, #order=c(3, 1, 1), 
                   xreg=cbind(dailyT))

Acf(diff(mod3blank$residuals), na.action = na.interp, lag.max = 100)
Pacf((mod3blank$residuals), na.action = na.interp, lag.max = 100)

plot(mod3blank$residuals)


mod3 <- arima(dailyInt_noWE/100, order=c(2, 1, 2), 
                   seasonal=list(order=c(0, 1, 1), period=7),
                   xreg=cbind(dailyT))


tsdiag(mod3)
qqnorm(mod3$residuals)
qqline(mod3$residuals)







# ressq <- mod1$residuals * mod1$residuals
# Acf(ressq)
# Pacf(ressq)
# 
# mod2spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
#                        distribution.model = "std")
# mod2 <- ugarchfit(spec = mod2spec, data = mod1$residuals, solver = "hybrid")
# tsdiag(mod2)
# 
# Pacf(mod2@fit$residuals)
# Acf(mod1$residuals)
# 
# qqnorm(mod2@fit$residuals)
