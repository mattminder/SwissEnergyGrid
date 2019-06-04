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

# Remove first value (incomplete week)
dailyInt = apply.daily(energyTS$TotConsumption, sum)[-1]

# Remove leap days
remove_leap_days <- function(ts){
  idx = index(ts)
  select <- (format(idx, "%d")=="29" & format(idx, "%m")=="02") | 
    as.numeric(format(idx, "%Y"))<2009 | as.numeric(format(idx, "%Y"))>=2019
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
sundayIx = ix$wday == 1 | holidayIx

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

mod1 <- arima(dailyInt/100, order=c(1, 0, 3), 
              seasonal = list(order=c(2, 1, 2), period=7),
              xreg=cbind(saturdayIx, sundayIx, dailyT, #sw1, cw1, #sw2, cw2,
                         sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))
tsdiag(mod1)


par(mfrow=c(1,1))
qqnorm(mod1$residuals)
qqline(mod1$residuals)


# --------------------- MODEL 2 ---------------------------------
# ------- regress out weekly, yearly and holidays, garch --------
ressq <- mod1$residuals * mod1$residuals
Acf(ressq)
Pacf(ressq)

mod2spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                       distribution.model = "std")
mod2 <- ugarchfit(spec = mod2spec, data = mod1$residuals, solver = "hybrid")
tsdiag(mod2)

Pacf(mod2@fit$residuals)
Acf(mod1$residuals)

qqnorm(mod2@fit$residuals)
