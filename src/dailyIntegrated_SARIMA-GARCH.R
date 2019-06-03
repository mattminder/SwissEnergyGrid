### IMPORTS ####
library(gdata)
library(xts)
library(dplyr)
library(forecast)
library(astsa)
library(tis)
library(rugarch)
library(tseries)
# ------------------------- load the data-------------------------------
# Set path
homePath = "/Users/yvesrychener/Studium/TimeSeries/SwissEnergyGrid/"
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
dailyInt = (dailyInt-mean(dailyInt))/(var(dailyInt)[1]) # does not have to be statistically 100% correct
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



# --------------------- MODEL 1 ---------------------------------
# ---------- regress out weekly, yearly and holidays ------------
weekPer = 7
yearPer = 365
t = 1:length(dailyInt)
reg1 <- lm(dailyInt ~ saturdayIx + sundayIx + dailyT)



res1 <- dailyInt - fitted(reg1)
plot(res1)
plot(diff(res1))

par(mfrow=c(3,1))
Acf(res1)
Acf(diff(diff(res1,7)))
Acf(diff(res1, 7))

Pacf(res1)
Pacf(diff(diff(res1,7)))
Pacf(diff(res1, 7))
par(mfrow=c(1,1))

res_decomp = stl(res1, "p", robust = T)
plot(res_decomp)


# --------------------------- sarima -------------------------------------
#regressor = cbind(matrix(holidayIx, ncol=1), matrix(sundayIx, ncol=1),  matrix(weekendIx, ncol=1), matrix(dailyT, ncol=1))
sw1 <- sin(2*pi*t/weekPer)
cw1 <- cos(2*pi*t/weekPer)
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
meanreg = cbind(matrix(saturdayIx, ncol = 1), matrix(sundayIx, ncol = 1), matrix(dailyT, ncol = 1))



#msarima <- Arima(dailyInt, order = c(2, 0, 2), seasonal = list(order = c(1, 1, 2), period = 7), xreg = meanreg)
msarima <- Arima(dailyInt, order = c(2, 0, 2), seasonal = list(order = c(0, 1, 2), period = 7), xreg = meanreg)
tsdiag(msarima)
cpgram(msarima$residuals)
std_res = msarima$residuals / sd(msarima$residuals)

qqnorm(std_res)
qqline(std_res)

ressqd = msarima$residuals*msarima$residuals
par(mfrow=c(2,1))
Acf(ressqd, lag.max = 50)
Pacf(ressqd, lag.max = 50)
par(mfrow=c(1,1))
# suggest taking garch(1,1) 


# without sar term does not work, gives singularity error
msarima <- Arima(dailyInt, order = c(2, 0, 2), seasonal = list(order = c(1, 1, 2), period = 7), xreg = meanreg)
lnames <- c(paste0("ar", which(sapply(msarima$model$phi, function(th) {
  isTRUE(all.equal(th, 0))
}))), paste0("ma", which(sapply(msarima$model$theta, function(th) {
  isTRUE(all.equal(th, 0))
}))))
constraints <- rep(list(0), length(lnames))
names(constraints) <- lnames
order <- c(length(msarima$model$phi), length(msarima$model$theta))

model <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),
                        external.regressors = meanreg), 
                        mean.model = list(armaOrder = order, include.mean = TRUE, 
                        external.regressors = meanreg), 
                        distribution.model = "norm", fixed.pars=constraints) 

fitmodel <- ugarchfit(spec = model, data = dailyInt, solver = "hybrid")

par(mfrow=c(2,1))
plot(fitmodel@fit$residuals, type='l')
acf(fitmodel@fit$residuals)

par(mfrow=c(1,1))
qqnorm(fitmodel@fit$residuals)
qqline(fitmodel@fit$residuals)
cpgram(fitmodel@fit$residuals)
