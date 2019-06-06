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
plotPath = paste(homePath, "res/plots/", sep="")
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

sundayIx = ix$wday == 0 | holidayIx
saturdayIx = ix$wday == 6 & !sundayIx


# -------------------- temperature ------------------------------
dailyT = apply.daily(energyTS$Temperature, mean)
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

mod1blank <- arima(dailyInt/10e3, #order=c(1, 0, 3), 
              xreg=cbind(saturdayIx, sundayIx, #dailyT,
                         sw1,cw1,sw2,cw2,#sw3,cw3,#sw4,cw4,sw5,cw5,
                         sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))

# pdf(paste(plotPath, "Daily_TempFit_before.pdf", sep=""), height=8, width=8,
#     useDingbats = F)
par(mfrow=c(2,2))
Acf(diff(mod1blank$residuals))#, lag.max = 100)
Pacf(diff(mod1blank$residuals))#, lag.max = 100)

Acf(diff(mod1blank$residuals), lag.max = 30)
Pacf(diff(mod1blank$residuals), lag.max = 30)
# dev.off()

mod1 <- arima(dailyInt/10e3, order=c(3, 1, 1), 
              xreg=cbind(saturdayIx, sundayIx, #dailyT))
                         sw1, cw1, sw2, cw2,#, sw3, cw3,
                         sy1,cy1))#,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))
tsdiag(mod1)
#abs(mod1$coef) - 2*diag(mod1$var.coef) < 0


# Still seasonality within residuals

par(mfrow=c(1,1))
qqnorm(mod1$residuals)
qqline(mod1$residuals)


# --------------------- MODEL 2 ---------------------------------
# ---------- regress out yearly and holidays, sarima ------------
mod2blank <- arima(dailyInt/10e3,  
                   xreg=cbind(saturdayIx, sundayIx, #dailyT, 
                              sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))
par(mfrow=c(2,2))
Acf(diff(mod2blank$residuals, 7))#, lag.max = 100)
Pacf(diff(mod2blank$residuals, 7))#, lag.max = 100)

Acf(diff(mod2blank$residuals, 7), lag.max = 30)
Pacf(diff(mod2blank$residuals, 7), lag.max = 30)
# dev.off()




mod2 <- arima(dailyInt/10e3, order=c(1, 0, 3),
              seasonal = list(order=c(2, 1, 2), period=7),
              xreg=cbind(saturdayIx, sundayIx, #dailyT,
                         sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))
# mod2 <- arima(dailyInt/10e3, order=c(1, 0, 3), 
#               seasonal = list(order=c(2, 1, 2), period=7),
#               xreg=cbind(saturdayIx, sundayIx, dailyT, 
#                          sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5))


pdf(paste(plotPath, "daily_weeksarima_diag.pdf", sep=""), height=8, width=4,
    useDingbats = F)
tsdiag(mod2)
dev.off()


pdf(paste(plotPath, "daily_weeksarima_diag2.pdf", sep=""), height=8, width=4,
    useDingbats = F)
par(mfrow=c(2,1))
cpgram(mod2$residuals)

qqnorm(mod2$residuals)
qqline(mod2$residuals)
dev.off()

pdf(paste(plotPath, "daily_weeksarima_acf.pdf", sep=""), height=6, width=8,
    useDingbats = F)
par(mfrow=c(1,1))
Acf(mod2$residuals)
dev.off()


# still structure in residuals, residuals non normal


# # --------------------- MODEL 2.1 ---------------------------------
# # -------------------- adding GARCH ------------------------------
# # Create regressor matrix, for Arima package 
# transf <- function(x) {
#   matrix(x, ncol=1)
# }
# reg21list = list(saturdayIx, sundayIx, dailyT, sw1, cw1, sw2, cw2,
#                sy1,cy1,sy2,cy2,sy3,cy3,sy4,cy4,sy5,cy5)
# reg21 = do.call(cbind, lapply(reglist, transf))
# 
# mod21_sar <- forecast::Arima(dailyInt/10e3, order = c(1, 0, 3), 
#                              seasonal = list(order = c(2, 1, 2),  period = 7),
#                              xreg = reg21)
# 
# # Set constraints to do GARCH-SARIMA, adapted from Practical 4
# lnames <- c(paste0("ar",which(sapply(mod21_sar$model$phi, function(th){ isTRUE(all.equal(th,0))}))),
#             paste0("ma",which(sapply(mod21_sar$model$theta, function(th){ isTRUE(all.equal(th,0))}))))
# constraints <- rep(list(0), length(lnames))
# names(constraints) <- lnames
# order <- c(length(mod21_sar$model$phi),length(mod21_sar$model$theta))
# 
# 
# 
# mod21_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),
#                                                external.regressors = reg21), 
#                          mean.model = list(armaOrder = order, include.mean = TRUE, 
#                                            external.regressors = reg21), 
#                          distribution.model = "norm", fixed.pars=constraints) 
# 
# mod21 <- ugarchfit(spec = mod21_spec, data = dailyInt/10e3, solver = "hybrid")
# 
# par(mfrow=c(2,1))
# plot(mod21@fit$residuals, type='l')
# acf(mod21@fit$residuals)
# 
# par(mfrow=c(1,1))
# qqnorm(mod21@fit$residuals)
# qqline(mod21@fit$residuals)
# cpgram(mod21@fit$residuals)
# 
# # absolutely terrible
# 
# 
# # -------------------- non-normal resiudals -------------------------
# 
# mod22_spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0, 0)),
#                          mean.model = list(armaOrder = order, include.mean = TRUE, 
#                                            external.regressors = reg21), 
#                          distribution.model = "std", fixed.pars=constraints) 
# 
# mod22 <- ugarchfit(spec = mod22_spec, data = dailyInt/10e6, solver = "hybrid")
# 
# par(mfrow=c(2,1))
# plot(mod22@fit$residuals, type='l')
# acf(mod22@fit$residuals)
# 
# par(mfrow=c(1,1))
# qqPlot(mod22@fit$residuals, distribution = "t", df=5)
# 
# qqnorm(mod22@fit$residuals)
# qqline(mod22@fit$residuals)
# cpgram(mod22@fit$residuals)

# not that much better



# --------------------- MODEL 3 ---------------------------------
# ------------- set weekends and holidays to NA  ----------------
dailyInt_noWE <- dailyInt
dailyInt_noWE[saturdayIx] <- NA
dailyInt_noWE[sundayIx] <- NA

mod3blank <- arima(dailyInt_noWE/10e3, #order=c(3, 1, 1), 
                   xreg=cbind(sy1, cy1, sy2, cy2, sy3, cy3, sy4, cy4, sy5, cy5))

Acf(diff(mod3blank$residuals), na.action = na.interp, lag.max = 100)
Pacf((mod3blank$residuals), na.action = na.interp, lag.max = 100)

plot(mod3blank$residuals)


mod3 <- arima(dailyInt_noWE/10e3, order=c(2, 1, 2), 
                   seasonal=list(order=c(0, 1, 1), period=7),
                   xreg=cbind(sy1, cy1, sy2, cy2, sy3, cy3, sy4, cy4, sy5, cy5))


pdf(paste(plotPath, "daily_noWE_diag.pdf", sep=""), height=8, width=4,
    useDingbats = F)
tsdiag(mod3)
dev.off()

pdf(paste(plotPath, "daily_noWE_diag2.pdf", sep=""), height=8, width=4,
    useDingbats = F)
par(mfrow=c(2,1))
qqnorm(mod3$residuals)
qqline(mod3$residuals)
cpgram(mod2$residuals)
dev.off()

pdf(paste(plotPath, "daily_noWE_acf.pdf", sep=""), height=6, width=8,
    useDingbats = F)
par(mfrow=c(1,1))
Acf(mod3$residuals, na.action = na.pass)
dev.off()
# looks pretty good





