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

holidayIx = (day >= 24 & month == 12) | (day <= 2 & month == 1)
weekendIx = ix$wday %in% c(0,6)

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

holidayIx = holidayIx | goodfridayIx | pentecostIx  | 
  ascensionIx | firstOfAugustIx | easterMondayIx | ascensionFridayIx

# -------------------- temperature ------------------------------
dailyT = apply.daily(energyTS$Temperature, mean)[-1]
dailyT = remove_leap_days(dailyT)



# --------------------- MODEL 1 ---------------------------------
# ---------- regress out weekly, yearly and holidays ------------
weekPer = 7
yearPer = 365
t = 1:length(dailyInt)
reg1 <- lm(dailyInt ~ holidayIx + weekendIx +
             # sin(2*pi*t/weekPer)+cos(2*pi*t/weekPer) +
             # sin(2*pi*t*2/weekPer) + cos(2*pi*t*2/weekPer) +
             # sin(2*pi*t*3/weekPer) + cos(2*pi*t*3/weekPer) +
             sin(2*pi*t/yearPer)+cos(2*pi*t/yearPer) +
             sin(2*pi*t*2/yearPer)+cos(2*pi*t*2/yearPer) +
             sin(2*pi*t*3/yearPer)+cos(2*pi*t*3/yearPer) +
             sin(2*pi*t*4/yearPer)+cos(2*pi*t*4/yearPer) +
             sin(2*pi*t*5/yearPer)+cos(2*pi*t*5/yearPer) +
             dailyT)



res1 <- dailyInt - fitted(reg1)
plot(res1)
plot(diff(res1))

par(mfrow=c(3,1))
Acf(res1)
Acf(diff(res1))
Acf(diff(res1, 7))

Pacf(res1)
Pacf(diff(res1))
Pacf(diff(res1, 7))
par(mfrow=c(1,1))

res_decomp = stl(res1, "p", robust = T)
plot(res_decomp)


# --------------------------- sarima -------------------------------------
mod1 <- sarima(res1, 1, 1, 5, 1, 1, 1, 7)


# ---------------------- sarima grid search ------------------------------
# aic_grid = list()
# aic_param = list()
# for (p in 1:5) {
#   for (q in 1:5) {
#     for (P in 1:3) {
#       for (Q in 1:3) {
#         fit <- sarima(res1, p, 1, q, P, 1, Q, 7)
#         aic_grid = c(aic_grid, fit$AIC)
#         aic_param = c(aic_param, c(p, q, P, Q))
#       }
#     }
#   }
# }
# 
# res1_sarima = sarima(res1, 1, 1, 0, 3, 1, 1, 7)

# ds isch fr nix


