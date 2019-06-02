### IMPORTS ####
library(gdata)
library(xts)
library(dplyr)
library(forecast)
library(SMPracticals)
library(tseries)
library(glue)
library(astsa)
# ------------------------- load the data---------------------------
# Set path
homePath = "/Users/myfiles/Documents/EPFL/M_II/TSE/SwissEnergyGrid/"
dataPath = paste(homePath, "data/", sep="")
outPath = paste(homePath, "res/Robjects/", sep="")

energyTS = readRDS(paste(outPath, "energyTS.RDS", sep=""))

# the bars tell us that there is no trend

# ------------------------- remove leap days -----------------------
# remove leap days
remove_leap_days <- function(ts){
  idx = index(ts)
  select <- (format(idx, "%d")=="29" & format(idx, "%m")=="02") | 
    as.numeric(format(idx, "%Y"))<2009 | as.numeric(format(idx, "%Y"))>=2019
  result = ts[!select]
}
consumption = remove_leap_days(energyTS$TotConsumption)


# ------------------------- perform stl ----------------------------
cons_stl = stl(consumption, "per", robust = TRUE)

# -------------------- investigation of residuals ------------------
cons_stl_rem = cons_stl$time.series[,"remainder"]
plot(cons_stl_rem)
cons_stl_diff = diff(cons_stl_rem)
plot(cons_stl_diff)

acf(cons_stl_diff)
pacf(cons_stl_diff)
# No conclusive model to fit, large amount of sign params

cons_stl_diffdiff = diff(cons_stl_diff)
acf(cons_stl_diffdiff)
pacf(cons_stl_diffdiff)
# Better, but pacf suggests large AR dependence structure


# ------------------- separate weekends -------------------------
get_weekends = function(ts) {
  ix = index(ts)
  weekend = c(0, 6)
  return(ts[ix$wday %in% weekend])
}

get_workdays = function(ts) {
  ix = index(ts)
  weekend = c(0, 6)
  return(ts[!ix$wday %in% weekend])
}

cons_we = get_weekends(consumption)
cons_wd = get_workdays(consumption)

plot(cons_we)
plot(cons_wd)

plot(cons_we[0:20000])
plot(cons_wd[0:20000])

we_stl = stl(cons_we, 36000, robust=T)
plot(we_stl)


# ------------------- fit yearly sinusoidal -------------------------
t <- 1:length(consumption)
per = 24*4*365

# fit with one harmonic
sinfit <- lm(consumption ~ sin(2*pi/per*t)+cos(2*pi/per*t))
plot(t, consumption, type="l")
lines(fitted(sinfit)~t,col=4,lty=2)

sinfit_res = consumption-fitted(sinfit)
plot(sinfit_res)


# fit with multiple harmonics
fit_periodicity = function(y, per, nharmonics=1) {
  t <- 1:length(y)
  form = "y ~ sin(2*pi/per*t)+cos(2*pi/per*t)"
  for (i in 2:nharmonics) {
    form = glue("{form} + sin(2*pi/per/{i}*t)+cos(2*pi/per/{i}*t)")
  }
  sinfit <- lm(as.formula(form))
  sinfit
}

sinfit2 <- fit_periodicity(consumption, per, 2)

plot(t, consumption, type="l")
lines(fitted(sinfit2)~t,col=4,lty=2)

# doesn't really add anything to the fit


# ------------- fit year sinus and weekday -------------------------
days = index(consumption)
t <- 1:length(consumption)
per = 24*4*365

weekend_days = c(0, 6)
workday = as.numeric(!days$wday%in%weekend_days)
weday = as.numeric(days$wday%in%weekend_days)
weekfit <- lm(consumption ~ sin(2*pi/per*t)+cos(2*pi/per*t) + workday)

week_res = consumption - fitted(weekfit)
par(mfrow=c(2,1))
plot(week_res)
plot(week_res[60000:68000])


# ------------- fit year and week sinus  -------------------------
days = index(consumption)
t <- 1:length(consumption)
perY = 24*4*365
perD = 24*4
perW = 24*4*7

weekfit2 <- lm(consumption ~ sin(2*pi/perY*t)+cos(2*pi/perY*t) + 
                sin(2*pi/perW*t)+cos(2*pi/perW*t))

week2_res = consumption - fitted(weekfit2)
par(mfrow=c(2,1))
plot(week2_res)
plot(week2_res[60000:68000])



# ------------- fit year, day sinus and weekday -------------------------
# days = index(consumption)
# t <- 1:length(consumption)
# perY = 24*4*365
# perD = 24*4
# 
# weekend_days = c(0, 6)
# workday = as.numeric(!days$wday%in%weekend_days)
# weekfit <- lm(consumption ~ sin(2*pi/perY*t)+cos(2*pi/perY*t) + workday + 
#                 sin(2*pi/perD*t)+cos(2*pi/perD*t))
# 
# week_res = consumption - fitted(weekfit)
# par(mfrow=c(2,1))
# plot(week_res)
# plot(week_res[60000:68000])



# -----------------------fit sarima  ---------------------------------
acf(week_res)
pacf(week_res)

diff_week_res = diff(week_res)
diff_week_res[1] = 0
acf(diff_week_res)
pacf(diff_week_res)

# Somehow doesn't work
#mod <- sarima(week_res, 1, 0, 1, 1, 1, 1, 24*4)
#plot(mod)


# ---------------------- regress out temperature ---------------------
temperature = remove_leap_days(energyTS$Temperature)

par(mfrow=c(2,1))
plot(consumption)
plot(-temperature)

tempfit = lm(consumption ~ temperature)
tempres = consumption - fitted(tempfit)

par(mfrow=c(2,1))
plot(tempres)
# there seems still to be some structure within the residuals

plot(diff(tempres))
# Looks garch-y

# ------------- regress out smoothed temperature ---------------------
smooth_temp = smooth.spline(temperature)
smotempfit = lm(consumption ~ smooth_temp$y + workday)

smotemp_res = consumption - fitted(smotempfit)
plot(smotemp_res)
plot(week_res)
