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

# the bars tell us that there is no trend

# ------------------------- remove leap days -------------------------------
# remove leap days
remove_leap_days <- function(ts){
  idx = index(ts)
  select <- (format(idx, "%d")=="29" & format(idx, "%m")=="02") | as.numeric(format(idx, "%Y"))<2009 | as.numeric(format(idx, "%Y"))>=2019
  result = ts[!select]
}
consumption = remove_leap_days(energyTS$TotConsumption)

# plot
plot(consumption)
# plot ~2 weeks of data
plot(consumption[2300:(4000)])

# add frequency attribute
attr(consumption, 'frequency') <- 35040 # 15min, yearly frequency

# trend seasonality residual decomposition
plot(stl(consumption, "per", robust = TRUE))
# strong seasonality

# ------------------------- periodigram and acf -------------------------------
# acf of ts does not make much sense since it is clearly not stationnary
acf(consumption)

# periodigram
per = spectrum(consumption)
plot(per)
# get frequency with highest amplitude
per$freq[which.max(per$spec)]

# ------------------------- fit sarima (000)(000) with yearly regressors -------------------------------
t = 1:length(consumption)
s_yearly = matrix(sin(2*pi*t/(24*4*365)), ncol = 1)
c_yearly = matrix(cos(2*pi*t/(24*4*365)), ncol = 1)

regressor = cbind(c_yearly, s_yearly)
fit_yearly = Arima(as.ts(consumption), order = c(0, 0, 0), seasonal = c(0, 0, 0), xreg = regressor)



# ------------------------- create regressor for weekly behavior -------------------------------
n_per_week = 24 * 4 * 7
n_weeks = floor(length(consumption)/n_per_week)
mat = scale(matrix(consumption[1:(n_weeks*n_per_week)], nrow = n_per_week))
typical_week = rowMeans(mat)
week_regressor = matrix(rep(typical_week, n_weeks+1)[1:length(consumption)], ncol=1)

regressor2 = cbind(c_yearly, s_yearly, week_regressor)
fit_year_weekly = Arima(as.ts(consumption), order = c(0, 0, 0), seasonal = c(0, 0, 0), xreg = regressor2)


# ------------------------- ACF analysis -------------------------------
Acf(fit_year_weekly$residuals, lag.max = 2000) # suggests differencing seasonal 672, seasonality of 96 is weaker
Pacf(diff(fit_year_weekly$residuals, lag = n_per_week), lag.max = 500) #-> suggests seasonality of 96->daily
Pacf(diff(fit_year_weekly$residuals, lag = n_per_week), lag.max = 2000) # -> suggests seasonality of 672 -> weekly


# fitting arima
fit_year_weekly = Arima(as.ts(consumption), order = c(1, 0, 1), seasonal = list(order=c(0,1,1),period=n_per_week), xreg = regressor)

