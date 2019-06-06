### IMPORTS ####
library(gdata)
library(xts)
library(dplyr)
library(forecast)
library(ggplot2)

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
consumption = consumption[as.numeric(format(index(consumption), "%Y")) < 2018]
f = remove_leap_days(energyTS$TotConsumption)
forward = f[as.numeric(format(index(f), "%Y")) >= 2018]

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
Acf(consumption)

# periodigram
per = spectrum(consumption)
plot(per)
# get frequency with highest amplitude
per$freq[which.max(per$spec)]


# convert to msts
consumption_msts = msts(consumption, seasonal.periods = c(4*24, 4*24*7, 4*24*365))
# fit tbats
start_time <- Sys.time()
fit_tbats = tbats(consumption_msts)
end_time <- Sys.time()
print(end_time - start_time)

# diagnostic plots
res = residuals(fit_tbats)
plot(res)

qqnorm(res)
qqline(res)

cpgram(res)

par(mfrow=c(2,1))
Acf(res)
Pacf(res)
par(mfrow=c(1,1))

# forecasting 20 days
predictions = forecast(fit_tbats, h = 20*4*24)
yl = 'Forecast'
predictions %>%
  autoplot(include = 104, ylab = yl) +
  geom_line(
    aes(
      x = as.numeric(time(predictions$mean)),
      y = as.numeric(forward[1:(20*4*24)])
    ),
    col = "red")

# forecasting only 1 day
predictions = forecast(fit_tbats, h = 4*24)
yl = 'Forecast'
predictions %>%
  autoplot(include = 104, ylab = yl) +
  geom_line(
    aes(
      x = as.numeric(time(predictions$mean)),
      y = as.numeric(forward[1:(4*24)])
    ),
    col = "red")
