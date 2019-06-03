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
dailyInt = apply.daily(energyTS$TotConsumption, sum)[-1]

# Remove value after 1. january 2019
dailyInt = dailyInt[as.numeric(format(index(dailyInt), "%Y")) < 2019]

plot(dailyInt)

# ------------------------- Analysis step 1-------------------------------
spectrum(dailyInt)
# peak correspond to yearly cycle plus harmonics
par(mfrow= c(2,1))
Acf(dailyInt)
Pacf(dailyInt)
par(mfrow=c(1,1))
# ACF clearly suggests taking seasonal difference, D=1, with Period 7

# ------------------------- Trying to Regress out yearly cycle-------------
per = 365
t = 1:length(dailyInt)
reg2 <- lm(dailyInt ~ sin(2*pi*t/per)+cos(2*pi*t/per) + 
             sin(2*pi*t*2/per)+cos(2*pi*t*3/per) + 
             sin(2*pi*t*3/per)+cos(2*pi*t*4/per) + 
             sin(2*pi*t*4/per)+cos(2*pi*t*5/per))
res2 <- dailyInt - fitted(reg2)


plot(res2)
spectrum(res2)
# peaks now correspond to roughly weekly cycle
par(mfrow= c(2,1))
Acf(res2)
Pacf(res2)
par(mfrow=c(1,1))

# ------------------------- Trying to Regress out yearly cycle & week-------------
sunday = as.numeric(as.POSIXlt(index(dailyInt))$wday==0)
saturday = as.numeric(as.POSIXlt(index(dailyInt))$wday==6)

per = 365
t = 1:length(dailyInt)
reg2 <- lm(dailyInt ~ saturday + sunday + sin(2*pi*t/per)+cos(2*pi*t/per) + 
             sin(2*pi*t*2/per)+cos(2*pi*t*3/per) + 
             sin(2*pi*t*3/per)+cos(2*pi*t*4/per) + 
             sin(2*pi*t*4/per)+cos(2*pi*t*5/per))
res2 <- dailyInt - fitted(reg2)


plot(res2)
spectrum(res2)
par(mfrow= c(2,1))
Acf(res2)
Pacf(res2)
par(mfrow=c(1,1))
