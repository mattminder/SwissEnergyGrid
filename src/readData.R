######
### IMPORTS ####
library(gdata)
library(xts)
library(dplyr)

#####
### READING ENERGY FILES ####
# Data availability:
# The data used for this project can be found under (as per 19.3.2019):
# https://www.swissgrid.ch/de/home/operation/grid-data.html

# Set path
homePath = "/Users/myfiles/Documents/EPFL/M_II/TSE/SwissEnergyGrid/"
dataPath = paste(homePath, "data/", sep="")
outPath = paste(homePath, "res/Robjects/", sep="")

# Read files
filenames = paste(dataPath,
                  "EnergieUebersichtCH-", 2009:2019,".xls",sep="")
tmp = list()
for (i in 1:length(filenames)) {
  # Only read first 21 columns, those present in all years
  tmp[[i]] = read.xls(filenames[i], sheet=3, encoding="latin1")[ , 1:21]
  tmp[[i]] = tmp[[i]][-1, ]
}
rawdata = do.call(rbind, tmp)
rm(tmp)
#saveRDS(rawdata, "rawdata.RDS")

# Redefine column names
names(rawdata) = c("Time", "TotEndUserConsumption", "TotProduction", 
                   "TotConsumption","NetOutflow", "GridFeedIn", "PosSecControl",
                   "NegSecControl", "PosTertControl", "NegTertControl", 
                   "XBorder_CH_AT", "XBorder_AT_CH", "XBorder_CH_DE",
                   "XBorder_DE_CH", "XBorder_CH_FR", "XBorder_FR_CH",
                   "XBorder_CH_IT", "XBorder_IT_CH", "Transit", "Import", 
                   "Export")

# Convert data to numeric
rawdataNum <- mutate_all(rawdata, function(x) as.numeric(as.character(x)))
rawdataNum[ , 1] = rawdata[ , 1] # Replace NA's created by dates

# Create time series
energyTS = xts(rawdataNum[ ,-1], order.by=as.POSIXlt(rawdataNum[ , 1]))

# Save rawdata
saveRDS(rawdataNum, paste(outPath, "energyData.RDS", sep=""))



######
### READING TEMPERATURE FILES ####
# source: http://www.agrometeo.ch/de/meteorology/datas

get_tempData = function(name) {
  tempFull = read.csv(paste(dataPath, "agrometeo-data-", name, ".csv", sep=""), sep=";", 
                      skip=0, stringsAsFactors = F)
  temp = as.numeric(tempFull[, 2])
  temp
}

temps = list()
temps$noflen = get_tempData("NOFLEN")
temps$pully = get_tempData("PULLY")
temps$thundorf = get_tempData("THUNDORF")
temps$kriessern = get_tempData("KRIESSERN")
temps$aesch = get_tempData("AESCH-BL")

allTemps = do.call(cbind, temps)

temp = rowMeans(allTemps, na.rm=TRUE)


# Convert 10min to 15min time series by linear interpolation
totTempVals = length(temp)
exactVals = temp[seq(1, totTempVals, 3)]
at10 = temp[seq(2, totTempVals, 3)] # values at 10/40 past
at20 = temp[seq(3, totTempVals, 3)] # values at 20/50 past
at15 = (at10 + at20)/2

converted = as.vector(rbind(exactVals, at15)) # Convert into 1 vector

energyTS$Temperature <- converted

plot(energyTS$Temperature)

#####
### SAVE ####
# Save timeseries
saveRDS(energyTS, paste(outPath, "energyTS.RDS", sep=""))

