library(gdata)
library(xts)


# ------------------------- load the data---------------------------
# Set path
homePath = "/Users/myfiles/Documents/EPFL/M_II/TSE/SwissEnergyGrid/"
dataPath = paste(homePath, "data/", sep="")
resPath = paste(homePath, "res/", sep="")
outPath = paste(homePath, "res/RObjects/", sep="")

energyTS = readRDS(paste(outPath, "energyTS.RDS", sep=""))

# ---------------------- remove leap days ---------------------------
remove_leap_days <- function(ts){
  idx = index(ts)
  select <- (format(idx, "%d")=="29" & format(idx, "%m")=="02") | 
    as.numeric(format(idx, "%Y"))<2009 | as.numeric(format(idx, "%Y"))>=2019
  result = ts[!select]
}
consumption = remove_leap_days(energyTS$TotConsumption)


pdf(paste(resPath, "plots/EntireSeries.pdf", sep=""), width=8, height=4)
plot(consumption)
dev.off()

pdf(paste(resPath, "plots/OneWeek.pdf", sep=""), width=8, height=4)
plot(consumption[60200:61535])
dev.off()

pdf(paste(resPath, "plots/OneDay.pdf", sep=""), width=8, height=4)
plot(consumption[60284:60480])
dev.off()


# add frequency attribute
attr(consumption, 'frequency') <- 35040 # 15min, yearly frequency

# trend seasonality residual decomposition
png(paste(resPath, "plots/STLdecomp.png", sep=""), width=1200, height=1200)
plot(stl(consumption, "per", robust = TRUE))
dev.off()

periodogram = spectrum(consumption)
par(mar = rep(0, 4))
png(paste(resPath, "plots/Periodogram.png", sep=""), width=1200, height=800)
plot(periodogram)
dev.off()


