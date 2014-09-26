## I will be testing the ET in BioCro
## Date: Nov 20 2013

data(weather05)

## Without water stress
res <- BioGro(weather05, day1 = 90, dayn=150,
              soilControl = soilParms(wsFun="none",
                soilDepth=20, soilLayers=50, hydrDist=TRUE))
plot(res, plot.kind='ET')
plot(res, plot.kind="SW")

resd <- as.data.frame(res[1:12])
xyplot(SoilWatCont ~ Hour | factor(DayofYear), data = resd,
       subset = DayofYear %in% c(183, 184, 185))

xyplot(CanopyTrans ~ Hour | factor(DayofYear), data = resd,
       subset = DayofYear %in% c(183, 184, 185), ylim = c(0,20))

## With water stress
res <- BioGro(weather05, soilControl = soilParms(wsFun="thresh", soilDepth=1.5))
plot(res, plot.kind='ET')
plot(res, plot.kind="SW")
plot(res)
plot(res$StomatalCondCoefs)

resd <- as.data.frame(res[1:12])
xyplot(SoilWatCont ~ Hour | factor(DayofYear), data = resd,
       subset = DayofYear %in% c(183, 184, 185))

xyplot(CanopyTrans ~ Hour | factor(DayofYear), data = resd,
       subset = DayofYear %in% c(183, 184, 185), ylim = c(0,20))

res <- BioGro(weather05, soilControl = soilParms(wsFun="exp", soilDepth=1.5))
plot(res, plot.kind='ET')
plot(res, plot.kind="SW")
plot(res)
plot(res$StomatalCondCoefs)

resd <- as.data.frame(res[1:12])
xyplot(SoilWatCont ~ Hour | factor(DayofYear), data = resd,
       subset = DayofYear %in% c(183, 184, 185))

xyplot(CanopyTrans ~ Hour | factor(DayofYear), data = resd,
       subset = DayofYear %in% c(183, 184, 185), ylim = c(0,20))

## Days 184 and 185 produce ridiculous high values of ET, why?

doy184 <- weather05[weather05$doy == 184,]
res.c.184 <- numeric(24)

for(i in 0:23){
  solar <- doy184[doy184$hour == i, "solarR"]
  temp <- doy184[doy184$hour == i, "DailyTemp.C"]
  rh <- doy184[doy184$hour == i, "RH"]
  ws <- doy184[doy184$hour == i, "WindSpeed"]
  prcp <- doy184[doy184$hour == i, "precip"]
  res.c.184[i] <- CanA(7, 184, i, solar, temp, rh, ws)$CanopyTrans
}
