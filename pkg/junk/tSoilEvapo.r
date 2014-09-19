## Testing soil evaporation

## Function for relationship between relative soil moisture and soil evaporation
rawc <- seq(0,1,0.01)
sef <- function(rawc, theta=5){
##  ans <- 1 - (1 + 1.3 * rawc)^theta
  ans <- exp(theta * rawc) / exp(theta)
  ans
}

plot(rawc, sef(rawc, 5))

data(weather05)

weather05$precip <- 0
weather05$precip[5000] <- 100

res <- BioGro(weather05, photoControl = photoParms(vmax=1e-6),
              iRhizome=1e-7, soilControl=soilParms(iWatCont=0, soilLayers=5,
                               FieldC=0.35, WiltP=0.2))

xyplot(SoilEvaporation ~ ThermalT, data = res)

plot(res, plot.kind="SW")

## OK so 1 Mg/ha/hr is equivalent to
## 1 m3 water / ha / hr
## 1e-4 m / hr
## 1e3 mm / hr

data(weather05)

res2 <- BioGro(weather05, day1=120, dayn=121, soilControl = soilParms(wsFun='none', soilLayers=1))

sum(res2$SoilEvaporation) * 1e-1 ## This is in mm now
sum(res2$CanopyTrans) * 1e-1 ## This is in mm now

plot(res2, plot.kind="SW")

res3 <- BioGro(weather05, soilControl = soilParms(wsFun='linear'))

sum(res3$SoilEvaporation) * 1e-1 ## This is in mm now
sum(res3$CanopyTrans) * 1e-1 ## This is in mm now it is 1247 in the old model which seems about right

## Isolating the soil evaporation function
##		/* Penman-Monteith */
PhiN <- 2
Temp <- 25
temp <- 25
SlopeFS = 0.338376068 +  0.011435897 * Temp +  0.001111111 * Temp^2;
Delta <- 4098 * (0.6108 * exp(17.27 * temp/(temp + 237.3)))/(temp + 237.3)^2

Evaporation = (SlopeFS * PhiN + LHV * PsycParam * SoilBoundaryLayer * DeltaPVa) / (LHV * (SlopeFS + PsycParam));
