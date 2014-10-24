## Testing some aspects of BioGro

## Testing the effect of different water stress models on BioCro
data(weather05)

soilP0 <- soilParms(wsFun = 'none', soilDepth=10)
## ans0 <- BioGro(weather05, day1=120, dayn=121, soilControl = soilP0)
ans0 <- BioGro(weather05, soilControl = soilP0)
sum(ans0$SoilEvaporation)
sum(ans0$CanopyTrans)

xyplot(ans0$SoilEvaporation + ans0$CanopyTrans ~ ans0$ThermalT)

soilP1 <- soilParms(wsFun = 'thresh', smthresh=0.3, soilDepth = 2)
ans1 <- BioGro(weather05, soilControl = soilP1)
plot(ans1, plot.kind="stress")


## If I aggregate by day
ctbd <- aggregate(ans0$CanopyTrans , by = list(ans0$DayofYear), FUN=sum)
names(ctbd) <- c("day","CT")
ctbd$CTmm <- ctbd$CT * 0.1
xyplot(CTmm ~ day , data = ctbd, ylab="CT in mm")

xyplot(ans0$LAI + ans1$LAI ~ ans0$ThermalT)

xyplot(ans0$CanopyTrans + ans1$CanopyTrans ~ ThermalT, data = ans0)

## Here the difference is because the leaf area index is affected
xyplot(ans0$LeafReductionCoefs + ans1$LeafReductionCoefs ~ ThermalT, data = ans0)

## Effect of the new sunML function on Canopy level transpiration

##
data(weather05)

soilP <- soilParms(soilLayers = 5, wsFun= 'logistic', hydrDist = TRUE)
res <- BioGro(weather05, soilControl = soilP)
plot(res, plot.kind="stress")

## png('./soilwater-new-sunML.png')
plot(res, plot.kind='SW')
## dev.off()
## png('./canopyTrans-new-sunML.png')
plot(res$CanopyTrans)
## dev.off()

## Looking at soil water potential
data(weather04)

soilP <- soilParms(soilLayers = 1, soilDepth=1.5 , wsFun = 'lwp')
res <- BioGro(weather05, soilControl = soilP)
plot(res, plot.kind="stress")

soilP <- soilParms(soilLayers = 5, soilDepth=1.5 , wsFun = 'linear', hydrDist=TRUE)
res <- BioGro(weather05, soilControl = soilP)
plot(res, plot.kind="SW")
plot(res, plot.kind="stress")

soilP <- soilParms(soilLayers = 1, soilDepth=1.5 , wsFun = 'linear')
res <- BioGro(weather05, soilControl = soilP)
plot(res, plot.kind="SW")
plot(res, plot.kind="stress")
plot(res)


plot(res$LeafPsimVec)


     aws <- seq(0,0.4,0.001)
     wats.P <- numeric(length(aws))
     wats.L <- numeric(length(aws))
     phi2 <- 1e-7
     for(i in 1:length(aws)){
     wats.P[i] <- wtrstr(1,1,aws[i],0.5,0.37,0.2,2e-2,phi2=phi2)$wsPhoto
     wats.L[i] <- wtrstr(1,1,aws[i],0.5,0.37,0.2,2e-2,phi2=phi2)$wsSpleaf
     }
     
     xyplot(wats.P + wats.L ~ aws,
            xlab="Soil Water",
            ylab="Stress Coefficient")


## Transpiration for maize
     data(weather05)
     res <- MaizeGro(weather05, plant.day = 110, emerge.day = 120, harvest.day=300,
                       MaizePhenoControl = MaizePhenoParms(R6 = 2300),
                     laiControl = laiParms(lai.method='ind-leaf-Lizaso', Aex=1500))

res2 <- BioGro(weather05, soilControl = soilParms(wsFun='none', soilLayers=1))
sum(res2$CanopyTrans + res2$SoilEvapo)
max((res2$Leaf + res2$Stem))

res2 <- BioGro(weather05, soilControl = soilParms(wsFun='linear', soilLayers=1, soilDepth=4))
sum(res2$CanopyTrans + res2$SoilEvapo)
max((res2$Leaf + res2$Stem))

res2 <- BioGro(weather05, soilControl = soilParms(wsFun='logistic',
                            soilLayers=1, phi1=1e-5, soilDepth = 3.3))
sum(res2$CanopyTrans + res2$SoilEvapo)
max((res2$Leaf + res2$Stem))

res2 <- BioGro(weather05, soilControl = soilParms(wsFun='exp', soilLayers=1, soilDepth=3.6))
sum(res2$CanopyTrans + res2$SoilEvapo)
max((res2$Leaf + res2$Stem))

res2 <- BioGro(weather05, soilControl = soilParms(wsFun='lwp', soilLayers=1, leafPotTh=-2500,
                            soilDepth = 3.2))
sum(res2$CanopyTrans + res2$SoilEvapo)
max((res2$Leaf + res2$Stem))
