## Testing some aspects of BioGro

## Testing the effect of different water stress models on BioCro
data(weather05)

soilP0 <- soilParms(wsFun = 'none', soilDepth=2)
## ans0 <- BioGro(weather05, day1=120, dayn=121, soilControl = soilP0)
ans0 <- BioGro(weather05, soilControl = soilP0)
plot(ans0)
plot(ans0, plot.kind="SW")
plot(ans0, plot.kind="ET")
plot(ans0, plot.kind="cumET")
plot(ans0, plot.kind="stress")

soilP1 <- soilParms(wsFun = 'thresh', smthresh=0.3, soilDepth = 2)
ans1 <- BioGro(weather05, soilControl = soilP1)
plot(ans1)
plot(ans1, plot.kind="SW")
plot(ans1, plot.kind="ET")
plot(ans1, plot.kind="cumET")
plot(ans1, plot.kind="stress")

## Effect of the new sunML function on Canopy level transpiration
soilP <- soilParms(soilLayers = 5, wsFun= 'logistic', hydrDist = TRUE)
ans2 <- BioGro(weather05, soilControl = soilP)
plot(ans2)
plot(ans2, plot.kind="SW")
plot(ans2, plot.kind="ET")
plot(ans2, plot.kind="cumET")
plot(ans2, plot.kind="stress")

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

soilP <- soilParms(soilLayers = 5, soilDepth=1.5 , wsFun = 'linear', hydrDist=TRUE, lrt=0.6, lrf=1e-3)
res <- BioGro(weather05, day1=120, dayn=270, soilControl = soilP)
plot(res)
plot(res, plot.kind="SW")
plot(res, plot.kind="ET")
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
soilP <- soilParms(soilLayers = 1, soilDepth=1.5 , wsFun = 'linear')
res <- MaizeGro(weather05, plant.day = 110, emerge.day = 120, harvest.day=300,
                MaizePhenoControl = MaizePhenoParms(R6 = 2300),
                soilControl = soilP,
                laiControl = laiParms(lai.method='ind-leaf-Lizaso', Aex=1500))
plot(res)

res2 <- BioGro(weather05, soilControl = soilParms(wsFun='none', soilLayers=1))
plot(res2)

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
