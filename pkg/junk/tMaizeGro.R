## Data from Illinois

data(weather05)

laiP <- laiParms(lai.method="ind-leaf-Lizaso", Aex=750, LLx=1500)
nitroP <- MaizeNitroParms(Vmax.b1=0.569, kLN=0.1)
phenoP <- MaizePhenoParms(base.temp = 8)
seneP <- seneParms(senStem=1000, senLeaf=1000, senRoot=1200)
photoP <- MaizePhotoParms(vmax = c(56), Rd = c(3,2))

#soilP <- soilParms(soilLayers = 1, wsFun= 'linear', soilDepth = 0.3, iWatCont=0.3)
#soilP <- soilParms(soilLayers = 5, wsFun= 'thresh', hydrDist=TRUE)
#soilP <- soilParms(soilDepth = 2.5, wsFun = 'none', soilLayers=5)
#soilP <- soilParms(wsFun='none', soilDepth=1.5, rsec=0.3,
#                   iWatCont=0.5,
#                   soilDepths = c(0.05,0.2,0.5,0.8,1,1.5),
#                   soilLayers=5, hydrDist=TRUE)

soilP <- soilParms(soilLayers=1, soilDepth=2.5, wsFun='thresh')

## Effect of N fertilization on maximum leaf size and longevity
###
res <- MaizeGro(weather05, plant.day = 110, emerge.day = 111, harvest.day=280,
                plant.density=9,
                laiControl = laiP,
                soilControl = soilP,
                photoControl = photoP,
                MaizeNitroControl= nitroP,
                MaizeSeneControl = seneP,
                MaizePhenoControl = phenoP)

plot(res)
plot(res, plot.kind="SW")
plot(res, plot.kind="ET")
plot(res, plot.kind="cumET")
plot(res, plot.kind="LAI")
plot(res, plot.kind="pheno")

res1 <- MaizeGro(weather04, plant.day = 110, emerge.day = 120, harvest.day=300,
                laiControl = laiParms(lai.method="ind-leaf-Lizaso", Aex=800, LLx=1000),
                MaizePhenoControl = MaizePhenoParms(base.temp = 8))

## Testing Phenology and growing degree days
xyplot(TTTc ~ DayofYear , data = res)

xyplot(LAI ~ DayofYear , data = res)

xyplot(CanopyAssim ~ DayofYear , data = res)

xyplot(CanopyTrans ~ DayofYear , data = res)

xyplot(LeafNVec ~ TTTc, data = res)

xyplot(VmaxVec ~ TTTc, data = res)

## Some testing with different weather data

