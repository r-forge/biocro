## Load in the data

## Biomass partitioning data for hybrid A
data(maizeBioA)

phenoP <-  MaizePhenoParms(R1=750, R6 = 1450, base.temp=8)

cfs.a <- idbpm(maizeBioA, MaizePhenoControl = phenoP)
cfs.a[8] <- 1e-4

## weather data
data(ausWea)

laiPa <- laiParms(lai.method="ind-leaf-Lizaso", Aex=880, a1=-4.51, a2=-0.089, LL = 850)
nitroPa <- MaizeNitroParms()
seneP <- MaizeSeneParms(senLeaf=700, senStem = 700)
canopyP <- canopyParms(Sp = 2.5)
photoP <- MaizePhotoParms(vmax = 56, alpha = 0.06)
phenoP <-  MaizePhenoParms(R1=750, R6 = 1450, base.temp=8)

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- cfs.a

res.a <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                  MaizeSeneControl = seneP,
                  canopyControl = canopyP,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.a, maizeBioA)

## Looking at hybrid B
data(maizeBioB)

phenoP <-  MaizePhenoParms(R1=750, R6 = 1450, base.temp=8)

cfs.b <- idbpm(maizeBioB, MaizePhenoControl = phenoP)

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- cfs.b

res.b <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                  canopyControl = canopyP,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.b, maizeBioB)

## Looking at hybrid C
data(maizeBioC)

phenoP <-  MaizePhenoParms(R1=750, R6 = 1450, base.temp=8)

cfs.c <- idbpm(maizeBioC, MaizePhenoControl = phenoP)

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- cfs.c

res.c <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                  canopyControl = canopyP,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.c, maizeBioC)

## Looking at hybrid D
data(maizeBioD)

phenoP <-  MaizePhenoParms(R1=750, R6 = 1450, base.temp=8)

cfs.d <- idbpm(maizeBioD, MaizePhenoControl = phenoP)

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- cfs.d

res.d <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                  canopyControl = canopyP,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.d, maizeBioD)

## Looking at hybrid E
data(maizeBioE)

phenoP <-  MaizePhenoParms(R1=750, R6 = 1450, base.temp=8)

cfs.e <- idbpm(maizeBioE, MaizePhenoControl = phenoP)

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- cfs.e

res.e <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                  canopyControl = canopyP,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.e, maizeBioE)


## Trying to optimize for Hybrid A

op <- OpMaizeGro(phen = 1,
                 iCoef = cfs.a,
                 cTT = res.a$TTTc,
                 WetDat=ausWea,
                 data = maizeBioA,
                 plant.day = 116,
                 emerge.day = 120,
                 harvest.day=258,
                 plant.density= 7,
                 lat=27,
                 MaizePhenoControl = phenoP, 
                 laiControl= laiPa,
                 MaizeCAllocControl= mCallocaP)

op2 <- OpMaizeGro(phen = 2,
                 iCoef = cfs.a,
                 cTT = res.a$TTTc,
                 WetDat=ausWea,
                 data = maizeBioA,
                 plant.day = 116,
                 emerge.day = 120,
                 harvest.day=258,
                 plant.density= 7,
                 lat=27,
                 MaizePhenoControl = phenoP, 
                 laiControl= laiPa,
                 MaizeCAllocControl= mCallocaP)

## Testing the results

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- op2$coefs

res.a <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.a, maizeBioA)

## optimizing pheno stage 3
op3 <- OpMaizeGro(phen = 3,
                 iCoef = op2$coefs,
                 cTT = res.a$TTTc,
                 WetDat=ausWea,
                 data = maizeBio1,
                 plant.day = 116,
                 emerge.day = 120,
                 harvest.day=258,
                 plant.density= 7,
                 lat=27,
                 MaizePhenoControl = phenoP, 
                 laiControl= laiPa,
                 MaizeCAllocControl= mCallocaP)

## Testing the results

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- op3$coefs
mCallocaP[8] <- op3$coefs[9]
mCallocaP[9] <- op3$coefs[8]

res.a <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.a, maizeBio1)

## optimizing pheno stage 3
op4 <- OpMaizeGro(phen = 3,
                 iCoef = op3$coefs,
                 cTT = res.a$TTTc,
                 WetDat=ausWea,
                 data = maizeBio1,
                 plant.day = 116,
                 emerge.day = 120,
                 harvest.day=258,
                 plant.density= 9,
                 lat=27,
                 MaizePhenoControl = phenoP, 
                 laiControl= laiPa,
                 MaizeCAllocControl= mCallocaP)

## Testing the results

mCallocaP <- MaizeCAllocParms()
mCallocaP[1:13] <- op3$coefs

res.a <- MaizeGro(ausWea, plant.day = 116, emerge.day = 120, harvest.day=258, plant.density= 7,
                  lat=27,
                  MaizePhenoControl = phenoP, 
                   laiControl= laiPa,
                   MaizeNitroControl = nitroPa,
                   MaizeCAllocControl= mCallocaP)

plot(res.a, maizeBio1)
