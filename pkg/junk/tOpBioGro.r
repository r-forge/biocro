## Testing the different optimization methods 

data(weather05)

## Some coefficients
pheno.ll <- phenoParms(kStem1=0.47, kLeaf1=0.48, kRoot1=0.05, kRhizome1=-1e-4,
                       kStem2=0.64, kLeaf2=0.15, kRoot2=0.21, kRhizome2=-1e-4,
                       kStem3=0.56, kLeaf3=0.01, kRoot3=0.13, kRhizome3=0.3, 
                       kStem4=0.56, kLeaf4=0.01, kRoot4=0.13, kRhizome4=0.3,
                       kStem5=0.56, kLeaf5=0.01, kRoot5=0.13, kRhizome5=0.3,
                       kStem6=0.56, kLeaf6=0.01, kRoot6=0.13, kRhizome6=0.3)

system.time(ans <- BioGro(weather05, phenoControl = pheno.ll))
system.time(ans0 <- BioGro(weather05))

ans.dat <- data.frame(ThermalT=ans$ThermalT, Stem=ans$Stem, Leaf=ans$Leaf,
                      Root=ans$Root, Rhizome=ans$Rhizome, Grain=ans$Grain,
                      LAI=ans$LAI)
sel.rows <- seq(1,nrow(ans.dat),400)
simDat <- ans.dat[sel.rows,]
plot(ans,simDat)

## Residual sum of squares before the optimization

ans0 <- BioGro(weather05)
RssBioGro(simDat,ans0)

## Limited testing suggests that this method works best
## phen = 0 optimizes all phenological stages sequentially
## Some may not converge and can be re-run
## Can also try other methods
system.time(op <- OpBioGro(phen=0, WetDat=weather05, data = simDat))
## or
## system.time(cop <- constrOpBioGro(phen = 0, WetDat = weather05, data = simDat))
            
phenoP <- phenoParms()
phenoP[7:30] <- op$coefs[7:30-6]

ans.op <- BioGro(WetDat = weather05, phenoControl = phenoP)
RssBioGro(simDat,ans.op)
plot(ans.op, simDat)




