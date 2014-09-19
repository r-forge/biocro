## The R version has not been fully tested but it was used to prototype the C version

## Limited test for the R version
## Water infiltration without evapotranspiration
## Starting with a dry soil
## Without hydraulic distribution
data(weather05)
ans.bg <- BioGro(weather05)

## First scenario
## We have both precipitation and canopy transpiration
layers <- 10
precip <- weather05$precip[2952:7224]
cws <- rep(0.33,layers)
mat <- matrix(nrow=length(precip),ncol=layers)
mat2 <- matrix(nrow=length(precip),ncol=layers)
mat3 <- matrix(nrow=length(precip),ncol=layers)
for(i in seq_along(precip)){
  
  res <-  soilML(precip[i],ans.bg$CanopyTrans[i],cws,soilDepth=1.2,0.33,0.13,rootDB=10,soilLayers=layers,LAI=3,k=0.68,AirTemp=0,IRad=0,winds=3,RelH=0.8, soilType=6)

  cws <- res[,1]

  mat[i,] <- res[,1]
  mat2[i,] <- res[,5]
  mat3[i,] <- res[,6]
}

## Soil Moisture by layer
matplot(1:nrow(mat),mat, type="l")
plot(1:nrow(mat),apply(mat,1,mean),type="l")
## Canopy Transpiration as distributed on the soil layers
matplot(1:nrow(mat2),mat2, type="l")



## Previous scenario plus hydraulic distribution

layers <- 10
precip <- weather05$precip[2952:7224]
cws <- rep(0.33,layers)
mat <- matrix(nrow=length(precip),ncol=layers)
mat2 <- matrix(nrow=length(precip),ncol=layers)
mat3 <- matrix(nrow=length(precip),ncol=layers)
for(i in seq_along(precip)){
  
  res <-  soilML(precip[i],ans.bg$CanopyTrans[i],cws,soilDepth=1.5,0.33,0.13,rootDB=10,soilLayers=layers,LAI=3,k=0.68,AirTemp=0,IRad=0,winds=3,RelH=0.8, hydrDist=1, soilType=6)

  cws <- res[,1]

  mat[i,] <- res[,1]
  mat2[i,] <- res[,5]
  mat3[i,] <- res[,6]
}

## Soil Moisture by layer
matplot(1:nrow(mat),mat, type="l")
plot(1:nrow(mat),apply(mat,1,mean),type="l")
## Canopy Transpiration as distributed on the soil layers
matplot(1:nrow(mat2),mat2, type="l")


## Testing the overall model 
data(weather05)
wetdat <- weather05

soil.ll0.2.1L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=1) ## optimized soil depth
soil.ll0.2.2L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=2) ## optimized soil depth
soil.ll0.2.3L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=3) ## optimized soil depth
soil.ll0.2.4L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=4) ## optimized soil depth
soil.ll0.2.5L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=5) ## optimized soil depth
soil.ll0.2.6L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=6) ## optimized soil depth
soil.ll0.2.10L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=10) ## optimized soil depth
soil.ll0.2.36L <- soilParms(wsFun="none", FieldC=0.34, WiltP=0.13 , soilDepth = 1.5, soilLayers=36) ## optimized soil depth


sene.ll <- seneParms(senLeaf=2700)

ans0.2.1L <- BioGro(wetdat, soilControl = soil.ll0.2.1L, seneControl=sene.ll)
ans0.2.2L <- BioGro(wetdat, soilControl = soil.ll0.2.2L, seneControl=sene.ll)
ans0.2.3L <- BioGro(wetdat, soilControl = soil.ll0.2.3L, seneControl=sene.ll)
ans0.2.4L <- BioGro(wetdat, soilControl = soil.ll0.2.4L, seneControl=sene.ll)
ans0.2.5L <- BioGro(wetdat, soilControl = soil.ll0.2.5L, seneControl=sene.ll)
ans0.2.6L <- BioGro(wetdat, soilControl = soil.ll0.2.6L, seneControl=sene.ll)
ans0.2.10L <- BioGro(wetdat, soilControl = soil.ll0.2.10L, seneControl=sene.ll)
ans0.2.36L <- BioGro(wetdat, soilControl = soil.ll0.2.36L, seneControl=sene.ll)

xyplot(ans0.2.1L$SoilWatCont +
       ans0.2.2L$SoilWatCont +
       ans0.2.3L$SoilWatCont +
       ans0.2.4L$SoilWatCont +
       ans0.2.5L$SoilWatCont +
       ans0.2.6L$SoilWatCont +
       ans0.2.10L$SoilWatCont +
       ans0.2.36L$SoilWatCont ~ ans0.2.1L$DayofYear, type="l", auto.key=TRUE,
       ylab = "Soil Moisture")

xyplot(ans0.2.1L$SoilEvaporation +
       ans0.2.2L$SoilEvaporation +
       ans0.2.3L$SoilEvaporation +
       ans0.2.4L$SoilEvaporation +
       ans0.2.5L$SoilEvaporation +
       ans0.2.6L$SoilEvaporation +
       ans0.2.10L$SoilEvaporation +
       ans0.2.36L$SoilEvaporation ~ ans0.2.1L$DayofYear, type="l", auto.key=TRUE,
       ylab = "Soil Evaporation")

xyplot(ans0.2.1L$CanopyT + ans0.2.2L$CanopyT + ans0.2.3L$CanopyT +
       ans0.2.4L$CanopyT + ans0.2.5L$CanopyT + ans0.2.6L$CanopyT +
       ans0.2.10L$CanopyT ~ ans0.2.1L$ThermalT, type="l", auto.key=TRUE,
       ylab = "Canopy Transpiration")

xyplot(ans0.2.1L$Stem + ans0.2.2L$Stem + ans0.2.3L$Stem +
       ans0.2.4L$Stem + ans0.2.5L$Stem + ans0.2.6L$Stem +
       ans0.2.10L$Stem ~ ans0.2.1L$ThermalT, type="l", auto.key=TRUE,
       ylab = "Stem Biomass")

## More in detail looking at each depth

## This section depends on msummarize a function on a personal package
cwsMat <- ans0.2.10L$cwsMat
colnames(cwsMat) <- factor(soil.ll0.2.10L$soilDepths)[-1]
cws.dat <- data.frame(cwsMat,DOY=ans0.2.10L$DayofYear)
cws.datS <- msummarize(DOY ~ X0.15 + X0.3 + X0.45 + X0.6 + X0.75 + X0.9 + X1.05 + X1.2 + X1.35 + X1.5 ,FUN=mean, data = cws.dat)
xyplot(cws.datS[,2] + cws.datS[,3] + cws.datS[,4] + cws.datS[,5] + cws.datS[,6] + cws.datS[,7]
       + cws.datS[,8] + cws.datS[,9] + cws.datS[,10] ~ DOY, data = cws.datS, type="l",ylab="SWC")
cws.datST <- data.frame(DOY=rep(cws.datS$DOY,each=10),Depth=rep(-c(soil.ll0.2.10L$soilDepths[-1]),179),SWC=c(t(cws.datS[,-1])))
cws.datST2 <- cws.datST[cws.datST$DOY %in% seq(90,300,30),]
xyplot(Depth ~ SWC ,groups =factor(DOY), data = cws.datST2,  type="o", auto.key=TRUE)

