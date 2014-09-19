## Test script for the Century Model

iSoilP0 <- iSoilP

iSoilP <- centuryParms(timestep="week")
iSoilP[1:9] <- iSoilP0[1:9]

#nsim <- 52*1000 ## weekly
nsim <- 52 * 1000 ## yearly

mat <- matrix(nrow=nsim,ncol=9)
mat2 <- matrix(nrow=nsim,ncol=4)
mat3 <- matrix(nrow=nsim,ncol=9)
respVec <- numeric(nsim)
MinNVec <- numeric(nsim)
  
iLeafL <- 100 * 0.1
iStemL <- 200 * 0.1
iRootL <- 200 * 0.1
iRhizL <- 100 * 0.1

add <- TRUE
for(i in 1:nsim){

  mat2[i,] <- c(iLeafL,iStemL,iRootL,iRhizL) 
  
  aLeafL <- iLeafL * 0.1
  aStemL <- iStemL * 0.1
  aRootL <- iRootL * 0.1
  aRhizL <- iRhizL * 0.1

##  res <- Century(aLeafL,aStemL,aRootL,aRhizL,0.23,15,10,50,centuryControl=iSoilP, verbose=FALSE)
  res <- CenturyC(aLeafL,aStemL,aRootL,aRhizL,0.23,15,10,50,centuryControl=iSoilP,soilType=5)

  iLeafL <- iLeafL - aLeafL
  iStemL <- iStemL - aStemL
  iRootL <- iRootL - aRootL
  iRhizL <- iRhizL - aRhizL
  
  for(j in 1:9) iSoilP[j] <- res$SCs[j] 
  iSoilP$iMinN <- res$MinN
  
  mat[i,] <- res$SCs
  mat3[i,] <- res$SNs
  
  respVec[i] <- res$Resp
  MinNVec[i] <- res$MinN

if(add){  
    if(i %% 52 == 0){
      iLeafL <- 100 * 0.1
      iStemL <- 200 * 0.1
      iRootL <- 200 * 0.1
      iRhizL <- 100 * 0.1
    }
  }
}

mat.9d <- mat

## Plotting the different soil carbon pools
labs <- paste("SC",1:9,sep="")
xyplot(mat[,1] + mat[,2] + mat[,3] +
       mat[,4] + mat[,5] + mat[,6] +
       mat[,7] + mat[,8] + mat[,9] ~ I(1:nsim/52),type="l",ylab="carbon",
       col=rainbow(9),
       lty=1, lwd=3,
       key=list(text=list(labs),lines=TRUE,col=rainbow(9)),
       main = paste("total annual biomass = 6 Mg/ha/yr"))

## plotting the organic matter accumulation
## Assume bulk density of 1.4 Mg/m2 up to 30cm = 420 kg
om.9 <- rowSums(mat.9[,5:8]) * 1e-3 * (1/420) * 2 * 100
om.6 <- rowSums(mat.6[,5:8]) * 1e-3 * (1/420) * 2 * 100
om.3 <- rowSums(mat.3[,5:8]) * 1e-3 * (1/420) * 2 * 100
pdf("./om-century.pdf")
xyplot(om.9 + om.6 + om.3
       ~ I(1:nsim/52),type="l",ylab="organic matter (%)",
       xlab = "years",
       col = c("purple","blue","green"),
       lty=1, lwd=3,
       key=list(text=list(c("9 Mg/ha/year","6 Mg/ha/year","3 Mg/ha/year")),
         lines = TRUE, col = c("purple","blue","green")))
dev.off()

## Decline in OM
om.9d <- rowSums(mat.9d[,5:8]) * 1e-3 * (1/420) * 2 * 100
pdf("./om-decline-century.pdf")
xyplot(om.9d
       ~ I(1:nsim/52),type="l",ylab="organic matter (%)",
       xlab = "years",
       col = "purple",
       lty=1, lwd=3,
       key=list(text=list(c("0.6 Mg/ha/year")),
         lines = TRUE, col = c("purple")))
dev.off()

## Plotting the plant litter
xyplot(mat2[,1] + mat2[,2] + mat2[,3] + mat2[,4] ~ 1:nsim,
       type="l",ylab="carbon",col=1:4,
       key=list(text=list(c("Leaf","Stem","Root","Rhiz")),lines=TRUE,col=1:4))

xyplot(respVec ~ 1:nsim , ylab="Respiration")

xyplot(MinNVec ~ 1:nsim , ylab="Mineralized N") 


## Plotting the different nitrogen carbon pools
labs <- paste("SN",1:9,sep="")
xyplot(mat3[,1] + mat3[,2] + mat3[,3] +
       mat3[,4] + mat3[,5] + mat3[,6] +
       mat3[,7] + mat3[,8] + mat3[,9] ~ 1:nsim,type="l",ylab="nitrogen",
       col=rainbow(9),
       lty=1, lwd=3,
       key=list(text=list(labs),lines=TRUE,col=rainbow(9)))


## Testing another aspect of Century
data(weather05)

litter <- c(0,0,0,0)
iC <- rep(0.5,9) # These are in Mg ha^-1
nyears <- 1e4
Cmat <- matrix(ncol=9,nrow=nyears)
litterMat <- matrix(ncol=4, nrow=nyears)

iRhiz <- 0.06
Rhiz <- numeric(nyears)
yield <- numeric(nyears)

sene.ll <- seneParms(senRoot = 3500, senRhizome = 3500)
photo.ll <- photoParms(vmax = 8, alpha = 0.025)

for(i in 1:nyears){

  century.ll <- centuryParms(SC1 = iC[1], SC2 = iC[2], SC3 = iC[3], SC4 = iC[4], SC5 = iC[5], SC6 = iC[6], SC7 = iC[7], SC8 = iC[8], SC9 = iC[9], Litter = litter)
  res <- BioGro(weather05, centuryControl = century.ll,
                iRhizome = iRhiz, seneControl = sene.ll,
                photoControl = photo.ll)
  last <- length(res$Stem)
  litter <- c(res$AboveLitter[last]*0.3, res$AboveLitter[last]*0.7, res$BelowLitter[last]*0.7, res$BelowLitter[last]*0.3) * 0.45 ## Biomass is 45% carbon
  litterMat[i,] <- litter
  iC <- res$SCpools
  Cmat[i,] <- iC
  iRhiz <- res$Rhizome[last]
  Rhiz[i] <- iRhiz
  yield[i] <- max(I(res$Stem + res$Leaf))
  
}

Cmat <- rbind(rep(0,9),Cmat)
litterMat <- rbind(rep(0,4),litterMat)

matplot(1:I(nyears+1),Cmat, type="l")

xyplot(Cmat[,5] + Cmat[,6] + Cmat[,7] + Cmat[,8] ~ 1:I(nyears+1),
       type = 'l', auto.key=TRUE)

xyplot(I(Cmat[,5] + Cmat[,6] + Cmat[,7] + Cmat[,8]) ~ 1:I(nyears+1),
       type = 'l', auto.key=TRUE,
       ylab = "Soil Carbon (Mg/ha)")

matplot(1:I(nyears+1),litterMat, type="l", ylab="litter")

xyplot(c(0.06,Rhiz) ~ c(0,1:nyears), type = "o")

xyplot(yield ~ 1:nyears , ylab = "Dry biomass", xlab = "years", type="o")
