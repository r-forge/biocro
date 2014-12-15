## First I need to figure out the units in Century
## Let's say that I want CenturyR to run on g/m2
## I know that the inputs might be about 1-200 g/m2
## of biomass and that a soil with 5% OM with 0.5 of 
## that being C results in 9750 g/m2 for the first 30cm
## With a distribution of 0.01, 0.19, 0.8 into f,s,p pools
## means 7800 g/m2 p, 1850 g/m2 s and 97.5 g/m2 f

centP <- centuryParms(SC6=97.5, SC7=1850, SC8=7800)
res <- Century(200,200,200,200, 0.3, 25, 2, 0, centuryControl = centP)

## The result of 1.38 g/m2 of soil respiration mean that
## in a day this was emitted from the soil. In one second
## 1.38 / (24*60*60) were emitted. This is in mol equivalents
## (44 g/mol) 1.38 / (24*60*60) / 44. Now in micro mols
## 1.38 / (24*60*60) * 1e6 result in 0.363 micro/mol/m2/s
## This is ok as I do not have all the pools stabilized here.
## What if I run the model for a long time?

centP <- centuryParms(SC6=97.5, SC7=1850, SC8=7800, timestep='week')

nsim <- 2e4

ans <- numeric(nsim)

LL <- 100
SL <- 100
RL <- 100
RhL <- 100

for(i in 1:nsim){

  res <- Century(LL,SL,RL,RhL, 0.2, 20, 2, 0, centuryControl = centP)
  ans[i] <- res$Resp  / (7*24*60*60) * (1/44) * 1e6
  centP[1:8] <- res$SCs[1:8]
  if(nsim %% 52 == 0){
      LL <- 100
      SL <- 100
      RL <- 100
      RhL <- 100    
  }else{
      LL <- 0
      SL <- 0
      RL <- 0
      RhL <- 0
  }
}

xyplot(ans ~ 1:nsim, ylab = "soil resp (micro mol/m2/s)", type='l')

## MaizeGro vs BioGro
data(weather05)
res <- BioGro(weather05, centuryControl = centuryParms(om=0.5))

## Convert Mg to g, ha to m^2, hr to sec, g to mol, mol to micro mol
plot(res$RespVec * 1e6 * 1e-4 * (1/3600) * (1/44) * 1e6,
     ylab = "CO2 efflux (micro mol/m2/s)") ## This is micro mol CO2 /m2/s

res <- MaizeGro(weather05,
                plant.day=120,
                emerge.day=130,
                harvest.day=280,
                laiControl = laiParms(lai.method="ind-leaf-Lizaso"),
                centuryControl = centuryParms(SC1=0,SC2=0,
                    SC3=0,SC4=0,SC5=0,om=5))


plot(res$mRespVec * 1e6 * 1e-4 * (1/3600) * (1/44) * 1e6,
     ylab = "CO2 efflux (micro mol/m2/s)") ## This is micro mol CO2 /m2/s)


## Test script for the Century Model

iSoilP0 <- iSoilP

iSoilP <- centuryParms(timestep="week")
iSoilP[1:9] <- iSoilP0[1:9]

#nsim <- 52*1000 ## weekly
nsim <- 52 * 10000 ## yearly

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
labs <- paste("SC",1:8,sep="")
xyplot(mat[,1] + mat[,2] + mat[,3] +
       mat[,4] + mat[,5] + mat[,6] +
       mat[,7] + mat[,8] ~ I(1:nsim/52),type="l",ylab="carbon",
       col=rainbow(8),
       lty=1, lwd=3,
       key=list(text=list(labs),lines=TRUE,col=rainbow(8)),
       main = paste("total annual biomass = 6 Mg/ha/yr"))

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
