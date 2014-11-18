## Simple example
nlay <- 50
res <- CanA(lai=3, doy=200, hr=12,
            solar=1500, temp=25, rh=0.7,
            windspeed=2, nlayers=nlay)

res$CanopyAssim

##  2 layers = 0.00347
##  3 layers = 0.00378
##  4 layers = 0.00394
## 10 layers = 0.00422
## 20 layers = 0.00431
## 40 layers = 0.00435
## 50 layers = 0.00435

## Example for a full day by layer
data(weather05)
doy200 <- weather05[weather05$doy == 200,]

lai <- 5
nlay <- 5
chi.l <- 1
lat <- 42
tmp2 <- NULL

for(i in 1:24){
  doy <- doy200[i,2]
  hr  <- doy200[i,3]
  solar <- doy200[i,4]
  temp <- doy200[i,5]
  rh <- 0.5
  ws <- 4

  tmp <- CanA(lai,doy,hr,solar,temp,rh,ws,
              nlayers=nlay,chi.l=chi.l,
              lat = lat)$LayMat

  tmp <- cbind(hour=hr, layers=1:nlay,tmp)
  tmp2 <- rbind(tmp2,tmp)
     
}

ttle <- paste("LAI = ",lai,
              "   layers = ",nlay,
              "   chi.l = ",chi.l,
              "   lat = ",lat, sep="")

xyplot(IDir + IDiff ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Irradiance (micro mol/m2/s)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))

## png("../figs/LeafSun-Shade.png")
xyplot(Leafsun + Leafshade ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Leaf Area (m2/m2)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))
## dev.off()

xyplot(AssimSun + AssimShade ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Assimilation (micro mol/m2/s)",
       main=ttle, 
       auto.key=TRUE, data = as.data.frame(tmp2))

## png("../figs/TransSun-Shade.png")
xyplot(I(TransSun*Leafsun) + I(TransShade*Leafshade) ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Transpiration (kg/m2/hr)",       
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))
## dev.off()

## png("../figs/AssimSun-Shade.png")
xyplot(I(AssimSun*Leafsun) + I(AssimShade*Leafshade) ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Assimilation (micro mol /m2 ground /s)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))
## dev.off()

## png("../figs/DeltaSun-Shade.png")
xyplot(DeltaSun + DeltaShade ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Delta temperature",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))
#dev.off()

## png("../figs/CondSun-Shade.png")
xyplot(CondSun + CondShade ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Conductance (mmol/m2/s)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))
## dev.off()

## png("../figs/CondSun-Shade.png")
xyplot(CondSun + CondShade ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Conductance (mmol/m2/s)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))
## dev.off()

xyplot(RH ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Relative Humidity",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(WindS ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Wind speed (m/s)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(CanHeight ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Canopy Height (m)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))




## Example for a full day
data(weather05)
doy200 <- weather05[weather05$doy == 200,]

lai <- 3
nlay <- 8
tmp2 <- NULL

for(i in 1:24){
  doy <- doy200[i,2]
  hr  <- doy200[i,3]
  solar <- doy200[i,4]
  temp <- doy200[i,5]
  rh <- 0.5
  ws <- 4

  tmp <- CanA(lai,doy,hr,solar,temp,rh,ws,nlayers=nlay,chi.l=1)$LayMat

  tmp <- cbind(hour=hr, layers=1:nlay,tmp)
  tmp2 <- rbind(tmp2,tmp)
     
}

## xyplot(TransSun + TransShade ~ hour | factor(layers), type='l',
##        xlab = "Layers", ylab="Transpiration",
##        auto.key=TRUE, data = as.data.frame(tmp2))

xyplot(Leafsun + Leafshade ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Leaf Area (m2/m2)",
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(AssimSun + AssimShade ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Assimilation (micro mol/m2/s)",
       auto.key=TRUE, data = as.data.frame(tmp2))

xyplot(I(TransSun*Leafsun) + I(TransShade*Leafshade) ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Transpiration (kg/m2/hr)",
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(I(AssimSun*Leafsun) + I(AssimShade*Leafshade) ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Assimilation (micro mol /m2 ground /s)",
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(DeltaSun + DeltaShade ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Delta temperature",
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(CondSun + CondShade ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Conductance (mmol/m2/s)",
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))







data(doy124)
     dat2 <- NULL
     tmp2 <- matrix(ncol=5,nrow=24)
     tmp3 <- matrix(ncol=4,nrow=24)
     layers <- 10
     for(i in 1:24){
         
        lai <- doy124[i,1]
        doy <- doy124[i,3]
        hr  <- doy124[i,4]
      solar <- doy124[i,5]
       temp <- doy124[i,6]
         rh <- doy124[i,7]
         ws <- doy124[i,8]

##       tmp1 <- CanA(lai,doy,hr,solar,temp,rh,ws, StomataWS=1, nlayers=layers, lnControl=lnParms(LeafN=2,lnFun="linear",kpLN=0))
       tmp1 <- CanA(lai,doy,hr,solar,temp,rh,ws, StomataWS=1, nlayers=layers, lnControl=lnParms(LeafN=3.7,lnFun="linear"), chi.l=0.5)

        tmp2[i,1] <- tmp1$CanopyAssim
       tmp2[i,2] <- tmp1$CanopyTrans
       tmp2[i,3] <- tmp1$TranEpen
       tmp2[i,4] <- tmp1$TranEpries
       tmp2[i,5] <- tmp1$CanopyCond
        
       dat1 <- data.frame(hour=i,layer=1:layers, as.data.frame(tmp1$LayMat))

       dat2 <- rbind(dat2,dat1)
     }

tmp2 <- as.data.frame(tmp2)

names(tmp2) <- c("CanopyAssim","CanopyTrans","TranPen","TranEpries","CanopyCond")

xyplot(CanopyTrans + TranPen + TranEpries ~ 1:24, data = tmp2,
       type='l', auto.key=TRUE)


## Plot of Irradiance for the 10 layers
##png("old-sunML.png")
xyplot(IDir + IDiff ~ hour | factor(layer),type="o",
     data = dat2, xlab="hour",layout=c(2,layers/2),col=c("blue","green"),lwd=1.5,
        ylab=expression(paste("Irradiance (",mu,"mol ",m^-2," ",s^-1,")")))
##dev.off()

## Plot of TempDiff for the 10 layers
#pdf("LeafTemp.pdf")
xyplot(DeltaSun + DeltaShade ~ hour | factor(layer),type="o",
      data = dat2, xlab="hour",layout=c(2,layers/2),col=c("blue","green"),lwd=1.5,
            ylab="Delta temperature (Celsius)")
#dev.off()

## Plot of Leaf area (sunlit and shaded) for the 10 layers
xyplot(Leafsun + Leafshade ~ hour | factor(layer),type="o",
      data = dat2, xlab="hour",layout=c(2,layers/2),col=c("blue","green"),lwd=1.5,
            ylab=expression(paste("Leaf Area (",m^2," ",m^-2,")")))


## Plot of Transpiration for the 10 layers
xyplot(TransSun + TransShade ~ hour | factor(layer),type="o",
     data = dat2, xlab="hour",layout=c(2,layers/2),col=c("blue","green"),lwd=1.5,
           ylab=expression(paste("Transpiration (mm ",H[2],"O ",m^-2," ",s^-1,")")))


## Plot of Assimilation for the 10 layers
   xyplot(AssimSun + AssimShade ~ hour | factor(layer),type="o",
      data = dat2, xlab="hour",layout=c(2,layers/2),col=c("blue","green"),lwd=1.5,
            ylab=expression(paste("Assimilation (",mu,"mol ",m^-2," ",s^-1,")")))

     ## Plot of Conductance for the 10 layers
     xyplot(CondSun + CondShade ~ hour | factor(layer),type="o",
      data = dat2, xlab="hour",layout=c(2,layers/2),col=c("blue","green"),lwd=1.5,
            ylab=expression(paste("Conductance (mmol ",m^-2," ",s^-1,")")))

     ## What does relative humidity look like?
     xyplot(RH ~ layer | factor(hour), data = dat2, auto.key=TRUE, type="l")

     ## What does the wind speed look like?
     xyplot(WindS ~ layer, data = dat2, auto.key=TRUE, type="a")

## Testing the effect of N distribution
dat2.no <- dat2

LeafN.no <- dat2.no$LeafN
LeafN.li <- dat2$LeafN
## Plot of Leaf Nitrogen
##pdf("LeafNitrogen.pdf")
     xyplot(LeafN.no + LeafN.li ~  layer,type="o", subset = hour == 12,
      data = dat2, xlab="layer",col=c("blue","green"),lwd=1.5,
            ylab=expression(paste("Leaf Nitrogen (g ",m^-2,")")))
##dev.off()



Vmax.no <- dat2.no$Vmax
Vmax.li <- dat2$Vmax
## Plot of Vmax
##pdf("Vmax.pdf")
     xyplot(Vmax.no + Vmax.li ~  layer,type="o", subset = hour == 12,
      data = dat2, xlab="layer",col=c("blue","green"),lwd=1.5,
            ylab=expression(paste("Vmax (",mu,"mol ",m^-2," ",s^-1,")")))
##dev.off()


## Let's do some math 2g m^-2 times 8 = 16 g total
## This gives 
(Atot.no <- sum(dat2.no$AssimSun + dat2.no$AssimShade)) ## 1002
## How do I distribute the same ammount of N more efficiently?
## Let us say I start with 3g m^-2
sum(subset(dat2,hour==12)$LeafN)*(8/10) ## 16 the same total ammount of N
(Atot.li <- sum(dat2$AssimSun + dat2$AssimShade)) ## 1160



## Plot of Assimilation for the 10 layers
##pdf("Assim.pdf")
xyplot(I(dat2$AssimSun + dat2$AssimShade) +
       I(dat2.no$AssimSun + dat2.no$AssimShade) ~ hour | factor(layer),type="o",
       data = dat2, xlab="hour",layout=c(2,layers/2),col=c("blue","green"),lwd=1.5,
       ylab=expression(paste("Assimilation (",mu,"mol ",m^-2," ",s^-1,")")),
       key=list(text=list(c("exp","const")),lines=TRUE,points=TRUE,
         col=c("blue","green"),type="o",pch=21))
##dev.off()

## Testing the effect of chi.l

     data(doy124)
     tmp <- numeric(24)
     tmp2 <- numeric(24)
     
     for(i in 1:24){
        lai <- doy124[i,1]
        doy <- doy124[i,3]
        hr  <- doy124[i,4]
      solar <- doy124[i,5]
       temp <- doy124[i,6]
         rh <- doy124[i,7]
         ws <- doy124[i,8]
     
       tmp[i] <- CanA(lai,doy,hr,solar,temp,rh,ws,chi.l=1)$CanopyAssim
       tmp2[i] <- CanA(lai,doy,hr,solar,temp,rh,ws,chi.l=0.01)$CanopyAssim
     
     }

png("new-canopy-photosynthesis-chi.l=0.01-chi=1.png")
     xyplot(tmp + tmp2 ~ c(0:23),
            auto.key=TRUE,
            type="l",lwd=2,
            xlab="Hour",
            ylab=expression(paste("Canopy assimilation (kg  ",
                m^-2," ",h^-1,")")))
dev.off()


## What is the value of evaporation for a hot summer day?
data(weather05)
doy200 <- weather05[weather05$doy == 200,]

lai <- 5
nlayers <- 45
trns <- numeric(24)
tpen <- numeric(24)
tpries <- numeric(24)


for(i in 1:24){
  doy <- doy200[i,2]
  hr  <- doy200[i,3]
  solar <- doy200[i,4]
  temp <- doy200[i,5]
##  rh <- doy200[i,6]
##  ws <- doy200[i,7]
  rh <- 0.5
  ws <- 4

  trns[i] <- CanA(lai,doy,hr,solar,temp,rh,ws,chi.l=1, nlayers=nlayers)$CanopyTrans
  tpen[i] <- CanA(lai,doy,hr,solar,temp,rh,ws,chi.l=1, nlayers=nlayers)$TranEpen
  tpries[i] <- CanA(lai,doy,hr,solar,temp,rh,ws,chi.l=1, nlayers=nlayers)$TranEpries
     
}

sum(trns)
sum(tpen)
sum(tpries)



LayerWindSpeed <- 2
kappa <- 0.41
WindSpeedHeight <- 2
CanopyHeight <- 1
ZetaCoef = 0.026;
ZetaMCoef = 0.13;
dCoef <- 0.77
Zeta = ZetaCoef * CanopyHeight;
	Zetam = ZetaMCoef * CanopyHeight;
	d = dCoef * CanopyHeight;

	ga0 = kappa^2 * LayerWindSpeed;
	ga1 = log((WindSpeedHeight + Zeta - d)/Zeta);
	ga2 = log((WindSpeedHeight + Zetam - d)/Zetam);
	(ga = ga0/(ga1*ga2))

## ws = 0.5, ga = 0.00923
## ws = 1, ga = 0.0185
## ws = 2, ga = 0.0369
