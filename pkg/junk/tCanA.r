## Simple example
nlay <- 10
res <- CanA(lai=3, doy=200, hr=12,
            solar=1500, temp=25, rh=0.7, ##photoControl=photoParms(b1=7),
            windspeed=2, nlayers=nlay)

res$CanopyAssim
res$CanopyTrans

##  2 layers = 0.00347
##  3 layers = 0.00378
##  4 layers = 0.00394
## 10 layers = 0.00422
## 20 layers = 0.00431
## 40 layers = 0.00435
## 50 layers = 0.00435

## Example for a full day by layer
data(boo14.200)

lai <- 5
nlay <- 10
chi.l <- 1
lat <- 42
tmp2 <- NULL

for(i in 1:24){
  doy <- boo14.200[i,2]
  hr  <- boo14.200[i,3]
  solar <- boo14.200[i,4]
  temp <- boo14.200[i,5]
  rh <- boo14.200[i,6]
  ws <- boo14.200[i,7]

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

xyplot(Leafsun + Leafshade ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Leaf Area (m2/m2)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(AssimSun + AssimShade ~ hour | factor(layers), type='l',
       xlab = "Layers", ylab="Assimilation (micro mol/m2/s)",
       main=ttle, 
       auto.key=TRUE, data = as.data.frame(tmp2))

xyplot(I(TransSun*Leafsun) + I(TransShade*Leafshade) ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Transpiration (kg/m2/hr)",       
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")),
       data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(I(AssimSun*Leafsun) + I(AssimShade*Leafshade) ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Assimilation (micro mol /m2 ground /s)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(DeltaSun + DeltaShade ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Delta temperature",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))

xyplot(CondSun + CondShade ~ hour | factor(layers), type='l',
       xlab = "hour", ylab="Conductance (mmol/m2/s)",
       main=ttle, 
       key=simpleKey(text=c("Sun","Shade")), data = as.data.frame(tmp2), layout=c(nlay,1))

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

## Code to test transpiration
data(boo14.200)
dat2 <- NULL
tmp2 <- matrix(ncol=5,nrow=24)
tmp3 <- matrix(ncol=4,nrow=24)
layers <- 50
lai <- 5
doy <- 200
photoP <- photoParms(b1=7)
## stress <- c(rep(1,13),rep(0.5,5),rep(1,6))

for(i in 1:24){
         
    hr  <- boo14.200[i,3]
    solar <- boo14.200[i,4]
    temp <- boo14.200[i,5]
    rh <- boo14.200[i,6]
    ws <- boo14.200[i,7]
    
    tmp1 <- CanA(lai,doy,hr,solar,temp,rh,ws, ## stress=stress[i],
                 nlayers=layers, photoControl=photoP)

    tmp2[i,1] <- tmp1$CanopyAssim
    tmp2[i,2] <- tmp1$CanopyTrans
    tmp2[i,3] <- tmp1$TranEpen
    tmp2[i,4] <- tmp1$TranEpries
    tmp2[i,5] <- tmp1$CanopyCond
        
    dat1 <- data.frame(hour=i,layer=1:layers, as.data.frame(tmp1$LayMat))

    dat2 <- rbind(dat2,dat1)
}

tmp2 <- as.data.frame(tmp2)

names(tmp2) <- c("CanopyAssim","CanopyTrans",
                 "TranPen","TranEpries","CanopyCond")

xyplot(CanopyTrans + TranPen + TranEpries ~ 1:24, data = tmp2,
       type='o',
       key= simpleKey(text=c("Penman-Monteith","Penman","Priestly"),
           lines=TRUE, points=FALSE),
           xlab='hour',
           ylab="Transpiration (mm/h)")

xyplot(CanopyCond ~ 1:24, data = tmp2,
       type='o',
       xlab='hour',
       ylab="Conductance (mmol/m2/s)")

apply(tmp2[,2:4], 2, sum)



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



