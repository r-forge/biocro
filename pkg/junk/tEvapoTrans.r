## Testing the EvapoTrans function through the CanA function

data(weather05)
doy200 <- weather05[weather05$doy == 200,]

dat2 <- NULL
tmp2 <- matrix(ncol=5,nrow=24)
tmp3 <- matrix(ncol=4,nrow=24)
layers <- 10
lai <- 5

for(i in 1:24){

  doy <- doy200[i,2]
  hr  <- doy200[i,3]
  solar <- doy200[i,4]
  temp <- doy200[i,5]
  rh <- doy200[i,6]
  ws <- doy200[i,7]

##       tmp1 <- CanA(lai,doy,hr,solar,temp,rh,ws, StomataWS=1, nlayers=layers, lnControl=lnParms(LeafN=2,lnFun="linear",kpLN=0))
       tmp1 <- CanA(lai = lai,doy,hr,solar,temp,rh,ws, StomataWS=1, nlayers=layers, lnControl=lnParms(LeafN=3.7,lnFun="linear"), chi.l=0.5)

  tmp2[i,1] <- tmp1$CanopyAssim
  tmp2[i,2] <- tmp1$CanopyTrans
  tmp2[i,3] <- tmp1$TranEpen
  tmp2[i,4] <- tmp1$TranEpries
  tmp2[i,5] <- tmp1$CanopyCond
  
  dat1 <- data.frame(hour=i,layer=1:layers, as.data.frame(tmp1$LayMat))
  
  dat2 <- rbind(dat2,dat1)
}

transp <- data.frame(tmp2)
names(transp) <- c("CA","CT","Epen","Epries","Ccond")
## The returned values are in kg/m2/hr how about converting this to mm?
## kg to Mg multiply by 1e-3
## m/hr to mm/hr multiply by 1e3

sum(transp$CT)
sum(transp$Epen)
sum(transp$Epries)

xyplot(CT + Epen + Epries ~ 0:23, data = transp, type = 'l', auto.key=TRUE)
