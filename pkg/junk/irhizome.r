## What is the effect of the initial rhizome?

data(weather04)


ans <- BioGro(weather04, iRhizome = 0.06, irtl = 1e-30)
plot(ans)

ans <- BioGro(weather04, iRhizome = 0.6, irtl = 1e-30)
plot(ans)

## Still with a very small term that controls the change from rhizome
## to leaf the amount of rhizome at the end of the growing season is
## very high. This is a result of the very high translocation
## coefficient later in the growing season.


## Effect of initial rhizome on time to reach ceiling yield

irhiz <- 0.6
lgth <- 20
rhiz <- numeric(20)
ag <- numeric(20)
bg <- numeric(20)

wet <- weather04

wet$precip <- 0.99 * wet$precip

for(i in 1:lgth){

  if(i == 1) day1 <- 140
  else day1 <- 102
  
  res <- BioGro(wet, day1 = day1, iRhizome = irhiz, irtl = 1e-30,
                soilControl = soilParms(iWatCont = 0.27))

  irhiz <- max(res$Rhizome)

  rhiz[i] <- irhiz
  bg[i] <- max(res$Rhizome + res$Root)
  ag[i] <- max(res$Stem + res$Leaf) * 0.67

}

xyplot(rhiz + ag ~ 1:lgth, type = 'o', auto.key=TRUE)
  
xyplot(I(ag/bg) ~ 1:lgth, main = "Shoot:root ratio",
       ylim=c(0,2),
       panel = function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(h = 1)
       })
