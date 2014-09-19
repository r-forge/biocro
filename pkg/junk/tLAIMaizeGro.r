## Testing the different models for caluclating LAI in the Maize model

data(weather05)

## The default method is based on thermal time alone
laiP0 <- laiParms(lai.method = "TT", max.lai = 10)
ans0 <- MaizeGro(weather05, plant.day = 110, emerge.day = 120, harvest.day=300,
                 laiControl = laiP0)

## The spla method is based on specific leaf area
laiP1 <- laiParms(lai.method = "spla", max.lai = 10)
ans1 <- MaizeGro(weather05, plant.day = 110, emerge.day = 120, harvest.day=300,
                 laiControl = laiP1)

## The spla method is based on specific leaf area
laiP2 <- laiParms(lai.method = "ind-leaf-Lizaso", max.lai = 10, Aex = 700, LLx = 1000)
ans2 <- MaizeGro(weather05, plant.day = 110, emerge.day = 120, harvest.day=300,
                 laiControl = laiP2)

## Birch-Discontinuos
laiP3 <- laiParms(lai.method = "Birch-Discontinuous", max.lai = 10)
ans3 <- MaizeGro(weather05, plant.day = 110, emerge.day = 120, harvest.day=300,
                 laiControl = laiP3)

## Birch-continuos
laiP4 <- laiParms(lai.method = "Birch-Continuous", max.lai = 10)
ans4 <- MaizeGro(weather05, plant.day = 110, emerge.day = 120, harvest.day=300,
                 laiControl = laiP4)

xyplot(ans0$LAI + ans1$LAI + ans2$LAI +
       ans3$LAI + ans4$LAI ~ ans0$TTTc,
       auto.key=TRUE, type = 'l', ylim = c(0,7),
       ylab = "LAI")

matplot(ans2$TTTc, ans2$LAImat, type = 'l')

dat <- data.frame(ThermalT = ans1$TTTc, Stem = ans1$Stem, Leaf = ans1$Leaf, Root = ans1$Root, Grain = ans1$Grain, LAI = ans1$LAI)

ss <- sample(length(dat$ThermalT[dat$ThermalT > 5]), 10)

dat2 <- dat[c(3500,3700,4000,5000,5500),]

plot(ans1, dat2)
