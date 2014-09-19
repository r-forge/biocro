## Testing the soil multilayer model

data(weather05)
soilP <- soilParms(soilLayers = 5, hydrDist=TRUE, soilDepth = 2)
ans <- BioGro(weather05, day1 = 90, dayn=350, soilControl = soilP)

plot(ans)
plot(ans, plot.kind="SW")
plot(ans, plot.kind="ET")
plot(ans, plot.kind="cumET")
