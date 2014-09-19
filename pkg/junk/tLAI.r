## Testing the leaf area functions

laifun <- function(time, A, k, tm){

  leaf.area <- A/(1 + exp(-k*(time - tm)))

  leaf.area

}

leafArea <- function(LN, Ax, LNx, a1, a2){

  .ans0 <- (LN - LNx)/(LNx - 1)
  ans <- Ax * exp(a1 * .ans0^2 + a2 * .ans0^3)
  ans

}

lnumb <- 1:20

a1 <- 20 * 0.25 - 10.61
a2 <- 20 * 0.27 - 5.99
leafA <- leafArea(lnumb, 650, 0.67*20, a1, a2)

xyplot(leafA ~ lnumb)
