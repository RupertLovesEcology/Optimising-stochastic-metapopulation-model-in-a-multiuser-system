##########################################

# create the density dependant survival inhibition loops

#####################################################################################################

## Remove everything
rm(list = ls())
library(deSolve)      # function solving for density dependence    

#create the egg density correcting function and variables
surv.mult.egg.up <- 1.03
surv.mult.egg.upmid <- 1.03
surv.mult.egg.mid <- 0.94
surv.mult.egg.lo <- 0.8
surv.mult.egg.lo.lo <- 0.2
surv.mult.egg.lo.lo.lo <- 0.05

K.egg.up <- 1
K.egg.upmid <- 0.98
K.egg.mid <- 0.90
K.egg.lo <- 0.7
K.egg.lo.lo <- 0.3
K.egg.lo.lo.lo <- 0.01

K.egg.vec <- c(K.egg.up,K.egg.upmid, K.egg.mid,K.egg.lo, K.egg.lo.lo, K.egg.lo.lo.lo)
surv.mult.egg.vec <- rev(c(surv.mult.egg.up, surv.mult.egg.upmid, surv.mult.egg.mid, surv.mult.egg.lo, surv.mult.egg.lo.lo, surv.mult.egg.lo.lo.lo))
DD.dat <- data.frame(K.egg.vec, surv.mult.egg.vec)

#the formula for the function
SS<-getInitial(surv.mult.egg.vec~SSlogis(K.egg.vec,alpha,xmid,scale),data=DD.dat)
fit.expd.egg <- nls(surv.mult.egg.vec ~ a/((exp((b-K.egg.vec)/c)) + 1),
                    data = DD.dat,
                    algorithm = "port",
                    start = c(a = as.numeric(SS["alpha"]), b = as.numeric(SS["xmid"]), c = as.numeric(SS["scale"])),
                    trace = TRUE,
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))

eggfunc <- function(x) (coef(fit.expd.egg)[1])/((exp((coef(fit.expd.egg)[2] - x)/coef(fit.expd.egg)[3]))+ 1)
intArea <- integrate(eggfunc,lower=0,upper=1.02)                    
area1 <- (as.numeric(intArea[1])/1.02)                    
s.mult.egg.iter <- as.numeric(coef(fit.expd.egg)[1])/((exp((coef(fit.expd.egg)[2] - 1.02)/coef(fit.expd.egg)[3]))+ 1) 

# density dependant inhibition of egg laying, use 
eggfunc <- function(x) (0.9991499)/((exp((0.9380068 - x)/-0.02672893))+ 1)
# and set 
s.mult.egg.iter <- 0.04442657 

## invoke a density-feedback function on tadpole survival to year 1
# density feedback survival multiplier for tadpoles hinges on the density of other tadpoles in the pond
#form of the curve from "Effect of Stocking Density on the Survival and Growth of 
# Hoplobatrachus occipitalis (Günther, 1858) (Amphibia: Dicroglossidae) 
# of Tadpoles Reared in Ponds from Benin, Godome, 2018, International Journal of Aquaculture"  &
# An Analysis of Density Effects and Predation in Bufo Americanus Tadpoles
# from Brockelman 1969
surv.mult.up <- 1.0
surv.mult.upmid <- 0.58
surv.mult.mid <- 0.19
surv.mult.lo <- 0.10

K.up <- 1
K.upmid <- 0.83
K.mid <- 0.45
K.lo <- 0.3

K.tad.vec <- c(K.up,K.upmid, K.mid,K.lo)
surv.mult.tad.vec <- rev(c(surv.mult.up, surv.mult.upmid, surv.mult.mid, surv.mult.lo))
plot(K.tad.vec, surv.mult.tad.vec, pch=19)

# Bleasdale
# y = (a + bx)^(-1/c)
DD.dat <- data.frame(K.tad.vec, surv.mult.tad.vec)
param.init <- c(-2.41e-01, 1.54, 1.17)
fit.expd.tad <- nls(surv.mult.tad.vec ~ (a + (b*K.tad.vec))^(-1/c), 
                    data = DD.dat,
                    algorithm = "port",
                    start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                    trace = TRUE,      
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
plot(K.tad.vec, surv.mult.tad.vec, pch=19, xlab="K", ylab="reduction Tadpole survival to 1 yr")
K.pred.tad.vec <- seq(K.lo,1,0.01)
pred.surv.tad.mult <- (as.numeric(coef(fit.expd.tad)[1]) + (K.pred.tad.vec * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
lines(K.pred.tad.vec, pred.surv.tad.mult, lty=2, lwd=1, col="red")

# for density dependant inhibition of tadpole survival, use 
(as.numeric(0.8227) + (K.pred.tad.vec * as.numeric(0.5914)))^(-1/as.numeric(0.1579))
# to calculate s.mult.iter.tad

## invoke a density-feedback function on Juvenile survival from year 1 to year 2
# density feedback survival multiplier for juveniles hinges on the density of themselves (but more strongly than adults)
### Make an approximation of the curve, see work by Berven on density dependence 
surv.mult.up <- 1.1
surv.mult.upmid <- 1
surv.mult.midmidup <- 0.5
surv.mult.mid <- 0.3
surv.mult.midlo <- 0.19
surv.mult.lo <- 0.10

K.up <- 1
K.upmid <- 0.89
K.mid <- 0.75
K.midlo <- 0.6
K.midlolo <- 0.3
K.lo <- 0.1

K.juv.vec <- c(K.up,K.upmid, K.mid,K.midlo,K.midlolo,K.lo)
surv.mult.juv.vec <- rev(c(surv.mult.up, surv.mult.upmid,surv.mult.midmidup,surv.mult.mid,surv.mult.midlo,surv.mult.lo))
plot(K.juv.vec, surv.mult.juv.vec, pch=19,main = "the curve I want to emulate")

D.dat <- data.frame(K.juv.vec, surv.mult.juv.vec)

# use the deSolve package to determine your starting parameters for the nls function
# see more here     https://datascienceplus.com/first-steps-with-non-linear-regression-in-r/
SS<-getInitial(surv.mult.juv.vec~SSlogis(K.juv.vec,alpha,xmid,scale),data=D.dat)

#the formula for the function
fit.expd.juv <- nls(surv.mult.juv.vec ~ a/((exp((b-K.juv.vec)/c)) + 1),
                    data = D.dat,
                    algorithm = "port",
                    start = c(a = as.numeric(SS["alpha"]), b = as.numeric(SS["xmid"]), c = as.numeric(SS["scale"])),
                    trace = TRUE,
                    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))

fit.expd.adult <- nls(surv.mult.juv.vec ~ a/((exp((b-K.juv.vec)/c)) + 1),
                      data = D.dat,
                      algorithm = "port",
                      start = c(a = as.numeric(SS["alpha"]), b = as.numeric(SS["xmid"]), c = as.numeric(SS["scale"])),
                      trace = TRUE,
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))


# for density dependant inhibition of survival in adult life stages use 
as.numeric(1.2131120)/((exp((0.5508745 - K.rel.juv)/-0.1825602))+ 1)
# to calculate s.mult.iter.juv and s.mult.iter.ad











K.rel.juv <- 0.8

K.rel.juv <- 1

K.rel.juv <- 1.2

aa <- as.numeric(1.2131120)/((exp((0.5508745 - K.rel.juv)/-0.1825602))+ 1)

bb <- as.numeric(coef(fit.expd.juv)[1])/((exp((coef(fit.expd.juv)[2] - K.rel.juv)/coef(fit.expd.juv)[3]))+ 1)
aa
bb



a/((exp((b-K.juv.vec)/c)) + 1),


coefAd1 <- 1.2131120  
coefAd2 <- 0.5508745 
coefAd3 <- -0.1825602  


as.numeric(coef(fit.expd.juv)[1])/((exp((coef(fit.expd.juv)[2] - K.rel.juv)/coef(fit.expd.juv)[3]))+ 1)

coefAd1 <- 1.2131120  
coefAd2 <- 0.5508745 
coefAd3 <- -0.1825602 



coefEgg1 <- 0.9991499 
coefEgg2 <- 0.9380068 
coefEgg3 <- -0.02672893 

# this should work to replace the existing
eggfunc <- function(x) 0.9991499/((exp((0.9380068  - x)/-0.02672893))+ 1)
intArea <- integrate(eggfunc,lower=0,upper=1.02)                    
area1 <- (as.numeric(intArea[1])/1.02)     
s.mult.egg.iter <- as.numeric(0.9991499)/((exp((0.9380068 - 1.02)/-0.02672893))+ 1) 

# ??
coefTad1 <- 0.8227
coefTad2 <- 0.5914 
coefTad3 <- 0.1579 

# sdjvnsadijnvl;p

coefAd1 <- 1.2131120  
coefAd2 <- 0.5508745 
coefAd3 <- -0.1825602 

