
rm(list = ls())

# adjust the movement distances between each of the wetlands
movementType <- read.csv("C:/Workspace/movementType.csv", header = FALSE, sep = ",", dec = ".")
wetlandMovement <- read.csv("C:/Workspace/wetlandMovement.csv", header = FALSE, sep = ",", dec = ".")
## resistance modifiers for dispersal categories A-G *note F is a reserved term for FALSE*
A <- 0.5 # DS along river
B <- 0.75 # DS cross river
C <- 1 # Normal overland
D <- 1.25 # US along river
E <- 1.5 # US cross river
G <- 1 # 50 km +

## And a temporary holder for the seasonally adjusted matrix 
bass1 <- bass2 <- bass <- seq(1,nrow(wetlandMovement),1)


for (i in 1:nrow(movementType)) {
  bass1[i] <- paste("From",bass[i])
}
rownames(wetlandMovement) <- bass1

for (i in 1:nrow(wetlandMovement)) {
  bass2[i] <- paste("To",bass[i])
}
colnames(wetlandMovement) <- bass2
adj.wetlandMovement <- wetlandMovement
adj.wetlandMovement[] <- NA

### Create an array of distance between wetlands where the distance is adjusted to reflect the type of journey ie downstream easier/shorter, upstream harder/further.
for (i in 1:nrow(wetlandMovement)) {
  for (j in 1:ncol(wetlandMovement)) {
    if (i == j) {
      adj.wetlandMovement[i,j] <- NA
    } else {
      adj.wetlandMovement[i,j] <- wetlandMovement[i,j]*get(movementType[i,j])
      if (adj.wetlandMovement[i,j] > 50) {  adj.wetlandMovement[i,j] <- 50 }
    }
  }
}

# convert to metres
adj.wetlandMovement <- adj.wetlandMovement * 1000
# save the final movement data
write.csv(adj.wetlandMovement,'G:/My Drive/University Milestones/Cost benefit analysis 2_genetic algorithm/csvs for code/GA Meta adjMovement.csv')



## Levy flight Setup code
## calculate the equation to generate a levy flight model distance generator
min.dist <- 50
max.dist <- 50000
lmax.dist <- log(max.dist)
max.prob.mov <- 0.3
nmoves <- 100000
u <- simm.levy(1:nmoves, l0=min.dist, mu = 2.3)
hist.out <- hist(log(u[[1]]$dist),main="", xlab="log distance (m)", ylab="frequency", br=50)
range(u[[1]]$dist, na.rm=T)

lcounts <- log(hist.out$counts)
lprob.scaled <- log(max.prob.mov * hist.out$counts/max(hist.out$counts))
ldist <- hist.out$mids
ldt <- na.omit(data.frame(ldist,lcounts,lprob.scaled))
ldat.clean <- ldt[which(ldt$ldist <= lmax.dist), ]

if (length(which(is.infinite(ldat.clean$lcounts)==T)) > 0) {
  ldat.rem <- which(is.infinite(ldat.clean$lcounts)==T)
  ldat.clean <- ldat.clean[-ldat.rem,]
}

plot(ldat.clean$ldist, ldat.clean$lprob.scaled, pch=19, ylab="log scaled prob", xlab="log dist (m)")
fit.levy <- lm(lprob.scaled ~ ldist, data=ldat.clean)
abline(fit.levy, lty=2, col="red")
summary(fit.levy)
exp(range(ldat.clean$ldist, na.rm=T))

levy.int <- as.numeric(coef(fit.levy)[1])
levy.sl <-  as.numeric(coef(fit.levy)[2])

dist.vec <- seq(min.dist, max.dist, 1)
pr.pred <- exp(levy.int + (log(dist.vec))*levy.sl)
plot(dist.vec, pr.pred, type="l")

dist.rand <- sample(dist.vec, 1, replace=T, prob=pr.pred)
dist.rand
#hist(dist.rand)

dist.vec <- seq(min.dist, max.dist, 1)
pr.pred <- exp(levy.int + (log(dist.vec))*levy.sl)

# the line equation working backwards from distance is then pr.pred <- exp(levy.int + (log(dist.vec))*levy.sl)
# so we can rearrange the equation algebraically to 
 distance <- exp((ln(PrDist) - levy.int)/levy.sl)
 # maximum distance is reached at PrP <- 4.69e-05
 # minimum distance is reached at PrP <- 0.3559407
 # use PrDist <- runif(1,4.69e-05,0.3559407)
# levy.sl <- -1.294724
 # levy.int <- 4.04117
 
 Pr <- seq(4.69e-05,0.3559407,0.000001)
 dis <- exp((log(Pr) - levy.int)/levy.sl)
 plot(dis, Pr, type="l")
 

# Thus to generate a distance we use 
 levy.sl <- -1.294724
 levy.int <- 4.04117
 PrDist <- runif(1,4.69e-05,0.3559407)
 distance <- exp((log(PrDist) - levy.int)/levy.sl)
 
 # or for a single line
distance <- round(exp((log(runif(1,4.69e-05,0.3559407)) - levy.int)/levy.sl),0)
distance 