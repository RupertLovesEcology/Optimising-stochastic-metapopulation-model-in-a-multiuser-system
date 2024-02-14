# this code analyses the combined outputs from desktop computers (LOFI GAMetapop V17.R) and high performance computer
# HPCGA_MetapopV19.r) 
# these outputs underpin the publication ""
# author Rupert Mathwin  rupert.mathwin.ecology@gmail.com

## Remove everything
rm(list = ls())


# load libraries
library(pls)
library(vip)
library(caret)
library(tidyverse)
library(tidyverse)
library(ggplot2)
library(vegan) 


# custom function to remove duplicates and sort by V40
TidyUp <- function(dat) {
  # remove duplicate rows
  dat <- dat[!duplicated(dat), ]
  
  # sort by V40
  dat <- dat[order(dat$V40, decreasing = TRUE),]
  return(dat)
}
## let's get some data!!
## OK Reminder, this is what I am working with:
## 86 numerics (criteria 1 - 16, svWetlands 1-23, svTotal 1, eWaterOther 1-23, eWaterFrog 1-23, eWaterTotal 1) 
dirNm <- 'G:/My Drive/University Milestones/Cost benefit analysis 2_genetic algorithm/csvs for analysis/'
solutionsAll <- read.csv(paste0(dirNm, 'iterSeed CombinedMASTER.csv'))
solutionsAll <- solutionsAll[ , -1]
solutionsAll <- TidyUp(solutionsAll)

names <- NULL  # then past in colnames 
for (i in 1:16) {
  names <- append(names, paste0('Weighting', i))
}
for (i in 1:23) {
  names <- append(names, paste0('svWL', i))
}
names <- append(names, 'svTotal')
for (i in 1:23) {
  names <- append(names, paste0('eW.other', i)) 
}
for (i in 1:23) {
  names <- append(names, paste0('eW.frog', i)) 
}
names <- append(names, 'eWTotal')
colnames(solutionsAll) <- names

# and add a column for all eWaterevents that weren't specifically for redgum or lignum (I'm calling them frogs) 
solutionsAll <- solutionsAll %>%   
  mutate(eW.frogTotal = rowSums(select(., 64:86)))

# and the total of eWaters for that wetland (other and frog)
for (i in 1:23) {
  solutionsAll[paste0("eW.both", i)] <- rowSums(solutionsAll[, c(paste0("eW.other", i), paste0("eW.frog", i))], na.rm = TRUE)
}

# let's transform the count of eWater events to the volume of water delivered
wLmeta <- read.csv(paste0(dirNm, 'GA wetlandMetaDataV4.csv'))

volTotal <- rep(0, nrow(solutionsAll))
volFrogTotal <- rep(0, nrow(solutionsAll))

for (i in 1:23) {  # this code can be improved if I have time
  for (j in 1:nrow(solutionsAll)) {
    volTotal[j] <- volTotal[j] + (solutionsAll[[paste0('eW.both', i)]][j] * wLmeta[i, 5])
    volFrogTotal[j] <- volFrogTotal[j] + (solutionsAll[[paste0('eW.frog', i)]][j] * wLmeta[i, 5])
  }
}

################################################################################
# some quick plots to visualise 
hist(solutionsAll$svTotal, breaks = 200, xlim = c(300, 800))


# eWater vs survival
plot(y = solutionsAll$eW.frogTotal, x = solutionsAll$svTotal, pch = 16, cex = 0.5, xlab = 'sum of surviving subpopulations', 
     ylab = 'sum of eWater events (frogs)')


plot(y = solutionsAll$eWTotal, x = solutionsAll$svTotal, pch = 16, cex = 0.5, xlab = 'sum of surviving subpopulations', 
     ylab = 'sum of all eWater events')

# and now the same with volume of water used
sV <- solutionsAll$svTotal
tempdf <- data.frame(sV, volTotal, volFrogTotal)

# some quick plots to visualise volume
plot(y = tempdf$volTotal, x = tempdf$sV, pch = 16, cex = 0.5, xlab = 'sum of surviving subpopulations', 
     ylab = 'vol of all eWater events')
lm_model <- lm(tempdf$volTotal ~ tempdf$sV)
abline(lm_model, col = "red", lwd = 2)  # Add the trendline to the plot

plot(y = tempdf$volFrogTotal, x = tempdf$sV, pch = 16, cex = 0.5, xlab = 'sum of surviving subpopulations', 
     ylab = 'vol of all eWater events (frogs)')
lm_model <- lm(tempdf$volFrogTotal ~ tempdf$sV)
abline(lm_model, col = "red", lwd = 2)  # Add the trendline to the plot

#####################################################################################################################
# lets do a quick check of how the variance changes with sample size
# subset just the parameters and the objective function 
solutions <- solutionsAll[, 1:16]
solutions$V40 <- solutionsAll$svTotal
colnames(solutions) <- c('Criteria1', 'Criteria2', 'Criteria3', 'Criteria4', 'Criteria5', 'Criteria6', 
                         'Criteria7', 'Criteria8', 'Criteria9', 'Criteria10', 'Criteria11', 'Criteria12', 
                         'Criteria13', 'Criteria14', 'Criteria15', 'Criteria16', 'V40')

# we'll start  in increments of 50 for the best and worst values
objBmin <- objWmin <- objBmax <- objWmax <- objBmn <- objWmn <- objBsd <- objWsd <- ranges <- seq(50, 5000, 50)

objBmin[] <- objWmin[] <- objBmax[] <- objWmax[] <- objBmn[] <- objWmn[] <- objBsd[] <- objWsd[] <- 0

for (i in 1:length(ranges)) {
  Range <- ranges[i] 
  best <- solutions[1:Range, 17]
  worst <- solutions[(nrow(solutions) - Range):nrow(solutions), 17] 
  objBmn[i] <- mean(best)
  objWmn[i] <- mean(worst)
  objBsd[i] <- sd(best)
  objWsd[i] <- sd(worst)
  objBmin[i] <- min(best)
  objWmin[i] <- min(worst)
  objBmax[i] <- max(best)
  objWmax[i] <- max(worst)
}

tempdf <- data.frame(ranges, objBmn, objWmn, objBsd, objWsd, objBmin, objBmax, objWmax, objWmin)

tempdf <- tempdf %>%
  mutate(Bup = objBmn + objBsd) %>% 
  mutate(Bdown = objBmn - objBsd) %>% 
  mutate(Wup = objWmn + objWsd) %>% 
  mutate(Wdown = objWmn - objWsd) 

plot(x = tempdf$ranges, y = tempdf$objBmn, main = 'best solutions', xlim = c(0, 900))
lines(x = tempdf$ranges, y = tempdf$Bup, col = 'red')
lines(x = tempdf$ranges, y = tempdf$Bdown, col = 'red')
points(x = tempdf$ranges, y = tempdf$objBmin, col = 'blue')
points(x = tempdf$ranges, y = tempdf$objBmax, col = 'blue')

plot(x = tempdf$ranges, y = tempdf$objWmn, main = 'worst solutions', xlim = c(0, 900))
lines(x = tempdf$ranges, y = tempdf$Wup, col = 'red')
lines(x = tempdf$ranges, y = tempdf$Wdown, col = 'red')
points(x = tempdf$ranges, y = tempdf$objWmin, col = 'blue')
points(x = tempdf$ranges, y = tempdf$objWmax, col = 'blue')

# based on this I could use the best and worst 600 solutions for comparisons
# 850 also looks OK. 1100 still produces strong patterns 



## all of the individual criteria weights max mean etc
## not sure if I want this, storing for now
meanW <- meanB <- matrix(data = 0, nrow = length(ranges), ncol = 16)
sdW <- sdB <- matrix(data = 0, nrow = length(ranges), ncol = 16)
maxW <- maxB <- matrix(data = 0, nrow = length(ranges), ncol = 16)
minW <- minB <- matrix(data = 0, nrow = length(ranges), ncol = 16)
absW <- absB <- matrix(data = 0, nrow = length(ranges), ncol = 16)

for (i in 1:length(ranges)) {
  Range <- ranges[i] 
  best <- solutionsAll[1:Range, 1:16]
  worst <- solutionsAll[(nrow(solutionsAll) - Range):nrow(solutionsAll), 1:16] 
 
  for (j in 1:16) {
    meanB[i, j] <- mean(best[ ,j])
    meanW[i, j] <- mean(worst[ ,j])
    sdB[i, j] <- sd(best[ ,j])
    sdW[i, j] <- sd(worst[ ,j])
    maxB[i, j] <- max(best[ ,j])
    maxW[i, j] <- max(worst[ ,j])
    minB[i, j] <- min(best[ ,j])
    minW[i, j] <- min(worst[ ,j])
    absB[i, j] <- abs(maxB[i,j] - minB[i,j]) 
    absW[i, j] <- abs(maxW[i,j] - minW[i,j]) 
  }
} 

################################################################################
## First we will compare the suitability of partial least squares regression and principal component regression to evaluate 
# the importance of the 16 Criteria with the sum of surviving wetlands produced in the model 
# it borrows from https://www.statology.org/partial-least-squares-in-r/ and https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
## and selects which of the models is the best fit 
###############################################
## OK Let's start with partial least squares regression selecting the number  of components for inclusion using both the 
# one sigma heuristic (the model with fewest components that is still less than one standard error away from the overall best model)
# and the randomisation approach

##  run the partial least squares regression model using cross validation
GA.plsrCV <- plsr(V40 ~ Criteria1 + Criteria2 + Criteria3 + Criteria4 + Criteria5 + Criteria6 + Criteria7 + Criteria8 +
                    Criteria9 + Criteria10 + Criteria11 + Criteria12 + Criteria13 + Criteria14 + Criteria15 + Criteria16,
                  data = solutions, scale = F, validation = "CV") 

# view summary of CV model 
summary(GA.plsrCV)

## view summary plots 
validationplot(GA.plsrCV)
validationplot(GA.plsrCV, val.type="MSEP")

# select the number of axes/principle components to include
GA.plsrCV.onesigma <- selectNcomp(GA.plsrCV, method = "onesigma", plot = TRUE)
GA.plsrCV.permut <- selectNcomp(GA.plsrCV, method = "randomization", plot = TRUE) 


## create the training dataset to evaluate the model fit
# Extract every 10th row into a new dataframe
indices <- seq(10, nrow(solutions), by = 10)
trainingSet <- solutions[indices, ]
testingSet <- solutions[-indices, 1:16]
targets <- solutions[-indices, 17]

# and run the new model with a seed
set.seed(22)
GA.plsrCVTEST <- plsr(V40 ~ Criteria1 + Criteria2 + Criteria3 + Criteria4 + Criteria5 + Criteria6 + Criteria7 + Criteria8 +
                        Criteria9 + Criteria10 + Criteria11 + Criteria12 + Criteria13 + Criteria14 + Criteria15 + Criteria16,
                      data = trainingSet, scale = F, validation = "CV") 


## use the recommended number of components to predict and get the root mean squared error (model fit)
set.seed(40)
# predict with the GA.plsrCVTEST model
plsr_predCV <- predict(GA.plsrCVTEST, testingSet, ncomp = 4) # For PLS one sigma use 4 components

#calculate RMSE for plsr using CV (one sigma)
plsrCV1sigmafit <- sqrt(mean((plsr_predCV - targets)^2))

## use the recommended number of components to predict and get the root mean squared error (model fit)
set.seed(40)
# predict with the GA.plsrCVTEST model
plsr_predCV <- predict(GA.plsrCVTEST, testingSet, ncomp = 5) # For PLS randomisation use 5 components

#calculate RMSE for plsr using CV (one sigma)
plsrCVrandfit <- sqrt(mean((plsr_predCV - targets)^2))

##########################################################################################
## Now repeat this process for principle component regression

##  run the principle component regression using 10-fold cross validation
GA.pcrCV <- pcr(V40 ~ Criteria1 + Criteria2 + Criteria3 + Criteria4 + Criteria5 + Criteria6 + Criteria7 + Criteria8 +
                  Criteria9 + Criteria10 + Criteria11 + Criteria12 + Criteria13 + Criteria14 + Criteria15 + Criteria16,
                data = solutions, scale = F, validation = "CV") 

# view summary of CV model 
summary(GA.pcrCV)

## view summary plots 
validationplot(GA.pcrCV)
validationplot(GA.pcrCV, val.type="MSEP")

GA.pcrCV.onesigma <- selectNcomp(GA.pcrCV, method = "onesigma", plot = TRUE)  
GA.pcrCV.permut <- selectNcomp(GA.pcrCV, method = "randomization", plot = TRUE) 


# test these models 
set.seed(22)
GA.pcrCVTEST <- pcr(V40 ~ Criteria1 + Criteria2 + Criteria3 + Criteria4 + Criteria5 + Criteria6 + Criteria7 + Criteria8 +
                      Criteria9 + Criteria10 + Criteria11 + Criteria12 + Criteria13 + Criteria14 + Criteria15 + Criteria16, 
                    data = trainingSet, scale = F, validation = "CV") 
set.seed(40)

# predict with the GA.pcrCVTEST model
pcr_predCV <- predict(GA.pcrCVTEST, testingSet, ncomp = 12) # 1 sigma

#calculate RMSE for PCR using CV and 1 sigma
pcrCV1sigmafit <- sqrt(mean((pcr_predCV - targets)^2))

set.seed(40)

# predict with the GA.pcrCVTEST model
pcr_predCV <- predict(GA.pcrCVTEST, testingSet, ncomp = 16) # randomisation

#calculate RMSE for PCR using CV and randomisation
pcrCVrandfit <- sqrt(mean((pcr_predCV - targets)^2))

########################################################################################
## now to finalise which model we like best
cat('the lineup is: \n plsrCV1sigmafit = ', plsrCV1sigmafit, '\n plsrCVrandfit = ', plsrCVrandfit, '\n pcrCV1sigmafit = ', pcrCV1sigmafit, '\n pcrCVrandfit = ', pcrCVrandfit)

## the best model fit is the plsr using 5 components 
## so now I have the best model what are the principal components and the contributions of each Criteria
bestModel <- GA.plsrCV
PCs <- 5

# OK shifting to vignette: https://cran.r-project.org/web/packages/pls/pls.pdf
plot(bestModel, ncomp = PCs, asp = 1, line = TRUE) # meh

## relative importance of each criteria on the first principle components
selectedPCcoef <- coef(bestModel, ncomp = bestModel$ncomp, comps = 1:PCs, intercept = FALSE)

# Plot the realtive importance of the criteria on the model
var_importance <- varImp(bestModel)  ## note to self: this is also plotted in 'Plotting GAMetapop V*.R'
var_importance$Criteria <- as.factor(rownames(var_importance))
var_importance$Criteria <- factor(var_importance$Criteria, levels = var_importance$Criteria[order(var_importance$Overall)])

varImpPlot <- ggplot(var_importance) +
  geom_col(aes(Overall, Criteria), fill = 'blue', width = 0.8) + 
  theme_bw() + 
  ylab("") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 8), breaks = seq(1, 8, by = 1))

varImpPlot ## finished in graphics software

########################

###################################################################################################
# now then, let's explore the PCA, PCoA and NMDS family
# just quickly, if I coose to select the treatments by standard deviation, how many are included at 2
# and 3 standard deviations from the mean
# to add both eWater instead (uses eW.frog and eW.other columns)
samp <- 1100
startCol <- 41
endCol <- 86
testPCA <- solutionsAll[1:samp, startCol:endCol]
testPCA <- rbind(testPCA, solutionsAll[(nrow(solutionsAll) - samp + 1):(nrow(solutionsAll)), startCol:endCol])

# lets name the solutions by treatment so we can label them and check if we need to 
namesA <- namesB <- NULL
for (i in 1:samp) {
  namesA <- append(namesA, paste0('H', i)) 
  namesB <- append(namesB, paste0('L', i)) 
}
names <- append(namesA, namesB)
rownames(testPCA) <- names

# make a vector for later for the 'treatment'                       
Treatment <- append(rep('high', samp), rep('low', samp))

# lets remove any wetland that don't vary at all (weren't watered) 
holder <- NULL
for (i in 1:ncol(testPCA)) {
  aa <- max(testPCA[ ,i ]) - min(testPCA[ ,i ])
  if (aa == 0) {
    cat(" remove column ", i, "\n")
    holder <- append(holder, i)   
  }
}
dim(testPCA)
testPCA <- testPCA[ , -holder]
dim(testPCA)

## the core nmds for this version of testPCA
##########################################################################
# run the fit 2 dimensions 
nmds_results <- metaMDS(comm = testPCA, distance = "euclidian", k = 2, try = 1000, zerodist = "add")  

nmds_results

# Shepards test/goodness of fit
goodness(nmds_results) # Produces a results of test statistics for goodness of fit for each point


stressplot(nmds_results) # Produces a Shepards diagram

#Now we can plot the NMDS
plot(nmds_results)
# colourblind colour 1 "#009E73"  green High
# colourblind colour 2 "#0072B2"  blue  arrows
# colourblind colour 3 "#D55E00"  red   Low

# nmds ordination for 50 high (green) and 50 low (red/brown) with characteristics in blue 

pl <- ordiplot(nmds_results, type = "n")
ordiellipse(nmds_results, groups = Treatment, scaling = "symmetric", kind = "ehull", col = c("#009E73", "#D55E00"), lwd=2, label = F)
# points(pl, "sites", pch=18, col="black")
points(pl, "sites", pch=18, col=c(rep("#009E73", samp),rep("#D55E00", samp)))
text(pl, "species", arrows = TRUE, length = 0.1, col = "blue", lwd = 1, cex = 0.6)

nmds_results

# which actions were significantly different between good and bad (simper analysis)
testPCA.Env <- rep('high', samp)
testPCA.Env <- append(testPCA.Env, rep('low', samp))
testPCA.Env <- as.matrix(testPCA.Env)
rownames(testPCA.Env) <- names
colnames(testPCA.Env) <- c('Outcome')
testPCA.Env <- as.data.frame(testPCA.Env)
testPCA.Env$Outcome <- as.factor(testPCA.Env$Outcome)

## start with a simper
sim <- with(testPCA.Env, simper(testPCA, Outcome, permutations = 99))
summary(sim)


## OK so everything is significant and the plot is hilarious
# lets do a summary plot

best <- testPCA[1:samp, ]
worst <- testPCA[(samp + 1):nrow(testPCA), ]
best1 <- best/100
worst1 <- worst/100

boxplot(best1, las = 2, xlim = c(1.5,(ncol(best1) - 0.5)), ylim = c(-0.1, 12), 
        border = "#009E73")
boxplot(worst1,add=TRUE,border = "#D55E00", xaxt = 'n', yaxt = 'n')
for (i in 1:ncol(best1)) {
  abline(v = i + 0.5, col ='#E6E6E6', lty = 2)
}




















# here be dragons


## are these influences positive or negative??
########################
# plot of  Criteria vs ObjFun
for (i in 1:16) {
  plot(x = solutionsAll$V40, y = solutionsAll[ , i], main = paste0('Criteria',i))
  # Fit a linear regression model
  lm_model <- lm(solutionsAll[, i] ~ solutionsAll$V40)
  abline(lm_model, col = "red", lwd = 2)  # Add the trendline to the plot
  
}

## to explore a tiny bit more for the best and worst values
Range <- 1100  #can play with this

best <- solutionsAll[1:Range, 1:16]
worst <- solutionsAll[(nrow(solutionsAll) - Range):nrow(solutionsAll), 1:16]

# make a wee table to hold ma summaries
Weightings <- matrix(data = NA, ncol = 2, nrow = 16)
colnames(Weightings) <- c(paste0('best',Range),  paste0('worst', Range))
rz <- NULL
for (i in 1:16) {
  a <- paste0('Criteria', i)
  rz <- append(rz, a)
}
rownames(Weightings) <- rz

# add the best
for (i in 1:16) {
  M <- round(mean(best[ , i]), 2)
  S <- round(sd(best[ , i]), 2)
  E <- paste0(M,' (', S, ')')
  Weightings[i, 1] <- E
}

# add the worst
for (i in 1:16) {
  M <- round(mean(worst[ , i]), 2)
  S <- round(sd(worst[ , i]), 2)
  E <- paste0(M,' (', S, ')')
  Weightings[i, 2] <- E
}

Weightings













## summarise the direction of the criteria weightings (sadly this looks shit) 
samp <- 200
best <- solutionsAll[1:samp, 1:16]
worst <- solutionsAll[(samp + 1):nrow(solutionsAll), 1:16]

####################################
lox <- c(4, 14, 8, 9, 7, 12, 6, 3, 2, 13, 15, 5, 16, 10, 11, 1)
edges <- c(1,(ncol(best)))
bestFill <- adjustcolor("#009E73", alpha.f = 0.25)
worstFill <- adjustcolor("#D55E00", alpha.f = 0.25)
lab <- weights <- paste("weight", 1:16)



boxplot(best, las = 2, outline = F, xlim = edges, ylim = c(-10.5, 10.5), at = (lox - 0.25),
        border = "white", col = 'white', boxwex = 0.5, labels = lab, xaxt = 'n')

axis(1, at = lox, labels = paste("weight", 1:16), las = 2)

for (i in 1:(ncol(best) -1)) {
  abline(v = i + 0.5, col ='#E6E6E6', lty = 2)
}
boxplot(best, las = 2, add = TRUE, outline = F, at = (lox - 0.25), xlim = edges, 
        ylim = c(-10.5, 10.5), border = "#009E73", boxwex = 0.5, xaxt = 'n', yaxt = 'n', col = bestFill)
boxplot(worst, add = TRUE, border = "#D55E00", xaxt = 'n', yaxt = 'n', xlim = edges, 
        outline = F, boxwex = 0.5, ylim = c(-10.5, 10.5), at = (lox + 0.25), col = worstFill)



library(tidyverse)
library(factoextra)




