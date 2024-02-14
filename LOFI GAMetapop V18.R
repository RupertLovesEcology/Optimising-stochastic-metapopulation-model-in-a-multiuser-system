# this code runs on a PC (a opposed to HPCGA_MetapopV19.r that runs on an HPC). 
# We here model environmental water delivery to 23 wetlands between Weir 2 and Weir 3 in the Murray Darling Basin (Australia)
# eWater is prioritised for redgums and lignum and we determine the remaining water deliveries based on 16 criterai. 
# The objective fitness function is the number of surviving Litoria raniformis (southern bell frog) subpopulations across 
# 100 25-year hydric scenarios
# these outputs underpin the publication ""
# Mathwin, R., Zecchin, A., Mason, K., Gibbs, M., Wassens, S., and Bradshaw, C.J.A.

## Remove everything
rm(list = ls())

# Bootstrap all Libraries on the 'Master'
library(data.table)
library(pracma)
library(EnvStats)
library(rlang)
library(parallel)
library(foreach)
library(future)
library(doFuture)

# for rng in parallel 
library(rngtools)
library(doRNG)

# load the csvs (mostly as data.tables)
dirNm <- c('G:/My Drive/University Milestones/Cost benefit analysis 2_genetic algorithm/')

# load the starting populations
startPops <- fread(paste0(dirNm, 'csvs for analysis/GA Metapop startPops.csv'), select = c(2:1001)) # the starting populations         

# Arid hydrology (as data.tables for speed)
StayWet <- fread(paste0(dirNm,'csvs for analysis/GA StayWet 200AridScenarios.csv'))
setnames(StayWet, 1, 'Rw')
setcolorder(StayWet, c(names(StayWet)[-1], 'Rw'))

WetDry <- fread(paste0(dirNm,'csvs for analysis/GA WetDry 200AridScenarios.csv'))
setnames(WetDry, 1, 'Rw')
setcolorder(WetDry, c(names(WetDry)[-1], 'Rw'))

# load the wetland data
wetlandMetadata <- fread(paste0(dirNm, 'csvs for analysis/GA wetlandMetaDataV4.csv'))

## Populate underlying movement data (site-to-site distance which incorporates site-to-site movement type)
adj.wetlandMovement <- fread(paste0(dirNm, 'csvs for analysis/GA Meta adjMovement.csv'), select=c(2:24)) # removed rownames to keep as data.table

source(paste0(dirNm, 'R Files/matrixOperators.r'))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# custom functions
# breed the 100 fittest solutions
breed <- function(DF) {
  DF2 <- NULL
  for (i in 1:(nrow(DF)/2)) {
    # Randomly select two rows
    random_rows <- sample(1:nrow(DF), 2, replace = FALSE)
    row1 <- DF[random_rows[1], ]
    row2 <- DF[random_rows[2], ]
    row1 <-  as.vector(row1)
    row2 <-  as.vector(row2)
    
    # Select a single random crossover point along their length
    cross_over_point <- sample(1:(ncol(DF)-1), 1)
    
    # Create two new rows by crossing over at the point
    aa <- c(row1[1:cross_over_point], row2[(cross_over_point + 1):ncol(DF)])
    DF2 <- rbind(DF2, aa)
    aa <- c(row2[1:cross_over_point], row1[(cross_over_point + 1):ncol(DF)])
    DF2 <- rbind(DF2, aa)
  }
  return(DF2)
}

# mutate the solutions (1/16 probability: approximately one genome will mutate in each child) 
# the extent of mutation is resampled from an equilateral triangular distribution which scales smaller as the GA progresses
# the mutate function for the GA (updated 20231112)
mutate <- function(DF, iVal) {
  bb <- as.numeric(6.565435)/((exp((31.54871 - iVal)/-239.7932))+ 1) 
 aa <- bb * -1 
  for (i in 1:nrow(DF)) {
    for (j in 1:ncol(DF)) {
      if (runif(1, 0, 1) <= 1/16) {
        k <- rtri(1, min = aa, max = bb, mode = 0)
        DF[i,j] <- as.numeric(DF[i, j]) + k
        if (DF[i,j] > 10) { # reflect
          x <- 10 - as.numeric(DF[i,j])
          DF[i,j] <- 10 - abs(x)
        }
        if (DF[i,j] < -10) { # reflect
          x <- 10 + as.numeric(DF[i,j])
          DF[i,j] <- -10 + abs(x)
        }
      } 
    }
  }
  return(DF)
}

#create the egg density correcting function and variables
eggfunc <- function(x) (0.9991499)/((exp((0.9380068 - x)/-0.02672893))+ 1)

## function to allocate each emigrant to its new location
## use # destination <- determine.destination(adj.wetlandmove[wetlands],emigrants)
determine.destination <- function(journey, emigrnts) { 
  # journey <- as.numeric(journey)
  moveDist <- c()
  DSLoss <- USLoss <- 0
  for (imiCol in 1:emigrnts) {
    aa <- round(exp((log(runif(1,4.69e-05,0.3559407)) - 4.04117)/-1.294724),0)
    if ((aa > journey[1,23] | is.na(journey[1,23])) && runif(1,0,1) >= 0.5) {
      DSLoss <- DSLoss + 1 
    } else if ((aa > journey[1,1] | is.na(journey[1,1])) && runif(1,0,1) >= 0.5) {
      USLoss <- USLoss + 1 
    } else {
      moveDist <- append(moveDist,which.min(abs(journey - aa)))
    }
  }
  # quick reminder c) is DS loss, d) is US loss
  moveDist <- append(moveDist, DSLoss)
  moveDist <- append(moveDist, USLoss)
  return(moveDist)
}

# to quickly determine ages and mortality for immigrants from outside the reach
DetermineAges <- function(immiNums,mort) {
  fromUS <- fromDS <- c()
  # these proportions match the stable age structure of the population
  bins <- c(0, 0.66, 0.93, 0.987, 1)
  
  while (sum(immiNums[1,1]) > 0) {
    if (runif(1, 0, 1) < mort) {
      fromDS <- append(fromDS, cut(runif(1,0,1), breaks = bins, labels = FALSE, right = FALSE))
    } 
    immiNums[1,1] <- immiNums[1,1] - 1
  }
  while (sum(immiNums[1,2]) > 0) {
    if (runif(1, 0, 1) < mort) {
      fromUS <- append(fromUS,cut(runif(1, 0, 1), breaks = bins, labels = FALSE, right = FALSE))
    } 
    immiNums[1,2] <- immiNums[1,2] - 1
  }
  if (is.null(fromDS)) { fromDS <- NA }
  if (is.null(fromUS)) { fromUS <- NA }
  length(fromDS) <- max(length(fromDS), length(fromUS))
  length(fromUS) <- max(length(fromDS), length(fromUS))
  immiAges <- rbind(fromDS, fromUS)
  return(immiAges)
}

# the main metapopulation simulation function
MetaSim <- function(x, input) {
  # should output  86 numerics (criteria 1 - 16, svWetlands 1-23, svTotal 1, eWaterOther 1-23, eWaterFrog 1-23, eWaterTotal 1) 
  # the components of the return vector (evaluation) 
  scoringCriteria <- unlist(input[1:16])
  svWetlands <- rep(0, 23)
  svTotal <- 0
  eWaterOther <- rep(0, 23)
  eWaterFrog <- rep(0, 23)
  eWaterTotal <- 0
  
  
  eWaterAnnual2 <-  eWaterAnnual <- 429750 # volume to be delivered in m3
  # priority sites for other taxa
  redgumList <- c(4, 8, 9, 11, 13, 14, 15)
  lignumList <- c(5, 7)
  frognumList <- c(3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 17, 22, 23)
  eWaterTally <- c(0,0,0)
  
 # run parameters
  iterations <- 100
  simYears <- 25
  wetlandNum <- nrow(wetlandMetadata)
  
  ## RV of 17 spp. of anuran papers mean dispersal rates 16% +/- 14%
  ## Base probability of movement (will be adjusted using site K)
  wetMove <- 0.16
  
  ## mortality probability during movement
  moveSurv <- 0.5
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  # Create population dynamics 
  longev <- 5
  age.vec <- seq(0, longev, 1)
  lage <- length(age.vec)
  sex.ratio <- 0.5
  stages <- lage
  ## set population storage matrices n.mat
  n.mat <- array(data = 0, dim = (stages * (simYears + 1) * wetlandNum))
  dim(n.mat) <- c(stages, (simYears + 1), wetlandNum)
  
  popmat <- matrix(0, nrow = stages,ncol = stages)
  colnames(popmat) <- age.vec[1:stages]
  rownames(popmat) <- c('Fertility','Survival to 1 year','Survival to 2 year','Survival to 3 year','Survival to 4 year','Survival to 5 year')
  ## fertility data 
  min.clutch <- 1487.05
  max.clutch <- 5021.871
  mean.clutch <- 3244
  prop.breeding <- c(0,rep(1, 5))
  
  #duration data (eggs and tadpoles)
  hatch.dur <- c(2, 4)
  tadpole.dur <- c(70, 80) # duration (days) L. raniformis 23 deg Cree 1984
  
  ## egg survival 
  hatch.pr <- c(0.92578, 0.06448) # Christy PhD L. aurea (probability of fertilisation as a proxy for hatch probability)
  
  # tadpole survival 
  tp.s.alpha <- 3.644738
  tp.s.beta <- 10.37349
  
  # Adult annual survival Turner et al 2022
  ad.s.yr <- c(0.23, 0.08)
  
  # We hereafter consider 4 population sizes, each represented with 4 density feedback values (one for each of egg, tad, juve, adult)
  # Small, medium, large, very large #NOTE small sites aren't used in this model
  # they can be invoked to represent small, intensively manaeged breeding sites (see L. aurea in Sydney Olympic Park)
  # Assign the number of spawning masses that each site holds and convert it to the number of female eggs (uses mean clutch size)
  K.rel.egg.dem.vec <- c(NA, 28, 80, 129)
  K.rel.egg.dem.vec <- K.rel.egg.dem.vec * mean.clutch * sex.ratio
  
  # tadpoles surviving to 1 year old will inhibit themselves NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO THE ACTUAL NUMBER WILL BE DOUBLE  
  K.rel.tad.dem.vec <- c(NA, 270, 800, 1417)
  
  # number of first years present will influence survival from 1 to 2 years old but more strongly than for older age brackets         
  K.rel.juv.dem.vec <- c(NA, 108, 290, 570)
  
  # number of first years present will influence adult survival NOTE THIS EQUATION ONLY CONSIDERS FEMALES SO ACTUAL NUMBER WILL BE DOUBLE 
  #These numbers are used for the probability of emigration (ie above K high movement, below K lower movement) 
  K.rel.adult.dem.vec <- c(NA, 108, 290, 570)
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  # for generating immigration from up and downstream
  fromUSdist <- adj.wetlandMovement[1,]
  fromUSdist[1,1] <- 0
  fromDSdist <- adj.wetlandMovement[23,]
  fromDSdist[1,23] <- 50
  
  # Create a holder for losses out the top and bottom of the reach during emmigration
  emiLoss <- array(data = 0, dim = 2)
  dim(emiLoss) <- c(1,2)
  colnames(emiLoss) <- c('DS Loss', 'US Loss')
  immigrantAges <- matrix(data = NA, nrow = 2, ncol = 1)
  rownames(immigrantAges) <- c('fromDS', 'fromUS') 
  
  ## create a data.table to prioritise eWater delivery
  score <- rep(0, 23)
  WLand <- seq(1:23)
  critScores <- data.table(WLand, score)
  critScores[,vol:=wetlandMetadata$wetlandVolume]
  
  ## Version 13.0
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  ## 23 populations of three possible sizes (M,L,XL)
  ## 25 simulated years 
  
  ## set up arrays to track emigration and immigration
  moveAgeNm <- c('1-2yr','2-3yr','3-4yr','4-5yr','5+yr') 
  ## emigrants is from wetland# to wetland # and third dimension is age class (only 1+ - 5+ age classes move)
  emigrants <- array(data = 0, dim = wetlandNum * wetlandNum * 5)
  dim(emigrants) <- c(wetlandNum, wetlandNum, 5) 
  
  emi.gen <- array(data = 0, dim = 5)
  dim(emi.gen) <- c(1*5)
  
  immigrants <- array(data = 0, dim = wetlandNum * 5)
  dim(immigrants) <- c(wetlandNum,5)
  
  # one to hold each years immigrants
  bass2 <- seq(1,wetlandNum,1)
  for (i in 1:wetlandNum) {
    bass2[i] <- paste('To', i)
  }
  rownames(immigrants) <- bass2
  colnames(immigrants) <- moveAgeNm
  immigrantDestinations <- 0
  destinations <- c()
  
  #reset the trackers
  #  cat('iteration is ', iter, '\n') 
 
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  ## And NOW for the actual model!
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  # Create deepcopies (and not just pointers)
  for (iter in 1:iterations) {
    #  set the seed for each iteration to make outputs consistant
    set.seed(iter)
    
    # select the hydrological scenario
    minR <- 1 + (iter*simYears)
    maxR <- simYears + (iter*simYears)
    # subset fillTS and stayWetTS by the min and max rows we want (inclusive) 
    fillTS.iter <- subset(WetDry, Rw %in% minR:maxR)
    stayWetTS.iter <- subset(StayWet, Rw %in% minR:maxR)
    startP <-   unlist(startPops[, ..iter]) 
    
    for (st in 1:wetlandNum) {
      BottomVal <- st * 6
      n.mat[1:6,1,st] <-  unlist(startP[(BottomVal - 5):BottomVal])
      class(n.mat) <- 'numeric'
    }
    
    fillTS <- duplicate(fillTS.iter, shallow = F)
    stayWetTS <- duplicate(stayWetTS.iter, shallow = F)
    
    ## The Second nested loop: run the current projection set up for 25 years (simYears) 
    for (yr in 1:(simYears)) {
      if (yr == 1) {
        n.mat[,2:26,] <- 0
      }
      
      # evaluate the watering priorities # YEAR 1 IS ALWAYS WET
      if (yr > 1 &&  sum(fillTS.iter[yr, .SD, .SDcols = 1:23]) < 23) { # evaluate the watering priorities # YEAR 1 IS ALWAYS WET 
        evalList <- which((fillTS.iter[yr] == F), arr.ind = F) # only consider wetlands which will be dry
        for (dry in 1:(length(evalList))) {
          candidate <- evalList[dry]
          
          #  Q1 has this wetland been dry for 1 year only (to speed everyting up I'm using data.table so the atypical syntax is a necesity)
          if (fillTS.iter[yr, ..candidate] == F) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[1]] }
          
          #  Q2 has this wetland been dry for 2 year only 
          if (yr > 2 && sum(fillTS.iter[(yr-1):yr, ..candidate] == 0)) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[2]] }
          
          #  Q3 has this wetland been dry for 3 year only
          if (yr > 3 && sum(fillTS.iter[(yr-2):yr, ..candidate] == 0)) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[3]] }
          
          #  Q4 has this wetland been dry for 4 year only
          if (yr > 4 && sum(fillTS.iter[(yr-3):yr, ..candidate] == 0)) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[4]] }
          
          #  Q5 has this wetland been dry for 5 or more years
          if (yr > 5 && sum(fillTS.iter[(yr-4):yr, ..candidate] == 0)) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[5]] }
          
          # is this site Medium
          if (wetlandMetadata$wetlandSize[candidate] == 2) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[6]] }
          
          # is this site Large
          if (wetlandMetadata$wetlandSize[candidate] == 3) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[7]] }
          
          # is the site very large
          if (wetlandMetadata$wetlandSize[candidate] == 4) { critScores[candidate, score := critScores[candidate, score] + scoringCriteria[8]] }
          
          # did this site remain wet last year
          if (yr >= 3 && fillTS.iter[(yr - 1), ..candidate] == T && stayWetTS.iter[yr, ..candidate] == T) { 
            critScores[candidate, score := critScores[candidate, score] + scoringCriteria[9]]
          }
          
          # how far is the nearest site inside the reach (in km) * Criteria 10
          # distance to closest site in km * criteria
          if (candidate == 1) {
            nearestWL <- 0.05  
          } else if (candidate == 23) {
            nearestWL <- 0.3
          } else {
            nearestWL <- (min(adj.wetlandMovement[candidate,(candidate-1):(candidate + 1)], na.rm = T)/1000)
            critScores[candidate, score := critScores[candidate, score] + (nearestWL*scoringCriteria[10])]
          }
          
          # is this a redgum site
          if (candidate %in% redgumList) { 
            critScores[candidate, score := critScores[candidate, score] + scoringCriteria[11]]
          }
          
          # is this a lignum site  # 12
          if (candidate %in% lignumList) { 
            critScores[candidate, score := critScores[candidate, score] + scoringCriteria[12]]
          }
          
          # Let's assume perfect detection. Were there adult bell frogs calling there in the previous year
          if (sum(n.mat[2:stages,(yr-1),candidate]) > 0 && fillTS.iter[(yr-1), ..candidate] == T) {
            critScores[candidate, score := critScores[candidate, score] + scoringCriteria[13]]
            
          }
          
          # Let's assume perfect detection. Were there adult bell frogs calling there 2 years ago but the site was dry last year
          if ((sum(n.mat[2:stages,(yr-2),candidate]) > 0 && fillTS.iter[(yr-2), ..candidate] == T) && (fillTS.iter[(yr-1), ..candidate] == F)) {
            critScores[candidate, score := critScores[candidate, score] + scoringCriteria[14]]
          }
          
          # how many time was this site wet in the last 5 years (as a proportion 0-1)
          if (yr >= 6) {
            propWet <- (sum(fillTS.iter[(yr-1):(yr-5), ..candidate])/5)
            critScores[candidate, score := critScores[candidate, score] + (propWet*scoringCriteria[15])]
          }
          
          # rank of wetland from upstream to downstream (upstream scores higher from 1:1/23)
          # I haven't implemented the inverse because it can reverse by criteria dipping to -10
          critScores[candidate, score := critScores[candidate, score] + (scoringCriteria[16] * ((24-candidate)/23))]
          
        }        
        # to track our water use for other taxa (overRide)
        overRide.WL <- c()
        overRide.vol <- 0
        eWaterTally[(yr%%3) + 1] <- 0
        maxWetlands <- 10 - sum(eWaterTally)
        #redgums priority sites are watered 1 in 3. I let the first year be wet so this decision making doesn't kick in till yr 4
        # making sure this doesn't blow out the hydro threshold by periodically watering all 7 sites
        for (i in 1:length(redgumList)) {
          site <- redgumList[i]
          if (yr >= 4 && sum(fillTS.iter[(yr-2):(yr-1), ..site]) == 0 && 
              eWaterAnnual2 >= wetlandMetadata$wetlandVolume[site] && maxWetlands > 0) {
            maxWetlands <- maxWetlands - 1
            eWaterAnnual2 <- eWaterAnnual2 - wetlandMetadata$wetlandVolume[site]
            overRide.WL <- append(overRide.WL, site)
            overRide.vol <- overRide.vol + wetlandMetadata$wetlandVolume[site]
            eWaterOther[site] <- eWaterOther[site] + 1
            eWaterTotal <- eWaterTotal + 1
          }
        }
        
        # lignum is watered every third year in this reach (on average), 2 dry is a threshold
        for (i in 1:length(lignumList)) {
          site <- lignumList[i]
          if (yr >= 4 && sum(fillTS.iter[(yr-2):(yr-1), ..site]) == 0 && 
              eWaterAnnual2 >= wetlandMetadata$wetlandVolume[site] && maxWetlands > 0) {
            maxWetlands <- maxWetlands - 1
            eWaterAnnual2 <- eWaterAnnual2 - wetlandMetadata$wetlandVolume[site]
            overRide.WL <- append(overRide.WL, site)
            overRide.vol <- overRide.vol + wetlandMetadata$wetlandVolume[site]
            eWaterOther[site] <- eWaterOther[site] + 1
            eWaterTotal <- eWaterTotal + 1
          }
        }
        
        # how much water remains for frogs?
        eWater.remaining <- eWaterAnnual - overRide.vol
        if (eWater.remaining < 0) { eWater.remaining <- 0 }
        
        #  sort the candidate wetlands by score
        critScores <- critScores[order(critScores$score, decreasing = T),]  
        # if we haven't added the next priority (overRide) and there is enough water left (eWater.remaining) then deliver eWater  
        for (i in 1:23) {
          if (critScores$vol[i] <= eWater.remaining && (critScores$WLand[i] %in% overRide.WL) == F && maxWetlands > 0) {
            maxWetlands <- maxWetlands - 1
            site <- critScores$WLand[i]
            eWater.remaining <- eWater.remaining - critScores$vol[i]
            overRide.WL <- append(overRide.WL, site)
            eWaterFrog[site] <- eWaterFrog[site] + 1
            eWaterTotal <- eWaterTotal + 1
          }
        }
        
        eWaterTally[((yr%%3) + 1)] <- length(overRide.WL)
        
        # now change our breeding outcomes by wetting the eWatered sites (holding on to overRideWL to help tracking)  
        if (length(overRide.WL) >= 1) {
          for (i in 1:(length(overRide.WL))) {
            s <- paste0('V',overRide.WL[i])
            fillTS.iter[yr, (s) := T]
          }
        }
        
        # reset our trackers
        overRide.WL <- c()
        overRide.vol <- 0
        critScores <- critScores[order(critScores$WLand, decreasing = F),]         
        critScores$score[] <- 0
        eWaterAnnual2 <- eWaterAnnual
      }
      ## Cycle through each of the wetlands, The Wetlands Loop from upstream to downstream (1 - 23)
      for (wetland in 1:wetlandNum) {
        # if no frogs are alive at this wetland this year then skip this step
        if (sum(n.mat[, yr, wetland]) == 0) { break }
        # stochastically resample the population matrix (popmat) at this wetland for the coming year
        # assign fertility
        popmat[1,] <- round(sex.ratio*runif(stages, min=min.clutch, max=max.clutch) * prop.breeding, 0)
        # assign survival to metamorphosis
        hatch <- pmin(rnorm(1, mean = hatch.pr[1], sd=hatch.pr[2]), 1)
        tad.sv <- pmax(pmin(rbeta(1, tp.s.alpha, tp.s.beta),1), 0) 
        # duration of tadpole and egg stages
        postmet.dur <- 365 - round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
        # sample a daily probability of survival
        toad.daily.s <- nthroot(pmax(rnorm(1,ad.s.yr[1], ad.s.yr[2]), 0), 365)
        # survival from metamorphosis to 1yr
        toad.s <- hatch * tad.sv * (toad.daily.s ^ postmet.dur)
        surv.iter <- c(toad.s, pmax(rnorm(5,ad.s.yr[1], ad.s.yr[2]), 0))
        diag(popmat[2:(stages), ]) <- surv.iter[-stages]
        
        #note these are placeholder values for popmat.current
        popmat.current <- popmat
        
        # Implement density dependence effects for each of three wetland size classes and feed into popmat[,,]
        # note density effect on egg laying is applied after matrix multiplication 
        s.mult.iter.ad <- s.mult.iter.juv <-s.mult.iter.tad <- 1
        
        # instil density dependence for tadpoles growing into year 1 adults
        K.rel.tad <- n.mat[1, yr, wetland]/K.rel.tad.dem.vec[wetlandMetadata$wetlandSize[wetland]]
        if (!is.nan(K.rel.tad)) {
          if (K.rel.tad > 1.4)  { K.rel.tad <- 1.4 }
          if (K.rel.tad <= 1.4) {
            s.mult.iter.tad <- (as.numeric(0.8227) + (K.rel.tad * as.numeric(0.5914)))^(-1/as.numeric(0.1579))
            popmat.current[2,1] <- popmat[2,1] * s.mult.iter.tad  
          } 
        }
        
        # instil density dependence for juveniles (1 - 2 years) is driven by the  number of yr 1 present/competing per Berven
        K.rel.juv <- n.mat[2,yr,wetland]/K.rel.juv.dem.vec[wetlandMetadata$wetlandSize[wetland]] 
        if (!is.nan(K.rel.juv)) {
          if (K.rel.juv > 0.65)  { K.rel.juv <- 0.65 }
          if (K.rel.juv <= 0.65) {
            s.mult.iter.juv <- as.numeric(1.2131120)/((exp((0.5508745 - K.rel.juv)/-0.1825602))+ 1)
            popmat.current[3,2] <- popmat[3,2] * s.mult.iter.juv    } 
        } 
        
        # instill density dependence for adults  driven by the  number of yr 1s emerging from Berven 2009
        K.rel.adult <- sum(n.mat[2:6,yr,wetland])/K.rel.adult.dem.vec[wetlandMetadata$wetlandSize[wetland]]
        if (!is.nan(K.rel.adult)) {
          if (K.rel.adult > 0.65)  { K.rel.adult <- 0.65 }
          if (K.rel.adult <= 0.65) {
            s.mult.iter.ad <- as.numeric(1.2131120)/((exp((0.5508745 - K.rel.juv)/-0.1825602))+ 1) 
            for (adgens in 4:6) {
              popmat.current[adgens,(adgens-1)] <- (popmat.current[adgens,(adgens - 1)] * s.mult.iter.ad)  
            }
          }
        }  
        
        # set popmat.fail for use during dry years
        popmat.fail <- popmat.current
        popmat.fail[1,] <- 0  
        
        # always allow the first year to be wet at all sites 
        # this removes the possibility of all populations going stochastically extinct before evaluation commences.
        if (yr == 1) {
          n.mat[,yr+1,wetland] <- popmat.current %*% n.mat[,yr,wetland] 
        } else {
          # if it's wetland 1 (Banrock1) water every second year (Banrock2 is not watered for frogs)
          if (wetland == 1 && fillTS.iter[ yr, wetland, with = FALSE] == F) {
            n.mat[,yr+1,wetland] <- popmat.current %*% n.mat[,yr,wetland] 
          } else { 
            
            # if this year is wet
            if (fillTS.iter[yr, wetland, with = FALSE] == T) {
              # was it wet last year, AND were the river level high enough to staywet
              if (fillTS.iter[(yr - 1),wetland, with = FALSE] == T && stayWetTS.iter[ yr, wetland, with = FALSE] == T) {
                # apply the accumulation penalty whereby tadpoles are reduced by predators
                popmat.current[1,] <- popmat.current[1,]/250
              } 
              n.mat[,yr+1,wetland] <- popmat.current %*% n.mat[,yr,wetland] 
            } else {
              n.mat[,yr+1,wetland] <- popmat.fail %*% n.mat[,yr,wetland]
            }
          }
        }
        
        # round to whole frogs  
        n.mat[(2:6),yr+1,wetland] <- floor(n.mat[(2:6),yr+1,wetland])
        
        # density dependant inhibition of survival for spawn
        K.rel.egg.dem <-  K.rel.egg.dem.vec[wetlandMetadata$wetlandSize[wetland]]
        K.rel.egg <- n.mat[1, yr+1,wetland]/K.rel.egg.dem
        if (!is.nan(K.rel.egg) && (K.rel.egg > 0)) {
          if (K.rel.egg <= 1.02) {  
            intArea <-  integrate(eggfunc,lower=0,upper=K.rel.egg)
            area <- as.numeric(intArea[1])/K.rel.egg
            if (area >= 1) { area <- 1 }
            n.mat[1, yr+1, wetland] <- (round(n.mat[1, yr+1, wetland] * area))
          } else if (K.rel.egg > 1.02) {
            # if  > than the limit for the pond then calculate what happens to the first 99% of the pond limit (egg1) then apply the 99th%ile inhibition on the remaining eggs(remain)
            egg1 <- (K.rel.egg.dem * 0.9176419)
            remain <- (n.mat[1, yr+1, wetland] - K.rel.egg.dem) 
            n.mat[ 1, yr+1, wetland] <- (round(egg1 + (0.04442657 * remain)))  
          } else { 
            stop('Crashed at line 1894: the egg conversion value K.rel.egg is misbehaving')
          }  
        }
        ##  Calculate how many frogs disperse from this wetland 
        ## We are inside the wetland loop, so we store all dispersal (age and destination) and then complete these dispersals 
        # when iterating between years (after we have iterated all 23 wetland  this year)
        totalPop <- sum(n.mat[2:6,yr + 1,wetland])
        emigrants <- 0
        if (totalPop > 0) {
          # we use carrying capacity to determine the probability of moving using Pr(move) = n(t)/K * wetMove 
          # note this is an individual probability of movement, not a proportion of movement
          PrMove <- (totalPop/(K.rel.adult.dem.vec[wetlandMetadata$wetlandSize[wetland]])) * wetMove
          ## In order to maintain a PrMove of mean 0.16 +/- 0.14 we constrain the possible values 
          if (PrMove > 0.3) { PrMove <- 0.3 }
          if (PrMove < 0.02) { PrMove <- 0.02 }
          for (potMov in 1:totalPop) { 
            if (runif(1,0,1) < PrMove) { emigrants <- emigrants + 1 }
          }
        }
        
        #############   
        if (emigrants >= 1) {
          # Use the determine destinations function to determine the destinations that each departing frog arrives at   
          destinations <- determine.destination(adj.wetlandMovement[wetland,],emigrants)
          ## To calculate the number of emigrants lost at the upper and lower ends of the reach inside of the determine.destination function
          ## I have hidden them in the returned array. I will now extract them and add them instead to emi.Loss 
          # nb mortality is NOT applied here
          for (simple in 2:1) {
            emiLoss[1,simple] <- destinations[length(destinations)]
            destinations <- destinations[-(length(destinations))]
          }
          
          ## Determine the age class of each emigrant and remove them from n.mat (this includes emiLoss)
          emi.gen[] <- 0
          for (sss in 5:1) {
            emi.gen[sss] <- n.mat[sss+1,yr+1,wetland]/totalPop
            emi.gen[sss] <- round(emi.gen[sss] * emigrants, digits = 0)
            if(is.nan(emi.gen[sss])) { emi.gen[sss] <- 0 }
          }
          n.mat[2:6,yr+1,wetland] <- n.mat[2:6,yr+1,wetland] - emi.gen
          
          # deal with any remaining emigrants 
          emigrants2 <- emigrants - sum(emi.gen)
          rem <- 2 
          while (emigrants2 > 0) {
            if (n.mat[rem,yr+1,wetland] >= emigrants2) { 
              n.mat[rem,yr+1,wetland] <- n.mat[rem,yr+1,wetland] - emigrants2
              emi.gen[rem] <- emi.gen[rem] + emigrants2
              emigrants2 <- 0
            } else {
              emi.gen[rem] <- emi.gen[rem] + n.mat[rem,yr+1,wetland]
              emigrants2 <- emigrants2 - n.mat[rem, yr+1, wetland]
              n.mat[rem, yr+1, wetland] <- 0
            }
            rem <- rem + 1 
            if (rem > 6) { break('line 772 the remainder isnt working properly') }
          }  
          
          # quick error check
          for (remx in 2:6) {
            if (n.mat[remx,yr+1,wetland] < 0 | isFALSE(all.equal(n.mat[remx,yr+1,wetland], as.integer(n.mat[remx,yr+1,wetland])))) {
              stop('Some issues with the removal of emigrants')
            }
          }
          
          # remove the frogs that left the reach (age selected at random)
          lostFrogs <- sum(emiLoss)
          while (lostFrogs > 0) {
            kik <- round(runif(1,1,5),0)
            if(emi.gen[kik] > 0) { 
              emi.gen[kik] <- emi.gen[kik] - 1
              lostFrogs <- lostFrogs - 1
            }
          }
          
          ## Assign immigrants from emi.gen into immigrants, calculating movement survival
          # from the remaining frogs    
          lenDes <-  length(destinations) 
          while (lenDes > 0) {
            kik <- round(runif(1,1,5),0)
            if(emi.gen[kik] > 0) { 
              emi.gen[kik] <- emi.gen[kik] - 1 
              # check if it survived 
              if (runif(1,0,1) <= moveSurv) {
                immigrants[destinations[lenDes],kik] <- immigrants[destinations[lenDes],kik] + 1
              }
              lenDes <- lenDes -1  
            }
          }
          
          # call the DetermineAges function to create an array of immigrant ages from outside the reach. Store as immigrantAges. 
          if (sum(emiLoss, na.rm = TRUE) > 0) {
            immigrantAges <- DetermineAges(emiLoss,moveSurv)
            if (ncol(immigrantAges) > 0) {
              # determine arrival location by DS row then US row of ImmigrantAges, pasting into immigrants matrix
              for (imiCol in 1:ncol(immigrantAges)) {
                if (!is.na(immigrantAges[1,imiCol])) {
                  aa <- round(exp((log(runif(1,4.69e-05,0.3559407)) - 4.04117)/-1.294724),0)
                  bb <- which.min(abs(fromDSdist - aa))
                  immigrants[bb,immigrantAges[1,imiCol]] <- immigrants[bb,immigrantAges[1,imiCol]]  + 1
                } 
                if (!is.na(immigrantAges[2,imiCol])) {
                  aa <- round(exp((log(runif(1,4.69e-05,0.3559407)) - 4.04117)/-1.294724),0)
                  bb <- which.min(abs(fromUSdist - aa))
                  immigrants[bb,immigrantAges[2,imiCol]] <- immigrants[bb,immigrantAges[2,imiCol]]  + 1
                }
              }
            }
          }
        }
        # reset the arrays 
        immigrantAges[] <- NA
        emiLoss[] <- 0
        destinations <- c()
        
        ##  Last line of the 'wetland' loop which cycles through each of our wetlands in turn
      }
      
      ##  add all the immigrants to the respective n.mats and then reset the immigrants array
      if (sum(immigrants > 0)) {  
        for (kk in 1:wetlandNum) {
          for (kkk in 1:5) {
            if (immigrants[kk,kkk] > 0) {
              n.mat[kkk+1,yr+1,kk] <- (n.mat[kkk+1,yr+1,kk] + immigrants[kk,kkk])
            }
          }  
        }
      } 
      # reset array
      immigrants[] <- 0 
      
      #  If this is the 25th (and final) year then store all of the necessary values     
      if (yr == simYears) {
        for (k in 1:wetlandNum) {
          if (sum(n.mat[1:6,yr+1,k] > 0)) {
            svWetlands[k] <- svWetlands[k] + 1
            svTotal <- svTotal + 1
          }
        }
      } 
      
      # Last line of the simYears loop which controls how many years we run for (25 years)
    }
    
    # Last line of the iteration Loop (1)  
  }
  evaluation <- c(unlist(unname(scoringCriteria)), unlist(svWetlands), unlist(svTotal), eWaterOther, eWaterFrog, eWaterTotal)
#  message("MSE: ", x)
  return(evaluation)
  # last line of the MetaSim function  
}


###############################################################################################################
# let's setup the parallel architecture
cores <- (detectCores() - 1)
registerDoFuture()
plan(multisession, workers = cores) # stipulate how many workers. Seems to work faster than letting the computer choose. 

  fitSolutions <- fread('G:/My Drive/University Milestones/Cost benefit analysis 2_genetic algorithm/AAA Lofi HPC/iterSeeds Home/iterSeed resolved.csv')
  fitSolutions <- fitSolutions[ , -1]
  
   as.data.frame(fitSolutions)
 
  unfitSolutions <- NULL
  evaluatedSolutions <- NULL
  
  for (opt in 1:5000) { 
    
     print(opt)  
    
    NewPop <- breed(fitSolutions[,1:16]) # gains rownames here?? Is a Matrix
    # NewPop <- as.data.frame(NewPop) 
    rownames(NewPop) <- c()
    NewPop <- mutate(NewPop, opt)
    NewPop <- matrix(unlist(NewPop), nrow = 100, ncol = 16) # this solves the non-numeric issue
    NewPop <- as.data.frame(NewPop) # this can feed into lapply
    # NewPop <- split(NewPop, seq(nrow(NewPop))) # don't need a list for 'foreach'
    
    # Resolve the seed solutions
    evaluatedSolutions <- foreach(i=1:nrow(NewPop), .combine='rbind', .inorder=FALSE) %dorng% MetaSim(i, NewPop[i, ])
    
    # combine the evaluated children with their parents 
    evaluatedSolutions <- rbind(fitSolutions, evaluatedSolutions)
    evaluatedSolutions <-  as.data.frame(evaluatedSolutions)
    evaluatedSolutions <- evaluatedSolutions[order(evaluatedSolutions$V40, decreasing=TRUE), ] # objective function = col 'V40'
    #  select the 100 fittest individuals (objective function = sum of surviving populations) 
    fitSolutions <- evaluatedSolutions[1:100, ]
    unfitSolutions <- evaluatedSolutions[101:nrow(evaluatedSolutions),  ]
    write.csv(fitSolutions, paste0('G:/My Drive/University Milestones/Cost benefit analysis 2_genetic algorithm/AAA Lofi HPC/iterSeeds Home/HomeLofi setSeed fit', opt, '.csv')) 
    write.csv(unfitSolutions, paste0('G:/My Drive/University Milestones/Cost benefit analysis 2_genetic algorithm/AAA Lofi HPC/iterSeeds Home/HomeLofi setSeed Unfit', opt, '.csv')) 
    evaluatedSolutions <- NULL
    unfitSolutions <- NULL
  }
  
  
  ######################################################################################################################
  
  
  
  
  