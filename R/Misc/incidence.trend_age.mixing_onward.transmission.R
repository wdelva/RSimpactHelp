rm(list = ls())

setwd("/home/david/Dropbox/Niyukuri/Stress_Testing_Cluster_Code")

library(RSimpactCyan)
library(RSimpactHelper)

library(devtools)
library(Rcpp)
library(ape)
library(expoTree)
library(data.table)
library(RSimpactCyan)
library(RSimpactHelper)
library(readr)
library(phangorn)
library(dplyr)
library(adephylo)
library(treedater)
library(geiger)
library(picante)
library(igraph)
library(lme4)


############################################################
# 1. Run Simpact and read the simulation output (datalist) #
############################################################

datalist <- get(load("/home/david/Dropbox/Analysis_for_Manuscripts_December2017_Sequence_Coverage/MasterModelSubOptimalSeqCovearge.datalistB.RData"))


######################
# 2. Incidence trend #
######################

incidence.df.15.24 <- as.data.frame(incidence.calculator(datalist = datalist,
                                                         agegroup = c(15, 25), timewindow = c(10, 40)))

incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
incidence.df.15.24.women <- incidence.df.15.24$incidence[2]

incidence.df.15.24.men.95.ll <- incidence.df.15.24$incidence.95.ll[1]
incidence.df.15.24.women.95.ll <- incidence.df.15.24$incidence.95.ll[2]

incidence.df.15.24.men.95.ul <- incidence.df.15.24$incidence.95.ul[1]
incidence.df.15.24.women.95.ul <- incidence.df.15.24$incidence.95.ul[2]


incidence.df.25.34 <- as.data.frame(incidence.calculator(datalist = datalist,
                                                         agegroup = c(25, 35), timewindow = c(10, 40)))

incidence.df.25.34.men <- incidence.df.25.34$incidence[1]
incidence.df.25.34.women <- incidence.df.15.24$incidence[2]

incidence.df.25.34.men.95.ll <- incidence.df.25.34$incidence.95.ll[1]
incidence.df.25.34.women.95.ll <- incidence.df.25.34$incidence.95.ll[2]

incidence.df.25.34.men.95.ul <- incidence.df.25.34$incidence.95.ul[1]
incidence.df.25.34.women.95.ul <- incidence.df.25.34$incidence.95.ul[2]


incidence.df.35.44 <- as.data.frame(incidence.calculator(datalist = datalist,
                                                         agegroup = c(35, 45), timewindow = c(10, 40)))

incidence.df.35.44.men <- incidence.df.35.44$incidence[1]
incidence.df.35.44.women <- incidence.df.35.44$incidence[2]

incidence.df.35.44.men.95.ll <- incidence.df.35.44$incidence.95.ll[1]
incidence.df.35.44.women.95.ll <- incidence.df.35.44$incidence.95.ll[2]

incidence.df.35.44.men.95.ul <- incidence.df.35.44$incidence.95.ul[1]
incidence.df.35.44.women.95.ul <- incidence.df.35.44$incidence.95.ul[2]


##################################
# 3. Age mixing in transmissions #
##################################

# Source transmission networks builder function

source("/home/david/Dropbox/Niyukuri/Stress_Testing_Cluster_Code/transmNetworkBuilder.diff2.R")


# Function of age mixing data: agemixing.trans.df(datalist = datalist, 
#                                                 trans.network = trans.network)

# require library(data.table)

#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param trans.network Transmission networks obtained by by \code{\link{transNetBuilder4()}}.
#' @return a data table of all transmissions 
#' ID1 : ID of men  
#' ID2 : ID of women    
#' TOBID1 : Time of Birth of men     
#' TOBID2 : Time of Birth of women   
#' AgeID1 : Age of men at transmission time 
#' AgeID2 : Age of women at transmission time
#' AgeGap : Age gap between man & woman
#' infecttime : Infection time
#' samptime : Removal time

# 3.1 Function to compute the data table of age mixing in transmission

# Due to much computation time of transmNetworkBuilder.diff2() we put it external to age-mixing function



trans.network <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = datalist$itable$population.simtime[1])



agemixing.trans.df <- function(datalist = datalist, 
                               trans.network = trans.network){
  
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
  
  pers.infec <- pers.infec.raw[which(pers.infec.raw$InfectTime <= datalist$itable$population.simtime[1]),]
  
  # person table of infected individuals by seed event
  pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)
  
  # id of people who got infection by seed event: seeds.id
  seeds.id <- pers.table.seed$ID # do
  
  infectionTable <- vector("list", length(seeds.id))
  
  for (i in 1: length(seeds.id)) {
    
    trans.network.i <- as.data.frame(trans.network[[i]])
    
    if(nrow(trans.network.i) >= 3){
      
      trans.network.i <- trans.network.i[-1,]
      
      rrtable <- as.data.frame(cbind(trans.network.i$RecId, trans.network.i$DonId, 
                                     trans.network.i$GenderRec, trans.network.i$GenderDon,
                                     trans.network.i$TOBRec, trans.network.i$TOBDon,
                                     trans.network.i$InfecTime, trans.network.i$SampTime))
      
      names(rrtable) <- c("RecId", "DonId", "GenderRec",
                          "GenderDon", "TOBRec", "TOBDon", "InfecTime", "SampTime")
      
      rrtable.men <- subset(rrtable, rrtable$GenderDon=="0")
      rrtable.women <- subset(rrtable, rrtable$GenderDon=="1")
      
      ids.men <- rrtable.men$DonId
      ids.men.part.w <- rrtable.men$RecId
      age.gap.ID1 <- abs(rrtable.men$TOBDon) - abs(rrtable.men$TOBRec) # men are donors
      tob.men.ID1 <- rrtable.men$TOBDon
      tob.women.ID1 <- rrtable.men$TOBRec
      
      age.men.ID1 <- abs(rrtable.men$TOBDon) + rrtable.men$InfecTime
      age.women.ID1 <- abs(rrtable.men$TOBRec) + rrtable.men$InfecTime
      
      infectime.m <- rrtable.men$InfecTime
      samptime.m <- rrtable.men$SampTime
      
      ID1.m <- ids.men
      ID2.m <- ids.men.part.w
      age.gap.m <- age.gap.ID1
      infectime.m <- infectime.m 
      
      infectable.m <- cbind(ID1.m, ID2.m, tob.men.ID1, tob.women.ID1, age.men.ID1, age.women.ID1, age.gap.m, infectime.m, samptime.m)
      
      ids.women <- rrtable.women$DonId
      ids.women.part.m <- rrtable.women$RecId
      age.gap.ID2 <- abs(rrtable.women$TOBRec) - abs(rrtable.women$TOBDon) # men are receptors
      tob.men.ID2 <- rrtable.women$TOBRec
      tob.women.ID2 <- rrtable.women$TOBDon
      
      age.men.ID2 <- abs(rrtable.women$TOBDon) + rrtable.women$InfecTime
      age.women.ID2 <- abs(rrtable.women$TOBRec) + rrtable.women$InfecTime
      
      infectime.w <- rrtable.women$InfecTime
      samptime.w <- rrtable.women$SampTime
      
      ID1.w <- ids.women.part.m
      ID2.w <- ids.women
      age.gap.w <- age.gap.ID2
      infectime.w <- infectime.w
      
      infectable.w <- cbind(ID1.w, ID2.w, tob.men.ID2, tob.women.ID2, age.men.ID2, age.women.ID2, age.gap.w, infectime.w, samptime.w)
      
      infectable.i <- as.data.frame(rbind(infectable.m, infectable.w))
      
      names(infectable.i) <- c("ID1", "ID2", "TOBID1", "TOBID2", "AgeID1", "AgeID2", "AgeGap", "infecttime", "samptime")
      infectionTable[[i]] <- infectable.i
    }
    
    
  }
  
  
  infecttable <- rbindlist(infectionTable) 
  
  return(infecttable)
  
}


agemix.df <- agemixing.trans.df(datalist = datalist, 
                                trans.network = trans.network)


# 3.2. Fit mixed effects model to age mixing data

# require  library(lme4)

fit.agemix.trans <- function(datatable = agemix.df){
  
  datatable <- agemix.df
  
  agemix.inter <- lmer(AgeID2 ~ AgeID1 + (1|ID1), data = datatable) # choice of age_woman by a man
  
  return(agemix.inter)
  
}


agemix.fit <- fit.agemix.trans(datatable = agemix.df)



# 3.3. Visualisation of age mixing in transmission

coef.inter <- fixef(agemix.fit)

age.mix.intercept <- coef.inter[1]

age.mix.slope <- coef.inter[2]


x=agemix.df$AgeID1
y=agemix.df$AgeID2
plot(x, y, lwd=1, col = "blue",
     xlab = "Man Age",
     ylab = "Woman Age")
abline(coef = c(coef.inter[1], coef.inter[2]))


###########################
# 4. Onward transmissions #
###########################

# Function to return a list of how many onward infections produce by donors (after cquiring the infection)

# requires transNetBuilder4() function

#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param trans.network Transmission networks obtained by by \code{\link{transNetBuilder4()}}.
#' @return a vector of number of transmissions produced by one individual after ecquiring the infection

onwardtransmissions.dat <- function(datalist = datalist, 
                                    trans.network = trans.network){
  
  
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
  
  pers.infec.raw.died <- pers.infec.raw[pers.infec.raw$TOD != "Inf", ] # consider only these who had full time to transmit before they die
  
  pers.infec <- pers.infec.raw.died[which(pers.infec.raw.died$InfectTime <= datalist$itable$population.simtime[1]),]
  
  # person table of infected individuals by seed event
  pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)
  
  # id of people who got infection by seed event: seeds.id
  seeds.id <- pers.table.seed$ID # do
  
  
  # Onward transmissions in each transmission network
  
  onwardtransm <- vector("list", length(seeds.id))
  
  for (j in 1: length(seeds.id)) {
    
    trans.network.j <- as.data.frame(trans.network[[j]])
    
    trans.network.j <- trans.network.j[-1,] # remove the universal infector
    
    if(nrow(trans.network.j) > 1){ # consider transmission networks with at least one onward transmission
      
      d.j <- table(trans.network.j$DonId) # in the transmission table, the number of times DonId appears is the number of Onward transmissions after acuiring the infection 
      num.j <- as.data.frame(as.numeric(d.j))
      names(num.j) <- c("TransCount")
      onwardtransm[[j]] <- num.j
      
    }
    
  }
  
  onwardtransmissions <- rbindlist(onwardtransm) 
  
  count.dat <- onwardtransmissions$TransCount
  
  return(count.dat) # count.dat = all infections - seeds which didn;t produce at least one transmission
  
}


transm.count <- onwardtransmissions.dat(datalist = datalist, 
                                        trans.network = trans.network)

# Advise for which distribution: power-law???

# Visualisation of onward distribution

hist(transm.count)

fit.pow <- power.law.fit(transm.count)




