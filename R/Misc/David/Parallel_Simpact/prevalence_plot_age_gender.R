library(tidyr)
library(ggplot2)

prevalence.df <- prevalence.calculator(datalist = datalist,
                                       agegroup = c(15, 30),timepoint = 30)

prevalence.df <-prevalence.plotter(datalist = datalist, agegroup = c(15, 50))


#Swaziland SHIMS 2011 incidence by gender
xmen <- c(0.8,6.6,21.3,36.6,47.0,45.5,42.5)
xwom <- c(14.3, 31.5, 46.7, 53.8, 49.1, 39.7, 31.6)
time <- 1:7

data.plot <- data.frame(cbind(time, xwom,xmen))
names(data.plot) <- c("time", "Women", "Men")
## wide to long format
data.plot <- data.plot %>%
  tidyr::gather(Gender, "Incidence", 2:3)

ggplot(data.plot,
       aes(x=time, y=Incidence, group=Gender, label = Incidence)) +
  geom_point() + geom_line() + aes(colour = Gender)  +
  xlim("18-19","20-24","25-29","30-34","35-39","40-44","45-49") +
  ylim(0,60) +
  geom_text(aes(label=Incidence), size=4, vjust=-.8, show.legend = FALSE) +
  labs(title = "2011 HIV Prevalence by Age Group in Swaziland",
       x = "Age Group", y = " HIV Prevalence") +  theme_bw() +
  theme(axis.text.x  = element_text(vjust=0.5, size=14),
        axis.title.x = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=14),
        axis.title.y = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5, size=16))



### Overall prevalence plotter


prevalence.plotter <- function(datalist = datalist,
                               agegroup = c(15, 50)){
  prev.agegroup <- agegroup
  timevect <- 0:datalist$itable$population.simtime[1]

  prevplot.df <- data.frame(timevect = timevect, prev = NA, ll = NA, ul = NA)
  rowindex <- 1
  for (timepoint in timevect){
    prev.tibble <- prevalence.calculator(datalist = datalist,
                                         agegroup = prev.agegroup, timepoint = timepoint)

    prevplot.df$prev[rowindex] <- as.numeric(prev.tibble$pointprevalence[3])
    prevplot.df$ll[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ll[3])
    prevplot.df$ul[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ul[3])
    rowindex <- rowindex + 1
  }
  prevplot.df <- prevplot.df %>% tidyr::gather(type, prevalence, prev:ul)
  spectrum.estimates <- data.frame(timevect = 3+(0:30), #1980:2010
                                   prevalence = 0.01 * c(0.01, 0.03, 0.07, 0.12, 0.19,
                                                         0.30, 0.44, 0.65, 0.93, 1.34, 1.91,
                                                         2.69, 3.78, 5.23, 7.11, 9.41, 12.00, 14.70,
                                                         17.32, 19.66, 21.61,23.12, 24.24, 25.01,
                                                         25.50, 25.80, 25.91, 25.86, 25.82, 25.80, 25.81),
                                   type = "spectrum.output")

  prevplot.df <- gtools::smartbind(prevplot.df, spectrum.estimates)

  prevplot.df$pointestimate <- prevplot.df$type=="prev"
  prevplot.df$spectrum.output <- prevplot.df$type=="spectrum.output"

  p <- ggplot(prevplot.df, aes(x=timevect+1977, y=prevalence, group=type)) +
    geom_line(aes(linetype = !(pointestimate | spectrum.output),
                  colour = spectrum.output)) +
    xlab("Simulation time") +
    ylab(" HIV prevalence") +
    #facet_wrap(~Gender, nrow=1)
    theme_bw() +
    theme(legend.position = "none")
  return(p)
}


## Sampling in a given time interval

# 1. Overa-all age-mixing patterns
###################################

survey.window <- c(20,25)

persontable <- dplyr::filter(datalist$ptable, TOD>=survey.window[1]) # ||TOD<=survey.window[2])

persontable.hiv <- dplyr::filter(persontable, InfectTime!="Inf")


indiv.samp <- persontable.hiv$ID

source("R/transmNetworkBuilder.diff3.R")

simpact.trans.net <- transmNetworkBuilder.diff3(datalist = datalist,
                                                endpoint = datalist$itable$population.simtime[1])

trans.network <- simpact.trans.net


# Define function for age-mixing in transmission networks

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

    if(nrow(trans.network.i) >= 2){

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
      age.gap.ID2 <- abs(rrtable.women$TOBRec) - abs(rrtable.women$TOBDon) # women are receptors
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


  infecttable <- as.data.frame(rbindlist(infectionTable))

  return(infecttable)

}


# Call the function

agemix.df <- agemixing.trans.df(datalist = datalist,
                                trans.network = trans.network)

agemix.df <- dplyr::filter(agemix.df, infecttime < survey.window[2])

indiv.w <- intersect(agemix.df$ID2, indiv.samp)
indiv.m <- intersect(agemix.df$ID1, indiv.samp)

rm.non.trans <- setdiff(indiv.samp, c(indiv.w,indiv.m)) # seed individuals who never transmit the infection

indiv.samp.trans <- intersect(indiv.samp, c(indiv.w, indiv.m))

length(indiv.w)+length(indiv.m) # 237
#  119         &   118
length(indiv.samp.trans) # 237

indiv.samp.men.agemix.df <- subset(agemix.df, agemix.df$ID1%in%indiv.samp.trans) # indiv.m)

indiv.samp.women.agemix.df <- subset(agemix.df, agemix.df$ID2%in%indiv.samp.trans) # indiv.w)

indiv.samp.agemix.df <- rbind(indiv.samp.men.agemix.df, indiv.samp.women.agemix.df)

indiv.samp.agemix.df <- indiv.samp.agemix.df[!duplicated(indiv.samp.agemix.df), ]

nrow(indiv.samp.agemix.df)

overall.age.mixing.df <- indiv.samp.agemix.df

# Fit mixed effects model to age mixing data

# require  library(lme4)

fit.agemix.trans <- function(datatable = agemix.df){

  datatable <- agemix.df

  agemix.inter <- lmer(AgeID2 ~ AgeID1 + (1|ID1), data = datatable) # choice of age_woman by a man

  return(agemix.inter)

}


overall.agemix.fit <- fit.agemix.trans(datatable = overall.age.mixing.df)



# Visualisation of age mixing in transmission

coef.inter <- fixef(overall.agemix.fit)

age.mix.intercept <- coef.inter[1] # Statistics to save

age.mix.slope <- coef.inter[2] # Statistics to save

#
x=agemix.df$AgeID1
y=agemix.df$AgeID2
plot(x, y, lwd=1, col = "blue",
     xlab = "Man Age",
     ylab = "Woman Age")
abline(coef = c(coef.inter[1], coef.inter[2]))


# 2. Age-mixing between younger women < 25 years and men of 25-40 years

indiv.samp.agemix.df$AgeID1 <- abs(indiv.samp.agemix.df$AgeID1)+indiv.samp.agemix.df$infecttime # absolute age

indiv.samp.agemix.df$AgeID2 <- abs(indiv.samp.agemix.df$AgeID2)+indiv.samp.agemix.df$infecttime # absolute age


w.25.m.25.40.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2<25, AgeID1>=25 & AgeID1<=40)


w.25.m.25.40.agemix.fit <- fit.agemix.trans(datatable = w.25.m.25.40.df)



# Visualisation of age mixing in transmission

coef.inter <- fixef(w.25.m.25.40.agemix.fit)

age.mix.intercept <- coef.inter[1] # Statistics to save

age.mix.slope <- coef.inter[2] # Statistics to save

#
x=agemix.df$AgeID1
y=agemix.df$AgeID2
plot(x, y, lwd=1, col = "blue",
     xlab = "Man Age",
     ylab = "Woman Age")
abline(coef = c(coef.inter[1], coef.inter[2]))


# 3. Women < 25 years and men 25 - 40 years

w.25.m.25.40.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2<25, AgeID1>=25 & AgeID1<=40)

# 4. Women < 25 years and men 41 - 49 years

w.25.m.41.49.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2<25, AgeID1>=41 & AgeID1<=49)

# 5. Women 25 - 40 years and men 25 - 40 years

w.25.40.m.25.40.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2>=25 & AgeID2<=40, AgeID1>=25 & AgeID1<=40)

# 6. Women 25 - 40 years and men 41 - 49 years

w.25.40.m.41.49.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2>=25 & AgeID2<=40, AgeID1>=41 & AgeID1<=49)


# 7. Women 41 - 49 years and men 25 - 40 years

w.41.49.m.25.40.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2>=41 & AgeID2<=49, AgeID1>=25 & AgeID1<=40)


# 8. Women 25 - 40 years and men 41 - 49 years

w.25.40.m.41.49.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2>=25 & AgeID2<=40, AgeID1>=41 & AgeID1<=49)


# 9. Young men < 25 years and women < 25 years !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

m.25.w.25.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2<25, AgeID1<25)


# 10. Young men < 25 years and women 25 - 40 years

m.25.w.25.40.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2>=25 & AgeID2<=40, AgeID1<25)

# 11. Young men < 25 years and women 41 - 49 years

m.25.w.41.49.df <- dplyr::filter(indiv.samp.agemix.df, AgeID2>=41 & AgeID2<=49, AgeID1<25)

