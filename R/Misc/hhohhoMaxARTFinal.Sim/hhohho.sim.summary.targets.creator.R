#This will be called through source()
#function to clean up the column or row names,
format.names <- function(x,replace = "X",...){
  gsub(replace,'',x)
}

#Recruitment
#maxart.study.pop <- recruit.study.clients(chunk.datalist.test)

#Dummy with no real recruitment
two.weeks <- 14/365.245
twelve.months <- 12/12
six.months <- 6/12
init.immediate <- 0
init.not.eligible <- -6/12
more.months <- twelve.months + six.months

hhohho.sim.summary.creator <- function(sim.datalist = chunk.datalist.test){

  #creator a range of summary statistics
  sim.datalist$ptable$primary <- 0

  #get the MaxART population group
  maxart.study.pop <- subset(sim.datalist$ptable, pfacility != "Not Hhohho")

  df.rw <- nrow(maxart.study.pop)

  #Dummy when eligible
  sample.eligible <- sample(c(two.weeks, twelve.months, six.months,
                                init.immediate, more.months,
                                init.not.eligible, NA), replace = TRUE, df.rw)

  maxart.study.pop$art.eligible <- maxart.study.pop$TreatTime + sample.eligible


  #dummy 40% as primary
  count.prim <- round(0.4 * nrow(maxart.study.pop), 0)
  maxart.study.pop$primary[which(maxart.study.pop$ID %in%  sample(maxart.study.pop$ID, count.prim))] <- 1

  all.prim.maxart <- sum(maxart.study.pop$primary==1)
  all.maxart <- nrow(maxart.study.pop)

  #Get the two primary and all clients distinction
  max.art.ret.study.pop.prim <- maxart.study.pop %>%
    filter(TreatTime !=Inf & primary ==1) %>%
    as.data.frame

  max.art.ret.study.pop.all <- maxart.study.pop %>%
    filter(TreatTime !=Inf)

  #alldf
  chunk.datalist.test.prim <- sim.datalist
  chunk.datalist.test.prim$ptable <- max.art.ret.study.pop.prim

  chunk.datalist.test.all <- sim.datalist
  chunk.datalist.test.all$ptable <- max.art.ret.study.pop.all

  ### ART coverage ####################

  hhohho.all.prim.art.coverage.var <- maxart.study.pop %>%
    filter(TreatTime!=Inf) %>%
    summarise(all.prim = sum(primary==1)/all.prim.maxart,
              all.all = n()/all.maxart,

              b4.elig.prim = sum(TreatTime < art.eligible & primary==1)/all.prim.maxart,
              b4.elig.all = sum(TreatTime < art.eligible )/all.maxart,

              within2wks.prim = sum(TreatTime>=art.eligible &
                                      TreatTime < (TreatTime+two.weeks) & primary==1)/all.prim.maxart,
              within2wks.all = sum(TreatTime>=art.eligible &
                                     TreatTime < (TreatTime+two.weeks) )/all.maxart,

              amongelig.prim = sum(TreatTime > art.eligible & primary==1)/all.prim.maxart,
              amongelig.all = sum(TreatTime > art.eligible )/all.maxart,

              within12m.prim = sum(TreatTime==art.eligible &
                                     TreatTime < (TreatTime+twelve.months) & primary==1)/all.prim.maxart,
              within12m.all = sum(TreatTime==art.eligible &
                                    TreatTime < (TreatTime+twelve.months) )/all.maxart,

              within6m.prim = sum(TreatTime==art.eligible &
                                    TreatTime < (TreatTime+six.months) & primary==1)/all.prim.maxart,
              within6m.all = sum(TreatTime==art.eligible &
                                   TreatTime < (TreatTime+six.months) )/all.maxart

    )

  prim.art.init <- dplyr::select(hhohho.all.prim.art.coverage.var, contains(".prim"))
  all.art.init <- dplyr::select(hhohho.all.prim.art.coverage.var, contains(".all"))

  hhohho.art.initiated.tar.values <- as.numeric(cbind(prim.art.init*100, all.art.init*100))
  #
  #when do we want the retention time
  maxart.ret.timepoint <- difftime(maxart.endtime, sim.start.full, units = "days")/365.242
  maxart.ret.timepoint <- round(as.numeric(maxart.ret.timepoint),0)
  maxart.starttime.ret <- round(as.numeric(difftime(maxart.starttime ,sim.start.full, units = "days")/365),0)

  #### Max Retention ####################################
  max.ret.tar <- length(row.names(max.art.retention.all))
  max.ret.tar.prim <- rep(NA, max.ret.tar)
  max.ret.tar.all <- rep(NA, max.ret.tar)
  max.ret.tar.list <- c(36, 6, 12)
  ret.age.group <- c(18, 150)

  for(ret in 1:max.ret.tar){

  ret.all.prim = ART.retention(datalist = chunk.datalist.test.prim,
                               agegroup = ret.age.group, #all ages
                               ARTtimewindow = c(maxart.starttime.ret, maxart.ret.timepoint),
                               retentiontimeMonths = max.ret.tar.list[ret], #6 months default
                               site="All")$percentage[1]

  ret.all.all = ART.retention(datalist = chunk.datalist.test.all,
                               agegroup = ret.age.group, #all ages
                               ARTtimewindow = c(maxart.starttime.ret, maxart.ret.timepoint),
                               retentiontimeMonths = max.ret.tar.list[ret], #6 months default
                               site="All")$percentage[1]

  max.ret.tar.prim[ret] <- ret.all.prim
  max.ret.tar.all[ret] <- ret.all.all

  }

  hhohho.art.retention.tar.values <- c(max.ret.tar.prim, max.ret.tar.all)

  # # ############### Viral Load none supresssion ############################################
  max.vl.sup.tar.dim <- length(row.names(max.vl.none.suppression.all))
  max.vl.sup.prim <- rep(NA, max.vl.sup.tar.dim)
  max.vl.sup.all <- rep(NA, max.vl.sup.tar.dim)
  max.vl.sup.list <- c(6, 12)

  for(vl in 1:max.vl.sup.tar.dim){

    vl.sup.prim = 100 - vl.suppressed(datalist = chunk.datalist.test.prim,
                                timepoint = maxart.ret.timepoint, vlcutoff = 1000,
                                lessmonths = max.vl.sup.list[vl], site="All")$percentage[[1]]

    vl.sup.all = 100 - vl.suppressed(datalist = chunk.datalist.test.all,
                               timepoint = maxart.ret.timepoint, vlcutoff = 1000,
                               lessmonths = max.vl.sup.list[vl], site="All")$percentage[[1]]

    max.vl.sup.prim[vl] <- ifelse(is.null(vl.sup.prim), NA, vl.sup.prim)
    max.vl.sup.all[vl] <- ifelse(is.null(vl.sup.all), NA, vl.sup.all)

  }

  hhohho.vl.none.suppression.tar.values <- c(max.vl.sup.prim, max.vl.sup.all)
  # #
  # #Mortality AIDS related and not
  max.mort.study.pop.prim <- maxart.study.pop %>%
    summarise(mort.prim = sum(primary == 1 & TOD != Inf & TOD > maxart.starttime & TOD < maxart.endtime) / n(),
            mort.all = sum(TOD != Inf & TOD > maxart.starttime & TOD < maxart.endtime) / n(),

            mort.aids.rel.prim = sum(primary == 1 & TOD != Inf & AIDSDeath == 1 &
                                       TOD > maxart.starttime & TOD < maxart.endtime) / n(),
            mort.aids.rel.all = sum(TOD != Inf & AIDSDeath == 1 &
                                      TOD > maxart.starttime & TOD < maxart.endtime) / n() )


  prim.mortality <- dplyr::select(max.mort.study.pop.prim, contains(".prim"))
  all.mortality <- dplyr::select(max.mort.study.pop.prim, contains(".all"))

  hhohho.mortality.tar.values <- as.numeric(cbind(prim.mortality * 100, all.mortality * 100))

  # # ######### Prevalence for multiple ages. End time 2007-March-31 ###################################
  prev.age2007.list <- format.names(row.names(hhohho.prev), replace = "A")

  swazi.sim.prev.2007.len <- length(prev.age2007.list)
  swazi.sim.prev.2007.f <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.m <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.fm <- rep(NA, swazi.sim.prev.2007.len)

  time.end.2007 <- as.numeric(difftime(as.Date("2007-03-31") ,sim.start.full, units = "days")/365.245)

  for(i in 1:swazi.sim.prev.2007.len){

    split.list.age2007 <- as.numeric(unlist(strsplit(prev.age2007.list[i], "[.]")))
    age.lower.2007 <- split.list.age2007[1]
    age.upper.2007 <- split.list.age2007[2]

    swazi.age.sim.prev <- prevalence.calculator(datalist = sim.datalist,
                                                agegroup = c(age.lower.2007, age.upper.2007),
                                                timepoint = time.end.2007)

    #Gender 0 <- male : 1 <- female
    swazi.sim.prev.2007.m[i] <- swazi.age.sim.prev$pointprevalence[1] * 100
    swazi.sim.prev.2007.f[i] <- swazi.age.sim.prev$pointprevalence[2] * 100
    swazi.sim.prev.2007.fm[i] <- swazi.age.sim.prev$pointprevalence[3] * 100

  }
  hhohho.prev.2007.tar.values <- c(swazi.sim.prev.2007.f, swazi.sim.prev.2007.m,
                                  swazi.sim.prev.2007.fm)

  # # #Swazi Age difference 2011-06-30
  hhohho.age.dif.time <- round(as.numeric(difftime(as.Date("2011-06-30") ,
                                                   sim.start.full, units = "days")/365.245),0)

  agemix.df.hhohho <- agemix.df.maker(sim.datalist)
  pattern.hhohho <- pattern.modeller(dataframe = agemix.df.hhohho, agegroup = c(18, 50),
                              timepoint = hhohho.age.dif.time, timewindow = 1, start = FALSE)

  pattern.hhohho <- as.data.frame(pattern.hhohho)

  pattern.hhohho.gender <- pattern.hhohho %>%
    group_by(Gender) %>%
    summarise(mean.AD = as.numeric(mean(AgeGap)),
              median.AD = as.numeric(median(AgeGap))
              ) %>%
    as.data.frame

  pattern.hhohho.all <- pattern.hhohho %>%
    summarise(mean.AD = as.numeric(mean(AgeGap)),
              median.AD = as.numeric(median(AgeGap))
              )%>%
    as.data.frame

  hhohho.AD.tar.values <- c(pattern.hhohho.gender$mean.AD[2], pattern.hhohho.gender$median.AD[2],
                            pattern.hhohho.gender$mean.AD[1], pattern.hhohho.gender$median.AD[1],
                            pattern.hhohho.all$mean.AD, pattern.hhohho.all$median.AD)


  #collect all the summary values
  sim.summary <- c(hhohho.art.initiated.tar.values, hhohho.art.retention.tar.values,
                   hhohho.vl.none.suppression.tar.values, hhohho.mortality.tar.values,
                   hhohho.prev.2007.tar.values, hhohho.AD.tar.values)

  # #Testing
  #sim.summary <- c(max.art.retention.tar.values, Sys.getpid())

  return(sim.summary)

}

#swazi.art.coverage.year.diff <- swazi.art.coverage[,1:2] - swazi.sim.art.coverage.year



