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

max.swazi.sim.summary.creator <- function(sim.datalist = chunk.datalist.test){

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

  max.all.prim.art.coverage.var <- maxart.study.pop %>%
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

  prim.art.init <- dplyr::select(max.all.prim.art.coverage.var, contains(".prim"))
  all.art.init <- dplyr::select(max.all.prim.art.coverage.var, contains(".all"))

  max.art.initiated.tar.values <- as.numeric(cbind(prim.art.init*100, all.art.init*100))

  # #max.art.coverage.tar <- data.frame(t(prim.art.init)*100, t(all.art.init)*100) %>%
  # #  setNames(c("PrimaryClients","AllClients"))
  # #row.names(max.art.coverage.tar) <- row.names(max.art.initiated.all)
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

  max.art.retention.tar.values <- c(max.ret.tar.prim, max.ret.tar.all)
  # # #max.art.ret.tar <- data.frame(PrimaryClients = max.ret.tar.prim,
  # # #                                 AllClients = max.ret.tar.all)
  # #
  # # #row.names(max.art.ret.tar) <- row.names(max.art.retention.all)
  # #
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

  max.vl.none.suppression.tar.values <- c(max.vl.sup.prim, max.vl.sup.all)
  # #
  # # #max.vl.sup.tar <- data.frame(PrimaryClients = max.vl.sup.prim,
  # # #                              AllClients = max.vl.sup.all)
  # # #row.names(max.vl.sup.tar) <- row.names(max.vl.none.suppression.all)
  # #
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

  max.mortality.tar.values <- as.numeric(cbind(prim.mortality * 100, all.mortality * 100))
  # #
  # # #max.mortality.tar <- data.frame(t(prim.mortality)*100, t(all.mortality)*100) %>%
  # # #  setNames(c("PrimaryClients","AllClients"))
  # # #row.names(max.mortality.tar) <- row.names(max.mortality.all)
  # #
  # ##### Growth rate ########################################
  gr.year.list <- as.numeric(format.names(names(swazi.growth.rate), replace = "X"))

  swazi.growth.rate.tar.values <- rep(NA, length(gr.year.list))

  #Growth rate is calculated per year.
  for(i in 1:length(gr.year.list)){
    from.time <- gr.year.list[i] - 1 - sim.start
    to.time <-  from.time + 1

    swazi.growth.rate.tar.values[i] <- pop.growth.calculator(datalist = sim.datalist,
                                                     timewindow = c(from.time, to.time))
  }
  # #
  # # # #Set names for the growth rate
  # # # #swazi.growth.rate.diff <- swazi.growth.rate - swazi.sim.growth.rate ###############################
  # # #
  # # ##### Incidence for multiple years, 15-49 yrs, time = March each year #############################
  inci.year.list <- as.numeric(format.names(names(swazi.inci.15.49), replace = "X"))
  swazi.inci.15.49.tar.values <- rep(NA, length(inci.year.list))

  for(i in 1:length(inci.year.list)){
    time.start <- inci.year.list[i] - 1 - sim.start
    time.end <- time.start + 1

    swazi.sim.inci <- incidence.calculator(datalist = sim.datalist,
                                              agegroup = c(15, 49), # <= hence not using 50.
                                              timewindow = c(time.start, time.end),
                                              only.active = "No")

    swazi.inci.15.49.tar.values[i] <- swazi.sim.inci$incidence[3] * 100

  }
  # #
  # # #swazi.inci.15.49.diff <- swazi.inci.15.49 - swazi.sim.inci.15.49  ######################
  # #
  # # #Incidence for multiple ages, start and end time set ###############################################
  inci.age2011.list <- format.names(row.names(swazi.inci.2011), replace = "A")
  swazi.sim.inci.2011.f <- rep(NA, length(inci.age2011.list))
  swazi.sim.inci.2011.m <- rep(NA, length(inci.age2011.list))

  time.start.2011 <- as.numeric(difftime(as.Date("2010-12-01") ,sim.start.full, units = "days")/365.245)
  time.end.2011 <- as.numeric(difftime(as.Date("2011-06-30") ,sim.start.full, units = "days")/365.245)

  for(i in 1:length(inci.age2011.list)){

    split.list.age2011 <- as.numeric(unlist(strsplit(inci.age2011.list[i], "[.]")))
    age.lower.2011 <- split.list.age2011[1]
    age.upper.2011 <- split.list.age2011[2]

    swazi.age.sim.inci <- incidence.calculator(datalist = sim.datalist,
                                           agegroup = c(age.lower.2011, age.upper.2011), # <= hence using exact.
                                           timewindow = c(time.start.2011, time.end.2011),
                                           only.active = "No")

    #Gender 0 <- male : 1 <- female
    swazi.sim.inci.2011.m[i] <- swazi.age.sim.inci$incidence[1] * 100
    swazi.sim.inci.2011.f[i] <- swazi.age.sim.inci$incidence[2] * 100

  }

  swazi.inci.2011.tar.values <- c(swazi.sim.inci.2011.f,swazi.sim.inci.2011.m)
  # #
  # # #swazi.sim.inci.2011.fm <- data.frame(swazi.sim.inci.2011.f, swazi.sim.inci.2011.m)
  # # #swazi.inci.age.2011.diff <- swazi.inci.2011[,1:2] - swazi.sim.inci.2011.fm
  # # ##################################################################################################
  # #
  # #
  # # ######### Prevalence for multiple ages. End time 2007-March-31 ###################################
  prev.age2007.list <- format.names(row.names(swazi.prev.2007), replace = "A")

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
  swazi.prev.2007.tar.values <- c(swazi.sim.prev.2007.f,
                                    swazi.sim.prev.2007.m,
                                    swazi.sim.prev.2007.fm)
  # #
  # # #swazi.prev.age.2007.diff <- swazi.prev.2007[,1:3] - swazi.sim.prev.2007
  # # ##################################################################################################
  # #
  # # #Prevalence for multiple ages and July 2006 and March 2007 [SDHS 2006-07] ########################
  swazi.hhohho.prev.age.list <- format.names(row.names(hhohho.prev), replace = "A")

  swazi.hhohho.sim.prev.age.yr.len <- length(swazi.hhohho.prev.age.list)
  swazi.hhohho.sim.prev.age.yr.f <- rep(NA, swazi.hhohho.sim.prev.age.yr.len)
  swazi.hhohho.sim.prev.age.yr.m <- rep(NA, swazi.hhohho.sim.prev.age.yr.len)
  swazi.hhohho.sim.prev.age.yr.fm <- rep(NA, swazi.hhohho.sim.prev.age.yr.len)

  swazi.hhohho.prev.df <- sim.datalist
  swazi.hhohho.prev.df.ptable <- subset(swazi.hhohho.prev.df$ptable, pfacility!="Not Hhohho")
  swazi.hhohho.prev.df$ptable <- swazi.hhohho.prev.df.ptable

  for(i in 1:swazi.hhohho.sim.prev.age.yr.len){

    split.list.id <- as.numeric(unlist(strsplit(swazi.hhohho.prev.age.list[i], "[.]")))

    time.end.age.yr <- round(as.numeric(difftime(as.Date(paste0("2007-03-31")),
                                           sim.start.full, units = "days")/365.245),0)

    age.lower.age.yr <- split.list.id[1]
    age.upper.age.yr <- split.list.id[2]


    swazi.hhohho.sim.age.yr.sim.prev <- prevalence.calculator(datalist = swazi.hhohho.prev.df,
                                                       agegroup = c(age.lower.age.yr, age.upper.age.yr),
                                                       timepoint = time.end.age.yr)

    #Gender 0 <- male : 1 <- female
    swazi.hhohho.sim.prev.age.yr.m[i] <- swazi.hhohho.sim.age.yr.sim.prev$pointprevalence[1] * 100
    swazi.hhohho.sim.prev.age.yr.f[i] <- swazi.hhohho.sim.age.yr.sim.prev$pointprevalence[2] * 100
    swazi.hhohho.sim.prev.age.yr.fm[i] <- swazi.hhohho.sim.age.yr.sim.prev$pointprevalence[3] * 100
  }
  swazi.hhohho.prev.tar.values <- c(swazi.hhohho.sim.prev.age.yr.f,
                                    swazi.hhohho.sim.prev.age.yr.m,
                                    swazi.hhohho.sim.prev.age.yr.fm)
  # #
  # #
  # #Prevalence for multiple ages and time Year ends 31 March Year  ##################################
  prev.age.year.list <- format.names(row.names(swazi.prev.age.year), replace = "A")

  swazi.sim.prev.age.yr.len <- length(prev.age.year.list)
  swazi.sim.prev.age.yr.f <- rep(NA, swazi.sim.prev.age.yr.len)
  swazi.sim.prev.age.yr.m <- rep(NA, swazi.sim.prev.age.yr.len)

  for(i in 1:swazi.sim.prev.age.yr.len){

    split.list.id <- as.numeric(unlist(strsplit(prev.age.year.list[i], "[.]")))

    time.end.age.yr <- as.numeric(difftime(as.Date(paste0(split.list.id[3],"-03-31")) ,
                                           sim.start.full, units = "days")/365.245)

    age.lower.age.yr <- split.list.id[1]
    age.upper.age.yr <- split.list.id[2]


    swazi.sim.age.yr.sim.prev <- prevalence.calculator(datalist = sim.datalist,
                                                agegroup = c(age.lower.age.yr, age.upper.age.yr),
                                                timepoint = time.end.age.yr)

    #Gender 0 <- male : 1 <- female
    swazi.sim.prev.age.yr.m[i] <- swazi.sim.age.yr.sim.prev$pointprevalence[1] * 100
    swazi.sim.prev.age.yr.f[i] <- swazi.sim.age.yr.sim.prev$pointprevalence[2] * 100

  }
  swazi.prev.age.year.tar.values <- c(swazi.sim.prev.age.yr.f, swazi.sim.prev.age.yr.m)
  # #
  # # #swazi.prev.age.age.yr.diff <- swazi.prev.age.year[,1:2] - swazi.sim.prev.age.yr
  # # ##################################################################################################
  # #
  # #
  # # #Swazi Age difference 2011-06-30
  hhohho.age.dif.time <- round(as.numeric(difftime(as.Date("2007-03-30") ,
                                                   sim.start.full, units = "days")/365.245),0)

  chunk.datalist.test.hhohho <- sim.datalist
  hhohho.df <- subset(sim.datalist$ptable, pfacility !="Hhohho")
  chunk.datalist.test.hhohho$ptable <- hhohho.df

  agemix.df.hhohho <- agemix.df.maker(chunk.datalist.test.hhohho)
  pattern.hhohho <- pattern.modeller(dataframe = agemix.df.hhohho, agegroup = c(18, 50),
                              timepoint = hhohho.age.dif.time, timewindow = 1, start = FALSE)

  pattern.hhohho <- as.data.frame(pattern.hhohho)

  pattern.hhohho.df.m <- as.data.frame(subset(pattern.hhohho, Gender == "male"))
  pattern.hhohho.df.f <- as.data.frame(subset(pattern.hhohho, Gender == "female"))

  hhohho.mean.m.18.50.AD <- as.numeric(mean(pattern.hhohho.df.m$AgeGap))
  hhohho.median.m.18.50.AD <- as.numeric(median(pattern.hhohho.df.m$AgeGap))

  hhohho.mean.w.18.50.AD <- as.numeric(mean(pattern.hhohho.df.f$AgeGap))
  hhohho.median.w.18.50.AD <- as.numeric(median(pattern.hhohho.df.f$AgeGap))

  hhohho.mean.fm.18.50.AD <- as.numeric(mean(pattern.hhohho$AgeGap))
  hhohho.median.fm.18.50.AD <- as.numeric(median(pattern.hhohho$AgeGap))


  swazi.hhohho.AD.tar.values <- c(hhohho.mean.w.18.50.AD, hhohho.median.w.18.50.AD,
                                  hhohho.mean.m.18.50.AD, hhohho.median.m.18.50.AD,
                                  hhohho.mean.fm.18.50.AD, hhohho.median.fm.18.50.AD)
  #
  #Swazi Age difference 2011-06-30
  swazi.age.dif.time <- as.numeric(difftime(as.Date("2011-06-30") ,sim.start.full, units = "days")/365.245)

  agemix.df <- agemix.df.maker(sim.datalist)
  swazi.pattern <- pattern.modeller(dataframe = agemix.df, agegroup = c(18, 50),
                              timepoint = swazi.age.dif.time, timewindow = 1, start = FALSE)

  swazi.pattern.df <- as.data.frame(swazi.pattern)

  swazi.mean.fm.18.50.AD <- as.numeric(mean(swazi.pattern.df$AgeGap))
  swazi.median.fm.18.50.AD <- as.numeric(median(swazi.pattern.df$AgeGap))

  swazi.AD.tar.values <- c(swazi.mean.fm.18.50.AD, swazi.median.fm.18.50.AD)


  # Swazi ART retention for multiple ages and two age groups
  ret.age.year.list <- as.numeric(format.names(names(swazi.art.retention), replace = "X"))

  swazi.sim.ret.age.yr.len <- length(ret.age.year.list)

  swazi.sim.ret.age.yr.l15.6 <- rep(NA, swazi.sim.ret.age.yr.len)
  swazi.sim.ret.age.yr.g15.6 <- rep(NA, swazi.sim.ret.age.yr.len)

  swazi.sim.ret.age.yr.l15.12 <- rep(NA, swazi.sim.ret.age.yr.len)
  swazi.sim.ret.age.yr.g15.12 <- rep(NA, swazi.sim.ret.age.yr.len)

  less15 <- c(0,14)
  greater15 <- c(15,150)


  for(i in 1:swazi.sim.ret.age.yr.len){

    uptime.lim <- ret.age.year.list[i] - sim.start

    #compute retention and get the overall value (ignore gender)
    swazi.sim.ret.age.yr.l15.6.var <- ART.retention(datalist = sim.datalist,
                                                agegroup = less15,
                                                ARTtimewindow = c(0, uptime.lim),
                                                retentiontimeMonths = 6, #6 months default
                                                site="All")

    swazi.sim.ret.age.yr.l15.6[i] <- swazi.sim.ret.age.yr.l15.6.var$percentage[
      is.na(swazi.sim.ret.age.yr.l15.6.var$Gender)]


    swazi.sim.ret.age.yr.g15.6.var <- ART.retention(datalist = sim.datalist,
                                                agegroup = greater15,
                                                ARTtimewindow = c(0, uptime.lim),
                                                retentiontimeMonths = 6, #6 months default
                                                site="All")

    swazi.sim.ret.age.yr.g15.6[i] <- swazi.sim.ret.age.yr.g15.6.var$percentage[
      is.na(swazi.sim.ret.age.yr.g15.6.var$Gender)]


    swazi.sim.ret.age.yr.l15.12.var <- ART.retention(datalist = sim.datalist,
                                                     agegroup = less15,
                                                     ARTtimewindow = c(0, uptime.lim),
                                                     retentiontimeMonths = 12, #6 months default
                                                     site="All")

    swazi.sim.ret.age.yr.l15.12[i] <- swazi.sim.ret.age.yr.l15.12.var$percentage[
      is.na(swazi.sim.ret.age.yr.l15.12.var$Gender)]


    swazi.sim.ret.age.yr.g15.12.var <- ART.retention(datalist = sim.datalist,
                                                     agegroup = greater15,
                                                     ARTtimewindow = c(0, uptime.lim),
                                                     retentiontimeMonths = 12, #6 months default
                                                     site="All")

    swazi.sim.ret.age.yr.g15.12[i] <- swazi.sim.ret.age.yr.g15.12.var$percentage[
      is.na(swazi.sim.ret.age.yr.g15.12.var$Gender)]
  }

  swazi.sim.ret.age.yr <- data.frame(matrix(NA, nrow = 0, ncol = swazi.sim.ret.age.yr.len))
  swazi.sim.ret.age.yr <- rbind(swazi.sim.ret.age.yr,
                                swazi.sim.ret.age.yr.l15.6, swazi.sim.ret.age.yr.g15.6,
                                swazi.sim.ret.age.yr.l15.12, swazi.sim.ret.age.yr.g15.12)


  names(swazi.sim.ret.age.yr) <- names(swazi.art.retention)
  #
  swazi.art.retention.tar.values <- as.numeric(unlist(swazi.sim.ret.age.yr))
  # #
  # # #####################################################################################################
  # #
  # # #ART coverage for multiple Years and 15+ yo. Year ends 31 March Year
  art.coverage.year.list <- as.numeric(format.names(row.names(swazi.art.coverage), replace = "15.over."))

  art.coverage.year.len <- length(art.coverage.year.list)
  art.coverage.year.f <- rep(NA, art.coverage.year.len)
  art.coverage.year.m <- rep(NA, art.coverage.year.len)
  art.coverage.year.fm <- rep(NA, art.coverage.year.len)

  for(i in 1:art.coverage.year.len){

    time.end.cov <- as.numeric(art.coverage.year.list[i] - sim.start)

    ARTcov <- ART.coverage.calculator(datalist = sim.datalist, agegroup = c(18, 50),
                                      timepoint = time.end.cov, site="All")

    #Gender 0 <- male : 1 <- female
    art.coverage.year.m[i] <- ARTcov$ART.coverage[1] * 100
    art.coverage.year.f[i] <- ARTcov$ART.coverage[2] * 100
    art.coverage.year.fm[i] <- ARTcov$ART.coverage[3] * 100

  }
  swazi.art.coverage.tar.values <- c(art.coverage.year.f, art.coverage.year.m, art.coverage.year.fm)
  #
  #collect all the summary values
  sim.summary <- c(max.art.initiated.tar.values, max.art.retention.tar.values,
                   max.vl.none.suppression.tar.values, max.mortality.tar.values,
                   swazi.growth.rate.tar.values, swazi.inci.15.49.tar.values,
                   swazi.inci.2011.tar.values, swazi.prev.2007.tar.values,
                   swazi.prev.age.year.tar.values, swazi.hhohho.prev.tar.values,
                   swazi.AD.tar.values, swazi.hhohho.AD.tar.values,
                   swazi.art.retention.tar.values, swazi.art.coverage.tar.values)

  # #Testing
  #sim.summary <- c(max.art.retention.tar.values, Sys.getpid())

  return(sim.summary)

}



#swazi.art.coverage.year.diff <- swazi.art.coverage[,1:2] - swazi.sim.art.coverage.year



