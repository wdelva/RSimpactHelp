#This will be called through source()
#function to clean up the column or row names,
format.names <- function(x,replace = "X",...){
  gsub(replace,'',x)
}

#Recruitment
#maxart.study.pop <- recruit.study.clients(chunk.datalist.test)

#Dummy with no real recruitment
days.in.yr <- 365.245
two.weeks <- 14/days.in.yr
twelve.months <- 12/12
six.months <- 6/12
init.immediate <- 0
init.not.eligible <- -6/12
more.months <- twelve.months + six.months

#Testing from saved data
#load("temp/chunk.datalist.WVDBdJaAev.rda")
#sim.datalist <- chunk.datalist.test

pre.hhohho.sim.summary.creator <- function(sim.datalist = chunk.datalist.test){

  # #when do we want the retention time
  maxart.ret.timepoint <- difftime(study.endtime, sim.start.full, units = "days")/days.in.yr
  maxart.ret.timepoint <- maxart.ret.timepoint + tail(sim.datalist$itable$maxart.stepinterval,1)
  maxart.starttime.ret <- as.numeric(difftime(maxart.starttime, sim.start.full, units = "days")/days.in.yr)
  ret.end <- maxart.starttime.ret + (as.numeric(difftime(maxart.endtime,maxart.starttime, units = "days"))/days.in.yr +
                                       tail(sim.datalist$itable$maxart.stepinterval,1))

  maxart.ret.timepoint <- as.numeric(maxart.ret.timepoint)
  #To cater for the simulation without study
  if(maxart.starttime < sim.start.full){
    #creator a range of summary statistics
    sim.datalist$ptable$age <- maxart.ret.timepoint - sim.datalist$ptable$TOB
  }else{
    sim.datalist$ptable$age <- 10  #Give a dummy age No one is in the study
  }
  #get the MaxART population group
  maxart.study.pop <- dplyr::filter(sim.datalist$ptable,  age >= 18 & InfectTime!=Inf &
                                      InfectTime < maxart.ret.timepoint)

  #alldf
  chunk.datalist.test.all <- sim.datalist
  chunk.datalist.test.all$ptable <- maxart.study.pop

  ### MaxART ART coverage/Initiation ####################
  hhohho.all.art.coverage.var <- ART.coverage.calculator(datalist = chunk.datalist.test.all,
                                                         agegroup = c(18,150),
                                                         timepoint = maxart.ret.timepoint,
                                                         site = "All")

  #Check Hhohho ART coverage in this simulation
  hhohho.art.initiated.tar.values <- ifelse(nrow(subset(hhohho.all.art.coverage.var, is.na(Gender)))==1,
                                            subset(hhohho.all.art.coverage.var, is.na(Gender))$ART.coverage*100,
                                            NA)

  #### Max Retention ####################################
  max.ret.tar <- length(row.names(max.art.retention.all))
  max.ret.tar.all <- rep(NA, max.ret.tar)
  max.ret.tar.list <- c(ret.end, 6, 12)
  ret.age.group <- c(18, 150)

  for(ret in 1:max.ret.tar){

  ret.all.all <- ART.retention(datalist = chunk.datalist.test.all,
                               agegroup = ret.age.group, #all ages
                               ARTtimewindow = c(maxart.starttime.ret, maxart.ret.timepoint),
                               retentiontimeMonths = max.ret.tar.list[ret], #6 months default
                               site="All")

  # This is for all gender
  max.ret.tar.all[ret] <- ifelse(nrow(subset(ret.all.all, is.na(Gender))) ==1,
                                 subset(ret.all.all, is.na(Gender))$percentage, NA)

  #for gender specific  1 <- Female  0 <- Male
  }

  hhohho.art.retention.tar.values <- max.ret.tar.all

  ################# Viral Load none supresssion ###########################################

  #Shims2 start and end date August 2016 and March 2017
  shims2.start <- as.Date("2016-08-31")
  shims2.end <- as.Date("2017-03-31")
  shims2.vl.time <- shims2.start + floor(difftime(shims2.end, shims2.start)/2)
  shims.vl.timepoint <- as.numeric(difftime(shims2.vl.time,sim.start.full, units = "days")/days.in.yr)

  max.vl.sup.tar.dim <- length(row.names(max.vl.none.suppression.all))
  max.vl.sup.all <- rep(NA, max.vl.sup.tar.dim)
  max.vl.sup.list <- c(6, 12, 12)

  for(vl in 1:max.vl.sup.tar.dim){

    if(vl==3){
      tp.var <- shims.vl.timepoint
      df.var <- sim.datalist
      min.age.vl.sup <- 15
      max.age.vl.sup <- 150
    }else{
      tp.var <- maxart.ret.timepoint
      df.var <- chunk.datalist.test.all
      min.age.vl.sup <- 18
      max.age.vl.sup <- 150
    }

    vl.sup.all <- vl.suppressed(datalist = df.var,
                                agegroup = c(min.age.vl.sup, max.age.vl.sup),
                                timepoint = tp.var, vlcutoff = 1000,
                                lessmonths = max.vl.sup.list[vl], site="All")

    vl.sup.all <- subset(vl.sup.all, is.na(Gender))$Percentage

    max.vl.sup.all[vl] <- ifelse(vl != 3, 100 - vl.sup.all, vl.sup.all)

  }

  hhohho.vl.none.suppression.tar.values <- max.vl.sup.all

  # ##### Growth rate ########################################
  gr.year.list <- as.numeric(format.names(names(hhohho.growth.rate), replace = "X"))

  hhohho.growth.rate.tar.values <- rep(NA, length(gr.year.list))

  #Growth rate is calculated per 10 year.
  for(i in 1:length(gr.year.list)){
    t.change <- 2 #ifelse(gr.year.list[i] == 2012, 5, 10)

    #The datalist population log is from sim.start = 1 : sim.duration

    from.time <- gr.year.list[i] - t.change - sim.start
    to.time <-  from.time + t.change

    hhohho.growth.rate.tar.values[i] <- pop.growth.calculator(datalist = sim.datalist,
                                                             timewindow = c(from.time, to.time))

  }


  ###Incidence for multiple ages, start and end time set ###############################################
  inci.age2016.list <- format.names(row.names(hhohho.inci.2016), replace = "A")
  #female + CI
  hho.inci.2016.f <- rep(NA, length(inci.age2016.list))
  hho.inci.2016.m <- rep(NA, length(inci.age2016.list))
  hho.inci.2016.fm <- rep(NA, length(inci.age2016.list))
  #male + CI
  hho.inci.2016.fl <- rep(NA, length(inci.age2016.list))
  hho.inci.2016.ml <- rep(NA, length(inci.age2016.list))
  hho.inci.2016.fml <- rep(NA, length(inci.age2016.list))
  #male and female  + CI
  hho.inci.2016.fu <- rep(NA, length(inci.age2016.list))
  hho.inci.2016.mu <- rep(NA, length(inci.age2016.list))
  hho.inci.2016.fmu <- rep(NA, length(inci.age2016.list))

  #Time period
  time.start.2016 <- as.numeric(difftime(as.Date("2014-09-01") ,sim.start.full, units = "days")/days.in.yr)
  time.end.2016 <- as.numeric(difftime(as.Date("2016-11-03") ,sim.start.full, units = "days")/days.in.yr)


  for(i in 1:length(inci.age2016.list)){

    split.list.age <- as.numeric(unlist(strsplit(inci.age2016.list[i], "[.]")))
    age.lower <- split.list.age[1]
    age.upper <- split.list.age[2]

    hho.sim.inci <- incidence.calculator(datalist = sim.datalist,
                                               agegroup = c(age.lower, age.upper), # <= hence using exact.
                                               timewindow = c(time.start.2016, time.end.2016),
                                               only.active = "No")

    #Gender 0 <- male : 1 <- female
    #Subset a df summary statistic to determine gender

    #male
    hho.inci.2016.m[i] <- ifelse(nrow(subset(hho.sim.inci, Gender == 0)) == 1,
                                 subset(hho.sim.inci, Gender == 0)$incidence * 100, NA)
    hho.inci.2016.ml[i] <- ifelse(nrow(subset(hho.sim.inci, Gender == 0)) == 1,
                                  subset(hho.sim.inci, Gender == 0)$incidence.95.ll * 100, NA)
    hho.inci.2016.mu[i] <- ifelse(nrow(subset(hho.sim.inci, Gender == 0)) == 1,
                                  subset(hho.sim.inci, Gender == 0)$incidence.95.ul * 100, NA)

    #female
    hho.inci.2016.f[i] <- ifelse(nrow(subset(hho.sim.inci, Gender == 1)) == 1,
                                 subset(hho.sim.inci, Gender == 1)$incidence * 100, NA)
    hho.inci.2016.fl[i] <- ifelse(nrow(subset(hho.sim.inci, Gender == 1)) == 1,
                                  subset(hho.sim.inci, Gender == 1)$incidence.95.ll * 100, NA)
    hho.inci.2016.fu[i] <- ifelse(nrow(subset(hho.sim.inci, Gender == 1)) == 1,
                                  subset(hho.sim.inci, Gender == 1)$incidence.95.ul * 100, NA)

    #All gender
    hho.inci.2016.fm[i] <- ifelse(nrow(subset(hho.sim.inci, is.na(Gender))) == 1,
                                  subset(hho.sim.inci, is.na(Gender))$incidence * 100, NA)
    hho.inci.2016.fml[i] <- ifelse(nrow(subset(hho.sim.inci, is.na(Gender))) == 1,
                                   subset(hho.sim.inci, is.na(Gender))$incidence.95.ll * 100, NA)
    hho.inci.2016.fmu[i] <- ifelse(nrow(subset(hho.sim.inci, is.na(Gender))) == 1,
                                   subset(hho.sim.inci, is.na(Gender))$incidence.95.ul * 100, NA)

  }

  hho.inci.2016.tar.values <- c(hho.inci.2016.f, hho.inci.2016.fu, hho.inci.2016.fl,
                                hho.inci.2016.m, hho.inci.2016.mu, hho.inci.2016.ml,
                                hho.inci.2016.fm, hho.inci.2016.fmu, hho.inci.2016.fml)


  #####  Incidence for multiple ages, year, start and end time set ###############
  inci.age.yr.list <- format.names(row.names(hhohho.inci.year), replace = "A")
  #female + CI
  hho.inci.yr.f <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  hho.inci.yr.m <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  hho.inci.yr.fm <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  #male + CI
  hho.inci.yr.fl <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  hho.inci.yr.ml <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  hho.inci.yr.fml <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  #male and female  + CI
  hho.inci.yr.fu <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  hho.inci.yr.mu <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))
  hho.inci.yr.fmu <- rep(NA, length(inci.age.yr.list) * length(hhohho.inci.year.list))

  jj <- 1
  for(dt in hhohho.inci.year.list){

    prev.dt <- dt - 1
    #Time period
    time.start.yr <- as.numeric(difftime(as.Date(paste0(prev.dt,"-03-01")) ,sim.start.full, units = "days")/days.in.yr)
    time.end.yr <- as.numeric(difftime(as.Date(paste0(dt,"-02-28")) ,sim.start.full, units = "days")/days.in.yr)

    for(i in 1:length(inci.age2016.list)){

      split.list.age <- as.numeric(unlist(strsplit(inci.age2016.list[i], "[.]")))
      age.lower <- split.list.age[1]
      age.upper <- split.list.age[2]

      hho.sim.inci.yr <- incidence.calculator(datalist = sim.datalist,
                                           agegroup = c(age.lower, age.upper), # <= hence using exact.
                                           timewindow = c(time.start.yr, time.end.yr),
                                           only.active = "No")

      #Gender 0 <- male : 1 <- female
      hho.inci.yr.m[i] <- ifelse(nrow(subset(hho.sim.inci.yr, Gender == 0)) == 1,
                                 subset(hho.sim.inci.yr, Gender == 0)$incidence * 100, NA)
      hho.inci.yr.ml[i] <- ifelse(nrow(subset(hho.sim.inci.yr, Gender == 0)) == 1,
                                  subset(hho.sim.inci.yr, Gender == 0)$incidence.95.ll * 100, NA)
      hho.inci.yr.mu[i] <- ifelse(nrow(subset(hho.sim.inci.yr, Gender == 0)) == 1,
                                  subset(hho.sim.inci.yr, Gender == 0)$incidence.95.ul * 100, NA)

      hho.inci.yr.f[i] <- ifelse(nrow(subset(hho.sim.inci.yr, Gender == 1)) == 1,
                                 subset(hho.sim.inci.yr, Gender == 1)$incidence * 100, NA)
      hho.inci.yr.fl[i] <- ifelse(nrow(subset(hho.sim.inci.yr, Gender == 1)) == 1,
                                  subset(hho.sim.inci.yr, Gender == 1)$incidence.95.ll * 100, NA)
      hho.inci.yr.fu[i] <- ifelse(nrow(subset(hho.sim.inci.yr, Gender == 1)) == 1,
                                  subset(hho.sim.inci.yr, Gender == 1)$incidence.95.ul * 100, NA)

      hho.inci.yr.fm[i] <- ifelse(nrow(subset(hho.sim.inci.yr, is.na(Gender))) == 1,
                                  subset(hho.sim.inci.yr, is.na(Gender))$incidence * 100, NA)
      hho.inci.yr.fml[i] <- ifelse(nrow(subset(hho.sim.inci.yr, is.na(Gender))) == 1,
                                   subset(hho.sim.inci.yr, is.na(Gender))$incidence.95.ll * 100, NA)
      hho.inci.yr.fmu[i] <- ifelse(nrow(subset(hho.sim.inci.yr, is.na(Gender))) == 1,
                                   subset(hho.sim.inci.yr, is.na(Gender))$incidence.95.ul * 100, NA)

      jj <- jj + 1

    }

  }

  hho.inci.yr.tar.values <- c(hho.inci.yr.f, hho.inci.yr.fu, hho.inci.yr.fl,
                                hho.inci.yr.m, hho.inci.yr.mu, hho.inci.yr.ml,
                                hho.inci.yr.fm, hho.inci.yr.fmu, hho.inci.yr.fml)


  ################ END INCI ################################

  # # ######### Prevalence for multiple ages. End time 2007-March-31 ####################
  prev.age2007.list <- format.names(row.names(hhohho.prev), replace = "A")
  swazi.sim.prev.2007.len <- length(prev.age2007.list)

  swazi.sim.prev.2007.f <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.fl <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.fu <- rep(NA, swazi.sim.prev.2007.len)

  swazi.sim.prev.2007.m <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.ml <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.mu <- rep(NA, swazi.sim.prev.2007.len)

  swazi.sim.prev.2007.fm <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.fml <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.fmu <- rep(NA, swazi.sim.prev.2007.len)

  s.start <- as.Date("2006-07-31")
  s.end <- as.Date("2007-03-31")
  prev.time <- s.start + floor(difftime(s.end, s.start)/2)
  time.end.2007 <- as.numeric(difftime(prev.time,sim.start.full, units = "days")/days.in.yr)

  for(i in 1:swazi.sim.prev.2007.len){

    split.list.age2007 <- as.numeric(unlist(strsplit(prev.age2007.list[i], "[.]")))
    age.lower.2007 <- split.list.age2007[1]
    age.upper.2007 <- split.list.age2007[2]

    swazi.age.sim.prev <- prevalence.calculator(datalist = sim.datalist,
                                                agegroup = c(age.lower.2007, age.upper.2007),
                                                timepoint = time.end.2007)

    #Gender 0 <- male : 1 <- female
    swazi.sim.prev.2007.m[i] <- ifelse(nrow(subset(swazi.age.sim.prev, Gender == 0)) == 1,
                                       subset(swazi.age.sim.prev, Gender == 0)$pointprevalence * 100, NA)
    swazi.sim.prev.2007.ml[i] <- ifelse(nrow(subset(swazi.age.sim.prev, Gender == 0)) == 1,
                                        subset(swazi.age.sim.prev, Gender == 0)$pointprevalence.95.ll * 100, NA)
    swazi.sim.prev.2007.mu[i] <- ifelse(nrow(subset(swazi.age.sim.prev, Gender == 0)) == 1,
                                        subset(swazi.age.sim.prev, Gender == 0)$pointprevalence.95.ul * 100, NA)

    swazi.sim.prev.2007.f[i] <- ifelse(nrow(subset(swazi.age.sim.prev, Gender == 1)) == 1,
                                       subset(swazi.age.sim.prev, Gender == 1)$pointprevalence * 100, NA)
    swazi.sim.prev.2007.fl[i] <- ifelse(nrow(subset(swazi.age.sim.prev, Gender == 1)) == 1,
                                        subset(swazi.age.sim.prev, Gender == 1)$pointprevalence.95.ll * 100, NA)
    swazi.sim.prev.2007.fu[i] <- ifelse(nrow(subset(swazi.age.sim.prev, Gender == 1)) == 1,
                                        subset(swazi.age.sim.prev, Gender == 1)$pointprevalence.95.ul * 100, NA)

    swazi.sim.prev.2007.fm[i] <- ifelse(nrow(subset(swazi.age.sim.prev, is.na(Gender))) == 1,
                                        subset(swazi.age.sim.prev, is.na(Gender))$pointprevalence * 100, NA)
    swazi.sim.prev.2007.fml[i] <- ifelse(nrow(subset(swazi.age.sim.prev, is.na(Gender))) == 1,
                                         subset(swazi.age.sim.prev, is.na(Gender))$pointprevalence.95.ll * 100, NA)
    swazi.sim.prev.2007.fmu[i] <- ifelse(nrow(subset(swazi.age.sim.prev, is.na(Gender))) == 1,
                                         subset(swazi.age.sim.prev, is.na(Gender))$pointprevalence.95.ul * 100, NA)

  }
  hhohho.prev.2007.tar.values <- c(swazi.sim.prev.2007.f, swazi.sim.prev.2007.fu, swazi.sim.prev.2007.fl,
                                   swazi.sim.prev.2007.m, swazi.sim.prev.2007.mu, swazi.sim.prev.2007.ml,
                                  swazi.sim.prev.2007.fm, swazi.sim.prev.2007.fmu, swazi.sim.prev.2007.fml)

  ########## Prevalence for multiple ages and years. ###################################
  prev.age.yr.list <- format.names(row.names(hhohho.prev.year), replace = "A")
  swazi.sim.prev.yr.len <- length(prev.age.yr.list) * length(hhohho.prev.year.list)

  swazi.sim.prev.yr.f <- rep(NA, swazi.sim.prev.yr.len)
  swazi.sim.prev.yr.fl <- rep(NA, swazi.sim.prev.yr.len)
  swazi.sim.prev.yr.fu <- rep(NA, swazi.sim.prev.yr.len)

  swazi.sim.prev.yr.m <- rep(NA, swazi.sim.prev.yr.len)
  swazi.sim.prev.yr.ml <- rep(NA, swazi.sim.prev.yr.len)
  swazi.sim.prev.yr.mu <- rep(NA, swazi.sim.prev.yr.len)

  swazi.sim.prev.yr.fm <- rep(NA, swazi.sim.prev.yr.len)
  swazi.sim.prev.yr.fml <- rep(NA, swazi.sim.prev.yr.len)
  swazi.sim.prev.yr.fmu <- rep(NA, swazi.sim.prev.yr.len)

  for(dt.p in hhohho.prev.year.list){
    #get through the yr list
    prev.dt.p <- dt.p - 1

    s.start.yr <- as.Date(paste0(prev.dt.p,"-03-01"))
    s.end.yr <- as.Date(paste0(dt.p, "-02-28"))
    prev.time.yr <- s.start.yr + floor(difftime(s.end.yr, s.start.yr)/2)
    time.end.yr <- as.numeric(difftime(prev.time.yr,sim.start.full, units = "days")/days.in.yr)

    for(i in 1:swazi.sim.prev.yr.len){

      split.list.age.yr <- as.numeric(unlist(strsplit(prev.age.yr.list[i], "[.]")))
      age.lower.yr <- split.list.age.yr[1]
      age.upper.yr <- split.list.age.yr[2]

      swazi.age.sim.prev.yr <- prevalence.calculator(datalist = sim.datalist,
                                                  agegroup = c(age.lower.yr, age.upper.yr),
                                                  timepoint = time.end.yr)

      #Gender 0 <- male : 1 <- female
      swazi.sim.prev.yr.m[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, Gender == 0)) == 1,
                                       subset(swazi.age.sim.prev.yr, Gender == 0)$pointprevalence * 100, NA)
      swazi.sim.prev.yr.ml[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, Gender == 0)) == 1,
                                         subset(swazi.age.sim.prev.yr, Gender == 0)$pointprevalence.95.ll * 100, NA)
      swazi.sim.prev.yr.mu[i] <-  ifelse(nrow(subset(swazi.age.sim.prev.yr, Gender == 0)) == 1,
                                         subset(swazi.age.sim.prev.yr, Gender == 0)$pointprevalence.95.ul * 100, NA)

      swazi.sim.prev.yr.f[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, Gender == 1)) == 1,
                                       subset(swazi.age.sim.prev.yr, Gender == 1)$pointprevalence * 100, NA)
      swazi.sim.prev.yr.fl[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, Gender == 1)) == 1,
                                        subset(swazi.age.sim.prev.yr, Gender == 1)$pointprevalence.95.ll * 100, NA)
      swazi.sim.prev.yr.fu[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, Gender == 1)) == 1,
                                        subset(swazi.age.sim.prev.yr, Gender == 1)$pointprevalence.95.ul * 100, NA)

      swazi.sim.prev.yr.fm[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, is.na(Gender))) == 1,
                                        subset(swazi.age.sim.prev.yr, is.na(Gender))$pointprevalence * 100, NA)
      swazi.sim.prev.yr.fml[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, is.na(Gender))) == 1,
                                         subset(swazi.age.sim.prev.yr, is.na(Gender))$pointprevalence.95.ll * 100, NA)
      swazi.sim.prev.yr.fmu[i] <- ifelse(nrow(subset(swazi.age.sim.prev.yr, is.na(Gender))) == 1,
                                         subset(swazi.age.sim.prev.yr, is.na(Gender))$pointprevalence.95.ul * 100, NA)

    }
  }

  hhohho.prev.yr.tar.values <- c(swazi.sim.prev.yr.f, swazi.sim.prev.yr.fu, swazi.sim.prev.yr.fl,
                                   swazi.sim.prev.yr.m, swazi.sim.prev.yr.mu, swazi.sim.prev.yr.ml,
                                   swazi.sim.prev.yr.fm, swazi.sim.prev.yr.fmu, swazi.sim.prev.yr.fml)

  ##########################  END PREV #######################################


  ########  SHIMS2 prevalence ########
  #use the vl time point already calculated: shims.vl.timepoint
  prev.age.shims2.list <- format.names(row.names(hhohho.shims2.prev), replace = "A")

  hho.sim.prev.2017.len <- length(prev.age.shims2.list)
  hho.sim.prev.2017.fm <- rep(NA, hho.sim.prev.2017.len)
  hho.sim.prev.2017.fml <- rep(NA, hho.sim.prev.2017.len)
  hho.sim.prev.2017.fmu <- rep(NA, hho.sim.prev.2017.len)

  for(i in 1:hho.sim.prev.2017.len){

    split.list.age2017 <- as.numeric(unlist(strsplit(prev.age.shims2.list[i], "[.]")))
    age.lower.2017 <- split.list.age2017[1]
    age.upper.2017 <- split.list.age2017[2]

    hho.shim2.age.sim.prev <- prevalence.calculator(datalist = sim.datalist,
                                                agegroup = c(age.lower.2017, age.upper.2017),
                                                timepoint = shims.vl.timepoint)

    #Gender 0 <- male : 1 <- female
    hho.sim.prev.2017.fm[i] <- ifelse(nrow(subset(hho.shim2.age.sim.prev, is.na(Gender))) == 1,
                                      subset(hho.shim2.age.sim.prev, is.na(Gender))$pointprevalence * 100, NA)
    hho.sim.prev.2017.fml[i] <- ifelse(nrow(subset(hho.shim2.age.sim.prev, is.na(Gender))) == 1,
                                       subset(hho.shim2.age.sim.prev, is.na(Gender))$pointprevalence.95.ll * 100, NA)
    hho.sim.prev.2017.fmu[i] <- ifelse(nrow(subset(hho.shim2.age.sim.prev, is.na(Gender))) == 1,
                                       subset(hho.shim2.age.sim.prev, is.na(Gender))$pointprevalence.95.ul * 100, NA)

  }
  hho.prev.2017.tar.values <- c(hho.sim.prev.2017.fm, hho.sim.prev.2017.fml, hho.sim.prev.2017.fmu)

  ###Swazi Age difference 2011-06-30

  #The survey was between December 2010 and June 2011
  age.d.start <- as.Date("2010-12-31")
  age.d.end <- as.Date("2011-06-30")
  age.diff.time <- age.d.start + floor(difftime(age.d.end, age.d.start)/2)

  hhohho.age.dif.time <- round(as.numeric(difftime(age.diff.time ,
                                                   sim.start.full, units = "days")/days.in.yr),0)

  agemix.df.hhohho <- agemix.df.maker(sim.datalist)
  pattern.hhohho <- pattern.modeller(dataframe = agemix.df.hhohho, agegroup = c(18, 50),
                              timepoint = hhohho.age.dif.time, timewindow = 1, start = FALSE)[[1]]

  pattern.hhohho <- as.data.frame(pattern.hhohho)

  pattern.hhohho.gender <- pattern.hhohho %>%
    group_by(Gender) %>%
    summarise(mean.AD = as.numeric(mean(AgeGap, na.rm = TRUE)),
              median.AD = as.numeric(median(AgeGap, na.rm = TRUE))
              ) %>%
    as.data.frame

  pattern.hhohho.all <- pattern.hhohho %>%
    summarise(mean.AD = as.numeric(mean(AgeGap, na.rm = TRUE)),
              median.AD = as.numeric(median(AgeGap, na.rm = TRUE))
              )%>%
    as.data.frame

  hhohho.AD.tar.values <- c(pattern.hhohho.gender$mean.AD[2], pattern.hhohho.gender$median.AD[2],
                            pattern.hhohho.gender$mean.AD[1], pattern.hhohho.gender$median.AD[1],
                            pattern.hhohho.all$mean.AD, pattern.hhohho.all$median.AD)


  #collect all the summary values
  sim.summary <- c(hhohho.growth.rate.tar.values,
                   hhohho.art.initiated.tar.values,
                   hhohho.art.retention.tar.values,
                   hho.inci.2016.tar.values,
                   hho.inci.yr.tar.values,
                   hhohho.vl.none.suppression.tar.values,
                   hhohho.prev.2007.tar.values,
                   hhohho.prev.yr.tar.values,
                   hho.prev.2017.tar.values,

                   hhohho.AD.tar.values)

  return(sim.summary)

}

