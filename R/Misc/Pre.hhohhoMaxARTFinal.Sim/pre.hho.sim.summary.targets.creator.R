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


  #when do we want the retention time
  maxart.ret.timepoint <- difftime(maxart.endtime, sim.start.full, units = "days")/days.in.yr
  maxart.ret.timepoint <- round(as.numeric(maxart.ret.timepoint),0)
  maxart.starttime.ret <- round(as.numeric(difftime(maxart.starttime, sim.start.full, units = "days")/days.in.yr),0)
  ret.end <- round(as.numeric(difftime(maxart.endtime,maxart.starttime, units = "days"))/days.in.yr,0) * 12

  #creator a range of summary statistics
  sim.datalist$ptable$age <- sim.datalist$itable$maxart.starttime[1] - sim.datalist$ptable$TOB
  sim.datalist$ptable$age.end.study <- maxart.ret.timepoint - sim.datalist$ptable$TOB

  #turn 18 during or before MaxART end and Got infected during or before MaxART end
  #get the MaxART population group
  maxart.study.pop <- dplyr::filter(sim.datalist$ptable,
                                      age.end.study >=18 & InfectTime!=Inf &
                                      InfectTime < maxart.ret.timepoint)

  #Total clients
  all.maxart <- nrow(maxart.study.pop)

  #alldf
  chunk.datalist.test.all <- sim.datalist
  chunk.datalist.test.all$ptable <- maxart.study.pop

  ### MaxART ART coverage ####################
  hhohho.all.art.coverage.var <- maxart.study.pop %>%
     dplyr::filter(TreatTime!=Inf) %>%
     summarise(all.all = n()/all.maxart)

  all.art.init <- dplyr::select(hhohho.all.art.coverage.var, contains(".all"))
  hhohho.art.initiated.tar.values <- as.numeric(all.art.init*100)


  #### Max Retention ####################################
  max.ret.tar <- length(row.names(max.art.retention.all))
  max.ret.tar.all <- rep(NA, max.ret.tar)
  max.ret.tar.list <- c(ret.end, 6, 12)
  ret.age.group <- c(18, 150)

  for(ret in 1:max.ret.tar){

  ret.all.all = ART.retention(datalist = chunk.datalist.test.all,
                               agegroup = ret.age.group, #all ages
                               ARTtimewindow = c(maxart.starttime.ret, maxart.ret.timepoint),
                               retentiontimeMonths = max.ret.tar.list[ret], #6 months default
                               site="All")$percentage[1]

  max.ret.tar.all[ret] <- ret.all.all

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
      maxart.ret.timepoint <- shims.vl.timepoint
      vl.datalist <- sim.datalist
      }else{vl.datalist <- chunk.datalist.test.all}

    vl.sup.all <- vl.suppressed(datalist = vl.datalist,
                                timepoint = maxart.ret.timepoint, vlcutoff = 1000,
                               lessmonths = max.vl.sup.list[vl], site="All")

    vl.sup.all <- subset(vl.sup.all, is.na(Gender))$Percentage[[1]]

    max.vl.sup.all[vl] <- ifelse(vl != 3, 100 - vl.sup.all, vl.sup.all)

  }

  hhohho.vl.none.suppression.tar.values <- max.vl.sup.all
  #########

  ###Mortality AIDS related and not
  maxart.starttime.mort <- round(as.numeric(difftime(maxart.starttime,sim.start.full, units = "days")/days.in.yr),0)
  maxart.endtime.mort <- round(as.numeric(difftime(maxart.endtime,sim.start.full, units = "days")/days.in.yr),0)

  max.mort.study.pop.prim <- maxart.study.pop %>%
    summarise(
            mort.all = sum(TOD != Inf & TOD > maxart.starttime.mort
                           & TOD < maxart.endtime.mort, na.rm = TRUE ) / n(),
            mort.aids.rel.all = sum(TOD != Inf & AIDSDeath == 1 &
                                      TOD > maxart.starttime.mort &
                                      TOD < maxart.endtime.mort, na.rm = TRUE) / n() )

  all.mortality <- dplyr::select(max.mort.study.pop.prim, contains(".all"))

  hhohho.mortality.tar.values <- as.numeric(all.mortality * 100)


  # ##### Growth rate ########################################
  gr.year.list <- as.numeric(format.names(names(hhohho.growth.rate), replace = "X"))

  hhohho.growth.rate.tar.values <- rep(NA, length(gr.year.list))

  #Growth rate is calculated per 10 year.
  for(i in 1:length(gr.year.list)){
    ifelse(gr.year.list[i] == 2012, t.change <- 5, t.change <- 10)

    from.time <- gr.year.list[i] - t.change - sim.start
    to.time <-  from.time + t.change

    hhohho.growth.rate.tar.values[i] <- pop.growth.calculator(datalist = sim.datalist,
                                                             timewindow = c(from.time, to.time))


  }

  # # ######### Prevalence for multiple ages. End time 2007-March-31 ###################################
  prev.age2007.list <- format.names(row.names(hhohho.prev), replace = "A")

  swazi.sim.prev.2007.len <- length(prev.age2007.list)
  swazi.sim.prev.2007.f <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.m <- rep(NA, swazi.sim.prev.2007.len)
  swazi.sim.prev.2007.fm <- rep(NA, swazi.sim.prev.2007.len)

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
    swazi.sim.prev.2007.m[i] <- swazi.age.sim.prev$pointprevalence[1] * 100
    swazi.sim.prev.2007.f[i] <- swazi.age.sim.prev$pointprevalence[2] * 100
    swazi.sim.prev.2007.fm[i] <- swazi.age.sim.prev$pointprevalence[3] * 100

  }
  hhohho.prev.2007.tar.values <- c(swazi.sim.prev.2007.f, swazi.sim.prev.2007.m,
                                  swazi.sim.prev.2007.fm)


  ########  SHIMS2 prevalence ########
  #use the vl time point already calculated: shims.vl.timepoint
  prev.age.shims2.list <- format.names(row.names(hhohho.shims2.prev), replace = "A")

  hho.sim.prev.2017.len <- length(prev.age.shims2.list)
  hho.sim.prev.2017.fm <- rep(NA, hho.sim.prev.2017.len)

  for(i in 1:hho.sim.prev.2017.len){

    split.list.age2017 <- as.numeric(unlist(strsplit(prev.age.shims2.list[i], "[.]")))
    age.lower.2017 <- split.list.age2017[1]
    age.upper.2017 <- split.list.age2017[2]

    hho.shim2.age.sim.prev <- prevalence.calculator(datalist = sim.datalist,
                                                agegroup = c(age.lower.2017, age.upper.2017),
                                                timepoint = shims.vl.timepoint)

    #Gender 0 <- male : 1 <- female
    hho.sim.prev.2017.fm[i] <- hho.shim2.age.sim.prev$pointprevalence[3] * 100

  }
  hho.prev.2017.tar.values <- c(hho.sim.prev.2017.fm)

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
                   hhohho.vl.none.suppression.tar.values,
                   hhohho.mortality.tar.values,
                   hhohho.prev.2007.tar.values,
                   hho.prev.2017.tar.values,
                   hhohho.AD.tar.values)


  return(sim.summary)

}

