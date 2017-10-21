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

#Testing from saved data
load("temp/chunk.datalist.WVDBdJaAev.rda")
sim.datalist <- chunk.datalist.test

pre.hhohho.sim.summary.creator <- function(sim.datalist = chunk.datalist.test){

  #creator a range of summary statistics
  sim.datalist$ptable$age <- sim.datalist$itable$maxart.starttime[1] - sim.datalist$ptable$TOB

  #get the MaxART population group
  maxart.study.pop <- filter(sim.datalist$ptable,  age >= 18 & TreatTime!=Inf)

  df.rw <- nrow(maxart.study.pop)

  #Dummy when eligible
  sample.eligible <- sample(c(two.weeks, twelve.months, six.months,
                                init.immediate, more.months,
                                init.not.eligible, NA), replace = TRUE, df.rw)

  maxart.study.pop$art.eligible <- maxart.study.pop$TreatTime + sample.eligible

  #interpolation of eligibility date
  #CD4atInfection (t1,c1)  CH4atDeath (t2, c2) ??? t2 not always known after simulation???
  #maxart.study.pop$art.eligible.int <- NA

  #interpolate use naive y = m * x + c where m = (c1 - c2)/(t1 - t2) and c = c1 - m * t1

  all.maxart <- nrow(maxart.study.pop)

  #Get all infected clients
  max.art.ret.study.pop.all <- maxart.study.pop %>%
    filter(TreatTime !=Inf)

  #alldf
  chunk.datalist.test.all <- sim.datalist
  chunk.datalist.test.all$ptable <- max.art.ret.study.pop.all

  ### ART coverage ####################

  hhohho.all.prim.art.coverage.var <- maxart.study.pop %>%
    filter(TreatTime!=Inf) %>%
    summarise(all.all = n()/all.maxart,

              b4.elig.all = sum(TreatTime < art.eligible, na.rm = TRUE )/all.maxart,

              within2wks.all = sum(TreatTime>=art.eligible &
                                     TreatTime < (TreatTime+two.weeks) , na.rm = TRUE)/all.maxart,

              amongelig.all = sum(TreatTime > art.eligible, na.rm = TRUE )/all.maxart,

              within12m.all = sum(TreatTime==art.eligible &
                                    TreatTime < (TreatTime+twelve.months), na.rm = TRUE )/all.maxart,

              within6m.all = sum(TreatTime==art.eligible &
                                   TreatTime < (TreatTime+six.months), na.rm = TRUE )/all.maxart

    )

  all.art.init <- dplyr::select(hhohho.all.prim.art.coverage.var, contains(".all"))

  hhohho.art.initiated.tar.values <- as.numeric(all.art.init*100)
  #
  #when do we want the retention time
  maxart.ret.timepoint <- difftime(maxart.endtime, sim.start.full, units = "days")/365.242
  maxart.ret.timepoint <- round(as.numeric(maxart.ret.timepoint),0)
  maxart.starttime.ret <- round(as.numeric(difftime(maxart.starttime ,sim.start.full, units = "days")/365),0)

  #### Max Retention ####################################
  max.ret.tar <- length(row.names(max.art.retention.all))
  max.ret.tar.all <- rep(NA, max.ret.tar)
  max.ret.tar.list <- c(36, 6, 12)
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

  # # ############### Viral Load none supresssion ############################################
  max.vl.sup.tar.dim <- length(row.names(max.vl.none.suppression.all))
  max.vl.sup.all <- rep(NA, max.vl.sup.tar.dim)
  max.vl.sup.list <- c(6, 12)

  for(vl in 1:max.vl.sup.tar.dim){

    vl.sup.all = 100 - vl.suppressed(datalist = chunk.datalist.test.all,
                               timepoint = maxart.ret.timepoint, vlcutoff = 1000,
                               lessmonths = max.vl.sup.list[vl], site="All")$percentage[[1]]

    max.vl.sup.all[vl] <- ifelse(is.null(vl.sup.all), NA, vl.sup.all)

  }

  hhohho.vl.none.suppression.tar.values <- max.vl.sup.all
  #########

  ###Mortality AIDS related and not
  max.mort.study.pop.prim <- maxart.study.pop %>%
    summarise(
            mort.all = sum(TOD != Inf & TOD > maxart.starttime & TOD < maxart.endtime, na.rm = TRUE ) / n(),
            mort.aids.rel.all = sum(TOD != Inf & AIDSDeath == 1 &
                                      TOD > maxart.starttime & TOD < maxart.endtime, na.rm = TRUE) / n() )

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
                                                             timewindow = c(from.time, to.time)) * 100
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
  time.end.2007 <- as.numeric(difftime(prev.time,sim.start.full, units = "days")/365.245)

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
                   #hhohho.art.initiated.tar.values,
                   hhohho.art.retention.tar.values,
                   hhohho.vl.none.suppression.tar.values,
                   hhohho.mortality.tar.values,
                   hhohho.prev.2007.tar.values,
                   hhohho.AD.tar.values)

  # #Testing
  #sim.summary <- c(max.art.retention.tar.values, Sys.getpid())

  return(sim.summary)

}

#swazi.art.coverage.year.diff <- swazi.art.coverage[,1:2] - swazi.sim.art.coverage.year



