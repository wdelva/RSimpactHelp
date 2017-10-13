#This will be called through source()

#function to clean up the column or row names,
format.names <- function(x,replace = "X",...){
  gsub(replace,'',x)
}

sim.summary.creator.all <- function(sim.datalist = chunk.datalist.test){

  ###### Growth rate ##############################################################
  gr.year.list <- as.numeric(format.names(names(growth.rate), replace = "X"))

  growth.rate.tar.values <- rep(NA, length(gr.year.list))

  #Growth rate is calculated per year.
  for(i in 1:length(gr.year.list)){
    from.time <- gr.year.list[i] - 1 - sim.start
    to.time <-  from.time + 1

    growth.rate.tar.values[i] <- pop.growth.calculator(datalist = sim.datalist,
                                                     timewindow = c(from.time, to.time)) * 100
  }

  ###Incidence for multiple ages, start and end time set ###############################################
  inci.age2011.list <- format.names(row.names(inci.2011), replace = "A")
  sim.inci.2011.f <- rep(NA, length(inci.age2011.list))
  sim.inci.2011.m <- rep(NA, length(inci.age2011.list))

  time.start.2011 <- as.numeric(difftime(as.Date("2010-12-01") ,sim.start.full, units = "days")/365.245)
  time.end.2011 <- as.numeric(difftime(as.Date("2011-06-30") ,sim.start.full, units = "days")/365.245)

  for(i in 1:length(inci.age2011.list)){

    split.list.age2011 <- as.numeric(unlist(strsplit(inci.age2011.list[i], "[.]")))
    age.lower.2011 <- split.list.age2011[1]
    age.upper.2011 <- split.list.age2011[2]

    age.sim.inci <- incidence.calculator(datalist = sim.datalist,
                                         agegroup = c(age.lower.2011, age.upper.2011),
                                         timewindow = c(time.start.2011, time.end.2011),
                                         only.active = "No")

    #Gender 0 <- male : 1 <- female
    sim.inci.2011.m[i] <- age.sim.inci$incidence[1] * 100
    sim.inci.2011.f[i] <- age.sim.inci$incidence[2] * 100

  }

  inci.2011.tar.values <- c(sim.inci.2011.f,sim.inci.2011.m)

  ########### Prevalence for multiple ages. End time 2007-March-31 ###################################
  prev.age2007.list <- format.names(row.names(prev.2007), replace = "A")

  sim.prev.2007.len <- length(prev.age2007.list)
  sim.prev.2007.f <- rep(NA, sim.prev.2007.len)
  sim.prev.2007.m <- rep(NA, sim.prev.2007.len)
  sim.prev.2007.fm <- rep(NA, sim.prev.2007.len)

  time.end.2007 <- as.numeric(difftime(as.Date("2007-03-31") ,sim.start.full, units = "days")/365.245)

  for(i in 1:sim.prev.2007.len){

    split.list.age2007 <- as.numeric(unlist(strsplit(prev.age2007.list[i], "[.]")))
    age.lower.2007 <- split.list.age2007[1]
    age.upper.2007 <- split.list.age2007[2]

    age.sim.prev <- prevalence.calculator(datalist = sim.datalist,
                                                agegroup = c(age.lower.2007, age.upper.2007),
                                                timepoint = time.end.2007)

    #Gender 0 <- male : 1 <- female
    sim.prev.2007.m[i] <- age.sim.prev$pointprevalence[1] * 100
    sim.prev.2007.f[i] <- age.sim.prev$pointprevalence[2] * 100
    sim.prev.2007.fm[i] <- age.sim.prev$pointprevalence[3] * 100

  }
  prev.2007.tar.values <- c(sim.prev.2007.f, sim.prev.2007.m, sim.prev.2007.fm)

  #collect all the summary values
  sim.summary <- c(growth.rate.tar.values, inci.2011.tar.values, prev.2007.tar.values)

  return(sim.summary)

}


