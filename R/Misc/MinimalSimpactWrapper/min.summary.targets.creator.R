#This will be called through source()

#function to clean up the column or row names,
format.names <- function(x,replace = "X",...){
  gsub(replace,'',x)
}

#Recruitment
#maxart.study.pop <- recruit.study.clients(chunk.datalist.test)

#Dummy with no real recruitment
days.in.yr <- 365.245


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

  ########### Prevalence for multiple ages. End time 2007-March-31 ###################################
  prev.age2007.list <- format.names(row.names(prev.2007), replace = "A")

  #create NA results
  sim.prev.2007.len <- length(prev.age2007.list)
  sim.prev.2007.f <- rep(NA, sim.prev.2007.len)
  sim.prev.2007.m <- rep(NA, sim.prev.2007.len)
  sim.prev.2007.fm <- rep(NA, sim.prev.2007.len)

  #prevalence calculated as at 31 March 2007
  time.end.2007 <- as.numeric(difftime(as.Date("2007-03-31") ,sim.start.full, units = "days")/days.in.yr)

  for(i in 1:sim.prev.2007.len){

    split.list.age2007 <- as.numeric(unlist(strsplit(prev.age2007.list[i], "[.]")))
    age.lower.2007 <- split.list.age2007[1]
    age.upper.2007 <- split.list.age2007[2]

    age.sim.prev <- prevalence.calculator(datalist = sim.datalist,
                                          agegroup = c(age.lower.2007, age.upper.2007),
                                          timepoint = time.end.2007)

    #Gender 0 <- male : 1 <- female : NA <- All
    sim.prev.2007.m[i] <- ifelse(nrow(subset(age.sim.prev, Gender == 0)) == 1,
                                 subset(age.sim.prev, Gender == 0)$pointprevalence * 100,
                                 NA)

    sim.prev.2007.f[i] <- ifelse(nrow(subset(age.sim.prev, Gender == 1)) == 1,
                                 subset(age.sim.prev, Gender == 1)$pointprevalence * 100,
                                 NA)

    sim.prev.2007.fm[i] <- ifelse(nrow(subset(age.sim.prev, is.na(Gender))) == 1,
                                  subset(age.sim.prev, is.na(Gender))$pointprevalence * 100,
                                  NA)

  }
  prev.2007.tar.values <- c(sim.prev.2007.f, sim.prev.2007.m, sim.prev.2007.fm)

  #collect all the summary values
  sim.summary <- c(growth.rate.tar.values, prev.2007.tar.values)

  return(sim.summary)

}


