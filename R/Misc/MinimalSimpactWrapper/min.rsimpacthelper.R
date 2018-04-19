#These are RSimpactHelper scripts for functions used in
#Preliminary hhohho
################# simpact.config.inputs

agedistr.creator <- function(shape = 5,
                             scale = 65){
  agebins <- seq(0.5, 99.5, 1)
  probofstillalive <- 1 - pweibull(agebins, shape = shape, scale = scale)
  fractionsinagebins <- 100 * probofstillalive/sum(probofstillalive)
  agedist.data.frame <- data.frame(Age = c(agebins, 100.5),
                                   Percent.Male = c(fractionsinagebins, 0),
                                   Percent.Female = c(fractionsinagebins, 0))
  return(agedist.data.frame)
}

simpact.config.inputs <- function(design.points = 10,
                                  resample.count = 1, ...){

  input.parameters <- list(...)

  if(length(input.parameters)==0){stop("No parameters have been specified.
                                       Use e.g input.varied.params
                                       (conception.alpha_base = c(-3.6, -1.2))")}

  # Creating the LHS over the 0-1 uniform parameter space for the parameters to be estimated
  variables <- length(input.parameters)
  set.seed(1976)

  if(design.points > 4){
    add.design.points <- design.points %/% 4
    aug.points <- design.points %% 4

    rlhs <- lhs::improvedLHS(add.design.points, variables)

    for(i in 1:3){

      rlhs <- lhs::augmentLHS(rlhs, add.design.points)

    }
    if(aug.points > 0){

      rlhs <- lhs::augmentLHS(rlhs, aug.points)
    }

  }else{
    rlhs <- lhs::improvedLHS(design.points, variables)
  }

  if (resample.count > 1){
    for(i in 1:resample.count){

      # Re-sample from the same grid as initially selected
      set.seed(1)
      rlhs <- lhs::augmentLHS(rlhs, design.points)
    }
    rlhs <- tail(rlhs, design.points)
  }

  #Select the config parameters that will be varied from the input config
  lhs.df <- data.frame(matrix(NA, nrow = 1, ncol = variables))
  names(lhs.df) <- names(input.parameters)

  #Repeat have to meet the design.points
  lhs.df <- as.data.frame(lapply(lhs.df, rep, design.points))

  #Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)
  x.index <- 0
  for (j in names(input.parameters)){
    x.index <- x.index + 1
    min.var <- input.parameters[j][[1]][1]
    max.var <- input.parameters[j][[1]][2]
    col.index <- which(names(lhs.df)==j)
    set.seed(1976)
    lhs.df[col.index] <- qunif(rlhs[, x.index], min = as.numeric(min.var), max =as.numeric(max.var))
  }

  #Name the xdesign dataframe
  rlhs <- data.frame(rlhs)
  names(rlhs) <- paste0("xdesign",1:length(rlhs))

  ##This will create the input file for the simpact
  simpactInputParams <- cbind.data.frame(sim.id = 1:design.points, rlhs, lhs.df)

  return(simpactInputParams)

}

intervention.introduced <- function(simulation.type = "maxart",
                                    simulation.start = "1970-03-31"){

  start.time.sim <- as.Date(simulation.start)

  if(simulation.type == "maxart"){
    ### Eligibility increased from <200 to <350 in 2010 and then to <500 in 2015
    # Simulation starts in 1970. After 33 years (in 2003), ART is introduced.

    art.time1 <- as.Date("2003-03-31")
    art.intro <- list()
    art.intro$facilities.outfile.facilityxypos <- "" #reset the writing of the facilities.xy position to none.
    art.intro["time"] <- round(as.numeric(difftime(art.time1,start.time.sim, units = "days")/365.242),0)
    art.intro["diagnosis.baseline"] <- 0 # Reset to 0, from its original @ sim start
    art.intro["monitoring.cd4.threshold.prestudy"] <- 200
    art.intro["monitoring.cd4.threshold.instudy.controlstage"] <- 350
    art.intro["monitoring.cd4.threshold.poststudy"] <- 20000
    art.intro["monitoring.interval.piecewise.cd4s"] <- c("100,250")
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2014:500
    art.time2 <- as.Date("2007-03-31")
    art.intro2 <- list()
    art.intro2["time"] <- round(as.numeric(difftime(art.time2,start.time.sim, units = "days")/365.242),0)
    art.intro2["monitoring.cd4.threshold.prestudy"] <- 200
    art.intro2["monitoring.cd4.threshold.instudy.controlstage"] <- 350
    art.intro2["monitoring.interval.piecewise.cd4s"] <- c("200,350")
    art.intro2["diagnosis.genderfactor"] <- 1.5

    art.time3 <- as.Date("2010-03-31")
    art.intro3 <- list()
    art.intro3["time"] <- round(as.numeric(difftime(art.time3,start.time.sim, units = "days")/365.242),0)
    art.intro3["monitoring.cd4.threshold.prestudy"] <- 350
    art.intro3["monitoring.cd4.threshold.instudy.controlstage"] <- 350
    art.intro3["monitoring.interval.piecewise.cd4s"] <- c("350,500")
    art.intro3["diagnosis.genderfactor"] <- 1

    art.time4 <- as.Date("2014-03-31")
    art.intro4 <- list()
    art.intro4["time"] <- round(as.numeric(difftime(art.time4,start.time.sim, units = "days")/365.242),0)
    art.intro4["monitoring.cd4.threshold.prestudy"] <- 500
    art.intro4["monitoring.cd4.threshold.instudy.controlstage"] <- 500
    art.intro4["monitoring.interval.piecewise.cd4s"] <- c("500,650")
    art.intro4["diagnosis.genderfactor"] <- 0.5

    art.time7 <- as.Date("2016-10-01")
    art.intro7 <- list()
    art.intro7["time"] <- round(as.numeric(difftime(art.time7,start.time.sim, units = "days")/365.242),0)
    art.intro7["monitoring.cd4.threshold.prestudy"] <- 20000
    art.intro7["monitoring.interval.piecewise.cd4s"] <- c("500,2000")
    art.intro7["diagnosis.genderfactor"] <- 0.5

    #collect the art specific changes
    iv <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro7)


  }else if(simulation.type == "simpact-cyan"){

    #Swaziland standard of care NO Maxart

    art.time1 <- as.Date("2003-03-31")
    art.intro <- list()
    #art.intro$facilities.outfile.facilityxypos <- "" #reset the writing of the facilities.xy position to none.
    art.intro["time"] <- round(as.numeric(difftime(art.time1,start.time.sim, units = "days")/365.242),0)
    art.intro["diagnosis.baseline"] <- 0 # Reset to 0, from its original @ sim start
    art.intro["monitoring.cd4.threshold"] <- 200
    art.intro["monitoring.interval.piecewise.cd4s"] <- c("100,250")
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2014:500
    art.time2 <- as.Date("2007-03-31")
    art.intro2 <- list()
    art.intro2["time"] <- round(as.numeric(difftime(art.time2,start.time.sim, units = "days")/365.242),0)
    art.intro2["monitoring.cd4.threshold"] <- 200
    art.intro2["monitoring.interval.piecewise.cd4s"] <- c("200,350")
    art.intro2["diagnosis.genderfactor"] <- 1.5

    art.time3 <- as.Date("2010-03-31")
    art.intro3 <- list()
    art.intro3["time"] <- round(as.numeric(difftime(art.time3,start.time.sim, units = "days")/365.242),0)
    art.intro3["monitoring.cd4.threshold"] <- 350
    art.intro3["monitoring.interval.piecewise.cd4s"] <- c("350,500")
    art.intro3["diagnosis.genderfactor"] <- 1

    art.time4 <- as.Date("2014-03-31")
    art.intro4 <- list()
    art.intro4["time"] <- round(as.numeric(difftime(art.time4,start.time.sim, units = "days")/365.242),0)
    art.intro4["monitoring.cd4.threshold"] <- 500
    art.intro4["monitoring.interval.piecewise.cd4s"] <- c("500,650")
    art.intro4["diagnosis.genderfactor"] <- 0.5

    art.time5 <- as.Date("2016-10-01")
    art.intro5 <- list()
    art.intro5["time"] <- round(as.numeric(difftime(art.time5,start.time.sim, units = "days")/365.242),0)
    art.intro5["monitoring.cd4.threshold"] <- 20000
    art.intro5["monitoring.interval.piecewise.cd4s"] <- c("500,2000")
    art.intro5["diagnosis.genderfactor"] <- 0.5

    iv <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5 )

  }else{
    iv <- list()
  }

  return(iv)
}

readthedata <- function(modeloutput){
  path <- as.character(modeloutput["outputfile"])
  outputID <- as.character(modeloutput["id"])
  DestDir <- sub(pattern = paste0(outputID, "output.txt"), replacement = "", x = path, fixed=T)
  personlogfilename <- paste0(DestDir, outputID, "personlog.csv")
  relationlogfilename <- paste0(DestDir, outputID, "relationlog.csv")
  eventlogfilename <- paste0(DestDir, outputID, "eventlog.csv")
  treatmentlogfilename <- paste0(DestDir, outputID, "treatmentlog.csv")
  inputparamlogfilename <- paste0(DestDir, outputID, "settingslog.csv")
  periodiclogfilename <- paste0(DestDir, outputID, "periodiclog.csv")
  viralloadlogfilename <- paste0(DestDir, outputID, "hivviralloadlog.csv")
  facilitiesxylogfilename <- paste0(DestDir, outputID, "facilitypositions.csv")

  ptable <- data.table::fread(personlogfilename, sep = ",", skip = 0)
  vltable <- data.table::fread(viralloadlogfilename, sep = ",", skip = 0)
  rtable <- data.table::fread(relationlogfilename, sep = ",", skip = 0)
  readetable <- readcsvcolumns::read.csv.columns(eventlogfilename, has.header = FALSE,
                                                 column.types = "rssiirsiir")
  etable <- data.table::setDT(readetable)
  etable.colnames <- colnames(etable)
  if (ncol(etable) == 10){
    data.table::setnames(etable, etable.colnames,
                         c("eventtime", "eventname", "p1name", "p1ID", "p1gender",
                           "p1age", "p2name", "p2ID", "p2gender", "p2age"))
  } else {
    data.table::setnames(etable, etable.colnames,
                         c("eventtime", "eventname", "p1name", "p1ID", "p1gender",
                           "p1age", "p2name", "p2ID", "p2gender", "p2age", "extradescript", "value"))

  }
  ttable <- data.table::fread(treatmentlogfilename, sep  = ",", skip = 0)
  itable <- data.table::fread(inputparamlogfilename, sep  = ",", skip = 0)

  if (file.exists(periodiclogfilename) && file.exists(facilitiesxylogfilename) ){
    ltable <- data.table::fread(periodiclogfilename, sep = ",", skip = 0)

    ftable <- data.table::fread(facilitiesxylogfilename, sep = ",", skip = 0)

    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable,
                         ttable = ttable, itable = itable, ltable = ltable,
                         vltable = vltable, ftable = ftable)
  }else if(file.exists(facilitiesxylogfilename) && !file.exists(periodiclogfilename)){
    ftable <- data.table::fread(facilitiesxylogfilename, sep = ",", skip = 0)

    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable,
                         ttable = ttable, itable = itable, vltable = vltable, ftable = ftable)

  }else if(!file.exists(facilitiesxylogfilename) && file.exists(periodiclogfilename) ){
    ltable <- data.table::fread(periodiclogfilename, sep = ",", skip = 0)

    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable,
                         ttable = ttable, itable = itable, ltable = ltable, vltable = vltable)
  }else{
    outputtables <- list(ptable = ptable, rtable = rtable,
                         etable = etable, ttable = ttable, itable = itable, vltable = vltable)
  }

  return(outputtables)
}

pop.growth.calculator <- function(datalist = datalist,
                                  timewindow = c(1, 20)){

  end.popsize <- datalist$ltable %>%
    subset(Time==timewindow[2]) %>%
    dplyr::select(PopSize) %>%
    as.numeric()

  start.popsize <- ifelse(timewindow[1]==0,
                          yes = datalist$itable$population.nummen +
                            datalist$itable$population.numwomen,
                          no = datalist$ltable %>% subset(Time==timewindow[1]) %>%
                            dplyr::select(PopSize) %>% as.numeric())

  growth.rate <- log(end.popsize/start.popsize) / diff(timewindow)

  return(growth.rate)
}

prevalence.calculator <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timepoint = 30){

  # First we only take the data of people who were alive at the timepoint
  DTalive.infected <- alive.infected(datalist = datalist, timepoint = timepoint)

  DTalive.infected <- DTalive.infected %>%
    mutate(age = timepoint - TOB)

  raw.df <- DTalive.infected %>%
    subset(age >= agegroup[1]  & age < agegroup[2]) %>%
    as.data.frame()

  if(nrow(raw.df) > 0){
    # Now we apply dplyr function to get sum of cases and sum of exposure.time per gender.
    prevalence.df <-  raw.df %>%
      group_by(Gender) %>%
      summarise(popsize = n(),
                sum.cases = sum(Infected),
                pointprevalence = sum.cases / popsize,
                pointprevalence.95.ll = as.numeric(
                  binom.test(x=sum.cases, n=popsize)$conf.int)[1],
                pointprevalence.95.ul = as.numeric(
                  binom.test(x=sum.cases,n=popsize)$conf.int)[2]
    )
    prevalence.all.df <- raw.df %>%
      summarise(popsize = n(),
                sum.cases = sum(Infected),
                pointprevalence = sum.cases / popsize,
                pointprevalence.95.ll = as.numeric(
                  binom.test(x=sum.cases, n=popsize)$conf.int)[1],
                pointprevalence.95.ul = as.numeric(
                  binom.test(x=sum.cases,n=popsize)$conf.int)[2]
    )

    prevalence.all.df <- cbind(Gender = NA, prevalence.all.df)

    prevalence.df <- rbind(prevalence.df, prevalence.all.df)
  } else {
    prevalence.df <- as.data.frame(matrix(NA, 3, 6))

  }

  names(prevalence.df) <- c("Gender", "popsize", "sum.cases",
                            "pointprevalence", "pointprevalence.95.ll",
                            "pointprevalence.95.ul")

  return(prevalence.df)
}

alive.infected <- function(datalist = datalist,
                           timepoint = 40,
                           site = "All") {
  # arguments are the personlog data.table and a point in time
  DT <- datalist$ptable

  if(nrow(DT) > 0){
    if (site == "All") {
      DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint)
    } else{

      DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint & pfacility == site)
    }

    # Now we allocate infection status to all alive people in our table
    DTalive <-  DTalive %>%
      dplyr::mutate(Infected = (timepoint >= InfectTime)) %>%
      as.data.frame()
  }else{
    DTalive <- datalist$ptable %>%
      as.data.frame()
  }
  return(DTalive)
}

input.params.creator <- function(mortality.normal.weibull.shape = 5,
                                 mortality.normal.weibull.scale = 65,
                                 mortality.normal.weibull.genderdiff = 5,
                                 periodiclogging.interval = 1,
                                 syncrefyear.interval = 1,
                                 formation.hazard.type = "agegapry",
                                 person.art.accept.threshold.dist.type = "fixed",
                                 person.art.accept.threshold.dist.fixed.value = 0.5,
                                 person.eagerness.man.type = "independent",
                                 person.eagerness.woman.type = "independent",
                                 person.eagerness.man.dist.type = "gamma",
                                 person.eagerness.woman.dist.type = "gamma",
                                 person.eagerness.man.dist.gamma.a = 0.231989836885,#0.15#0.85#0.1
                                 person.eagerness.man.dist.gamma.b = 45,#70#100#3.5#5#10#20 #170
                                 person.eagerness.woman.dist.gamma.a = 1.49800270590903,#0.15#0.1
                                 person.eagerness.woman.dist.gamma.b = 45,#70#100#3.5#5#10#20#170
                                 person.agegap.man.dist.type = "normal",
                                 person.agegap.woman.dist.type = "normal",
                                 person.agegap.man.dist.normal.mu = 4.84302562640607, #-5
                                 person.agegap.woman.dist.normal.mu = 0, #2.5
                                 person.agegap.man.dist.normal.sigma = 3.33997313140938,
                                 person.agegap.woman.dist.normal.sigma = 1,
                                 formation.hazard.agegapry.numrel_man = -0.5,
                                 formation.hazard.agegapry.numrel_woman = -1.69901499932697,
                                 formation.hazard.agegapry.gap_factor_man_exp = -0.35,#-0.15#-0.5
                                 formation.hazard.agegapry.gap_factor_woman_exp = -0.35,#-0.15
                                 formation.hazard.agegapry.gap_factor_man_age = 0.05,
                                 formation.hazard.agegapry.gap_factor_woman_age = 0.05,
                                 formation.hazard.agegapry.meanage = -0.1,
                                 formation.hazard.agegapry.numrel_diff = -0.1,
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 formation.hazard.agegapry.gap_agescale_man = 0.23,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.1,#0.23
                                 formation.hazard.agegapry.eagerness_sum = 0.1,
                                 person.vsp.tofinalaids.x = 12,
                                 person.vsp.toaids.x = 7,
                                 formation.hazard.agegapry.eagerness_diff = -0.019990725233685,#-0.110975
                                 dissolution.alpha_0 = -0.52,#-0.1 # baseline
                                 dissolution.alpha_4 = -0.05,
                                 debut.debutage = 15,
                                 population.simtime = 40,
                                 population.nummen = 500, #1000 #2000
                                 population.numwomen = 500, #1000 #2000
                                 population.eyecap.fraction = 0.2,
                                 hivseed.type = "amount",
                                 hivseed.amount = 20,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 30,
                                 hivseed.time = 10,
                                 hivtransmission.param.a = -1.0352239,
                                 hivtransmission.param.b = -89.339994,
                                 hivtransmission.param.c = 0.4948478,
                                 hivtransmission.param.f1 = 1.9981126929339, #~1.6 =>hazard is x5 in 15yo
                                 hivtransmission.param.f2 = log(log(2.5) / log(5)) / 5,
                                 mortality.aids.survtime.C = 62,
                                 mortality.aids.survtime.k = -0.2,
                                 conception.alpha_base = -2.44785979535887, #-3
                                 diagnosis.baseline = -100,
                                 monitoring.cd4.threshold = 0.1,#Treatment will not start before schedule
                                 monitoring.fraction.log_viralload = 0.3,
                                 population.msm = "no",
                                 person.eagerness.man.msm.dist.type = "fixed",
                                 person.eagerness.man.msm.dist.fixed.value = 0,
                                 formationmsm.hazard.type = "simple",
                                 formationmsm.hazard.simple.alpha_0 = 2,
                                 formationmsm.hazard.simple.alpha_12 = -0.4,
                                 formationmsm.hazard.simple.alpha_5 = -0.2, # The factor α5 controls the relative importance of the age gap between the partners.
                                 formationmsm.hazard.simple.alpha_6 = 0,  # weight for sum of eagerness parameters
                                 formationmsm.hazard.simple.alpha_7 = 0,
                                 birth.pregnancyduration.dist.type = "fixed",
                                 birth.pregnancyduration.dist.fixed.value = 268/365, # just over 38 weeks
                                 birth.boygirlratio = 1.0/2.01, #0.5024876, #101:100
                                 simulation.type = "simpact-cyan",
                                 monitoring.cd4.threshold.prestudy = 1,   #it has to be set to zero initially (does not exist)
                                 monitoring.cd4.threshold.instudy.controlstage = 1,  #set to default to barely minimal (does not exist)
                                 monitoring.cd4.threshold.poststudy = 20000 # setting the initial
                                 #monitoring.cd4.threshold.instudy.transitionstage = Inf,
                                 #monitoring.cd4.threshold.instudy.interventionstage = Inf
){
  input.params <- list()
  input.params$mortality.normal.weibull.shape <- mortality.normal.weibull.shape
  input.params$mortality.normal.weibull.scale <- mortality.normal.weibull.scale
  input.params$mortality.normal.weibull.genderdiff <- mortality.normal.weibull.genderdiff
  input.params$periodiclogging.interval <- periodiclogging.interval
  input.params$syncrefyear.interval <- syncrefyear.interval
  input.params$formation.hazard.type <- formation.hazard.type
  input.params$person.art.accept.threshold.dist.type <- person.art.accept.threshold.dist.type
  input.params$person.art.accept.threshold.dist.fixed.value <- person.art.accept.threshold.dist.fixed.value
  input.params$person.eagerness.man.type <- person.eagerness.man.type
  input.params$person.eagerness.woman.type <- person.eagerness.woman.type
  input.params$person.eagerness.man.dist.type <- person.eagerness.man.dist.type
  input.params$person.eagerness.woman.dist.type <- person.eagerness.woman.dist.type
  input.params$person.eagerness.man.dist.gamma.a <- person.eagerness.man.dist.gamma.a
  input.params$person.eagerness.man.dist.gamma.b <- person.eagerness.man.dist.gamma.b
  input.params$person.eagerness.woman.dist.gamma.a <- person.eagerness.woman.dist.gamma.a
  input.params$person.eagerness.woman.dist.gamma.b <- person.eagerness.woman.dist.gamma.b
  input.params$person.agegap.man.dist.type <- person.agegap.man.dist.type
  input.params$person.agegap.woman.dist.type <- person.agegap.woman.dist.type
  input.params$person.agegap.man.dist.normal.mu <- person.agegap.man.dist.normal.mu
  input.params$person.agegap.woman.dist.normal.mu <- person.agegap.woman.dist.normal.mu
  input.params$person.agegap.man.dist.normal.sigma <- person.agegap.man.dist.normal.sigma
  input.params$person.agegap.woman.dist.normal.sigma <- person.agegap.woman.dist.normal.sigma
  input.params$formation.hazard.agegapry.numrel_man <- formation.hazard.agegapry.numrel_man
  input.params$formation.hazard.agegapry.numrel_diff <- formation.hazard.agegapry.numrel_diff
  input.params$formation.hazard.agegapry.numrel_woman <- formation.hazard.agegapry.numrel_woman
  input.params$formation.hazard.agegapry.gap_factor_man_age <- formation.hazard.agegapry.gap_factor_man_age
  input.params$formation.hazard.agegapry.gap_factor_woman_age <- formation.hazard.agegapry.gap_factor_woman_age
  input.params$formation.hazard.agegapry.meanage <- formation.hazard.agegapry.meanage
  input.params$formation.hazard.agegapry.gap_factor_man_exp <- formation.hazard.agegapry.gap_factor_man_exp
  input.params$formation.hazard.agegapry.gap_factor_woman_exp <- formation.hazard.agegapry.gap_factor_woman_exp
  input.params$formation.hazard.agegapry.gap_factor_man_const <- formation.hazard.agegapry.gap_factor_man_const
  input.params$formation.hazard.agegapry.gap_factor_woman_const <- formation.hazard.agegapry.gap_factor_woman_const
  input.params$formation.hazard.agegapry.gap_agescale_man <- formation.hazard.agegapry.gap_agescale_man
  input.params$formation.hazard.agegapry.gap_agescale_woman <- formation.hazard.agegapry.gap_agescale_woman
  input.params$formation.hazard.agegapry.eagerness_sum <- formation.hazard.agegapry.eagerness_sum
  input.params$formation.hazard.agegapry.eagerness_diff <- formation.hazard.agegapry.eagerness_diff
  input.params$dissolution.alpha_0 <- dissolution.alpha_0
  input.params$dissolution.alpha_4 <- dissolution.alpha_4
  input.params$debut.debutage <- debut.debutage
  input.params$population.simtime <- population.simtime
  input.params$population.nummen <- population.nummen
  input.params$population.numwomen <- population.numwomen
  input.params$population.maxevents <-  population.simtime * population.nummen * 5 #4events/person
  input.params$population.eyecap.fraction <- population.eyecap.fraction
  input.params$hivseed.type <- hivseed.type
  input.params$hivseed.amount <- hivseed.amount
  input.params$hivseed.age.min <- hivseed.age.min
  input.params$hivseed.age.max <- hivseed.age.max
  input.params$hivseed.time <- hivseed.time
  input.params$person.vsp.tofinalaids.x <- person.vsp.tofinalaids.x
  input.params$person.vsp.toaids.x <- person.vsp.toaids.x
  input.params$hivtransmission.param.a <- hivtransmission.param.a
  input.params$hivtransmission.param.b <- hivtransmission.param.b
  input.params$hivtransmission.param.c <- hivtransmission.param.c
  input.params$hivtransmission.param.f1 <- hivtransmission.param.f1
  input.params$hivtransmission.param.f2 <- hivtransmission.param.f2
  input.params$mortality.aids.survtime.k <- mortality.aids.survtime.k
  input.params$mortality.aids.survtime.C <- mortality.aids.survtime.C
  input.params$conception.alpha_base <- conception.alpha_base
  input.params$monitoring.fraction.log_viralload <- monitoring.fraction.log_viralload
  input.params$diagnosis.baseline <- diagnosis.baseline
  input.params$birth.boygirlratio <- birth.boygirlratio
  input.params$population.msm <- population.msm
  input.params$person.eagerness.man.msm.dist.type <- person.eagerness.man.msm.dist.type
  input.params$person.eagerness.man.msm.dist.fixed.value <- person.eagerness.man.msm.dist.fixed.value
  input.params$formationmsm.hazard.type <- formationmsm.hazard.type
  input.params$formationmsm.hazard.simple.alpha_0 <- formationmsm.hazard.simple.alpha_0
  input.params$formationmsm.hazard.simple.alpha_12 <- formationmsm.hazard.simple.alpha_12
  input.params$formationmsm.hazard.simple.alpha_5 <- formationmsm.hazard.simple.alpha_5 # The factor α5 controls the relative importance of the age gap between the partners.
  input.params$formationmsm.hazard.simple.alpha_6 <- formationmsm.hazard.simple.alpha_6  # weight for sum of eagerness parameters
  input.params$formationmsm.hazard.simple.alpha_7 <- formationmsm.hazard.simple.alpha_7
  input.params$birth.pregnancyduration.dist.type <- birth.pregnancyduration.dist.type
  input.params$birth.pregnancyduration.dist.fixed.value <- birth.pregnancyduration.dist.fixed.value

  if(simulation.type == "simpact-cyan"){
    input.params$monitoring.cd4.threshold <- monitoring.cd4.threshold
  }else{
    input.params$facilities.outfile.facilityxypos <- "${SIMPACT_OUTPUT_PREFIX}facilitypositions.csv"
    input.params$monitoring.cd4.threshold.prestudy <- monitoring.cd4.threshold.prestudy
    input.params$monitoring.cd4.threshold.instudy.controlstage <- monitoring.cd4.threshold.instudy.controlstage
    input.params$monitoring.cd4.threshold.poststudy <- monitoring.cd4.threshold.poststudy #post study everyone was now on test and study
  }

  return(input.params)
}

