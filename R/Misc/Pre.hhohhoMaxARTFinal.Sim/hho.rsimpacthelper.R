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
  set.seed(1)

  rlhs <- randomLHS(design.points, variables)
  if (resample.count > 1){
    for(i in 1:resample.count){

      # Re-sample from the same grid as initially selected
      rlhs <- augmentLHS(rlhs, design.points)
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
    set.seed(1)
    lhs.df[col.index] <- qunif(rlhs[, x.index], min = as.numeric(min.var), max =as.numeric(max.var))
  }

  #Name the xdesign dataframe
  rlhs <- data.frame(rlhs)
  names(rlhs) <- paste0("xdesign",1:length(rlhs))

  ##This will create the input file for the simpact
  simpactInputParams <- cbind.data.frame(sim.id = 1:design.points, rlhs, lhs.df)


  return(simpactInputParams)

}

intervention.introduced <- function(simulation.type = "simpact-cyan",
                                    simulation.start = "1970-03-31"){

  start.time.sim <- as.Date(simulation.start)

  if(simulation.type != "simpact-cyan"){
    ### Eligibility increased from <200 to <350 in 2010 and then to <500 in 2015
    # Simulation starts in 1970. After 33 years (in 2003), ART is introduced.

    art.time1 <- as.Date("2003-03-31")
    art.intro <- list()
    art.intro$facilities.outfile.facilityxypos <- "" #reset the writing of the facilities.xy position to none.
    art.intro["time"] <- round(as.numeric(difftime(art.time1,start.time.sim, units = "days")/365.242),0)
    art.intro["diagnosis.baseline"] <- 0 # Reset to 0, from its original @ sim start
    art.intro["monitoring.cd4.threshold.prestudy"] <- 200
    art.intro["monitoring.cd4.threshold.poststudy"] <- 20000
    art.intro["monitoring.interval.piecewise.cd4s"] <- "100,250"
    art.intro["diagnosis.genderfactor"] <- 2

    # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2014:500
    art.time2 <- as.Date("2007-03-31")
    art.intro2 <- list()
    art.intro2["time"] <- round(as.numeric(difftime(art.time2,start.time.sim, units = "days")/365.242),0)
    art.intro2["monitoring.cd4.threshold.prestudy"] <- 200
    art.intro2["monitoring.interval.piecewise.cd4s"] <- "200,350"
    art.intro2["diagnosis.genderfactor"] <- 1.5

    art.time3 <- as.Date("2010-03-31")
    art.intro3 <- list()
    art.intro3["time"] <- round(as.numeric(difftime(art.time3,start.time.sim, units = "days")/365.242),0)
    art.intro3["monitoring.cd4.threshold.prestudy"] <- 350
    art.intro3["monitoring.interval.piecewise.cd4s"] <- "350,500"
    art.intro3["diagnosis.genderfactor"] <- 1

    art.time4 <- as.Date("2014-03-31")
    art.intro4 <- list()
    art.intro4["time"] <- round(as.numeric(difftime(art.time4,start.time.sim, units = "days")/365.242),0)
    art.intro4["monitoring.cd4.threshold.prestudy"] <- 500
    art.intro4["monitoring.interval.piecewise.cd4s"] <- "500,650"
    art.intro4["diagnosis.genderfactor"] <- 0.5

    art.time7 <- as.Date("2016-10-01")
    art.intro7 <- list()
    art.intro7["time"] <- round(as.numeric(difftime(art.time7,start.time.sim, units = "days")/365.242),0)
    art.intro7["monitoring.cd4.threshold.prestudy"] <- 20000
    art.intro7["monitoring.interval.piecewise.cd4s"] <- "500,2000"
    art.intro7["diagnosis.genderfactor"] <- 0.5

    #reset stepinterval to last facility and end of study
    art.time5 <- as.Date("2016-09-02")
    stepinterval <- list()
    stepinterval["time"] <- round(as.numeric(difftime(art.time5,start.time.sim, units = "days")/365.242),0)
    stepinterval["maxart.stepinterval"] <- 1/12

    #reset stepinterval to the end of study
    art.time6 <- as.Date("2016-10-02")
    stepinterval1 <- list()
    stepinterval1["time"] <- round(as.numeric(difftime(art.time6,start.time.sim, units = "days")/365.242),0)
    stepinterval1["maxart.stepinterval"] <- 7/12

    iv <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro7, stepinterval, stepinterval1 )
  }else{
    iv <-list()
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


ART.retention <- function(datalist = datalist,
                          agegroup = c(15,150),
                          ARTtimewindow = c(20,34),
                          retentiontimeMonths = 12, #6 months default
                          site="All"){

  #get the treatment table
  DT.treatment <- datalist$ttable

  #Get people who have started treatment within timepoint
  if (site == "All") {
    DT.treatment <- subset(DT.treatment,
                           (TStart >= ARTtimewindow[1] &
                              TStart < ARTtimewindow[2])
    )
  } else{
    DT.treatment <- subset(DT.treatment,
                           (TStart >= ARTtimewindow[1] &
                              TStart < ARTtimewindow[2]) &
                             pfacility == site )
  }

  #Get the right age within the time window
  age.group.ART <- age.group.time.window(datalist = datalist,
                                         agegroup = agegroup,
                                         timewindow = ARTtimewindow, site = site)


  if(nrow(DT.treatment)>0){
    #get everyone last known treatment status
    DT.treatment <- DT.treatment %>%
      group_by(ID) %>%
      dplyr::filter(row_number()==n()) %>%
      as.data.frame

    DT.treatment$retention.time <- DT.treatment$TStart + retentiontimeMonths/12

    #restrict the list to the right agegroup
    DT.treatment <- subset(DT.treatment, ID %in% age.group.ART$ID)



    DT.treatment$retention <- DT.treatment$TEnd > DT.treatment$retention.time

    DT.treatment.retention <- DT.treatment %>%
      group_by(Gender) %>%
      summarise(TotalCases = n(),
                ART.retention = sum(retention),
                percentage = ART.retention / TotalCases * 100 ) %>%
      as.data.frame

  }else{ DT.treatment.retention <- as.data.frame(matrix(NA, 0, 4))}


  if(nrow(DT.treatment.retention)==0){
    DT.treatment.retention <- as.data.frame(matrix(NA, 3, 4))
  }else{
    DT.treatment.retention <- rbind(DT.treatment.retention,
                                    c(NA, nrow(DT.treatment), sum(DT.treatment$retention),
                                      sum(DT.treatment$retention)/nrow(DT.treatment)*100) )
  }

  names(DT.treatment.retention) <- c("Gender","TotalCases","ART.retention", "percentage")

  return(DT.treatment.retention)
}

age.group.time.window <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timewindow = c(15, 30),
                                  site="All"){

  DT <- datalist$ptable

  if(site=="All"){
    DTexists.timewindow <- DT
  }else{
    DTexists.timewindow <- subset(DT, pfacility==site)
  }

  #Convert lower age of interest into time
  DTexists.timewindow <- DTexists.timewindow %>%
    dplyr::mutate(LowerTimeAgeGroup = TOB + agegroup[1],
                  LowerTimeWindow = timewindow[1],
                  exposure.start = pmax(LowerTimeAgeGroup, LowerTimeWindow),
                  UpperTimeAgeGroup = TOB + agegroup[2], #Convert upper age of interest into time
                  UpperTimeWindow = timewindow[2],
                  exposure.end = pmin(TOD, pmin(UpperTimeAgeGroup, UpperTimeWindow)),
                  exposure.time = exposure.end - exposure.start,
                  #Exposure time (Everyone with exposure time greater than 0)
                  real.exposure.time = exposure.time > 0
    ) %>%
    subset(real.exposure.time == TRUE) %>%
    #remove those that did not have exposure (died before timewindow or outside min age group)
    as.data.frame

  return(DTexists.timewindow)
}

vl.suppressed <- function(datalist = datalist,
                          timepoint = 35, vlcutoff = 1000,
                          lessmonths = 6, site="All"){

  #Get the individual who started treatment before the time timepoint
  if (site == "All") {
    DT.infected <- subset(datalist$ptable,
                          TreatTime < timepoint)
  } else{

    DT.infected <- subset(datalist$ptable,
                          TreatTime < timepoint,
                          pfacility == site)
  }

  #log10 vlcutoff
  vl.cutoff <- log10(vlcutoff)

  #cut off in months
  .months <- timepoint - lessmonths/12


  #look at the treatment table
  #get the treatment started before the timepoint
  #was on treatment within the time window to check VL

  last.treatment <- datalist$ttable %>%
    subset(ID %in% DT.infected$ID) %>%
    group_by(ID) %>%
    dplyr::filter(row_number()==n()) %>%
    as.data.frame()

  #Get the VL on value due to starting treatment
  vl.table <-  datalist$vltable %>%
    group_by(ID) %>%
    dplyr::filter(Desc =="Started ART",
                  row_number()==n()) %>%
    as.data.frame()

  last.treatment.vl <- dplyr::left_join(x = last.treatment,
                                        y = vl.table,
                                        by = c("ID") )

  #Assuming that we know when Tend.
  #As long as Treatment was started before timepoint and AFTER "lessmonths"
  #clients still on still treatment, we keep the record

  last.treatment.vl <- last.treatment.vl %>%
    mutate(on.treatment.lessmonths = (.months < TEnd),
           vl.suppressed = (Log10VL < vl.cutoff),
           vl.treat.suppressed = on.treatment.lessmonths * vl.suppressed
    )


  if(nrow(last.treatment.vl)==0){

    vlSuppressed.TP <- as.data.frame(matrix(NA, 3, 4))

  }else{
    vlSuppressed.TP <- last.treatment.vl %>%
      group_by(Gender) %>%
      summarise(TotalCases = n(),
                VLSuppressed = sum(vl.treat.suppressed, na.rm = TRUE),
                Percentage = sum(vl.treat.suppressed, na.rm = TRUE)/n() *100) %>%
      as.data.frame()


    vlSuppressed.TP <- rbind(vlSuppressed.TP,
                             c(NA, nrow(last.treatment.vl),
                               sum(last.treatment.vl$vl.treat.suppressed, na.rm = TRUE),
                               sum(last.treatment.vl$vl.treat.suppressed, na.rm = TRUE)/nrow(last.treatment.vl)*100
                             ))

  }

  names(vlSuppressed.TP) <- c("Gender", "TotalCases", "VLSuppressed", "Percentage")


  return(vlSuppressed.TP)
}

pop.growth.calculator <- function(datalist = datalist,
                                  timewindow = c(0, 20)){

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

incidence.calculator <- function(datalist = datalist,
                                 agegroup = c(15, 30),
                                 timewindow = c(20, 30),
                                 only.active = "No"){
  time.of.lowerbound.agegroup <- datalist$ptable$TOB + agegroup[1]
  time.of.lowerbound.timewind <- timewindow[1]
  exposure.start <- pmax(time.of.lowerbound.agegroup, time.of.lowerbound.timewind)
  time.of.upperbound.agegroup <- datalist$ptable$TOB + agegroup[2]
  time.of.upperbound.timewind <- timewindow[2]
  time.of.HIV.infection <- datalist$ptable$InfectTime

  exposure.end <- pmin(time.of.HIV.infection,
                       pmin(time.of.upperbound.agegroup, time.of.upperbound.timewind))
  exposure.time <- exposure.end - exposure.start #This is naive exposure time, before tidying up
  real.exposure.time <- exposure.time > 0 #We create a vector to see who REALLY had exposure time
  exposure.time[real.exposure.time == FALSE] <- 0

  # Now we check who of the people with the real exposure time had the event
  # Their InfectTime must be after their exposure.time started, and before or at exposure.end
  infection.after.exposure.start <- datalist$ptable$InfectTime > exposure.start
  infection.before.or.at.exposure.end <- datalist$ptable$InfectTime <= exposure.end
  infection.in.timewindow <- infection.after.exposure.start & infection.before.or.at.exposure.end


  datalist$ptable$incident.case <- infection.in.timewindow
  datalist$ptable$exposure.times <- exposure.time

  raw.df <- data.frame(datalist$ptable)
  raw.plus.df <- raw.df

  if (only.active != "No"){
    norels.timespent.df <- timespentsingle.calculator(datalist = datalist,
                                                      agegroup = agegroup,
                                                      timewindow = timewindow,
                                                      type = only.active)
    raw.plus.df <- dplyr::left_join(x = raw.df,
                                    y = norels.timespent.df,
                                    by = c("ID" = "woman.ID"))
    raw.plus.df$sum.norels.timespent[is.na(raw.plus.df$sum.norels.timespent)] <- 0
    raw.plus.df$exposure.times <- raw.plus.df$exposure.times - raw.plus.df$sum.norels.timespent
  }

  raw.plus.filtered.df <- dplyr::filter(raw.plus.df, exposure.times >= 0)

  if(nrow(raw.plus.filtered.df) > 0){

    # Now we apply some dplyr function to get the sum of cases and sum of exposure.time per gender.
    incidence.df <- dplyr::summarise(dplyr::group_by(raw.plus.filtered.df, Gender),
                                     sum.exposure.time = sum(exposure.times),
                                     sum.incident.cases = sum(incident.case),
                                     incidence = sum(incident.case) / sum(exposure.times),
                                     incidence.95.ll = as.numeric(
                                       exactci::poisson.exact(x=sum(incident.case),
                                                              T=sum(exposure.times))$conf.int)[1],
                                     incidence.95.ul = as.numeric(
                                       exactci::poisson.exact(x=sum(incident.case),
                                                              T=sum(exposure.times))$conf.int)[2]
    )

    # Now we add the overall incidence to this dataframe
    incidence.all.df <- dplyr::summarise(raw.plus.filtered.df,
                                         sum.exposure.time = sum(exposure.times),
                                         sum.incident.cases = sum(incident.case),
                                         incidence = sum(incident.case) / sum(exposure.times),
                                         incidence.95.ll = as.numeric(
                                           exactci::poisson.exact(x = sum(incident.case),
                                                                  T=sum(exposure.times))$conf.int)[1],
                                         incidence.95.ul = as.numeric(
                                           exactci::poisson.exact(x = sum(incident.case),
                                                                  T=sum(exposure.times))$conf.int)[2]
    )


    incidence.all.df <- cbind(Gender = NA,incidence.all.df)

    incidence.df <- rbind(incidence.df, incidence.all.df)
  }else{

    incidence.df <- as.data.frame(matrix(NA, 3, 6))
    names(incidence.df) <- c("Gender", "sum.exposure.time", "sum.incident.cases",
                              "incidence", "incidence.95.ll", "incidence.95.ul")

  }
  return(incidence.df)
}

timespentsingle.calculator <- function(datalist = datalist,
                                       agegroup = c(15, 30),
                                       timewindow = c(20, 30),
                                       type = "Strict"){
  time.of.lowerbound.agegroup <- datalist$ptable$TOB + agegroup[1]
  time.of.lowerbound.timewind <- timewindow[1]
  exposure.start <- pmax(time.of.lowerbound.agegroup, time.of.lowerbound.timewind)
  time.of.upperbound.agegroup <- datalist$ptable$TOB + agegroup[2]
  time.of.upperbound.timewind <- timewindow[2]
  time.of.HIV.infection <- datalist$ptable$InfectTime

  exposure.end <- pmin(time.of.HIV.infection,
                       pmin(time.of.upperbound.agegroup, time.of.upperbound.timewind))

  #Thus the naive exposure time, before tidying up
  exposure.time <- exposure.end - exposure.start

  # We create a vector to see who REALLY had exposure time
  real.exposure.time <- exposure.time > 0
  exposure.time[real.exposure.time == FALSE] <- 0

  # Now we check who of the people with the real exposure time had the event
  # Their InfectTime must be after their exposure.time started, and before or at exposure.end
  infection.after.exposure.start <- datalist$ptable$InfectTime > exposure.start
  infection.before.or.at.exposure.end <- datalist$ptable$InfectTime <= exposure.end
  infection.in.timewindow <- infection.after.exposure.start & infection.before.or.at.exposure.end

  datalist$ptable$incident.case <- infection.in.timewindow
  datalist$ptable$exposure.times <- exposure.time
  datalist$ptable$exposure.start <- exposure.start
  datalist$ptable$exposure.end <- exposure.end

  raw.df <- data.frame(datalist$ptable)
  women.filtered.df <- dplyr::filter(raw.df,
                                     exposure.times > 0,
                                     Gender == 1)

  women.rel.status.change <- dplyr::filter(datalist$etable,
                                           eventname == "dissolution" | eventname == "formation",
                                           p2gender == 1,
                                           p2ID %in% women.filtered.df$ID)
  women.rel.status.change$woman.ID <- women.rel.status.change$p2ID
  women.debut <- dplyr::filter(datalist$etable,
                               eventname == "debut",
                               p1gender == 1,
                               p1ID %in% women.filtered.df$ID)
  women.debut$woman.ID <- women.debut$p1ID
  women.rel.status.change <- rbind(women.rel.status.change, women.debut)
  # Merely to facilitate visual inspection:
  women.rel.status.change <- women.rel.status.change[order(women.rel.status.change$woman.ID,
                                                           women.rel.status.change$eventtime), ]
  women.rel.status.change$change <- NULL  #in case the woman.rel.status.change is an empty df.
  women.rel.status.change$change[women.rel.status.change$eventname == "debut"] <- 0
  women.rel.status.change$change[women.rel.status.change$eventname == "formation"] <- 1
  women.rel.status.change <- women.rel.status.change %>% dplyr::group_by(woman.ID) %>%
    dplyr::mutate(numrels = cumsum(change))
  women.rel.status.change <- women.rel.status.change %>%
    dplyr::group_by(woman.ID) %>%
    dplyr::mutate(timespent = diff(c(eventtime, as.numeric(datalist$itable$population.simtime[1]))))
  # We only to keep rows for the events (periods) that the woman spent not in any relationships
  women.no.rels <- dplyr::filter(women.rel.status.change,
                                 numrels == 0)
  women.no.rels$norels.start <- women.no.rels$eventtime
  women.no.rels$norels.end <- women.no.rels$eventtime + women.no.rels$timespent

  # Now we augment the women.no.rels dataset with the exposure.start and
  #exposure.end data from the women.filtered.df dataset
  women.no.rels.plus <- dplyr::left_join(x = women.no.rels,
                                         y = women.filtered.df,
                                         by = c("woman.ID" = "ID"))

  # The real start of the time period of being single is pmax(norels.start, exposure.start)
  # The real end of the time period of being single is pmin(norels.end, exposure.end)
  women.no.rels.plus$norels.start.real <- pmax(women.no.rels.plus$norels.start,
                                               women.no.rels.plus$exposure.start)
  women.no.rels.plus$norels.end.real <- pmin(women.no.rels.plus$norels.end,
                                             women.no.rels.plus$exposure.end)

  # This is the naive time period spent single, before tidying up
  norels.time <- women.no.rels.plus$norels.end.real - women.no.rels.plus$norels.start.real

  # a vector to see which episodes were REALLY spent single during the exposure time window
  really.norels.time <- norels.time > 0
  women.no.rels.plus.clean <- women.no.rels.plus[really.norels.time, ]
  women.no.rels.plus.clean$norels.timespent <- women.no.rels.plus.clean$norels.end.real -
    women.no.rels.plus.clean$norels.start.real

  # For the Harling version
  # How far in the calendar year did the interview with the woman take place
  #(baseline at time of exposure.start, and then each time 1 year thereafter)
  fract.survey <- women.no.rels.plus.clean$exposure.start %% 1

  # How far in the calendar year did the period of being single start?
  fract.norels.start.real <- women.no.rels.plus.clean$norels.start.real %% 1

  move.to.last.survey <- ceiling(fract.survey - fract.norels.start.real)

  women.no.rels.plus.clean$norels.start.surveytime<-ceiling(women.no.rels.plus.clean$norels.start.real)+
    fract.survey - move.to.last.survey

  women.no.rels.plus.clean$norels.end.surveytime<-floor(women.no.rels.plus.clean$norels.end.real)+
    fract.survey - move.to.last.survey

  # Only counting time spent single in whole calendar year blocks
  women.no.rels.plus.clean$norels.timespent.Harling <- pmax(0,
                                                            women.no.rels.plus.clean$norels.end.surveytime -
                                                              women.no.rels.plus.clean$norels.start.surveytime)

  # Only counting time spent single in whole calendar year blocks
  women.no.rels.plus.clean %>% dplyr::mutate(norels.timespent.Harling = pmax(0, norels.end.surveytime -
                                                                               norels.start.surveytime))


  if (type == "Strict"){
    # Now we sum the true time spent being single in the exposure time window, per woman
    norels.timespent.df <- dplyr::summarise(dplyr::group_by(women.no.rels.plus.clean, woman.ID),
                                            sum.norels.timespent = sum(norels.timespent))
    norels.timespent.tobereturned.df <- norels.timespent.df
  } else { # That means the type is "Harling". Now we only sum norels.timespent in blocks of one year
    norels.timespent.df <- dplyr::summarise(dplyr::group_by(women.no.rels.plus.clean, woman.ID),
                                            sum.norels.timespent = sum(norels.timespent.Harling))
  }
  return(norels.timespent.df)
}

prevalence.calculator <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timepoint = 30){

  # First we only take the data of people who were alive at the timepoint
  DTalive.infected <- alive.infected(datalist = datalist, timepoint = timepoint)


  DTalive.infected <- DTalive.infected %>%
    mutate(age = timepoint - TOB)

  DTalive.infected.agegroup <- subset(DTalive.infected, age >= agegroup[1]  &
                                        age < agegroup[2])

  raw.df <- data.frame(DTalive.infected.agegroup)
  if(nrow(raw.df)>0){
    # Now we apply dplyr function to get sum of cases and sum of exposure.time per gender.
    prevalence.df <- dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                      popsize = length(Gender),
                                      sum.cases = sum(Infected),
                                      pointprevalence = sum(Infected) / length(Gender),
                                      pointprevalence.95.ll = as.numeric(
                                        binom.test(x=sum(Infected),n=length(Gender))$conf.int)[1],
                                      pointprevalence.95.ul = as.numeric(
                                        binom.test(x=sum(Infected),n=length(Gender))$conf.int)[2]
    )
    prevalence.all.df <- dplyr::summarise(raw.df,
                                          popsize = length(Gender),
                                          sum.cases = sum(Infected),
                                          pointprevalence = sum(Infected) / length(Gender),
                                          pointprevalence.95.ll = as.numeric(
                                            binom.test(x=sum(Infected),n=length(Gender))$conf.int)[1],
                                          pointprevalence.95.ul = as.numeric(
                                            binom.test(x=sum(Infected),n=length(Gender))$conf.int)[2]
    )


    prevalence.all.df <- cbind(Gender = NA, prevalence.all.df)

    prevalence.df <- rbind(prevalence.df, prevalence.all.df)
  } else {

    prevalence.df <- as.data.frame(matrix(NA, 3, 6))
    names(prevalence.df) <- c("Gender", "popsize", "sum.cases",
                              "pointprevalence", "pointprevalence.95.ll",
                              "pointprevalence.95.ul")
  }

  return(prevalence.df)
}

alive.infected <- function(datalist = datalist,
                           timepoint = 40,
                           site = "All") {
  # arguments are the personlog data.table and a point in time
  DT <- datalist$ptable

  if (site == "All") {
    DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint)
  } else{

    DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint & pfacility == site)
  }

  # Now we allocate infection status to all alive people in our table
  DTalive <-  DTalive %>%
    dplyr::mutate(Infected = (timepoint >= InfectTime))

  return(DTalive)
}

agemix.df.maker <- function(datalist){

  if(!is.list(datalist) |
     names(datalist)[[1]] != "ptable" |
     names(datalist)[[2]] != "rtable") {
    stop("Datalist wrong type")
  }

  dfrmale <- datalist$rtable %>%
    data.frame() %>%
    dplyr::rename(ID = ID1) %>%
    dplyr::mutate(relid = paste0(ID, ID2))

  dfrfemale <- datalist$rtable %>%
    data.frame() %>%
    dplyr::rename(ID = ID2) %>%
    dplyr::mutate(relid = paste0(ID1, ID))

  dfmale <- datalist$ptable %>%
    data.frame() %>%
    dplyr::filter(Gender == 0) %>%
    dplyr::left_join(dfrmale, by = "ID") %>%
    dplyr::select(-ID2)

  dffemale <- datalist$ptable %>%
    data.frame() %>%
    dplyr::filter(Gender == 1) %>%
    dplyr::left_join(dfrfemale, by = "ID") %>%
    dplyr::select(-ID1)

  df <- dplyr::bind_rows(dfmale, dffemale) %>%
    dplyr::arrange(Gender, ID, relid, FormTime) %>%
    dplyr::group_by(Gender, ID, relid) %>%
    dplyr::mutate(episodeorder = row_number(),
                  agerelform = FormTime - TOB,
                  agerelform = dplyr::first(agerelform)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Gender = factor(Gender, labels = c("male", "female")),
                  pagerelform = ifelse(Gender == "male",
                                       agerelform - AgeGap,
                                       agerelform + AgeGap))

  return(df)
}

pattern.modeller <- function(dataframe, agegroup,
                             timepoint, timewindow,
                             start = FALSE) {
  #Warnings
  if (!is.data.frame(dataframe)) {
    stop("dataframe wrong type")
  }

  if (length(agegroup) != 2) {
    stop("Need an upper and lower bound")
  }

  if (timepoint < 1) {
    stop("Time must be at least 1 year")
  }

  if (timewindow < 0) {
    stop("Window must be a whole number")
  }

  #Constants
  time <- timepoint
  window <- timepoint - timewindow
  lwrage <- agegroup[1]
  uprage <- agegroup[2]


  if (start == TRUE) {
    # This only includes relationships that started
    # in the time window

    men <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(episodeorder == 1 &
                      (FormTime <= time & FormTime >= window) &
                      age >= lwrage &
                      age < uprage &
                      Gender == "male" &
                      TOD > time) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

    women <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(episodeorder == 1 &
                      (FormTime <= time & FormTime >= window) &
                      age >= lwrage &
                      age < uprage &
                      Gender == "female" &
                      TOD > time) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

  } else {
    # This includes all relationships that were ongoing
    # at somepoint during the time window, but may have
    # started long before the time window.

    men <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(FormTime <= time &
                      DisTime > window &
                      age >= lwrage &
                      age < uprage &
                      Gender == "male" &
                      TOD > time) %>%
      dplyr::distinct(ID, relid, .keep_all = TRUE) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

    women <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(FormTime <= time &
                      DisTime > window &
                      age >= lwrage &
                      age < uprage &
                      Gender == "female" &
                      TOD > time) %>%
      dplyr::distinct(ID, relid, .keep_all = TRUE) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

  }

  comb <- bind_rows(men, women)

  modoutput <- matrix(nrow = 1, ncol = 14) %>%
    as.data.frame()
  colnames(modoutput) <- c("slopem", "slopew", "interceptm", "interceptw",
                           "powerm",  "lowerpowerm", "upperpowerm", "powerw",
                           "lowerpowerw", "upperpowerw", "bvarm", "bvarw",
                           "wvarm", "wvarw")

  agemix.pieces <- list(comb, modoutput)

  return(agemix.pieces)

}

input.params.creator <- function(mortality.normal.weibull.shape = 5,
                                 mortality.normal.weibull.scale = 65,
                                 mortality.normal.weibull.genderdiff = 0,
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
                                 person.eagerness.woman.dist.gamma.a = 0.231989836885,#0.15#0.1
                                 person.eagerness.woman.dist.gamma.b = 45,#70#100#3.5#5#10#20#170
                                 person.agegap.man.dist.type = "normal",
                                 person.agegap.woman.dist.type = "normal",
                                 person.agegap.man.dist.normal.mu = 0, #-5
                                 person.agegap.woman.dist.normal.mu = 0, #2.5
                                 person.agegap.man.dist.normal.sigma = 1,
                                 person.agegap.woman.dist.normal.sigma = 1,
                                 formation.hazard.agegapry.numrel_man = -0.5,
                                 formation.hazard.agegapry.numrel_woman = -0.5,
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
                                 formation.hazard.agegapry.eagerness_diff = -0.048,#-0.110975
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
                                 hivtransmission.param.f1 = log(5), #~1.6 =>hazard is x5 in 15yo
                                 hivtransmission.param.f2 = log(log(2.5) / log(5)) / 5,
                                 mortality.aids.survtime.C = 62,
                                 mortality.aids.survtime.k = -0.2,
                                 conception.alpha_base = -2.35, #-3
                                 diagnosis.baseline = -100,
                                 monitoring.cd4.threshold = 0.1,#Treatment will not start before schedule
                                 monitoring.fraction.log_viralload = 0.3,
                                 population.msm = "yes",
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
                                 monitoring.cd4.threshold.prestudy = 350,
                                 monitoring.cd4.threshold.instudy.controlstage = 350
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
    input.params$monitoring.cd4.threshold.prestudy <- 350
    input.params$monitoring.cd4.threshold.instudy.controlstage <- 350
    input.params$monitoring.cd4.threshold.poststudy <- 20000 #post study everyone was now on test and study
    #input.params$monitoring.cd4.threshold.instudy.transitionstage <- Inf
    #input.params$monitoring.cd4.threshold.instudy.interventionstage <- Inf
  }

  return(input.params)
}




