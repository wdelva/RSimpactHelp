#' Calculate the HIV incidence in a time window and for specific age groups and gender strata
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound) that should be retained, e.g. c(20, 30)
#' @return a dataframe with cases, exposure time, incidence estimate and surrounding confidence bounds,
#' for the specified time window, stratified by gender and the specified age groups
#' @examples
#' incidence.df <- incidence.calculator(datalist = datalist, agegroup = c(15, 30), timewindow = c(20, 30))


# 1. baseline.df: Everybody who ever existed in the population after HIV was introduced, minus the hiv seeds themselves.
baseline.datamaker <- function(datalist = datalist, hivseed.time = hivseed.time){
  DTP <- datalist$ptable
  DTP.selected.rows <- DTP$TOD > hivseed.time & DTP$InfectType != 0
  baseline.df <- DTP[DTP.selected.rows, ] # This has more variables than we strictly need, but we'll leave this for now.
  baseline.df$ID.Gender <- paste(baseline.df$ID, baseline.df$Gender, sep=".")
  return(baseline.df)
}
# 2. timeline.df: one row per period of exposure time. We can use 1-year intervals as periods of exposure
# except when there is a non-AIDS death or a transmission event, the interval would be shorter.
# m.timeline.df should also include an indicator variable for the (Calendar) time of observation.
# IDm, age0, age1, time0, time1.

timeline.maker <- function(datalist = datalist,
                           endfollowup = endfollowup,
                           hivseed.time = hivseed.time){
  DTP <- datalist$ptable
  DTP.selected.rows <- DTP$TOD > hivseed.time & DTP$InfectType != 0
  baseline.df <- DTP[DTP.selected.rows, ] # This has more variables than we strictly need, but we'll leave this for now.
  baseline.df$ID.Gender <- paste(baseline.df$ID, baseline.df$Gender, sep=".")
  exposure.end <- pmin(pmin(baseline.df$TOD, endfollowup), baseline.df$InfectTime) #ceiling(baseline.df$InfectTime)) # End of exposure time: given by death, end of simulation, or (when) infection (may be picked up by an annual HIV test).
  exposure.start <- pmin(baseline.df$TODebut, exposure.end)
  #timeline.df <- data.frame(ID.Gender=baseline.df$ID.Gender, TOB=baseline.df$TOB, exposure.start, exposure.end)
  timeline.df <- baseline.df
  timeline.df$exposure.start <- exposure.start
  timeline.df$exposure.end <- exposure.end
  someexposure.indicator <- exposure.start != exposure.end # No exposure time is not allowed
  timeline2.df <- timeline.df[someexposure.indicator, ] # No exposure time is not allowed
  assign("timeline2.df", value=timeline2.df)
  return(timeline2.df)
}


timeline.datamaker <- function(datalist = datalist,
                               endfollowup = endfollowup,
                               hivseed.time = hivseed.time){
  baseline.df <- baseline.datamaker(datalist, hivseed.time)
  exposure.end <- pmin(pmin(baseline.df$TOD, endfollowup), ceiling(baseline.df$InfectTime)) # End of exposure time: given by death, end of simulation, or (when) infection (may be picked up by an annual HIV test).
  exposure.start <- pmin(baseline.df$TODebut, exposure.end)
  timeline2.df <- timeline.maker(datalist, endfollowup, hivseed.time)
  someexposure.indicator <- exposure.start != exposure.end # No exposure time is not allowed
  baseline.df <- baseline.df[someexposure.indicator, ] # No exposure time is not allowed
  #return(baseline.timeline.df)
  return(list(baseline.df, timeline2.df))
}
# 3. outcome.df
# m.outcome.hiv.postest.df: one row with the time of first HIV+ test ##newname = event(y,x) Mark an event at time y, with value x. We can choose x=3
# InfectTime, status (==3)
# m.outcome.hiv.infect.df: one row with the exact time of HIV acquisition ##newname = event(y,x) Mark an event at time y, with value x. We can choose x=1
# InfectTime, status (==1)
# m.outcome.survey.ageMRP.df: one row with the age of the MRP at the time of the survey ##newname = tdc(y, x)
# The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP).
# SurveyTime, survey.ageMRP
# m.outcome.infect.ageMRP.df: one row with the age of the MRP at the time of the infection ##newname = tdc(y, x)
# The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP).
# InfectTime, infect.ageMRP


# When did each of the people have a "fiveyear" birthday of turning 20, 25, 30, ... 70 in the course of the simulation?

fiveyearbirthdays <- function(start.end.age){
  start.age <- start.end.age[1]
  end.age <- start.end.age[2]
  first.fiveyearBD <- 5 * ceiling(start.age/5)
  last.fiveyearBD <- 5 * floor(end.age/5)
  fiveyearBDs <- seq(first.fiveyearBD, last.fiveyearBD, by = 5)
  return(fiveyearBDs)
}

fiveyearBT.dataframe.maker <- function(df = baseline.timeline.df,
                                       endfollowup = endfollowup,
                                       hivseed.time = hivseed.time){
  ID.Gender <- df$ID.Gender #= baseline.timeline.df$ID.Gender
  start.age <- -df$TOB + hivseed.time
  end.age <- -df$TOB + endfollowup
  start.end.age.df <- data.frame(start.age = start.age, end.age = end.age)
  fiveyearBD.list.vect <- unlist(list(apply(X = start.end.age.df, MARGIN = 1, FUN = fiveyearbirthdays))) # List of fiveyear birthdays per person
  #nBDs.vect <- unlist(lapply(fiveyearBD.list.vect, length)) # Number of fiveyear birthdays per person
  n.expand <- round((endfollowup - hivseed.time) / 5)
  ID.Gender.expanded <- rep(ID.Gender, each = n.expand)
  TOB.expanded <- rep(df$TOB, each = n.expand)
  fiveyearBT.df <- data.frame(ID.Gender = ID.Gender.expanded,
                              fiveyear.birthtimes = fiveyearBD.list.vect + TOB.expanded,
                              fiveyear.birthdays = fiveyearBD.list.vect)
  return(fiveyearBT.df)
}

## April 2016 Update:
## We create time blocks of 5 years, and age blocks of 5 years.


outcome.datamaker <- function(datalist = datalist,
                              endfollowup = endfollowup,
                              hivseed.time = hivseed.time){
  timeline.list <- timeline.datamaker(datalist, endfollowup, hivseed.time)
  baseline.timeline.df <- tmerge(data1 = timeline.list[[1]], data2 = timeline.list[[2]], id = ID.Gender, tstart = timeline.list[[2]]$exposure.start, tstop = timeline.list[[2]]$exposure.end)
  DTL <- datalist$ltable
  DTR <- datalist$rtable
  DTE <- datalist$etable
  fiveyear.intervals.end <- 5*unique(ceiling(oneyear.intervals.end/5))

  outcome.fiveyearblock <- expand.grid(ID.Gender = baseline.timeline.df$ID.Gender, fiveyear.intervals.end = fiveyear.intervals.end)

  Surv.df <- tmerge(data1 = baseline.timeline.df, data2 = outcome.fiveyearblock, id = ID.Gender, survey.number = event(y = fiveyear.intervals.end, x = fiveyear.intervals.end)) # adding survey rounds

  # Adding the times at which the individuals had "fiveyear" birthdays
  outcome.fiveyearAGEblock <- fiveyearBT.dataframe.maker(df = baseline.timeline.df,
                                                         endfollowup = endfollowup,
                                                         hivseed.time = hivseed.time)

  Surv2.df <- tmerge(data1 = Surv.df, data2 = outcome.fiveyearAGEblock, id = ID.Gender, ageblock.number = event(y = fiveyear.birthtimes, x = fiveyear.birthdays)) # adding survey rounds

  Surv3.df <- tmerge(data1 = Surv2.df, data2 = Surv.df, id = ID.Gender, hiv.infect = event(y = InfectTime)) # adding HIV infection events

  Surv3.df$timeblock.cat <- cut(x = Surv3.df$tstop,
                                breaks = seq(from = 0, to = endfollowup, by = 5),
                                include.lowest = TRUE,
                                right = TRUE)
  max.age <- max(Surv3.df$ageblock.number) + 5
  Surv3.df$ageblock.cat <- cut(x = Surv3.df$tstop - Surv3.df$TOB,
                               breaks = seq(from = 0, to = max.age, by = 5),
                               include.lowest = TRUE,
                               right = TRUE)
  Surv3.df$exposuretime <- Surv3.df$tstop - Surv3.df$tstart
  # Now we sum exposure time, by timeblock.cat (T), gender (G) and ageblock.cat (A)
  Surv3.dt <- data.table(Surv3.df)
  Incidence.dt <- Surv3.dt[ ,py.T.G.A := sum(exposuretime),by = "timeblock.cat,Gender,ageblock.cat"] # And we calculate py
  Incidence.dt <- Surv3.dt[ ,cases.T.G.A := sum(hiv.infect),by = "timeblock.cat,Gender,ageblock.cat"] # And we calculate cases
  Incidence.dt <- Surv3.dt[ ,inc.T.G.A := sum(hiv.infect) / sum(exposuretime),by = "timeblock.cat,Gender,ageblock.cat"] # And we calculate incidence


  #oneyear.intervals.end <- tail(datalist$ltable$Time, -1) # c(1, 2, 3, ..., population.simtime) It does not make sense to have a survey at time=0 because nothing happened before that time.
  #outcome.survey <- expand.grid(ID.Gender = baseline.timeline.df$ID.Gender, oneyear.intervals.end = oneyear.intervals.end)


  #Surv.df <- tmerge(data1 = baseline.timeline.df, data2 = outcome.survey, id = ID.Gender, survey.number = event(y = oneyear.intervals.end, x = oneyear.intervals.end)) # adding survey rounds
  #Surv.df <- tmerge(data1 = Surv.df, data2 = Surv.df, id = ID.Gender, hiv.infect = event(y = InfectTime)) # adding HIV infection events
  #Surv.df$HIVposTime <- ceiling(Surv.df$InfectTime)
  #Surv.df <- tmerge(data1 = Surv.df, data2 = Surv.df, id = ID.Gender, hiv.postest = event(y = HIVposTime, x = rep(3, nrow(Surv.df)))) # adding first HIV positive test events

  #outcome.survey.ageMRP.df: one row with the age of the MRP at the time of the survey ##newname = tdc(y, x)
  # The argument y is assumed to be on the scale of the start and end time, if x is present the count is set to the value of x (e.g. age of MRP).
  #DTR.m <- DTR.w <- DTR
  #DTR.m$ID.Gender <- paste(DTR$IDm, 0, sep = ".")
  #DTR.w$ID.Gender <- paste(DTR$IDw, 1, sep = ".")
  # index of MRP
  #indx.m <- neardate(id1 = Surv.df$ID.Gender, id2 = DTR.m$ID.Gender, y1 = Surv.df$survey.number, y2 = DTR.m$FormTime, best = "prior")
  # WE MUST NOT FORGET TO CHECK IF THE PARTICIPANT WAS IN A RELATIONSHIP WITH THE MRP AT THE TIME OF THE SURVEY.
  # Harling et al. did not do this check.
  #indx.m <- ifelse((DTR.m$DisTime[indx.m] - Surv.df$survey.number) > 0, indx.m, NA) # Now we know the relationship with the MRP was ongoing.
  #indx.w <- neardate(id1 = Surv.df$ID.Gender, id2 = DTR.w$ID.Gender, y1 = Surv.df$survey.number, y2 = DTR.w$FormTime, best = "prior")
  # WE MUST NOT FORGET TO CHECK IF THE PARTICIPANT WAS IN A RELATIONSHIP WITH THE MRP AT THE TIME OF THE SURVEY.
  # Harling et al. did not do this check.
  #indx.w <- ifelse((DTR.w$DisTime[indx.w] - Surv.df$survey.number) > 0, indx.w, NA) # Now we know the relationship with the MRP was ongoing.
  #AgeGaps.m <- DTR.m[indx.m, AgeGap]
  #AgeGaps.w <- DTR.w[indx.w, AgeGap]
  #AgeGaps.m.and.w <- AgeGaps.m
  #AgeGaps.m.and.w[is.na(AgeGaps.m.and.w) & !is.na(AgeGaps.w)] <- AgeGaps.w[is.na(AgeGaps.m.and.w) & !is.na(AgeGaps.w)]
  #outcome.survey.agegapMRP.df <- data.table(ID.Gender= Surv.df$ID.Gender, survey.number=Surv.df$survey.number, AgeGapMRP=AgeGaps.m.and.w)
  #TestSurv.df <- tmerge(data1 = Surv.df, data2 = outcome.survey.ageMRP.df, id = ID.Gender,
  #                  AgeGapMRP = tdc(y = outcome.survey.agegapMRP.df$survey.number, x = outcome.survey.agegapMRP.df$AgeGapMRP))

  # Preparing data for "MIICD" coxph model with interval censored data
  #   The data must contain at last two columns: left and right. For interval censored data, the left
  #   and the right columns indicates lower and upper bounds of intervals respectively. Inf in the right
  #   column stands for right censored observations.

  # First, we forget about the hiv.infect event
  #TestSurv.df <- subset(TestSurv.df, hiv.infect != 1)

  # we treat each observation as a right censored one, except the ones that led to hiv.postest===3
  #TestSurv.df$left <- TestSurv.df$tstop - TestSurv.df$tstart
  #TestSurv.df$right <- Inf
  #TestSurv.df$right[TestSurv.df$hiv.postest==3] <- 1


  # tmerge.dfs <- list(baseline.timeline.df=baseline.timeline.df,
  #                    outcome.survey=outcome.survey)
  #return(TestSurv.df)
  return(Incidence.dt)
}
