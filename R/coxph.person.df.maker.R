#' Create a dataframe for Cox Proportional Hazard Regression analysis at the person-level.
#'
#' Create a dataframe with one row per woman who is at least for some time during the time window HIV-negative and between 15 and 30 years old.
#'
#' We can treat the data as left-truncated (women must be HIV-negative at the time of entering the cohort of 15>= year old women).
#' Because we assume that the study cohort's observation time starts at some arbitrary point in the simulation
#' (e.g. 20 years after HIV is introduced), some women will enter the cohort at an age > 15 (but still < 30).
#' time is the left-truncated age of entering the cohort
#' time2 is the age at event time or age at the end of the time window or agegroup if the event wasn't observed by then.
#' The event parameter indicates whether or not the HIV transmission event was observed or not.
#' type indicates that we have right-censored data.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound) that should be retained, e.g. c(30, 40)
#' @return a dataframe of individuals, with age at the time of entering the cohort (age.start),
#' age at the time of event or censoring (i.e. reaching agegroup boundary or time window boundary) time (age.end),
#' and an indicator of whether HIV infection happened.
#' @examples
#' coxph.person.df <- coxph.person.df.maker(datalist = datalist, agegroup = c(15, 30), timewindow = c(30, 40))



coxph.person.df.maker <- function(datalist = datalist.test,
                                 agegroup = c(15, 30),
                                 timewindow = c(30, 40)){
  time.of.lowerbound.agegroup <- datalist$ptable$TOB + agegroup[1]
  time.of.lowerbound.timewind <- timewindow[1]
  exposure.start <- pmax(time.of.lowerbound.agegroup, time.of.lowerbound.timewind)
  datalist$ptable$age.start <- exposure.start - datalist$ptable$TOB
  # But people who reach the lower boundary of the age group in the course of the simulation must
  # wait until the next survey round to join the cohort.
  age.at.next.survey <- ceiling(exposure.start)-datalist$ptable$TOB
  datalist$ptable$age.start[exposure.start==time.of.lowerbound.agegroup] <- age.at.next.survey[exposure.start==time.of.lowerbound.agegroup]

  time.of.upperbound.agegroup <- datalist$ptable$TOB + agegroup[2]
  time.of.upperbound.timewind <- timewindow[2]
  time.of.HIV.infection <- datalist$ptable$InfectTime
  exposure.end <- pmin(time.of.HIV.infection, pmin(time.of.upperbound.agegroup, time.of.upperbound.timewind))
  datalist$ptable$age.end <- exposure.end - datalist$ptable$TOB



  # We add the indicator "transmissionevent" to indicate whether HIV transmission happened or not.
  datalist$ptable$age.HIV.infection <- datalist$ptable$InfectTime - datalist$ptable$TOB
  # Their age.HIV.infection must be larger than age.start and smaller or equal to age.end
  infection.after.age.start <- datalist$ptable$age.HIV.infection > datalist$ptable$age.start
  infection.before.or.at.age.end <- datalist$ptable$age.HIV.infection <= datalist$ptable$age.end
  datalist$ptable$transmissionevent <- infection.after.age.start & infection.before.or.at.age.end

  coxph.person.df <- data.frame(datalist$ptable)
  coxph.person.df$ID.Gender <- paste0(coxph.person.df$ID, ".", coxph.person.df$Gender)

  # Now we need to create a "timeline dataframe" that we can tmerge with coxph.person.df
  # There must be one row for each one-year time interval per person, simulating surveys during the cohort study
  coxph.person.df$time.start <- coxph.person.df$age.start + coxph.person.df$TOB
  coxph.person.df$time.stop <- coxph.person.df$age.end + coxph.person.df$TOB
  exposure.time <- coxph.person.df$time.stop - coxph.person.df$time.start # This is the naive exposure time, before tidying up
  real.exposure.time <- exposure.time > 0 # We create a vector to see who REALLY had exposure time
  coxph.person.df <- subset(coxph.person.df, real.exposure.time)
  Surv.df <- tmerge(data1 = coxph.person.df, data2 = coxph.person.df, id = ID.Gender,
                    tstart = time.start, tstop = time.stop)
  survey.times <- expand.grid(ID.Gender = coxph.person.df$ID.Gender, surveys = seq(timewindow[1], timewindow[2]))
  Surv2.df <- tmerge(data1 = Surv.df, data2 = survey.times, id = ID.Gender, survey.event = event(surveys))
  Surv2.df <- Surv2.df[Surv2.df$Gender==1, ]#subset(Surv2.df, Gender == 1)
  Surv2.df$age.interval.start <- Surv2.df$tstart - Surv2.df$TOB
  Surv2.df$age.interval.stop <- Surv2.df$tstop - Surv2.df$TOB
  Surv2.df$infection.event <- as.numeric(Surv2.df$InfectTime==Surv2.df$tstop)

  # Next we must get the age difference with "the most recent partner" at the start of each interval.
  # How is this partner defined?
  # 1. The relationship must be ongoing, as defined by a start time of the RELATIONSHIP (not the episode) before the survey time
  # and an end time of the RELATIONSHIP (not the episode) after the survey. In other words, if the survey falls in between two episodes, the
  # relationship can still be deemed "ongoing", but only if:
  # 2. There was at least one episode that overlapped with the year before the survey
  # (i.e. its start date was before the survey and its end date was after the "survey - 1" year mark)
  # 3. There was at least one episode that overlapped with the year after the survey
  # (i.e. its start date was before the "survey + 1" year mark and its end date was after the survey)
  # 4. There may be multiple relationship that match this definition of ongoing relationships at the time of the survey.
  # Of those, we will choose the relationship with the start time of the RELATIONSHIP closest to the survey.

  episode.df <- agemix.df.maker(datalist = datalist.test)
  # We have to loop through the values of survey.event (the time points of the surveys)
  Surv3.df <- data.frame()
  for (survey.i in seq(timewindow[1], timewindow[2])){
    recent.rels <- most.recent.partner(episode.df = episode.df, survey.time = survey.i)
    recent.rels$ID <- as.numeric(recent.rels$ID)
    Surv.df.expanded <- left_join(x = subset(Surv2.df, tstart == survey.i),
                                  y = recent.rels,
                                  by = c("ID", "Gender", "time.start"))
    Surv3.df <- rbind(Surv3.df, Surv.df.expanded)
  }
  #Surv3.df$age.interval.start <- Surv3.df$tstart - Surv3.df$TOB
  #Surv3.df$age.interval.stop <- Surv3.df$tstop - Surv3.df$TOB
  Surv3.df <- arrange(Surv3.df, ID, age.interval.start)
  Surv3.df$age.difference.at.survey <- Surv3.df$page.at.survey - Surv3.df$age.interval.start
  #Surv3.df$seroconversion <- as.numeric(Surv3.df$transmissionevent)

  # Let's now go back to the original Surv2.df object, and augment it with the new variables of Surv3.df
  #Surv4.df <- tmerge(data1 = Surv2.df, data2 = Surv3.df, id = ID.Gender, age.interval.start = tdc(tstart, age.interval.start))
  #Surv5.df <- tmerge(data1 = Surv4.df, data2 = Surv3.df, id = ID.Gender, age.interval.stop = tdc(tstop, age.interval.stop))

  # Adding age difference at survey time dependent covariate
  Surv4.df <- tmerge(data1 = Surv2.df, data2 = Surv3.df, id = ID.Gender, age.difference.at.survey = tdc(age.interval.stop, age.difference.at.survey))


  #coxph(Surv(age.interval.stop, infection.event==1) ~ age.difference.at.survey, Surv4.df)

  return(Surv4.df)
}


most.recent.partner <- function(episode.df = agemix.df, survey.time){
  # 1. The relationship must be ongoing, as defined by a start time of the RELATIONSHIP (not the episode) before the survey time
  # and an end time of the RELATIONSHIP (not the episode) after the survey. In other words, if the survey falls in between two episodes, the
  # relationship can still be deemed "ongoing", but only if:
  # 2. There was at least one episode that overlapped with the year before the survey
  # (i.e. its start date was before the survey and its end date was after the "survey - 1" year mark)
  # 3. There was at least one episode that overlapped with the year after the survey
  # (i.e. its start date was before the "survey + 1" year mark and its end date was after the survey)
  # 4. There may be multiple relationship that match this definition of ongoing relationships at the time of the survey.
  # Of those, we will choose the relationship with the start time of the RELATIONSHIP closest to the survey.
  episode.df$survey.between.start.and.end.rel <- episode.df$rel.start < survey.time & episode.df$rel.end > survey.time
  episode.df$before.survey <- episode.df$FormTime < survey.time & (episode.df$DisTime + 1) > survey.time
  episode.df$after.survey <- episode.df$FormTime < (survey.time + 1) & episode.df$DisTime > survey.time

  # Let's now only keep the episodes where Gender == female
  episode.female.df <- subset(episode.df, Gender == "female")
  episode.female.df$Gender <- as.integer(1)
  # Per relid, the sum of survey.between.start.and.end.rel >= 1 AND the sum of before.survey > 1 AND the sume of after.survey > 1
  rel.ongoing <- dplyr::summarise(group_by(episode.female.df, relid),
                                   ID = first(ID),
                                   Gender = first(Gender),
                                   rel.start = first(rel.start),
                                   rel.end = first(rel.end),
                                   pagerelform = first(pagerelform),
                                   cond1 = sum(survey.between.start.and.end.rel) >= 1,
                                   cond2 = sum(before.survey) >= 1,
                                   cond3 = sum(after.survey) >= 1,
                                   all.cond = sum(survey.between.start.and.end.rel) >= 1 & sum(before.survey) >= 1 & sum(after.survey) >= 1)

  rel.ongoing <- subset(rel.ongoing, all.cond == TRUE)
  rel.ongoing$recency <- survey.time - rel.ongoing$rel.start
  rel.ongoing <- arrange(rel.ongoing, ID, recency)
  # And now we only retain the MOST recent relid per ID
  recent.rels <- dplyr::summarise(group_by(rel.ongoing, ID),
                                  relid = first(relid),
                                  Gender = first(Gender),
                                  rel.start = first(rel.start),
                                  rel.end = first(rel.end),
                                  pagerelform = first(pagerelform),
                                  recency = first(recency))

  # And now the final step: the age of the partner at the time of the survey:
  recent.rels$page.at.survey <- recent.rels$pagerelform - recent.rels$rel.start + survey.time
  recent.rels$time.start <- survey.time

  return(recent.rels)
}


survey.ages <- function(age.start.end){
  start.age <- age.start.end[1]
  end.age <- age.start.end[2]
  survey.ages.vect <- seq(start.age, end.age)
  return(survey.ages.vect)
}
