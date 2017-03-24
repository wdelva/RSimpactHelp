#' Calculate time spent being single.
#'
#' Calculate the time spent not being in any relationships.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound)
#' that should be retained, e.g. c(15, 30)
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound)
#' that should be retained, e.g. c(20, 30)
#' @param type If type is "Strict", all time spent being not being in any relationship
#' will be excluded from exposure time.
#' If type is "Harling", time will be excluded from exposure time in blocks of one year,
#' if the person spent that entire block not in any relationship.
#' @return a vector with time spent not being in any relationships, per person.
#' @examples
#' data(datalist)
#' timespentsingle <- timespentsingle.calculator(datalist = datalist,
#' agegroup = c(15, 30), timewindow = c(20, 30), type = "Strict")
#' timespentsingle
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export


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
