#' Calculate time spent being single.
#'
#' Calculates the time spent not being in any relationships, among the subset of
#' women-only who are in the specified age groups during the period of interest.
#' This function allows users to specify different definitions of "being
#' single". Under a "Strict" definition, all time spent not being in a
#' relationship will be calculated per person. Under the definition called "Harling", time
#' is excluded from being single if the person spent part of the year in a
#' relationship. The Harling definition is derived from Harling et al. 2014. 
#' This function is called by the \code{\link{incidence.calculator}} function
#' and used to optionally exclude exposure time in incidence calculations.
#'
#' @param datalist The list object that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper 
#'   bound) that should be retained. Should be expressed as a vector of two 
#'   integers. e.g. c(15, 30)
#' @param timewindow Boundaries of the time window (lower bound < time <= upper 
#'   bound) that should be retained. Should be expressed as a vector of two 
#'   integers. e.g. c(20, 30)
#' @param type If type is "Strict", all time spent being not being in any 
#'   relationship will be excluded from exposure time. If type is "Harling", 
#'   time will be excluded from exposure time in blocks of one year, if the 
#'   person spent that entire block not in any relationship.
#' @return A dataframe with a column of women indentifiers (woman.ID) and a
#'   column for the time they spent being single (sum.norels.timespent).
#' @examples
#' cfg <- list() 
#' modeloutput <- RSimpactCyan::simpact.run(configParams = cfg, destDir = "temp") 
#' dl <- readthedata(modeloutput) 
#' timespentsingle <- timespentsingle.calculator(datalist = dl, agegroup = c(15, 30), timewindow = c(20, 30), type = "Strict")
#' timespentsingle
#' 
#' @importFrom magrittr %>%
#' @import dplyr
#' @export


timespentsingle.calculator <- function(datalist = datalist,
                                       agegroup = c(15, 30),
                                       timewindow = c(20, 30),
                                       type = "Strict"){
  
  # Scalars
  lwr.agegroup <- agegroup[1]
  upr.agegroup <- agegroup[2]
  lwr.timewindow <- timewindow[1]
  upr.timewindow <- timewindow[2]
  pop.simtime <- datalist$itable$population.simtime[1]
  
  # Data frames from the list
  df <- datalist$ptable
  df.event <- datalist$etable
  
  # Calculate naive exposure time
  # 1. Determine the time person is 15yo
  # 2. Determine the time person is 30yo
  # 3. Their exposure starts whichever comes last: lwr bound age, or lwr window
  # 4. Their exposure ends whichever comes 1st: upr bound age, upr window, or
  # infected with HIV
  # 5. Subtract exposure start from exposure end
  df1 <- df %>%
    dplyr::mutate(time.lwr.agegroup = TOB + lwr.agegroup,
                  time.upr.agegroup = TOB + upr.agegroup,
                  exposure.start = pmax(time.lwr.agegroup, lwr.timewindow),        
                  exposure.end = pmin(InfectTime, 
                                      time.upr.agegroup, 
                                      upr.timewindow),
                  exposure.time = exposure.end - exposure.start)
  
  # If the exposure time is negative, then they had no exposure during the
  # timewindow or agegroup of interest
  df2 <- df1 %>%
    dplyr::mutate(had.exposure.time = exposure.time > 0,
                  exposure.time = ifelse(exposure.time < 0, 
                                         0, 
                                         exposure.time))
  
  # Now we check who of the people with the real exposure time had the event
  # Their InfectTime must be after their exposure.time started, 
  # and before or at exposure.end
  df3 <- df2 %>%
    dplyr::mutate(incident.case = (InfectTime > exposure.start & 
                                     InfectTime <= exposure.end))
  
  # Create dataset with only women
  women <- df3 %>%
    dplyr::filter(exposure.time > 0 & Gender == 1)
  
  # Create dataset of dissolution and formation events for women
  # Keep only events for the women that are in the women dataset
  women.rel.change <- df.event %>%
    dplyr::filter((eventname == "dissolution" | eventname == "formation"),
                  p2gender == 1,
                  p2ID %in% women$ID) %>%
    dplyr::rename(woman.ID = p2ID)
  
  # Create dataset of sexual debut events for women
  women.deb <- df.event %>%
    dplyr::filter(eventname == "debut",
                  p1gender == 1,
                  p1ID %in% women$ID) %>%
    dplyr::rename(woman.ID = p1ID)
  
  # Combine both events datasets and create change variable
  women.rel.change2 <- bind_rows(women.rel.change, women.deb) %>%
    dplyr::arrange(woman.ID, eventtime) 
  
  # Create variables for:
  # 1. Indicator of when a relationship was added or subtracted at event time
  # 2. The total number of relationships woman had at a given event time
  # 3. The difference between each subsequent event time, within a woman
  women.rel.change3 <- women.rel.change2 %>%
    dplyr::mutate(change = ifelse(eventname == "debut", 0,
                                  ifelse(eventname == "formation", 1, -1))) %>%
    dplyr::group_by(woman.ID) %>%
    dplyr::mutate(numrels = cumsum(change),
                  timespent = diff(c(eventtime, pop.simtime))) # Where only 1 event
                                                               # use simulation time
  
  # We only to keep rows for the events (periods) that the woman spent not in 
  # any relationships
  # Create vars for start and end of not being in a relationship
  women.zero.rels <- women.rel.change3 %>%
    dplyr::filter(numrels == 0) %>%
    dplyr::mutate(norels.start = eventtime,
                  norels.end = eventtime + timespent)
  
  # Now we augment the women.no.rels dataset with the exposure.start and
  # exposure.end data from the women
  women.zero.rels.plus <- women.zero.rels %>%
    dplyr::left_join(women, by = c("woman.ID" = "ID"))
  
  # Create indicator of real start and end time of being single 
  # 1. Real start is whichever comes last: start of being single, or exposure start
  # 2. Real end is whichever comes 1st: end of being single, or exposure end
  # 3. Time spent not being in a relationship
  # Then remove rows where no time was spent single
  women.zero.rels.plus2 <- women.zero.rels.plus %>%
    dplyr::mutate(norels.start.real = pmax(norels.start, exposure.start),
                  norels.end.real = pmin(norels.end, exposure.end),
                  norels.time = norels.end.real - norels.start.real) %>%
    dplyr::filter(norels.time > 0)
  
  # For the Harling version, need indicators for:
  # 1. How far in calendar year did interview with the woman take place
  # The baseline survey is at start of exposure, and then 1 year thereafter
  # 2. How far in calendar year did period of being single start
  # 3. Move to last survey
  # 4. In what annual survey did being single start
  # 5. In what annual survey did being single end
  women.zero.rels.plus3 <- women.zero.rels.plus2 %>%
    dplyr::mutate(frac.survey = exposure.start %% 1, # 1 
                  frac.norels.start.real = norels.start.real %% 1, # 2 
                  move.to.last.survey = ceiling(frac.survey - frac.norels.start.real), # 3 
                  norels.start.surveytime = 
                    ceiling(norels.start.real) + frac.survey - move.to.last.survey, # 4
                  norels.end.surveytime = 
                    floor(norels.end.real) + frac.survey - move.to.last.survey) # 5
  
  # Only counting time spent single in whole calendar year blocks
  women.zero.rels.plus4 <- women.zero.rels.plus3 %>%
    dplyr::mutate(norels.time.Harling = pmax(0, (norels.end.surveytime - 
                                                   norels.start.surveytime)))
  
  if (type == "Strict"){
    
    # Now we sum the true time spent being single in the exposure time window, 
    # per woman
    norels.time.df <- women.zero.rels.plus4 %>%
      dplyr::group_by(woman.ID) %>%
      dplyr::summarise(sum.norels.timespent = sum(norels.time))
    
  } else { 
    
    # That means the type is "Harling". 
    # Now we only sum norels.timespent in blocks of one year
    norels.time.df <- women.zero.rels.plus4 %>%
      dplyr::group_by(woman.ID) %>%
      dplyr::summarise(sum.norels.timespent = sum(norels.time.Harling))
    
  }
  
  return(norels.time.df)
  
}
