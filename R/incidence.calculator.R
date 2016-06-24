#' Calculate HIV incidence, overall and stratified.
#'
#' Calculate the HIV incidence in a time window and for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound) that should be retained, e.g. c(20, 30)
#' @return a dataframe with cases, exposure time, incidence estimate and surrounding confidence bounds,
#' for the specified time window and age group, overall, and stratified by gender
#' @examples
#' incidence.df <- incidence.calculator(datalist = datalist, agegroup = c(15, 30), timewindow = c(20, 30))

incidence.calculator <- function(datalist = datalist,
                                 agegroup = c(15, 30),
                                 timewindow = c(20, 30)){
  time.of.lowerbound.agegroup <- datalist$ptable$TOB + agegroup[1]
  time.of.lowerbound.timewind <- timewindow[1]
  exposure.start <- pmax(time.of.lowerbound.agegroup, time.of.lowerbound.timewind)
  time.of.upperbound.agegroup <- datalist$ptable$TOB + agegroup[2]
  time.of.upperbound.timewind <- timewindow[2]
  time.of.HIV.infection <- datalist$ptable$InfectTime
  exposure.end <- pmin(time.of.HIV.infection, pmin(time.of.upperbound.agegroup, time.of.upperbound.timewind))
  exposure.time <- exposure.end - exposure.start # This is the naive exposure time, before tidying up
  real.exposure.time <- exposure.time > 0 # We create a vector to see who REALLY had exposure time
  exposure.time[real.exposure.time == FALSE] <- 0

  # Now we check who of the people with the real exposure time had the event
  # Their InfectTime must be after their exposure.time started, and before or at exposure.end
  infection.after.exposure.start <- datalist$ptable$InfectTime > exposure.start
  infection.before.or.at.exposure.end <- datalist$ptable$InfectTime <= exposure.end
  infection.in.timewindow <- infection.after.exposure.start & infection.before.or.at.exposure.end


  datalist$ptable$incident.case <- infection.in.timewindow
  datalist$ptable$exposure.times <- exposure.time

  raw.df <- data.frame(datalist$ptable)

  # Now we apply some dplyr function to get the sum of cases and sum of exposure.time per gender.
  incidence.df <- dplyr::summarise(group_by(raw.df, Gender),
                                   sum.exposure.time = sum(exposure.times),
                                   sum.incident.cases = sum(incident.case),
                                   incidence = sum(incident.case) / sum(exposure.time),
                                   incidence.95.ll = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.time))$conf.int)[1],
                                   incidence.95.ul = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.time))$conf.int)[2]
                                   )

  # Now we add the overall incidence to this dataframe
  incidence.all.df <- cbind(Gender = NA,
                            dplyr::summarise(raw.df,
                                             sum.exposure.time = sum(exposure.times),
                                             sum.incident.cases = sum(incident.case),
                                             incidence = sum(incident.case) / sum(exposure.time),
                                             incidence.95.ll = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.time))$conf.int)[1],
                                             incidence.95.ul = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.time))$conf.int)[2]
                            ))
  incidence.df <- rbind(incidence.df, incidence.all.df)
  return(incidence.df)
}

# Incidence is calculated in 3 steps.
# 1. Calculate PY of exposure per person
# 2. Calculate whether the person had the event or not
# 3. Divide events by sum of PY.

# 1. Calculation of PY of exposure per person. Exposure time starts at the max
# of the lower bound of the timewindow and the time at which the person reaches
# the lower bound of the age group. Exposure time ends at the min of the upper
# bound of the timewindow, the time at which the person reaches the upper bound
# of the age group, and the time at which the person gets infected with HIV.
# Any negative exposure time will be reset to zero.