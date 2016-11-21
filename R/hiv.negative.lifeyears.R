#' A function that returns the total number of HIV-negative life years lived
#' between two points in simulation time, for a particular age group and gender.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow HIV negative time window e.g timewindow = 30.
#' @param gender HIV negative by this gender.
#' @param agegroup HIV ngative within this agegroup.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of HIV-negative life years lived.
#' @examples
#' hiv.negative.lifeyears <- hiv.negative.lifeyears(datalist = datalist, agegroup=c(15,40), timewindow=c(15,30), gender="Male", site="All")

hiv.negative.lifeyears <- function(datalist = datalist, agegroup = c(15, 30),
                                   timewindow = c(15, 30), gender = "Male", site = "All"){
  gender.id <-1
  if(gender!="Male"){gender.id = 0}

  person.alive.negative <- age.group.time.window(datalist = datalist,
                                                 agegroup = agegroup, timewindow = timewindow ,site="All")


  person.alive.negative <- person.alive.negative %>% mutate(exposure.end.Neg = pmin(exposure.end,InfectTime))
  person.alive.negative <- person.alive.negative %>% mutate(exposure.time.Neg = exposure.end.Neg - exposure.start)


  person.alive.negative <- subset(person.alive.negative, exposure.time.Neg > 0 & Gender == gender.id) #HIV negative individuals

  ## Sum all negative people's contibution ages
  hiv.neg.lifeyears <- with(person.alive.negative, sum(exposure.time.Neg))

  return(hiv.neg.lifeyears)

}
