#' A function that returns the total number of life-years lived between two points in simulation
#' time for a particular age group and gender
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow alive people within this simulation time e.g timewindow = c(15, 30).
#' @param agegroup alive people within this agegroup.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of life-years lived aggregated by gender.
#' @examples
#' life.years.lived <- life.years.lived(datalist = datalist, agegroup=c(15, 30), timewindow=c(15, 30), site="All")
#'
#' @importFrom magrittr %>%
#' @import dplyr

life.years.lived <- function(datalist = datalist, agegroup = c(15, 30), timewindow = c(15, 30), site="All"){

  person.alive.timewindow <- age.group.time.window(datalist = datalist,
                                                   agegroup = agegroup, timewindow = timewindow, site="All")


  person.alive.timewindow <- person.alive.timewindow %>% mutate(exposure.time.Alive = exposure.end - exposure.start)
  person.alive.timewindow <- subset(person.alive.timewindow, exposure.time.Alive > 0)

  ## Sum all people's contibution years

  lifeyears.lived <- data.frame(dplyr::summarise(dplyr::group_by(person.alive.timewindow, Gender),
                                                   TotalExposed = n(),
                                                   LifeYearsLived = sum(exposure.time.Alive) ))

  return(lifeyears.lived)

}

