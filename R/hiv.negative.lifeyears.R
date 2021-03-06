#' A function that returns the total number of HIV-negative life years lived
#' between two points in simulation time, for a particular age group and gender.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup HIV ngative within this agegroup.
#' @param timewindow HIV negative time window e.g timewindow = c(30, 40).
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of HIV-negative life years lived aggregated by gender.
#' @examples
#' data(datalist)
#' hiv.negative.lifeyears <- hiv.negative.lifeyears(datalist = datalist,
#' agegroup=c(15, 40), timewindow=c(15, 30), site="All")
#' hiv.negative.lifeyears
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

hiv.negative.lifeyears <- function(datalist = datalist, agegroup = c(15, 30),
                                   timewindow = c(15, 30), site = "All"){

  person.alive.negative <- age.group.time.window(datalist = datalist,
                                                 agegroup = agegroup,
                                                 timewindow = timewindow, site="All")


  person.alive.negative <- person.alive.negative %>%
    dplyr::mutate(exposure.end.Neg = pmin(exposure.end, InfectTime))
  person.alive.negative <- person.alive.negative %>%
    dplyr::mutate(exposure.time.Neg = exposure.end.Neg - exposure.start)

  #HIV negative individuals
  person.alive.negative <- subset(person.alive.negative, exposure.time.Neg > 0 )

  ## Sum all negative people's contibution ages
  hiv.neg.lifeyears <- data.frame(dplyr::summarise(dplyr::group_by(person.alive.negative, Gender),
                                                    TotalExposed = n(),
                                                    LifeYears = sum(exposure.time.Neg) ))

  return(hiv.neg.lifeyears)

}
