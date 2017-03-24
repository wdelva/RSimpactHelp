#' Calculate the number of HIV infections that occured between two time points.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow HIV transmission that took place between this simulation time e.g timewindow = c(20, 30)
#' @param agegroup HIV infection within this agegroup.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return data.frame with number of infection aggregated by Gender
#' @examples
#' data(datalist)
#' hiv.infected <- hiv.infected.group(datalist = datalist,
#' agegroup=c(15, 30), timewindow=c(20, 40), site="All")
#' hiv.infected
#'
#' @import dplyr
#' @export

hiv.infected.group <- function(datalist = datalist,
                               agegroup = c(15, 30),
                               timewindow = c(20, 30), site="All"){

  person.alive.infected <- age.group.time.window(datalist = datalist,
                                                 agegroup = agegroup,
                                                 timewindow = timewindow, site="All")

  person.alive.infected <- subset(person.alive.infected, InfectTime!=Inf)

  #exclude TOD but this should be excluded already
  person.alive.infected <- subset(person.alive.infected, InfectTime > exposure.start &
                                    InfectTime <= exposure.end)
  #aggregate by gender
  hiv.infect.count <- data.frame(dplyr::summarise(dplyr::group_by(person.alive.infected, Gender),
                                                  InfectedCount = n()))

  return(hiv.infect.count)
}
