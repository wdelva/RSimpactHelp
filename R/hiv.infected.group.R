#' Calculate the number of HIV infections that occured between two time points.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow HIV transmission that took place between this simulation time e.g timewindow = c(20,30)
#' @param gender HIV infected by this gender.
#' @param agegroup HIV infection within this agegroup.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a list with the transmission network data, as required by the epi2tree function.
#' @examples
#' hiv.infected <- hiv.infected.group(datalist = datalist, agegroup=c(15,30), timewindow=c(20,40), gender="Male", site="All")

hiv.infected.group <- function(datalist = datalist,
                                agegroup = c(15,30),
                                timewindow = c(20,30),
                                gender = "Male", site="All"){
  gender.id <-1
  if(gender!="Male"){gender.id = 0}
  #person.alive.infected <- alive.infected(datalist, timewindow, site="All") #Remove HIV negative individuals till upper of timewindow

  person.alive.infected <- age.group.time.window(datalist = datalist,
                                                 agegroup = agegroup, timewindow = timewindow,site="All")

  person.alive.infected <- subset(person.alive.infected, InfectTime!=Inf)

  person.alive.infected <- subset(person.alive.infected, InfectTime>exposure.start &
                                    InfectTime<=exposure.end & Gender== gender.id) #include TOD but this should be excluded already

  hiv.infect.count <- nrow(person.alive.infected)

  return(hiv.infect.count)
}
