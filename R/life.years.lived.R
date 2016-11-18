#' A function that returns the total number of life-years lived between two points in simulation
#' time for a particular age group and gender
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow alive people within this simulation time e.g timewindow = 30.
#' @param gender alive gender.
#' @param agegroup alive people within this agegroup.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of life-years lived.
#' @examples
#' life.years.lived <- life.years.lived(datalist = datalist, agegroup=c(15,30), timewindow=c(15,30), gender="Male", site="All")

life.years.lived <- function(datalist = datalist, agegroup = c(15,30),
                             timewindow = c(15,30), gender = "Male", site="All"){
  gender.id <-1
  if(gender!="Male"){gender.id = 0}

  person.alive.timewindow <- age.group.time.window(datalist = datalist, 
                                                   agegroup = agegroup, timewindow = timewindow, site="All")

  person.alive.timewindow <- subset(person.alive.timewindow, Gender == gender.id) #HIV negative individuals
  
  person.alive.timewindow <- person.alive.negative %>% mutate(lifeyears = pmin(timewindow[2],MaxAgeTW2) - MinAgeTW1)
  
  ## Sum time lived within this time window
  lifeyears.lived <- with(person.alive.timewindow, sum(lifeyears))

  return(lifeyears.lived)

}


### Disability Adjusted Life Years
## DALY = YLL + YLD (YLL - Years of Life Lost, YLD - Years Lost due to Disability)
#YLL = N*L (N number of deaths and L stardard life expectancy at age of death in years)
#YLD = P * DW (P number of prevalent cases and DW is the Disability Weight ) (healthmetricsandevaluation.org/gbd)
#HIV/AIDS    Average disability weight    Range               Source
#HIV cases                 0.135          0.123 - 0.136         GBD 1990(c) varies with age
#AIDS cases not on ART     0.505                                GBD 1990
#AIDS case on ART         0.167           0.145 - 0.469         GBD 2004






