#' Subset the datalist ptable to include only people who existed within the time window and of
#' agegroup.
#' Add the eposure start and end time.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Individuals within this agegroup
#' @param timewindow The window of interest within the simulation time e.g timewindow = c(15, 30)
#' @param site Select only the particular site from the study, if "all", use all sites.
#' @return a data.table that only includes people who were alive at the timewindow and are within
#' the agegroup
#' @examples
#' data(datalist)
#' agegroup.timewindow.df <- age.group.time.window(datalist = datalist, agegroup = c(15, 30),
#' timewindow = c(10, 40), site="All")
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

age.group.time.window <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timewindow = c(15, 30),
                                  site = "All"){

  DT <- datalist$ptable

  if(site=="All"){
    DTexists.timewindow <- DT
  }else{
    DTexists.timewindow <- subset(DT, pfacility==site)
  }

  #Convert lower age of interest into time
  DTexists.timewindow <- DTexists.timewindow %>%
    dplyr::mutate(LowerTimeAgeGroup = TOB + agegroup[1],
                  LowerTimeWindow = timewindow[1],
                  exposure.start = pmax(LowerTimeAgeGroup, LowerTimeWindow),
                  UpperTimeAgeGroup = TOB + agegroup[2], #Convert upper age of interest into time
                  UpperTimeWindow = timewindow[2],
                  exposure.end = pmin(TOD, pmin(UpperTimeAgeGroup, UpperTimeWindow)),
                  exposure.time = exposure.end - exposure.start,
                  #Exposure time (Everyone with exposure time greater than 0)
                  real.exposure.time = exposure.time > 0
                  ) %>%
    subset(real.exposure.time == TRUE) %>%
    #remove those that did not have exposure (died before timewindow or outside min age group)
    as.data.frame

  return(DTexists.timewindow)
}
