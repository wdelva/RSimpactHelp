#' Subset the datalist ptable to include only people who existed within the time window and of age
#' group.
#' Add the
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow Point in time at which the subset should be created and HIV status should be evaluated.
#' @param agegroup Individuals who have this agegroup within this timewindow
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a data.table that only includes people who were alive at the timepoint and that records their HIV status.
#' @examples
#' agegroup.timewindow.df <- age.group.time.window(datalist = datalist, agegroup = c(15,30),timewindow = c(15,30),site="All")

age.group.time.window <- function(datalist = datalist,
                                  agegroup = c(15,30),
                                  timewindow = c(15,30),
                                  site="All"){
  DT <- datalist$ptable
  if(site=="All"){
    DTexists.timewindow <- DT
  }else{
    facilities.df <- read.cv(datalist$itable$facilities.geo.coords)
    facilities.df <- filter(facilities.df, Facility = site[1])
    DTexists.timewindow <- subset(DT, XCoord==facilities.df$Longitude
                      & YCoord==facilities.df$Latitude)
  }

  #remove columns that we do not need
  DTexists.timewindow <- subset(DTexists.timewindow, select = -c(IDF, IDM,TODebut,FormEag,FormEagMSM,HSV2InfectOriginID,
                                                                 HSV2InfectTime) )
  #Convert lower age of interest into time
  DTexists.timewindow <- DTexists.timewindow %>% mutate(LowerTimeAgeGroup =TOB + agegroup[1])
  DTexists.timewindow$LowerTimeWindow <- timewindow[1]
  DTexists.timewindow <- DTexists.timewindow %>% mutate(exposure.start = pmax(LowerTimeAgeGroup, LowerTimeWindow))

  #Convert upper age of interest into time
  DTexists.timewindow <- DTexists.timewindow %>% mutate(UpperTimeAgeGroup = TOB + agegroup[2])
  DTexists.timewindow$UpperTimeWindow <- timewindow[2]
  DTexists.timewindow <- DTexists.timewindow %>% mutate(exposure.end = pmin(UpperTimeAgeGroup, UpperTimeWindow))

  #Exposure time (Everyone with exposure time greater than 0)
  DTexists.timewindow <- DTexists.timewindow %>% mutate(exposure.time = exposure.end - exposure.start)
  DTexists.timewindow <- DTexists.timewindow %>% mutate(real.exposure.time = exposure.time > 0)

  #remove those that did not have exposure (died before timewindow or outside min age group)
  DTexists.timewindow <- subset(DTexists.timewindow, real.exposure.time == TRUE )


  return(DTexists.timewindow)
}
