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
                                  timewindow = c(10,40),
                                  site="All"){
  DT <- datalist$ptable
  if(site=="All"){
    DTexists.timewindow <- subset(DT, TOB < (timewindow[2]- agegroup[1]) & TOD > timewindow[1]) # you need to be at least agegroup[1] to be able to be in this df
  }else{
    facilities.df <- read.cv(datalist$itable$facilities.geo.coords)
    facilities.df <- filter(facilities.df, Facility = site[1])
    DTexists.timewindow <- subset(DT, TOB < (timewindow[2]- agegroup[1]) & TOD > timewindow[1]
                      & XCoord==facilities.df$Longitude
                      & YCoord==facilities.df$Latitude)
  }

  # Case1
  # TOB > tw- and TOD < tw+ keep [age- <= TOD - TB <= age+]
  case1 <- subset(DTexists.timewindow, TOB > timewindow[1] & TOD < timewindow[2])
  case1 <- subset(case1, agegroup[1] <= (TOD-TOB))
  case1 <- case1 %>% mutate(MinAgeTW1 = 0, MaxAgeTW2 = TOD - TOB)

  # Case2
  # TOB > tw- and alive (TOD > tw+ or TOD=Inf) keep [age- <= tw+ - TB <= age+]
  case2 <- subset(DTexists.timewindow, TOB > timewindow[1] & (TOD < timewindow[2] | TOD==Inf))
  case2 <- subset(case2, agegroup[1] <= (timewindow[2] - TOB) & (timewindow[2] - TOB) <= agegroup[2])
  case2 <- case2 %>% mutate(MinAgeTW1 = 0, MaxAgeTW2 = timewindow[2] - TOB)

  # Case3
  # TOB > 0 and TOB < tw- and TOD < tw+ keep [age- <= TOD - TB <= age+]
  case3 <- subset(DTexists.timewindow, TOB > 0 & TOB < timewindow[1] & TOD < timewindow[2])
  case3 <- subset(case3, agegroup[1] <= (TOD - TOB) & (TOD - TOB) <= agegroup[2])
  case3 <- case3 %>% mutate(MinAgeTW1 = timewindow[1]-TOB, MaxAgeTW2 = TOD - TOB)

  # Case4
  # TOB > 0 and TOB < tw- and alive (TOD > tw+ or TOD=Inf) keep [age- <= tw+ - TB <= age+]
  case4 <- subset(DTexists.timewindow, TOB > 0 & TOB < timewindow[1] & (TOD > timewindow[2] | TOD ==Inf))
  case4 <- subset(case4, agegroup[1] <= (timewindow[2] - TOB) & (timewindow[2] - TOB) <= agegroup[2])
  case4 <- case4 %>% mutate(MinAgeTW1 = timewindow[1]-TOB, MaxAgeTW2 = timewindow[2] - TOB)

  # Case5
  # TOB < 0 and TOD < tw+ keep [age- <= TOD - TB <= age+]
  case5 <- subset(DTexists.timewindow, TOB < 0 & TOD < timewindow[2])
  case5 <- subset(case5, (timewindow[1] - TOB) <= agegroup[2])
  case5 <- case5 %>% mutate(MinAgeTW1 = timewindow[1]-TOB, MaxAgeTW2 = TOD - TOB)

  # Case6
  # TOB < 0 and alive (TOD > tw+ or TOD=Inf) keep [age- <= tw+ - TB <= age+]
  case6 <- subset(DTexists.timewindow, TOB < 0 & (TOD > timewindow[2] | TOD ==Inf))
  case6 <- subset(case6, (timewindow[1] - TOB) <= agegroup[2])
  case6 <- case6 %>% mutate(MinAgeTW1 = timewindow[1]-TOB, MaxAgeTW2 = timewindow[2] - TOB)

  DTexists.timewindow <- rbind(case1, case2, case3, case4, case5, case6)

  #check here (You should get TRUE)
  #(nrow(DTexists.timewindow) == length(unique(DTexists.timewindow$ID)))
  #Else dplyr::intersect(case1$ID, case5$ID) #check which case is a problem

  return(DTexists.timewindow)
}
