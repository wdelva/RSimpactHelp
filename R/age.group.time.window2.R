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
  DT$pfacility <- "NA"
  pf.index <- which(colnames(DT)=="pfacility")

  if(site=="All"){
    DTexists.timewindow <- DT
  }else{
    facilities.df <- datalist.test$ftable
    colnames(facilities.df) <- c("facility.xy","XCoord","YCoord")

    for(i in nrow(DT)){
      fc.id <- which.min(sqrt((DT[i,XCoord] - facilities.df$XCoord)^2 + (DT[i,YCoord] - facilities.df$YCoord)^2))
      DT[i, pf.index] <- facilities.df[fc.id, facility.xy]
    }

    DTexists.timewindow <- subset(DT, pfacility==site)
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


  ##time.of.HIV.infection <- datalist$ptable$InfectTime
  #exposure.end <- pmin(time.of.HIV.infection, pmin(time.of.upperbound.agegroup, time.of.upperbound.timewind))

  #Exposure time (Everyone with exposure time greater than 0)
  DTexists.timewindow <- DTexists.timewindow %>% mutate(exposure.time = exposure.end - exposure.start)
  DTexists.timewindow <- DTexists.timewindow %>% mutate(real.exposure.time = exposure.time > 0)
  DTexists.timewindow <- subset(DTexists.timewindow, real.exposure.time == TRUE )


  #exposure.time <- exposure.end - exposure.start
  #real.exposure.time <- exposure.time > 0 # We create a vector to see who REALLY had exposure time
  #DTexists.timewindow$exposure.time[real.exposure.time == FALSE] <- 0

  # Now we check who of the people with the real exposure time had the event

  # Their InfectTime must be after their exposure.time started, and before or at exposure.end
  ##infection.after.exposure.start <- DTexists.timewindow$InfectTime > exposure.start
  ##infection.before.or.at.exposure.end <- DTexists.timewindow$InfectTime <= exposure.end
  ##infection.in.timewindow <- infection.after.exposure.start & infection.before.or.at.exposure.end


  ##DTexists.timewindow$incident.case <- infection.in.timewindow
  ##DTexists.timewindow$exposure.times <- exposure.time

  # For code testing: we create a dataset with only the people who contribute exposure time:
  # You must be HIV negative at the time that you qualify the minumum for age range and time range.
  ##DTexists.timewindow <- filter(DTexists.timewindow,
  ##                              exposure.times > 0)



  return(DTexists.timewindow)
}
