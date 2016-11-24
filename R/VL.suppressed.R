#' A function that returns the percentage of clients that are virally suppressed 6 or more months
#' after ART initiation, at a paticular point in the simulation time
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timepoint alive people within this simulation time e.g timewindow = 30.
#' @param gender alive gender.
#' @param vlcutoff viral load below this threshold e.g <200 copies/mL is defined as suppressed
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of of people whose VL at ART initiation are within a given threshold
#' @examples
#' vl.suppressed <- vl.suppressed(datalist = datalist, agegr, timewindow=c(15,30), viralload=c(3,4), gender="Male", site="All")

vl.suppressed <- function(datalist = datalist, timepoint = 40, vlcutoff = 200, gender = "Male", site="All"){

  gender.id <-1
  if(gender!="Male"){gender.id = 0}

  DTalive.infected <- datalist$ptable
  if(site=="All"){
    DTalive.infected <- subset(DTalive.infected, TOB <= timepoint & TOD > timepoint & InfectTime!=Inf & Gender == gender.id)
  }else{
    facilities.df <- read.cv(datalist$itable$facilities.geo.coords)
    facilities.df <- filter(facilities.df, Facility = site)
    DTalive.infected <- subset(DTalive.infected, TOB <= timepoint & TOD > timepoint & InfectTime!=Inf & Gender == gender.id
                      & XCoord==facilities.df$Longitude
                      & YCoord==facilities.df$Latitude)
  }

  vlcutoff <- log10(200)
  six.months <- timepoint - 0.5

  DTalive.infected <- DTalive.infected %>%  mutate(ARTLess6mnths = (TreatTime <= six.months))
  raw.df <- data.frame(DTalive.infected)

  #calculate the VL from time of treatment start to now
  #log10[VL(sp,new)] = k * log10[VL(sp)]  (where is the time factor)?

  #Check if any of these individuals dropped out of treatment during treatment time
  event.df <- subset(datalist$etable, p1ID %in% raw.df$ID & eventname %in% c("dropout","chronicstage", "aidsstage","finalaidsstage")
                     & eventtime <= six.months )

  #acute stage - 2-4weeks
  #chornic stage period when the virus is living
  #aids stage - immune syste is badly damaged


  #What if they dropped out? Is it after the time point or before? Did they start treatment again?
  art.df <- subset(datalist$ttable, ID %in% raw.df$ID & TStart <= six.months)

  #get

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))



  # Should we correct for dropout and restart [Simpact simply uses the Set-Point Viral Load]
  # log10SPVL is not affected by treatment

  vl.atARTinit <- nrow(vl.atARTinit)

  return(vl.atARTinit)
}

# C <- 1325.05
# k <- -0.49
# logVL <- seq(2,6,by=2)
# logTS <- log10(C)+k*logVL
# filteredPersons <- subset(checkinf,InfectType == 1 & TOD < Inf & TreatTime == Inf)
# survTime <- log10(filteredPersons$TOD - filteredPersons$InfectTime)
# SPVL <- filteredPersons$log10SPVL
# plot(SPVL, survTime, type = "p", col = "blue",xlim=c(2,6.5),ylim=c(0,2.5))
# lines(logVL,logTS,col="green")


















