#' A function that returns the total number of person-years on ART provided between two points in a simulation time
#' age group and gender
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Individuals who have this agegroup within this timewindow
#' @param gender gender.
#' @param timewindow Point in time at which the subset should be created and HIV status should be evaluated.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a total number of people on ART for the given gender.
#' @examples
#' person.years.art <- person.years(datalist = datalist, agegroup=c(15,30), gender="Male", timewindow =c(15,30), site="All")

person.years <- function(datalist = datalist, agegroup =c(15,30), gender="Male", timewindow = c(30,30), site="All"){
  
  ## get the correct population within this agegroup and timewindow
  gender.id <-1
  if(gender!="Male"){gender.id = 0}
  
  person.onart <- age.group.time.window(datalist = datalist, 
                                                   agegroup = agegroup, timewindow = timewindow, site="All")

  
  person.onart <- subset(person.onart, TreatTime !=Inf)
  
  raw.df <- data.frame(person.onart)
  art.df <- subset(datalist$ttable, ID %in% person.onart$ID & TStart > timewindow[1] & TStart < timewindow[2])

  art.df$onARTYears <- NA
  art.df$onARTYears[art.df$TEnd =="Inf"] <- timewindow[2] - art.df$TStart[art.df$TEnd =="Inf"]
  art.df$onARTYears[art.df$TEnd !="Inf"] <- art.df$TEnd[art.df$TEnd !="Inf"] - art.df$TStart[art.df$TEnd !="Inf"]

  ##What if the person dropped out and come back again?
  art.df <- data.frame(dplyr::summarise(dplyr::group_by(art.df, ID), Gender = max(Gender),
                                        TStart = mean(TStart), TEnd=mean(TEnd), DiedNow = max(DiedNow), #placeholders
                                        onARTYears=sum(onARTYears)))

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  if (nrow(raw.df) > 0) {
    raw.df$ARTcoverage <- !is.na(raw.df$TStart)
    raw.df$Infected <- 1
    #Now we apply some dplyr function to get the sum of cases and population size per gender.
    raw.df$Gender <- as.character(raw.df$Gender)

    ART.coverage.df <- dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                  sum.cases = sum(Infected),
                                  sum.EverOnART = sum(ARTcoverage, na.rm = TRUE),
                                  sum.EveronARTyears = sum(onARTYears, na.rm = TRUE)
    )
  } else {
    ART.coverage.df <- data.frame(Gender = NA, sum.cases = NA,
                                  sum.EverOnART = NA, sum.EveronARTyears = NA)
  }
  return(ART.coverage.df)
}
