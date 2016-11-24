#' A function that returns the total number of person-years on ART provided between two points in a simulation time
#' age group and gender
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Individuals who have this agegroup within this timewindow
#' @param timewindow alive people within this simulation time e.g timewindow = c(15, 30).
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a total number of people on ART aggregated by Gender.
#' @examples
#' person.years.art <- art.person.years(datalist = datalist, agegroup=c(15, 30), timewindow =c(15, 30), site="All")
#'
#' @importFrom magrittr %>%
#' @import dplyr

art.person.years <- function(datalist = datalist, agegroup =c(15, 30), timewindow = c(15, 40), site="All"){

  person.onart <- age.group.time.window(datalist = datalist,
                                                   agegroup = agegroup, timewindow = timewindow, site="All")

  person.onart <- subset(person.onart, TreatTime !=Inf)

  raw.df <- data.frame(person.onart)
  art.df <- subset(datalist$ttable, ID %in% raw.df$ID) #treatment episodes for the person.onart
  raw.df.column.part <- raw.df[,c("ID","exposure.start", "exposure.end")]

  art.df <- left_join(x = art.df, y = raw.df.column.part, by = c("ID"))
  art.df <- subset(art.df, TStart < exposure.end & TEnd > exposure.start) #treatment episodes within the timewindow


  art.df <- art.df %>% mutate(exposure.start.art = pmax(exposure.start, TStart))
  art.df <- art.df %>% mutate(exposure.end.art = pmin(exposure.end, TEnd))
  art.df <- art.df %>% mutate(onARTYears = exposure.end.art - exposure.start.art)

  ##What if the person dropped out and come back again?
  art.df <- data.frame(dplyr::summarise(dplyr::group_by(art.df, ID), onARTYears=sum(onARTYears), ARTcoverage = TRUE ))

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- left_join(x = raw.df, y = art.df, by = c("ID"))

  if (nrow(raw.df) > 0) {
    raw.df$Infected <- 1
    #Now we apply some dplyr function to get the sum of cases and population size per gender.
    raw.df$Gender <- as.character(raw.df$Gender)

    ART.coverage.df <-data.frame( dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                                   TotalExposed = n(),
                                                   sum.cases = sum(Infected),
                                                   sum.EverOnART = sum(ARTcoverage, na.rm = TRUE),
                                                   sum.EverOnARTyears = sum(onARTYears, na.rm = TRUE)
    ))
  } else {
    ART.coverage.df <- data.frame(Gender = c(0, 1), sum.cases = c(NA, NA),
                                  sum.EverOnART = c(NA, NA), sum.EverOnARTyears = c(NA,NA))
  }

  ART.coverage.df$Gender[ART.coverage.df$Gender==0] <- "Woman"
  ART.coverage.df$Gender[ART.coverage.df$Gender==1] <- "Man"

  return(ART.coverage.df)
}
