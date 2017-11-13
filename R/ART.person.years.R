#' A function that returns the total number of person-years on ART provided between two
#' points in a simulation time age group and gender
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Individuals who have this agegroup within this timewindow
#' @param timewindow alive people within this simulation time e.g timewindow = c(15, 30).
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return a total number of people on ART aggregated by Gender.
#' @examples
#' data(datalist)
#' person.years.art <- art.person.years(datalist = datalist, agegroup=c(15, 30),
#' timewindow =c(15, 30), site="All")
#' person.years.art
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

art.person.years <- function(datalist = datalist,
                             agegroup =c(15, 30),
                             timewindow = c(15, 40),
                             site="All"){

  person.onart <- age.group.time.window(datalist = datalist,
                                        agegroup = agegroup,
                                        timewindow = timewindow, site="All")

  person.onart <- subset(person.onart, TreatTime !=Inf)

  raw.df <- data.frame(person.onart)

  #treatment episodes for the person.onart
  art.df <- subset(datalist$ttable, ID %in% raw.df$ID)
  raw.df.column.part <- raw.df[,c("ID","exposure.start", "exposure.end")]

  art.df <- dplyr::left_join(x = art.df, y = raw.df.column.part, by = c("ID"))

  #treatment episodes within the timewindow
  art.df <- subset(art.df, TStart < exposure.end & TEnd > exposure.start)


  art.df <- art.df %>% dplyr::mutate(exposure.start.art = pmax(exposure.start, TStart, na.rm = TRUE))
  art.df <- art.df %>% dplyr::mutate(exposure.end.art = pmin(exposure.end, TEnd, na.rm = TRUE))
  art.df <- art.df %>% dplyr::mutate(onARTYears = exposure.end.art - exposure.start.art)

  ##What if the person dropped out and come back again?
  art.df <- art.df %>%
    group_by(ID) %>%
    summarise(onARTYears=sum(onARTYears, na.rm = TRUE),
              ARTcoverage = TRUE ) %>%
    as.data.frame()

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- dplyr::left_join(x = raw.df, y = art.df, by = c("ID"))

  if (nrow(raw.df) > 0) {
    raw.df$Infected <- 1
    #Now we apply some dplyr function to get the sum of cases and population size per gender.
    raw.df$Gender <- as.factor(raw.df$Gender)

    ART.coverage.df <- raw.df  %>%
      group_by(Gender) %>%
      summarise(
               TotalExposed = n(),
               sum.cases = sum(Infected),
               sum.EverOnART = sum(ARTcoverage, na.rm = TRUE),
               sum.EverOnARTyears = sum(onARTYears, na.rm = TRUE))

    ART.coverage.df.all <- raw.df %>%
      summarise(
        Gender = NA,
        TotalExposed = n(),
        sum.cases = sum(Infected, na.rm = TRUE),
        sum.EverOnART = sum(ARTcoverage, na.rm = TRUE),
        sum.EverOnARTyears = sum(onARTYears, na.rm = TRUE))

    ART.coverage.df <- rbind(ART.coverage.df, ART.coverage.df.all)

  }else {
    ART.coverage.df <- as.data.frame(matrix(NA, 3, 5))
    names(ART.coverage.df) <- c("Gender", "TotalExposed", "sum.cases",
                                "sum.EverOnART", "sum.EverOnARTyears")
  }

  return(ART.coverage.df)
}
