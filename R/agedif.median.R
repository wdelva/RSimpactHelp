#' Age difference median
#'
#' Calculate the median age difference at a specified time point,
#' for the supplied age group, stratfied by gender
#'
#' @param df The dataframe that is produced by \code{\link{agedif.df.maker()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the age-mixing metrics should be calculated.
#' @param timewindow The length of time before the timepoint for which relationships should be included,
#' e.g. 1, representing one year before the timepoint. This should be a whole number.
#' @param start This is a logical indicating that only relationships starting after the beginning of the window
#' should be used. If start = FALSE relationships could start before the time window.
#'
#' @return a dataframe with median age difference and inter-quartile range(IQR)
#' for the specified time point and age group, overall, and stratified by gender
#'
#' @examples
#' agedifmedtable <- agedif.median(df = dataframe, agegroup = c(15, 30), timewindow = 1, timepoint = 30)

# dplyr, magrittr

agedif.median <- function(df,
                          agegroup,
                          timepoint,
                          timewindow,
                          start = FALSE) {

  #Warnings
  if (!is.data.frame(df)) {
    stop("Dataframe wrong type")
  }

  if (length(agegroup) != 2) {
    stop("Need an upper and lower bound")
  }

  if (timepoint < 1) {
    stop("Time must be at least 1 year")
  }

  if (timewindow < 0) {
    stop("Window must be a whole number")
  }

  #Constants
  time <- timepoint
  window <- timepoint - timewindow
  lwrage <- agegroup[1]
  uprage <- agegroup[2]

  # Create subset df of relationships that meet input criteria
  # Also remove duplicated relationships (from multiple episodes) by gender only
  # The same relationship will still be represented twice —once by the female
  # partner and once by the male partner.
  # However, there will only be one episode for each.

  if (start == TRUE) {

    subdf <- df %>%
      mutate(age = TOB + time) %>%
      filter((FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage) %>%
      distinct(Gender, relid)

  } else {

    subdf <- df %>%
      mutate(age = TOB + time) %>%
      filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage) %>%
      distinct(Gender, relid)
  }

  #Remove duplicated relationship episodes overall
  subdf2 <- subdf %>%
    distinct(relid)

  # Calculate median age difference by gender
  adbygender <- subdf %>%
    group_by(Gender) %>%
    summarise(n = n(),
              median = median(AgeGap),
              Q1 = as.numeric(summary(AgeGap)["1st Qu."]),
              Q3 = as.numeric(summary(AgeGap)["3rd Qu."]))

  # Calculate median age difference overall
  adoverall <- subdf2 %>%
    summarise(n = n(),
              median = median(AgeGap),
              Q1 = as.numeric(summary(AgeGap)["1st Qu."]),
              Q3 = as.numeric(summary(AgeGap)["3rd Qu."]))


  #Combine tables
  adtable <- bind_rows(adbygender, adoverall) %>%
    mutate(Gender = as.factor(ifelse(is.na(Gender),
                                     "Overall",
                                     ifelse(Gender == "male",
                                            "Men",
                                            "Women"))))

  return(adtable)
}
