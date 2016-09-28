#' Age difference median.
#'
#' Calculate the median and inter-quartile range (IQR) of age differences for
#' the population simulated in Simpact. The user specifies a time point, time
#' window, and age group for which they would like to obtain the median age
#' difference.
#'
#' The user can optionally choose whether or not they would like to examine only
#' relationships that began during their specified time window.
#'
#' The table that is returned contains the median stratified by gender, as well
#' as the overall median in the population. The user should be aware that the
#' number of relationships for men and women does not sum to the same number as
#' the overall number. This is because Simpact produces a complete closed
#' population and therefore some of the relationships used to calculate the
#' median for men, may also be included in the calculation for women. The
#' overall calculation, examines only unique relationships in the entire
#' population.
#'
#' @param df The dataframe that is produced by \code{\link{agemix.df.maker()}}
#' @param agegroup Boundaries of the age group that should be retained, e.g.
#'   c(15, 30). The interval is closed on the left and open on the right.
#' @param timepoint Point in time during the simulation to be used in the
#'   calculation.
#' @param timewindow The length of time before the timepoint for which
#'   relationships should be included, e.g. 1, representing one year before the
#'   timepoint. This should be a whole number.
#' @param start This is a logical indicating that only relationships starting
#'   after the beginning of the window should be used. If start = FALSE
#'   relationships could start before the time window. This is the default.
#'
#' @return A dataframe of class "tbl_df" "tbl" "data.frame"
#'
#' @examples
#' # In our  population we want to look at the median age difference for men
#' # and women 40 years into the simulation. We want to pretend we are
#' # replicating a survey where participants are included only if they are
#' # between the ages of 18 and 49. Furthermore, the are asked only about
#' # the relationships that happened in the previous 2 years.
#'
#' load(dataframe)
#' agedifmedtable <- agedif.median(df = dataframe, agegroup = c(18, 50), timewindow = 2, timepoint = 40)
#'
#' @importFrom magrittr %>%
#' @import dplyr


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

  if (start == TRUE) {

    subdf <- df %>%
      mutate(age = time - TOB) %>%
      filter((FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage &
               TOD > time) %>%
      distinct(Gender, relid, .keep_all = TRUE)

  } else {

    subdf <- df %>%
      mutate(age = time - TOB) %>%
      filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage &
               TOD > time) %>%
      distinct(Gender, relid, .keep_all = TRUE)
  }

  #Remove duplicated relationship episodes overall
  subdf2 <- subdf %>%
    distinct(relid, .keep_all = TRUE)

  # Calculate median age difference by gender
  adbygender <- subdf2 %>%
    group_by(Gender) %>%
    dplyr::summarise(n = n(),
              median = median(AgeGap),
              Q1 = as.numeric(summary(AgeGap)["1st Qu."]),
              Q3 = as.numeric(summary(AgeGap)["3rd Qu."]))

  # Calculate median age difference overall
  adoverall <- subdf2 %>%
    dplyr::summarise(n = n(),
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
