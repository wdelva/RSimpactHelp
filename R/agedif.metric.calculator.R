#' Calculate different age-mixing metrics, overall and stratified.
#'
#' Calculate the median age difference at a specified time point,
#' for the supplied age group, stratfied strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the age-mixing metrics should be calculated.
#' @param timewindow The length of time before the timepoint for which relationships should be included,
#' e.g. 1, representing one year before the timepoint. This should be a whole number.
#' @param  start This is a logical indicating that only relationships starting after the beginning of the window
#' should be used. If start = FALSE relationships could start before the time window.
#'
#' @return a dataframe with median age difference and inter-quartile range(IQR)
#' for the specified time point and age group, overall, and stratified by gender
#'
#' @examples
#' agemix.df <- agedif.metric.calculator(datalist = datalist, agegroup = c(15, 30), timepoint = 30)


agedif.metric.calculator <- function(datalist,
                                  agegroup,
                                  timepoint,
                                  timewindow,
                                  start = FALSE) {

  #Constants
  time <- timepoint
  window <- timepoint - timewindow
  lwrage <- agegroup[1]
  uprage <- agegroup[2]

  # Create dataframe that merges all the necessary variables from different tables
  # Each row represents a relationship episode
  dfrmale <- datalist$rtable %>%
    data.frame() %>%
    rename(ID = IDm) %>%
    mutate(relid = paste0(ID, IDw))

  dfrfemale <- datalist$rtable %>%
    data.frame() %>%
    rename(ID = IDw) %>%
    mutate(relid = paste0(IDm, ID))

  dfmale <- datalist$ptable %>%
    data.frame() %>%
    select(ID, Gender, TOB) %>%
    filter(Gender == 0) %>%
    left_join(dfrmale, by = "ID") %>%
    select(-IDw)

  dffemale <- datalist$ptable %>%
    data.frame() %>%
    select(ID, Gender, TOB) %>%
    filter(Gender == 1) %>%
    left_join(dfrfemale, by = "ID") %>%
    select(-IDm)

  df <- bind_rows(dfmale, dffemale)

  # Create subset df of relationships that meet input criteria
  # Also remove duplicated relationships (from multiple episodes) by gender only
  if(start == TRUE) {

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
              iqr = IQR(AgeGap),
              Q1 = median - (iqr / 2),
              Q3 = median + (iqr / 2)) %>%
    select(-iqr)

  # Calculate median age difference overall
  adoverall <- subdf2 %>%
    summarise(n = n(),
              median = median(AgeGap),
              iqr = IQR(AgeGap),
              Q1 = median - (iqr / 2),
              Q3 = median + (iqr / 2)) %>%
    select(-iqr)


  #Combine tables
  adtable <- bind_rows(adbygender, adoverall)
  return(adtable)

}
