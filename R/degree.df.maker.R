#' Degree distribution of recent sexual relationships.
#'
#' Creates a table of two columns, "ID" and "degree" from the df dataframe,
#' with specific degrees of new relationships or ongoing relationships among women within
#' specified age group. The degree distribution shows the total number of relationships that
#' were ongoing (if only.new == FALSE) at some time in the survey window
#' or that were newly formed (if only.new == TRUE) during the survey window.
#'
#' @param dataframe.df The dataframe that is produced by \code{\link{agemix.df.maker}}
#' @param survey.time Time point of the cross-sectional survey.
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that
#' should be retained, e.g. c(15, 30)
#' @param hivstatus HIV status at the time of the survey. Options are 2, means all; 0
#' means only HIV-negative, 1 means only HIV-positive.
#' @param window.width Time period before the survey e.g 1 year before the survey.
#' @param gender.degree The gender of the results
#' @param only.new Logical indicator. If TRUE, only relationships that were newly started
#' during window.width are counted
#' (i.e. the individual was NEVER in a relationship with these partners before
#'  the start of the window).
#' If FALSE, all relationships that were ongoing at some point during the
#' window.width are counted.
#'
#' @return A dataframe, degree.df that has the columns "ID", "degree", and a fraction of women
#' with >1 sexual partner in the last year.
#@examples
#data(dataframe.df)
#degree.df <- degree.df.maker(dataframe.df, agegroup = c(15, 30), hivstatus = 0,
#survey.time = 10, window.width = 1, only.new = TRUE)

#' @importFrom magrittr %>%
#' @importFrom stats aggregate
#' @import dplyr
#' @export

degree.df.maker <- function(df, agegroup = c(15, 30), hivstatus = 0,
                         survey.time = 40, window.width = 1, 
                         gender.degree = "female", only.new = TRUE)
{
  time.start.window <- survey.time - window.width
  
  reldf <- df %>%
    dplyr::group_by(Gender, relid) %>%
    dplyr::mutate(RelFormTime = min(FormTime),
                  RelDisTime = max(DisTime)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gender, relid, .keep_all = TRUE) %>%
    dplyr::filter(TOD > survey.time)
  
  {
    if (hivstatus == 0) {
      dfnew <- dfnew %>%
        dplyr::filter(InfectTime > survey.time)
    }
    else if (hivstatus == 1) {
      dfnew <- dfnew %>%
        dplyr::filter(InfectTime <= survey.time)
    }
    }
  
  {
    if (only.new) {
      dfnew <- dfnew %>%
        dplyr::filter(RelFormTime >= time.start.window, 
                      RelFormTime < survey.time, 
                      RelDisTime > time.start.window, 
                      Gender == gender.degree,
                      (survey.time - TOB) >= agegroup[1], 
                      (survey.time - TOB) < agegroup[2])
    }
    else {
      dfnew <- dfnew %>%
        dplyr::filter(RelFormTime < survey.time,
                      RelDisTime > time.start.window, 
                      Gender == gender.degree, 
                      (survey.time - TOB) >= agegroup[1],
                      (survey.time - TOB) < agegroup[2])
    }
  }
  
  degree.df <- dfnew %>%
    group_by(ID) %>%
    summarise(Degree = n())
  
  return(degree.df)
}


