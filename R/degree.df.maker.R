#' Degree distribution of persons in recent sexual relationships.
#'
#' Creates a dataframe for the degree distribution of persons specified from the
#' dataset produced by \code{\link{agemix.df.maker}}. This function allows users
#' to specify whether they want new relationships or ongoing relationships among
#' a specified age group, gender, or HIV status of persons in the dataset. 
#' 
#' @param df The dataframe that is produced by \code{\link{agemix.df.maker}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper
#'   bound) that should be retained. Should be expressed as a vector of two
#'   integers. e.g. c(15, 30)
#' @param hivstatus HIV status at the time of the survey. Options are: 0 = only
#'   HIV-negative, 1 = only HIV-positive, 2 = all
#' @param survey.time Numeric integer representing the time point of the
#'   cross-sectional glance at the data
#' @param window.width Numeric integer representing how long before the
#'   cross-sectional glance should you look at the data. E.g 1 year before the
#'   survey
#' @param gender.degree The gender or the persons in the data. c("male",
#'   "female")
#' @param only.new Logical indicator. If TRUE, only relationships that were
#'   newly started during window.width are counted (i.e. the individual was
#'   NEVER in a relationship with these partners before the start of the
#'   window). If FALSE, all relationships that were ongoing at some point during
#'   the window.width are counted.
#'   
#' @return A dataframe that has the columns "ID" and "Degree"
#' 
#' @examples
#' cfg <- list()
#' modeloutput <- RSimpactCyan::simpact.run(configParams = cfg, destDir = "temp")
#' dl <- readthedata(modeloutput)
#' rels.df <- agemix.df.maker(datalist = dl)
#' degree.df <- degree.df.maker(rels.df, agegroup = c(20, 30), hivstatus = 1)

#' @importFrom magrittr %>%
#' @importFrom stats aggregate
#' @import dplyr
#' @export

degree.df.maker <- function(df, agegroup = c(15, 30), hivstatus = 0,
                         survey.time = 40, window.width = 1, 
                         gender.degree = "female", only.new = TRUE){

  # The time the window starts
  time.start.window <- survey.time - window.width
  
  dfnew <- df %>%
    dplyr::group_by(Gender, relid) %>% # Since each episode is duplicated for each gender
    dplyr::mutate(RelFormTime = min(FormTime), # Var with the start of rel
                  RelDisTime = max(DisTime)) %>% # Var with end of rel
    dplyr::ungroup() %>%
    dplyr::distinct(Gender, relid, .keep_all = TRUE) %>% # Keep only first row per relationship
    dplyr::filter(TOD > survey.time) # Remove rels from people who died before survey
  
  # Removing relationships from either HIV pos or neg people
  if (hivstatus == 0) {
    dfnew <- dfnew %>%
      dplyr::filter(InfectTime > survey.time)
  } else if (hivstatus == 1) {
    dfnew <- dfnew %>%
      dplyr::filter(InfectTime <= survey.time)
  }

  if (only.new) {
    dfnew <- dfnew %>%
      dplyr::filter(RelFormTime >= time.start.window, # Rel starts after start of window
                    RelFormTime < survey.time, # Rel starts before the survey
                    RelDisTime > time.start.window, # Rel ends after start of window
                    Gender == gender.degree, # Keep specified gender only
                    (survey.time - TOB) >= agegroup[1], # Lower bound of age group
                    (survey.time - TOB) < agegroup[2]) # Upper bound of age group
  } else {
    dfnew <- dfnew %>%
      dplyr::filter(RelFormTime < survey.time, # Must start before the survey
                    RelDisTime > time.start.window, # Must end after the start of window
                    Gender == gender.degree, 
                    (survey.time - TOB) >= agegroup[1],
                    (survey.time - TOB) < agegroup[2])
  }

  # Now create table that has the degree for each person
  degree.df <- dfnew %>%
    dplyr::group_by(ID) %>% # By person
    dplyr::summarise(Degree = n()) # Count how many rows
  
  return(degree.df)
}


