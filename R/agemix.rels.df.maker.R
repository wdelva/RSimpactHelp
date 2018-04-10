#' Create dataframe of relationships for agemixing analysis.
#'
#' Produces a dataframe of relationships from the dataframe of relationship
#' episodes created by \code{\link{agemix.episodes.df.maker}}. The user
#' specifies a time point, time window, and age group for which they would like
#' to obtain this dataframe. This dataframe can then be used as the input for
#' fitting linear mixed effect models that describe the agemixing pattern with
#' \code{\link{amp.modeller}}.
#'
#'
#' @param dataframe The dataframe that is produced by \code{\link{agemix.episodes.df.maker}}
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
#' @return returns a dataframe of relationships
#'
#' @examples
#' data(persreldf)
#' agemix.rels.df <- agemix.rels.df.maker(dataframe = persreldf, agegroup = c(15, 30),
#' timewindow = 1, timepoint = 30, start = FALSE)
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export


agemix.rels.df.maker <- function(dataframe,
                            agegroup,
                            timepoint,
                            timewindow,
                            start = FALSE) {
  #Warnings
  if (!is.data.frame(dataframe)) {
    stop("dataframe wrong type")
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
    # This only includes relationships that started
    # in the time window

    men <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(episodeorder == 1 &
               (FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage &
               Gender == "male" &
               TOD > time) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

    women <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(episodeorder == 1 &
               (FormTime <= time & FormTime >= window) &
               age >= lwrage &
               age < uprage &
               Gender == "female" &
               TOD > time) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

  } else {
    # This includes all relationships that were ongoing
    # at somepoint during the time window, but may have
    # started long before the time window.

    men <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage &
               Gender == "male" &
               TOD > time) %>%
      dplyr::distinct(ID, relid, .keep_all = TRUE) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)

    women <- dataframe %>%
      dplyr::mutate(age = time - TOB) %>%
      dplyr::filter(FormTime <= time &
               DisTime > window &
               age >= lwrage &
               age < uprage &
               Gender == "female" &
               TOD > time) %>%
      dplyr::distinct(ID, relid, .keep_all = TRUE) %>%
      dplyr::mutate(agerelform0 = agerelform - lwrage)
  }
  comb <- bind_rows(men, women)
  return(comb)
}
