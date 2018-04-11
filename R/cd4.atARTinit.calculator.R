#' CD4 at ART initiation
#'
#' Calculate the number of people who started ART below a given cd4 count
#' threshold, within a given time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper
#'   bound) that should be retained, e.g. c(15, 30).
#' @param timewindow The window in simulation time during which ART initiations are considered.
#' @param cd4count The cd4 cell count threshold below which ART initiations are considered.
#' @param site Select only the particular site from the study, if "All", use all
#'   sites.
#' @return A dataframe with the number of people who initiated ART below the
#'   threshold and the total number of people who initiated ART at some point in the time window.
#' @examples
#' data(datalist)
#' cd4count.atARTInit <- cd4.atARTinit.calculator(datalist = datalist, agegroup = c(15,50))
#' cd4count.atARTInit
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @export

cd4.atARTinit.calculator <- function(datalist = datalist,
                          agegroup = c(15, 30),
                          timewindow = c(10, 15),
                          cd4count = 350,
                          site = "All"){

  cd4.atARTinit <- age.group.time.window(datalist = datalist,
                                         agegroup = agegroup,
                                         timewindow = timewindow,
                                         site = site)
  # This function filters out people who belonged to the age group at some point
  # during the time window

  cd4.atARTinit <- dplyr::filter(cd4.atARTinit,
                                 TreatTime != Inf) # People who started ART

  raw.df <- data.frame(cd4.atARTinit)
  art.df <- dplyr::filter(datalist$ttable,
                          ID %in% cd4.atARTinit$ID,
                          TStart > timewindow[1],
                          TStart < timewindow[2])

  ##What if the person dropped out and came back again? We only want the first ART initiation
  art.df <- dplyr::distinct(art.df,
                            ID,
                            .keep_all = TRUE) %>%
    dplyr::mutate(ART.start.CD4 = CD4atARTstart < cd4count)
  # ART.start.CD4 indicates those who started ART when the CD4 count was below the threshold

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- dplyr::left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  #provide a summary of those that are on treatment and those that started below a threshold
  cd4count.atARTInit <- data.frame(dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                                    TotalInits = n(),
                                                    InitsBelowThreshold = sum(ART.start.CD4)))

  return(cd4count.atARTInit)
}
