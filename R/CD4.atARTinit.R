#' A function that returns the total number of people who started their treatment when
#' their CD4count was below a specified threshold segregated by gender
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @param timewindow alive people within this simulation time e.g timewindow = c(20,30).
#' @param agegroup alive people within this agegroup.
#' @param cd4count threshold at the time of ART initiation.
#' @param site select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of people whose CD4 count at ART initiation is within a was below this threshold
#' @examples
#' data(datalist)
#' cd4.atARTinit <- cd4.atARTinit(datalist = datalist, agegroup=c(15, 30),
#'  timewindow=c(10, 40), cd4count=350, site="All")
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

cd4.atARTinit <- function(datalist = datalist,
                          agegroup = c(15, 30),
                          timewindow = c(15, 30),
                          cd4count = 350, site="All"){

  cd4.atARTinit <- age.group.time.window(datalist = datalist,
                                        agegroup = agegroup,
                                         timewindow = timewindow, site="All")

  cd4.atARTinit <- subset(cd4.atARTinit, TreatTime !=Inf) #HIV positive individuals

  raw.df <- data.frame(cd4.atARTinit)
  art.df <- subset(datalist$ttable, ID %in% cd4.atARTinit$ID &
                     TStart > timewindow[1] & TStart < timewindow[2])

  ##What if the person dropped out and come back again?
  art.df <- art.df %>%
    group_by(ID, Gender) %>%
    summarise(CD4atARTstart = min(CD4atARTstart, na.rm = TRUE)) %>%
    as.data.frame()

  #indicate those who started their treatment when their CD4 count was below a given threshold
  art.df <- art.df %>% dplyr::mutate(ART.start.CD4 = CD4atARTstart < cd4count)

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- dplyr::left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  #provide a summary of those that are on treatment and those that started below a threshold
  cd4count.atARTInit <- raw.df %>%
    group_by(Gender) %>%
    summarise(TotalCases = n(),
              LessCD4initThreshold = sum(ART.start.CD4, na.rm = TRUE)) %>%
    as.data.frame()

  return(cd4count.atARTInit)
}

