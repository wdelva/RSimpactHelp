#' A function that returns the total number of people between two points in simulation
#' time for a particular age group and gender whose VL and CD4count are as indicated
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow alive people within this simulation time e.g timewindow = 30.
#' @param gender alive gender.
#' @param agegroup alive people within this agegroup.
#' @param cd4count threshold at the time of ART initiation.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of of people whose CD4 count at ART initiation is within a given threshold
#' @examples
#' cd4.atARTinit <- cd4.atARTinit(datalist = datalist, agegroup=c(15,30), timewindow=c(10,40), cd4count=c(50,350), gender="Male", site="All")

cd4.atARTinit <- function(datalist = datalist, agegroup = c(15,30),
                             timewindow = c(15,30), cd4count=c(50,350), gender = "Male", site="All"){
  gender.id <-1
  if(gender!="Male"){gender.id = 0}

  cd4.atARTinit <- age.group.time.window(datalist = datalist,
                                                   agegroup = agegroup, timewindow = timewindow, site="All")

  cd4.atARTinit <- subset(cd4.atARTinit, TreatTime !=Inf & Gender == gender.id) #HIV positive Gender individuals

  raw.df <- data.frame(cd4.atARTinit)
  art.df <- subset(datalist$ttable, ID %in% cd4.atARTinit$ID & TStart > timewindow[1] & TStart < timewindow[2])

  ##What if the person dropped out and come back again?
  art.df <- data.frame(dplyr::summarise(dplyr::group_by(art.df, ID), Gender = max(Gender),
                                       CD4atARTstart=min(CD4atARTstart)))

  ##LastCD4 <- which.min(CD4atARTstart)

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  #select those individuals who started their treatment when their CD4 count was between the given threshold
  cd4count.atARTInit <- subset(raw.df, cd4count[1] <= CD4atARTstart &  CD4atARTstart <= cd4count[2])

  cd4count.atARTInit <- nrow(cd4count.atARTInit)

  return(cd4count.atARTInit)
}

