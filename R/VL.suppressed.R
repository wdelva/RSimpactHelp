#' A function that returns the total number of people between two points in simulation
#' time for a particular age group and gender whose VL and CD4count are as indicated
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow alive people within this simulation time e.g timewindow = 30.
#' @param gender alive gender.
#' @param agegroup alive people within this agegroup.
#' @param viralload log10 set point viral load threshold at the time of ART initiation.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of of people whose VL at ART initiation are within a given threshold
#' @examples
#' vl.atARTinit <- vl.atARTinit(datalist = datalist, agegroup=c(15,30), timewindow=c(15,30), viralload=c(3,4), gender="Male", site="All")

vl.suppressed <- function(datalist = datalist, agegroup = c(15,30),
                             timewindow = c(15,30), vlcutoff = 2, gender = "Male", site="All"){
  gender.id <-1
  if(gender!="Male"){gender.id = 0}

  vl.atARTinit- age.group.time.window(datalist = datalist,
                                                   agegroup = agegroup, timewindow = timewindow, site="All")

  vl.atARTinit <- subset(vl.atARTinit, TreatTime !=Inf & Gender == gender.id) #HIV positive Gender individuals


  vl.atARTinit <- subset(raw.df, viralload[1] <= log10SPVL &  log10SPVL <= viralload[2])

  # Should we correct for dropout and restart [Simpact simply uses the Set-Point Viral Load]
  # log10SPVL is not affected by treatment

  vl.atARTinit <- nrow(vl.atARTinit)

  return(vl.atARTinit)
}

