#' A function that returns the total number of people between two points in simulation
#' time for a particular age group and gender whose VL and CD4count are as indicated
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow alive people within this simulation time e.g timewindow = 30.
#' @param agegroup alive people within this agegroup.
#' @param viralload at the time of ART initiation.
#' @param cd4count threshold at the time of ART initiation.
#' @param site Select only the particular site from the study, if all ignore site/use all sites.
#' @return the total number of of people whose VL and CD4 count at ART initiation are given
#' @examples
#' cd4.vl.atARTinit <- cd4.vl.atARTinit(datalist = datalist.test, agegroup=c(15,40), timewindow=c(15,40), viralload=c(3,4), cd4count=c(350,500), site="All")

cd4.vl.atARTinit <- function(datalist = datalist, agegroup = c(15,30),
                             timewindow = c(15,30), viralload=c(3,4), cd4count=c(350,500), site="All"){

  cd4.vl.atARTinit <- age.group.time.window(datalist = datalist,
                                                   agegroup = agegroup, timewindow = timewindow, site="All")

  cd4.vl.atARTinit <- subset(cd4.vl.atARTinit, TreatTime !=Inf) #HIV positive individuals

  raw.df <- data.frame(cd4.vl.atARTinit)
  art.df <- subset(datalist$ttable, ID %in% cd4.vl.atARTinit$ID & TStart > timewindow[1] & TStart < timewindow[2])

  # Now we apply the left_join dplyr function to add the ART status to raw.df.
  raw.df <- left_join(x = raw.df, y = art.df, by = c("ID", "Gender"))

  #select those individuals who started their treatment when their CD4 count was between the given threshold
  raw.df <- raw.df%>% mutate(cd4count.atARTInit = (cd4count[1] <= CD4atARTstart &  CD4atARTstart <= cd4count[2]))

  #select those individuals who started their treatment when their vl was between the given threshold
  raw.df <- raw.df%>% mutate(vl.atARTinit = (viralload[1] <= log10SPVL &  log10SPVL <= viralload[2]))

  #provide a summary of those that are on treatment and those that started below a threshold
  cd4.vl.atARTinit <- data.frame(dplyr::summarise(dplyr::group_by(raw.df, Gender),
                                                    TotalCases = n(),
                                                    cd4count.atARTInit = sum(cd4count.atARTInit),
                                                    vl.atARTinit = sum(vl.atARTinit)))

  cd4.vl.atARTinit$Gender[cd4.vl.atARTinit$Gender==0] <- "Woman"
  cd4.vl.atARTinit$Gender[cd4.vl.atARTinit$Gender==1] <- "Man"


  #cd4.vl.atARTinit <- c(cd4count.atARTInit, vl.atARTinit)

  return(cd4.vl.atARTinit)
}

