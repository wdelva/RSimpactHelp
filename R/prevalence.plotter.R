#' Plot the overall HIV prevalence for the duration of the simulation.
#'
#' Plot the overall HIV prevalence, with confidence interval, for the duration of the simulation, for a specific age group.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 50)
#' @return a ggplot2 object
#' @examples
#' prevalence.df <- prevalence.calculator(datalist = datalist, agegroup = c(15, 30), timepoint = 30)
#'
#' @import gtools
#' @import ggplot2

prevalence.plotter <- function(datalist = datalist,
                                  agegroup = c(15, 50)){
prev.agegroup <- agegroup
timevect <- 0:datalist.test$itable$population.simtime

prevplot.df <- data.frame(timevect = timevect, prev = NA, ll = NA, ul = NA)
rowindex <- 1
for (timepoint in timevect){
  prev.tibble <- prevalence.calculator(datalist = datalist.test, agegroup = prev.agegroup, timepoint = timepoint)
  prevplot.df$prev[rowindex] <- as.numeric(prev.tibble$pointprevalence[3])
  prevplot.df$ll[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ll[3])
  prevplot.df$ul[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ul[3])
  rowindex <- rowindex + 1
}
prevplot.df <- prevplot.df %>% tidyr::gather(type, prevalence, prev:ul)
spectrum.estimates <- data.frame(timevect = 3+(0:30), #1980:2010
                                 prevalence = 0.01 * c(0.01, 0.03, 0.07, 0.12, 0.19, 0.30, 0.44, 0.65, 0.93, 1.34, 1.91,
                                              2.69, 3.78, 5.23, 7.11, 9.41, 12.00, 14.70, 17.32, 19.66, 21.61,
                                              23.12, 24.24, 25.01, 25.50, 25.80, 25.91, 25.86, 25.82, 25.80, 25.81),
                                 type = "spectrum.output")

prevplot.df <- gtools::smartbind(prevplot.df, spectrum.estimates)

prevplot.df$pointestimate <- prevplot.df$type=="prev"
prevplot.df$spectrum.output <- prevplot.df$type=="spectrum.output"

p <- ggplot(prevplot.df, aes(x=timevect+1977, y=prevalence, group=type)) +
  geom_line(aes(linetype = !(pointestimate | spectrum.output),
                colour = spectrum.output)) +
  xlab("Simulation time") +
  ylab(" HIV prevalence") +
  #facet_wrap(~Gender, nrow=1)
  theme_bw() +
  theme(legend.position = "none")
return(p)
}
