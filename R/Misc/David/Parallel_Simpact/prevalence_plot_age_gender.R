library(tidyr)
library(ggplot2)

prevalence.df <- prevalence.calculator(datalist = datalist,
                                       agegroup = c(15, 30),timepoint = 30)

prevalence.df <-prevalence.plotter(datalist = datalist, agegroup = c(15, 50))


#Swaziland SHIMS 2011 incidence by gender
xmen <- c(0.8,6.6,21.3,36.6,47.0,45.5,42.5)
xwom <- c(14.3, 31.5, 46.7, 53.8, 49.1, 39.7, 31.6)
time <- 1:7

data.plot <- data.frame(cbind(time, xwom,xmen))
names(data.plot) <- c("time", "Women", "Men")
## wide to long format
data.plot <- data.plot %>%
  tidyr::gather(Gender, "Incidence", 2:3)

ggplot(data.plot,
       aes(x=time, y=Incidence, group=Gender, label = Incidence)) +
  geom_point() + geom_line() + aes(colour = Gender)  +
  xlim("18-19","20-24","25-29","30-34","35-39","40-44","45-49") +
  ylim(0,60) +
  geom_text(aes(label=Incidence), size=4, vjust=-.8, show.legend = FALSE) +
  labs(title = "2011 HIV Prevalence by Age Group in Swaziland",
       x = "Age Group", y = " HIV Prevalence") +  theme_bw() +
  theme(axis.text.x  = element_text(vjust=0.5, size=14),
        axis.title.x = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=14),
        axis.title.y = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5, size=16))



### Overall prevalence plotter


prevalence.plotter <- function(datalist = datalist,
                               agegroup = c(15, 50)){
  prev.agegroup <- agegroup
  timevect <- 0:datalist$itable$population.simtime[1]

  prevplot.df <- data.frame(timevect = timevect, prev = NA, ll = NA, ul = NA)
  rowindex <- 1
  for (timepoint in timevect){
    prev.tibble <- prevalence.calculator(datalist = datalist,
                                         agegroup = prev.agegroup, timepoint = timepoint)

    prevplot.df$prev[rowindex] <- as.numeric(prev.tibble$pointprevalence[3])
    prevplot.df$ll[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ll[3])
    prevplot.df$ul[rowindex] <- as.numeric(prev.tibble$pointprevalence.95.ul[3])
    rowindex <- rowindex + 1
  }
  prevplot.df <- prevplot.df %>% tidyr::gather(type, prevalence, prev:ul)
  spectrum.estimates <- data.frame(timevect = 3+(0:30), #1980:2010
                                   prevalence = 0.01 * c(0.01, 0.03, 0.07, 0.12, 0.19,
                                                         0.30, 0.44, 0.65, 0.93, 1.34, 1.91,
                                                         2.69, 3.78, 5.23, 7.11, 9.41, 12.00, 14.70,
                                                         17.32, 19.66, 21.61,23.12, 24.24, 25.01,
                                                         25.50, 25.80, 25.91, 25.86, 25.82, 25.80, 25.81),
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


