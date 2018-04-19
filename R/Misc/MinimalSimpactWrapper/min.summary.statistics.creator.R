#Parameters  creation

#Simulation to compute 1. growth rate 2. Inci 3. Prev

#Set up simulation parameters and initiate simulation
#just the year when simulation started.
simulation.type <- "maxart"  #"simpact-cyn"
sim.start.full <- as.Date("1970-03-31")

simulation.name <- "maxart-int" #"maxart-noint"

if(simulation.name == "maxart-int"){
  maxart.starttime <- as.Date("2014-09-01")
  maxart.endtime <- as.Date("2017-08-31")
  final.type <- "YesStudy" #There is a study
}else{
  maxart.starttime <- as.Date("2027-09-01")
  maxart.endtime <- as.Date("2032-11-03")
  final.type <- "NoStudy" #There is no study
}
#year simulation start
sim.start <- as.numeric(substr(sim.start.full,1,4))

#Simulation-end and HIV seed year
sim.end.full <- as.Date("2018-07-30")
seed.hiv.date <- as.Date("1986-03-31")

#preserve study yrtime
study.starttime <- maxart.starttime
study.endtime <- maxart.endtime

#initial population
init.population.total <- 2000
women.frac <- 0.5253

#GrowthRate for all ages and gender
growth.rate <- read.table(
    text="X2015   X2016
  one.year   1.83    1.81", header=TRUE, stringsAsFactors = FALSE)

# Data was collected between July 2006 and March 2007. TimeTo: 2007-03
prev.2007 <- read.table(
      text="F.value    M.value  FM.value
  A15.49    31.1       19.7     25.9
  A50.150   11.7       17.9     14.2",
                              header=TRUE, stringsAsFactors = FALSE)

###################### target names ######################

## merge row names and col names
tar.name <- function(df, tar.type = "name"){
  apply(expand.grid(rownames(df), names(df), ".",tar.type), 1, paste0, collapse="" )
}

### get the real target values
tar.value <- function(df){
  return(as.numeric(unlist(df, use.names = FALSE)))
}

#Creating target names
target.variables <- c(tar.name(growth.rate, "growth.rate"),
                      tar.name(prev.2007, "prev.2007"))

#if you will need to calibrate SET THE target values correctly
target.values <- c(tar.value(growth.rate),
                   tar.value(prev.2007) )

