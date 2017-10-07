#I presume this will be called through source()
#Idea is to give date and statistics. Calculated using start of the simulation - that time.

#Parameters  creation

#Set up simulation parameters and initiate simulation
sim.start <- 1970 #just the year when simulation started.
sim.start.full <- as.Date("1970-03-31")
maxart.starttime <- as.Date("2014-09-01")
maxart.endtime <- as.Date("2017-08-31")
sim.end.full <- as.Date("2019-03-31")
seed.hiv.date <- as.Date("1986-03-31")

#initial population
init.population.total <- 3000
women.frac <- 0.5253

#GrowthRate for all ages and gender
growth.rate <- read.table(
    text="X2015   X2016
  value   1.83    1.81", header=TRUE, stringsAsFactors = FALSE)

growth.rate.tar.names <- paste(names(growth.rate),"growth.rate", sep = ".")

############# 2011 Incidence ages specific: SHIMS 1 #############################
#consider for validation:  [all       3.14   1.65   2.38  Date range(2010-12    2011-06)
inci.2011 <- read.table(
      text="F.value    M.value
  A25.29    2.25       2.36
  A45.49    1.20       0.02", header=TRUE, stringsAsFactors = FALSE)

# Data was collected between July 2006 and March 2007. TimeTo: 2007-03
prev.2007 <- read.table(
      text="F.value    M.value  FM.value
  A15.49    31.1       19.7     25.9
  A50.150   11.7       17.9     14.2",
                              header=TRUE, stringsAsFactors = FALSE)

###################### target names

tar.name <- function(df, tar.type = "name"){
  apply(expand.grid(rownames(df), names(df), ".",tar.type), 1, paste0,collapse="" )
}

#Creating target names
target.variables <- c(growth.rate.tar.names, tar.name(inci.2011, "inci.2011"),
                      tar.name(prev.2007, "prev.2007"))

#Testing
#target.variables <- c(max.art.retention.tar.names,"node.id")


