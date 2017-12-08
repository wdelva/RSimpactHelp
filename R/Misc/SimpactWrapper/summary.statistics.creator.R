#I presume this will be called through source()
#Idea is to give date and statistics. Calculated using start of the simulation - that time.

#Parameters  creation

#Set up simulation parameters and initiate simulation
#just the year when simulation started.
sim.start.full <- as.Date("1970-03-31")
sim.start <- as.numeric(substr(sim.start.full,1,4))
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
  one.year   1.83    1.81", header=TRUE, stringsAsFactors = FALSE)

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

## merge row names and col names
tar.name <- function(df, tar.type = "name"){
  apply(expand.grid(rownames(df), names(df), ".",tar.type), 1, paste0,collapse="" )
}

### get the real target values
tar.value <- function(df){
  return(as.numeric(unlist(df, use.names = FALSE)))
}

tar.value(prev.2007)

#Creating target names
target.variables <- c(tar.name(growth.rate, "growth.rate"), tar.name(inci.2011, "inci.2011"),
                      tar.name(prev.2007, "prev.2007"))

#if you will need to calibrate SET THE target values correctly
target.values <- c(tar.value(growth.rate), tar.value(inci.2011), tar.value(prev.2007) )

