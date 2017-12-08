#Parameters  creation

#Set up simulation parameters and initiate simulation
#just the year when simulation started.
simulation.type <- "maxart"  #"simpact-cyn"
sim.start.full <- as.Date("1970-03-31")
sim.start <- as.numeric(substr(sim.start.full,1,4))
maxart.starttime <- as.Date("2014-09-01")
maxart.endtime <- as.Date("2016-11-03")
sim.end.full <- as.Date("2019-03-31")
seed.hiv.date <- as.Date("1986-03-31")

#initial population
init.population.total <- 2000
women.frac <- 0.5253

################### ART initiation ##############################
#removing this, need pfacility, TOD for eligibility from interpolation

max.art.initiated.all <- read.table(
         text="AllClients
   all              84.6", header=TRUE, stringsAsFactors = FALSE)
#   before.eligible  2.2
#   within2.weeks    64
#   among.eligible   90
#   within12.months  91
#   within6months    87
#TimeTo: 2016-11

##########################  ART Retention  #############################
# expressed as rate totalRetained/TotalInitiated * 100
# within the study time period (study.start - study.end)

max.art.retention.all <- read.table(
  text="AllClients
  all             73
  at6.months      87
  at12.months     79", header=TRUE, stringsAsFactors = FALSE)
#timeTo: 2016-11

###########################  ViralLoadSuppression ##################
# expressed as rate i.e 100 - suppressed/TotalInitiated * 100
# within the study time period (study.start - study.end)
max.vl.none.suppression.all <- read.table(
  text="AllClients
  at6.mths.ns           6
  at12.mths.ns          11
  at12.mths.s.shim2   72.1", header=TRUE, stringsAsFactors = FALSE)

#Check validation
#atleast6.months  8       2016-11

#############  Mortality  ###########################
# expressed as rate Number of death/ All clients * 100
# within the study time period (study.start - study.end)
max.mortality.all <- read.table(
    text="AllClients
  all            1.87
  aids.related   1.49", header=TRUE, stringsAsFactors = FALSE)

#timeTo 2016-11

################# Growth rate ##########################################
#The average annual continuos growth rate Iog(Nt/No)/t * 100
#"population growth rate gives the average annual percent change in the population,"
## http://www.citypopulation.de/Swaziland.html (%/yr)
#change 1976 - 1986 | #1986 - 1997
#1997 - 2007  |  #2007 - 2012

hhohho.growth.rate <- read.table(
  text = "X1986    X1997      X2007      X2012
year10  0.0625     0.03324    0.01015    0.0022", header=TRUE, stringsAsFactor = FALSE)

#year10    6.45    3.38    1.02    0.22" log(1+6.45/100)

#Keep these as growth rates equivalent


############  Hhohho prevalence ##############################################
#Data was collected between July 2006 and March 2007 [SDHS 2006-07]:
#TimeTo = midpoint [2007-03, 2006-07]
hhohho.prev <- read.table(
      text="F.value    M.value    FM.value
  A0.150    24.2       17.4       21.0
  A15.49    33.8       23.1       28.9
  A50.150   11.2       20.7       15.4",
                          header=TRUE, stringsAsFactors = FALSE)

############  Hhohho prevalence SHIMS2 #######################################
#Data was collected between August 2016 and March 2017 [SDHS2 2016-17]:
#TimeTo = midpoint [2017-03, 2016-08]
hhohho.shims2.prev <- read.table(
  text="FM.value
  A15.150   25.7", header=TRUE, stringsAsFactors = FALSE)

################### Hhohho age diff ##########################################
# Data was collected between December 2010 and June 2011 [SHIMS 2011] TimeTo: 2011-06
hhohho.age.diff <- read.table(
      text="F.value    M.value  FM.value
  mean      7.51       5.53     6.76
  meadian   6.82       5.34     5.19  ",
                              header=TRUE, stringsAsFactors = FALSE)

########################## TARGET NAMES ############################################
tar.name <- function(df, tar.type = "name"){
  apply(expand.grid(rownames(df),".", names(df), ".",tar.type), 1, paste0,collapse="" )
}


### get the real target values
tar.value <- function(df){
  return(as.numeric(unlist(df, use.names = FALSE)))
}


#Creating target names
target.variables <- c(tar.name(hhohho.growth.rate, "gr"),
                      tar.name(max.art.initiated.all, "init"),
                      tar.name(max.art.retention.all, "ret"),
                      tar.name(max.vl.none.suppression.all, "vlsup"),
                      tar.name(max.mortality.all, "mort"),
                      tar.name(hhohho.prev, "prev"),
                      tar.name(hhohho.shims2.prev, "prev.s2"),
                      tar.name(hhohho.age.diff, "agediff")  )


#if you will need to calibrate SET THE target values correctly
target.values <- c(tar.value(hhohho.growth.rate), tar.value(max.art.initiated.all),
                   tar.value(max.art.retention.all),tar.value(max.vl.none.suppression.all),
                   tar.value(max.mortality.all), tar.value(hhohho.prev),
                   tar.value(hhohho.shims2.prev), tar.value(hhohho.age.diff))

