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
init.population.total <- 1000
women.frac <- 0.5253


################### ART initiation ##############################
max.art.initiated.all <- read.table(
        text="PrimaryClients   AllClients
  all              77.6         84.6
  before.eligible  4.6          2.2
  within2.weeks    61           64
  among.eligible   85           90
  within12.months  85           91
  within6months    79           87", header=TRUE, stringsAsFactors = FALSE)

#TimeTo: 2016-11

##########################  ART Retention  #############################
max.art.retention.all <- read.table(
  text="PrimaryClients AllClients
  all             72      73
  at6.months      89      87
  at12.months     82      79", header=TRUE, stringsAsFactors = FALSE)
#timeTo: 2016-11

###########################  ViralLoadSuppression ##################
max.vl.none.suppression.all <- read.table(
  text="PrimaryClients AllClients
  at6.months       4      6
  at12.months      6      11", header=TRUE, stringsAsFactors = FALSE)

#Check validation
#atleast6.months  4      8       2016-11

#############  ViralLoadSuppression ###########################
max.mortality.all <- read.table(
    text="PrimaryClients AllClients
  all            1.01      1.87
  aids.related   0.51      1.49", header=TRUE, stringsAsFactors = FALSE)

#timeTo 2016-11

############  Hhohho prevalence ##############################################
#Data was collected between July 2006 and March 2007 [SDHS 2006-07]: TimeTo = 2007-03
hhohho.prev <- read.table(
      text="F.value    M.value    FM.value
  A0.150    24.2       17.4       21.0
  A15.49    33.8       23.1       28.9
  A50.150   11.2       20.7       15.4",
                          header=TRUE, stringsAsFactors = FALSE)

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

#Creating target names
target.variables <- c(tar.name(max.art.initiated.all, "init"),
                      tar.name(max.art.retention.all, "ret"),
                      tar.name(max.vl.none.suppression.all, "vlsup"),
                      tar.name(max.mortality.all, "mort"),
                      tar.name(hhohho.prev, "prev"),
                      tar.name(hhohho.age.diff, "agediff")  )

#Testing
#target.variables <- c(max.art.retention.tar.names,"node.id")

