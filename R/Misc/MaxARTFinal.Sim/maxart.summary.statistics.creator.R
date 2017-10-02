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

####################  Swaziland specific statitics ##########################.
#GrowthRate for all ages and gender
swazi.growth.rate <- read.table(
    text="X1985   X1990   X1995  X2000   X2005  X2010   X2015   X2016
  value   3.63    3.39    2.00   1.50    0.98   1.86    1.83    1.81",
                                header=TRUE, stringsAsFactors = FALSE)

swazi.growth.rate.tar.names <- paste(names(swazi.growth.rate),"growth.rate", sep = ".")

############ Incidence UNAIDS estimates (End of march each year) #############
swazi.inci.15.49 <- read.table(
  text="X1990  X1992  X1995  X2000  X2005  X2010  X2015
Value   1.24   2.9   4.66    3.7   3.08     3.09   2.36",
                               header=TRUE, stringsAsFactors = FALSE)

swazi.inci.15.49.tar.names <- paste(names(swazi.inci.15.49),"inci.15.49", sep = ".")

############# 2011 Incidence ages specific: SHIMS 1 #############################
#consider for validation:  [all       3.14   1.65   2.38  2010-12    2011-06]
swazi.inci.2011 <- read.table(
      text="F.value    M.value
  A18.19    3.84       0.84
  A20.24    4.17       1.66
  A25.29    2.25       2.36
  A30.34    2.78       3.12
  A35.39    4.10       0.44
  A40.44    2.07       1.24
  A45.49    1.20       0.02",
                              header=TRUE, stringsAsFactors = FALSE)

########### 2006-07 prevalence ages specific: Swaziland demographics and health survey SHDS ##
# Data was collected between July 2006 and March 2007. TimeTo: 2007-03
swazi.prev.2007 <- read.table(
      text="F.value    M.value  FM.value
  A15.19    10.1       1.9      5.8
  A20.24    38.4       12.4     26.5
  A25.29    49.2       27.8     39.3
  A30.34    45.2       43.7     44.6
  A35.39    37.7       44.9     40.8
  A40.44    27.9       40.7     32.7
  A45.49    21.4       27.9     24.0
  A50.54    24.3       28.3     26.0
  A55.59    9.6        17.4     12.7
  A60.150   7.0        13.3     9.5
  A15.49    31.1       19.7     25.9
  A50.150   11.7       17.9     14.2",
                              header=TRUE, stringsAsFactors = FALSE)

#A2.4      4.8        5.5      5.1       2007-03
#A5.9      3.6        4.8      4.2       2007-03
#A10.14    3.3        1.9      2.6       2007-03
#A2.150    22.1       14.9     18.8

#Note it is not possibloe to compute the young age incidence 2- 12 yr olds since these will
#not be sexually active more of which we would expect the infection to be through MTCT



############# Prevalence ages specific UNAIDS (End of March each year) ############
swazi.prev.age.year <- read.table(
         text="F.value    M.value
  A15.24.1990    2.6       1.5
  A15.49.1990    2.3       2.2
  A15.24.1995    17.4      7.5
  A15.49.1995    14.4      16.9
  A15.24.2000    21.9      8.1
  A15.49.2000    21.8      27.6
  A15.24.2005    18.1      6.7
  A15.49.2005    21.8      29.4
  A15.24.2010    17.4      7.2
  A15.49.2010    22.6      31.9
  A15.24.2015    16.7      7.3",
                              header=TRUE, stringsAsFactors = FALSE)
#A15.49.2015    NA    34.2  28.8 ### Maybe for validation
#A15.49.1990    2.3       2.2       2.2
#A15.49.1995    14.4      16.9      15.7
#A15.49.2000    21.8      27.6      24.8
#A15.49.2005    21.8      29.4      25.7
#A15.49.2010    22.6      31.9      27.2

############  Hhohho prevalence ##############################################
#Data was collected between July 2006 and March 2007 [SDHS 2006-07]: TimeTo = 2007-03
hhohho.prev <- read.table(
      text="F.value    M.value    FM.value
  A0.150    24.2       17.4       21.0
  A15.49    33.8       23.1       28.9
  A50.150   11.2       20.7       15.4  ",
                          header=TRUE, stringsAsFactors = FALSE)

#A2.14     4.1        3.6        3.9           2007-03

########## Swazi Age difference ######################################################
# Data was collected between December 2010 and June 2011 [SHIMS 2011] TimeTo: 2011-06
swazi.age.diff <- read.table(
     text="FM.value
  mean      5.19
  meadian   6.64  ", header=TRUE, stringsAsFactors = FALSE)

swazi.AD.tar.names <- paste("FM.value", row.names(swazi.age.diff),"AD", sep = ".")

################### Hhohho age diff ##########################################
# Data was collected between December 2010 and June 2011 [SHIMS 2011] TimeTo: 2011-06
hhohho.age.diff <- read.table(
      text="F.value    M.value  FM.value
  mean      7.51       5.53     6.76
  meadian   6.82       5.34     5.19  ",
                              header=TRUE, stringsAsFactors = FALSE)

##################### ART Retention UNAIDS (End of March each year) ############
#Consider using 2015 98 and 97 % for validation
swazi.art.retention <- read.table(
                      text="X2007  X2008  X2009  X2010  X2011  X2012  X2013  X2014
  less.15.at6.months        81     88     87     87     89     84     85     97
  greateq.15.at6.months     82     85     85     86     88     87     86     96
  less.15.at12.months       73     82     81     80     89     76     91     93
  greateq.15.at12.months    76     78     78     80     80     77     88     92",
                                  header=TRUE, stringsAsFactors = FALSE)
#consider validation
#less.15.at24.months       67    75    74   73   69    87    87   NA
#greateq.15.at24.months    70    71    71   71   69    82    87   NA

################ ART coverage UNAIDS (End of March each year) ####################
swazi.art.coverage <- read.table(
             text="F.value    M.value  FM.value
  15.over.2012     94         77       87
  15.over.2013     91         76       85
  15.over.2014     85         89       93
  15.over.2015     90         73       83", header=TRUE, stringsAsFactors = FALSE)
#14.less.2015     NA  NA  72 validation

##################################### END #######################################################

## target names
swazi.art.retention.tar.names <- c()
for(k in names(swazi.art.retention)){
  swazi.art.retention.tar.names <- c(swazi.art.retention.tar.names,
                                     paste(k, row.names(swazi.art.retention),"swazi.ret", sep = "."))
}

swazi.art.coverage.tar.names <- swazi.hhohho.AD.tar.names <- c()
swazi.hhohho.prev.tar.names <- swazi.prev.2007.tar.names <- c()

for(fmt in names(swazi.art.coverage)){
  swazi.art.coverage.tar.names <- c(swazi.art.coverage.tar.names,
                                     paste(fmt, row.names(swazi.art.coverage),"swazi.art.cov", sep = "."))

  swazi.hhohho.AD.tar.names <- c(swazi.hhohho.AD.tar.names,
                                 paste(fmt, row.names(hhohho.age.diff),"Hho.AD",sep = "."))

  swazi.hhohho.prev.tar.names <- c(swazi.hhohho.prev.tar.names,
                                   paste(fmt, row.names(hhohho.prev),"Hho.prev",sep = "."))

  swazi.prev.2007.tar.names <- c(swazi.prev.2007.tar.names,
                                 paste(fmt, row.names(swazi.prev.2007),"2007.swazi.prev",sep = "."))
}


max.art.initiated.tar.names <- max.art.retention.tar.names <- c()
max.vl.none.suppression.tar.names <- max.mortality.tar.names <- c()


for(max.tar in names(max.art.initiated.all[,1:2])){

  max.art.initiated.tar.names <- c(max.art.initiated.tar.names,
                                   paste(max.tar, row.names(max.art.initiated.all),"max.ART.init", sep = "."))

  max.art.retention.tar.names <- c(max.art.retention.tar.names,
                                   paste(max.tar, row.names(max.art.retention.all),"max.ret", sep = "."))

  max.vl.none.suppression.tar.names <- c(max.vl.none.suppression.tar.names,
                                         paste(max.tar, row.names(max.vl.none.suppression.all),"max.val",sep = "."))

  max.mortality.tar.names <- c(max.mortality.tar.names,
                               paste(max.tar, row.names(max.mortality.all),"max.mort",sep = "."))

}


swazi.prev.age.year.tar.names <- swazi.inci.2011.tar.names <- c()

for(fm in names(swazi.prev.age.year[,1:2])){

  swazi.prev.age.year.tar.names <- c(swazi.prev.age.year.tar.names,
                                   paste(fm, row.names(swazi.prev.age.year),"swazi.prev",sep = "."))

  swazi.inci.2011.tar.names <- c(swazi.inci.2011.tar.names,
                                 paste(fm, row.names(swazi.inci.2011),"swazi.inci.2011",sep = "."))

}

#Creating target names
target.variables <- c(max.art.initiated.tar.names, max.art.retention.tar.names,
                      max.vl.none.suppression.tar.names, max.mortality.tar.names,
                      swazi.growth.rate.tar.names, swazi.inci.15.49.tar.names,
                      swazi.inci.2011.tar.names, swazi.prev.2007.tar.names,
                      swazi.prev.age.year.tar.names, swazi.hhohho.prev.tar.names,
                      swazi.AD.tar.names, swazi.hhohho.AD.tar.names,
                      swazi.art.retention.tar.names, swazi.art.coverage.tar.names)

#Testing
#target.variables <- c(max.art.retention.tar.names,"node.id")


