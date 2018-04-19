#Parameters  creation

#Set up simulation parameters and initiate simulation
#just the year when simulation started.
simulation.type <- "maxart" #"maxart"  #"simpact-cyan"
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

#Simulation end and HIV seed year
sim.end.full <- as.Date("2025-03-31")
seed.hiv.date <- as.Date("1986-03-31")

#preserve study yrtime
study.starttime <- maxart.starttime
study.endtime <- maxart.endtime

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
# within the simulation time period (study.start - study.end)
# Retention will be measured as the time window of maxart.

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
# max.mortality.all <- read.table(
#     text="AllClients
#   all            1.87
#   aids.related   1.49", header=TRUE, stringsAsFactors = FALSE)

#timeTo 2016-11

################# Growth rate ##########################################
#The average annual continuos growth rate Iog(Nt/No)/t * 100
#"population growth rate gives the average annual percent change in the population,"
## http://www.citypopulation.de/Swaziland.html (%/yr)
#change 1976 - 1986 | #1986 - 1997
#1997 - 2007  |  #2007 - 2012

hhohho.growth.rate <- read.table(
  text = "X1986  X1988  X1990  X1992  X1994  X1996  X1998  X2000  X2002  X2004  X2006  X2008  X2010  X2012  X2014  X2016  X2018  X2020  X2022  X2024
year02  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X", header=TRUE, stringsAsFactor = FALSE)
#year10    6.45    3.38    1.02    0.22" log(1+6.45/100)

#Keep these as growth rates equivalent

############  Hhohho incidence ##############################################
#Data was collected between September 2014 and November 2017 [MaxART study]:
#TimeTo = between (2014-09-01 - 2017-11-03)
hhohho.inci.2016 <- read.table(
      text="F.val   Fu.val   Fl.val  M.val   Mu.val  Ml.val   FM.val   FMu.val  FMl.val
  A15.25    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A15.49    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A15.150   NA  NA  NA  NA  NA  NA  NA  NA  NA
  A25.34    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A35.44    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A45.49    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A50.150   NA  NA  NA  NA  NA  NA  NA  NA  NA", header=TRUE, stringsAsFactors = FALSE)

#zero-150 can not be accurately measured the model does not have MTCT


############  Hhohho incidence YOUNG Woman [2014 - 2025] ###############################

#TimeTo = between (2014 = 2013-03-01 - 2004-02-28 ; 2015 = 2014-03-01 - 2015-02-28)
hhohho.inci.year.list <- c(2014:2025)
hhohho.inci.year <- read.table(
  text="F.val   Fu.val   Fl.val  M.val   Mu.val  Ml.val   FM.val   FMu.val  FMl.val
  A15.18    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A18.25    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A15.25    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A25.50    NA  NA  NA  NA  NA  NA  NA  NA  NA", header=TRUE, stringsAsFactors = FALSE)


############  Hhohho prevalence SHIMS2 #######################################
#Data was collected between August 2016 and March 2017 [SDHS2 2016-17]:
#TimeTo = midpoint [2017-03, 2016-08]
hhohho.shims2.prev <- read.table(
  text="FM.val   FMl.val  FMu.val
  A15.150   25.7  NA  NA", header=TRUE, stringsAsFactors = FALSE)

############  Hhohho prevalence ##############################################
#Data was collected between July 2006 and March 2007 [SDHS 2006-07]:
#TimeTo = midpoint [2007-03, 2006-07]
hhohho.prev <- read.table(
  text="F.val   Fu.val   Fl.val  M.val   Mu.val  Ml.val   FM.val   FMu.val  FMl.val
  A15.25    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A15.49    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A15.150   NA  NA  NA  NA  NA  NA  NA  NA  NA
  A25.34    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A35.44    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A45.49    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A50.150   NA  NA  NA  NA  NA  NA  NA  NA  NA", header=TRUE, stringsAsFactors = FALSE)

# Get the year prev (2014 - 2025) where 2014 = midpoint(2013-03, 2014-03)
hhohho.prev.year.list <- c(2014:2025)
hhohho.prev.year <- read.table(
  text="F.val   Fu.val   Fl.val  M.val   Mu.val  Ml.val   FM.val   FMu.val  FMl.val
  A15.25    NA  NA  NA  NA  NA  NA  NA  NA  NA
  A25.50    NA  NA  NA  NA  NA  NA  NA  NA  NA", header=TRUE, stringsAsFactors = FALSE)

################### Hhohho age diff ##########################################
# Data was collected between December 2010 and June 2011 [SHIMS 2011] TimeTo: 2011-06
hhohho.age.diff <- read.table(
      text="F.value    M.value  FM.value
  mean      7.51       5.53     6.76
  meadian   6.82       5.34     5.19  ",
                              header=TRUE, stringsAsFactors = FALSE)

########################## TARGET NAMES ############################################
#DF has all the information in rows and col names to calculate the statistic
tar.name <- function(df, tar.type = "name"){
  apply(expand.grid(rownames(df),".", names(df), ".",tar.type), 1, paste0,collapse="" )
}

#Use the different year for variable
tar.yr.name <- function(df, df.yr, tar.type = "name"){
  df.names <- apply(expand.grid(rownames(df),".", names(df), ".",tar.type), 1, paste0,collapse="" )
  df.names <- apply(expand.grid(df.names, ".", df.yr), 1, paste0, collapse="")
  return(df.names)
}

#Creating target names
target.variables <- c(tar.name(hhohho.growth.rate, paste0("gr",".",final.type)),
                      tar.name(max.art.initiated.all, paste0("art.cov",".",final.type)),
                      tar.name(max.art.retention.all, paste0("ret",".",final.type)),
                      tar.name(hhohho.inci.2016, paste0("inci",".",final.type)),
                      tar.yr.name(hhohho.inci.year, hhohho.inci.year.list, paste0("inci.yr",".",final.type)),
                      tar.name(max.vl.none.suppression.all, paste0("vl.sup",".",final.type)),
                      tar.name(hhohho.prev, paste0("prev",".",final.type)),
                      tar.yr.name(hhohho.prev.year, hhohho.prev.year.list, paste0("prev.yr",".",final.type)),
                      tar.name(hhohho.shims2.prev, paste0("prev.s2",".",final.type)),
                      tar.name(hhohho.age.diff, paste0("agediff",".",final.type)  )
                      )


