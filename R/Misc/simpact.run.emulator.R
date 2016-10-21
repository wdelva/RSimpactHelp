##neatly load all the packages listed in p_load()
pacman::p_load(RSimpactCyan, RSimpactHelper, data.table, magrittr, dplyr,exactci,nlme,
               ggplot2, readcsvcolumns, survival, KMsurv, tidyr, lhs)

## First set up the runs from Simpact with the parameters that you want to vary.
simpact.set.simulation("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)

#### Set input params
##Specifying the initially chosen values for the simulation.
cfg <- input.params.creator(population.simtime = 40, population.numwomen = 500, population.nummen = 500)

#number of simulations repeats.
design.points <- 8000
simulation.number.count <- 0
#intervention introduced
# Simulation starts in 1977. After 27 years (in 2004), ART is introduced.
art.intro <- list()
art.intro["time"] <- 27
art.intro["diagnosis.baseline"] <- 0 # Reset to zero, from its original value of -100 at the start of the simulatio
art.intro["monitoring.cd4.threshold"] <- 100
art.intro["diagnosis.genderfactor"] <- 2

# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

art.intro2 <- list()
art.intro2["time"] <- 30
art.intro2["monitoring.cd4.threshold"] <- 200
art.intro2["diagnosis.genderfactor"] <- 1.5

art.intro3 <- list()
art.intro3["time"] <- 33
art.intro3["monitoring.cd4.threshold"] <- 350
art.intro3["diagnosis.genderfactor"] <- 1

art.intro4 <- list()
art.intro4["time"] <- 36
art.intro4["monitoring.cd4.threshold"] <- 500
art.intro4["diagnosis.genderfactor"] <- 0.5

# person.art.accept.threshold.dist.fixed.value

iv <- list(art.intro, art.intro2, art.intro3, art.intro4)


summary.statistics <- function(datalist.test){
  #summary statistics #################
  # 1. Population growth
  growth.rate <- pop.growth.calculator(datalist = datalist.test,
                                       timewindow = c(0, unique(datalist.test$itable$population.simtime)))
  # 2. Median age difference in most recent relationships (target 2-5 years)
  agemix.df <- agemix.df.maker(datalist.test)
  pattern <- pattern.modeller(dataframe = agemix.df,
                              agegroup = c(15, 30),
                              timepoint = 30,
                              timewindow = 1,
                              start = FALSE)
  # agedifmedtable$median[2] # Median age difference as reported by women
  median.AD <- as.numeric(median(pattern[[1]]$AgeGap[pattern[[1]]$Gender == "female"]))

  # 3. IQR age difference in most recent relationships (target IQR 2 - 6)
  Q1.AD <- as.numeric(summary(pattern[[1]]$AgeGap[pattern[[1]]$Gender == "female"])[2])
  Q3.AD <- as.numeric(summary(pattern[[1]]$AgeGap[pattern[[1]]$Gender == "female"])[5])
  #IQR.AD <- agedifmedtable$Q3[2] - agedifmedtable$Q1[2] # IQR as reported by women

  # 4. HIV prevalence in the window among men 15-25 (target 5% - 10%)
  prev <- prevalence.calculator(datalist = datalist.test, agegroup = c(15, 25), timepoint = 35)
  prev.men.15.25 <- prev$pointprevalence[1] # among men

  # 5. HIV prevalence in the window among men 25-50 (target 10% - 30%)
  prev <- prevalence.calculator(datalist = datalist.test, agegroup = c(25, 50), timepoint = 35)
  prev.men.25.50 <- prev$pointprevalence[1] # among men

  # 6. Point prevalence of concurrency. postsim function to be converted to RSimpactHelper function

  # 7. ART coverage among all HIV+ people aged 15-50 (target 15% - 40% in 2011)
  ARTcov <- ART.coverage.calculator(datalist = datalist.test, agegroup = c(15, 50), timepoint = 34) # 2011 is 34 years after 1977
  ART.cov.15.50 <- ARTcov$ART.coverage[1] # among men and women

  # 8. HIV incidence among women 15 <= x < 30 in the window 20-30 years after HIV introduction
  inc <- incidence.calculator(datalist = datalist.test, agegroup = c(15, 30), timewindow = c(30, 40), only.active = "Harling")
  incid.wom.15.30 <- inc$incidence[2] # among women

  # NOTE: We may want to also calculate a different type of HIV incidence, more in line with the Harling paper:
  # Where we only accumulate exposure time for the periods that a woman was in (a) relationship(s).
  # For now, this (correct) HIV incidence measure suffices

  # 9. % of women who had > 1 partner in the past 12 months
  degree.df <- degree.df.maker(agemix.df, agegroup = c(15, 30), hivstatus = 0, survey.time = 30, window.width = 1,
                               only.new = FALSE)
  frac.degreeGT1.wom.15.30 <- mean(degree.df$Degree > 1)

  out.test <- matrix(c(growth.rate, median.AD, Q1.AD, Q3.AD, prev.men.15.25, prev.men.25.50, ART.cov.15.50, incid.wom.15.30, frac.degreeGT1.wom.15.30),
                     nrow = 1, dimnames = list(NULL, c("growth.rate", "median.AD", "Q1.AD", "Q3.AD", "prev.men.15.25", "prev.men.25.50", "ART.cov.15.50", "incid.wom.15.30", "frac.degreeGT1.wom.15.30"))) # 9 summary statistics

  return(out.test)
}

# Creating an error function to catch the case when population.maxevents is reached before population.simtime is reached
errFunction <- function(e){
  if (length(grep("MAXEVENTS",e$message)) != 0)
    return("Not Complete")
  if (length(grep("internal event time",e$message)) != 0)
    return("Not Complete")
  # Een andere foutmelding dan MAXEVENTS, zorg dat we toch stoppen
  stop(e)
}


input.varied.params <- c("person.eagerness.dist.gamma.a", "person.eagerness.dist.gamma.b", "conception.alpha_base",
                         "formation.hazard.agegapry.numrel_man","formation.hazard.agegapry.numrel_woman",
                         "formation.hazard.agegapry.eagerness_diff","formation.hazard.agegapry.gap_factor_man_exp",
                         "formation.hazard.agegapry.gap_factor_woman_exp", "person.agegap.man.dist.normal.mu",
                         "person.agegap.woman.dist.normal.mu","person.agegap.man.dist.normal.sigma",
                         "person.agegap.woman.dist.normal.sigma")


## run one simulation to collect the parameters (Remove if not needed alternative way)
##testoutput.headers <- simpact.run(configParams = cfg ,destDir = "temp", intervention = iv, agedist = agedist.data.frame, seed = 1)
##datalist.headers <- readthedata(testoutput.headers)

## Create the list of parameters to be stimated with their min and max

# Creating the LHS over the 0-1 uniform parameter space for the parameters to be estimated
variables <- length(input.varied.params)

set.seed(1)
rlhs <- randomLHS(design.points, variables)

# Creating the LHS dataframe
#parameter list with input + default values

testoutput.headers <- simpact.run(configParams = cfg, destDir = "temp", agedist = agedist.data.frame, intervention = iv, seed = 1)

datalist.headers <- readthedata(testoutput.headers)

lhs.df <- as.data.frame(head(datalist.headers$itable,1))
unlink(paste(getwd(),"/temp/*", sep=""))


#lhs.df <- as.data.frame(simpact.getconfig(cfg))

#rm all the columns that have text + t
lhs.df <- within(lhs.df, rm(t, population.agedistfile, periodiclogging.outfile.logperiodic, logsystem.outfile.logevents,
                            logsystem.outfile.loglocation, logsystem.outfile.logpersons, logsystem.outfile.logrelations,
                            logsystem.outfile.logsettings, logsystem.outfile.logtreatments))

lhs.df <- as.data.frame(lapply(lhs.df, rep, design.points))

#Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)

lhs.df$person.eagerness.dist.gamma.a <- qunif(rlhs[ , 1], min = 0.1, max = 2)
lhs.df$person.eagerness.dist.gamma.b <- qunif(rlhs[ , 2], min = 5, max = 60)
lhs.df$conception.alpha_base <- qunif(rlhs[ , 3], min = -3.6, max = -1.2)
lhs.df$formation.hazard.agegapry.numrel_man <- qunif(rlhs[ , 4], min = -1.5, max = -0.1)
lhs.df$formation.hazard.agegapry.numrel_woman <- lhs.df$formation.hazard.agegapry.numrel_man
lhs.df$formation.hazard.agegapry.eagerness_diff <- qunif(rlhs[ , 5], min = -0.1, max = 0)
lhs.df$formation.hazard.agegapry.gap_factor_woman_exp <- qunif(rlhs[ , 6], min = -1.5, max = -0.4)
lhs.df$formation.hazard.agegapry.gap_factor_man_exp <- lhs.df$formation.hazard.agegapry.gap_factor_woman_exp
lhs.df$person.agegap.man.dist.normal.mu <- qunif(rlhs[ , 7], min = 0, max = 4)
lhs.df$person.agegap.woman.dist.normal.mu <- qunif(rlhs[ , 8], min = 0, max = 4)
lhs.df$person.agegap.man.dist.normal.sigma <- qunif(rlhs[ , 9], min = 0.5, max = 2)
lhs.df$person.agegap.woman.dist.normal.sigma <- qunif(rlhs[ , 10], min = 0.5, max = 2)

# Creating a dataframe for input AND output
inANDout.df <- cbind.data.frame(sim.id = 1:design.points,
                                rlhs,
                                lhs.df,
                                growth.rate = rep(NA, design.points),
                                median.AD = rep(NA, design.points),
                                Q1.AD = rep(NA, design.points),
                                Q3.AD = rep(NA, design.points),
                                prev.15.25 = rep(NA, design.points),
                                prev.25.50 = rep(NA, design.points),
                                ART.cov.15.50 = rep(NA, design.points),
                                incid.wom.15.30 = rep(NA, design.points),
                                frac.degreeGT1.wom.15.30 = rep(NA, design.points)
                                )

##inANDout.df.bk <- read.csv("C:/Users/tchibawara/Documents/MaxART/RSimpactHelp/RowUpdate-8000Points12Par_Partial2016-10-17.csv")
##inANDout.df <- inANDout.df.bk



# Creating a new Simpact4emulation function
# Generating New function from here
simpact4emulation <- function(sim.id, lhs.df, cfg, agedist.data.frame, iv){

  cfg$formation.hazard.agegapry.numrel_man <- lhs.df$formation.hazard.agegapry.numrel_man[sim.id]
  cfg$formation.hazard.agegapry.numrel_woman <- lhs.df$formation.hazard.agegapry.numrel_man[sim.id]
  cfg$formation.hazard.agegapry.eagerness_diff <- lhs.df$formation.hazard.agegapry.eagerness_diff[sim.id]
  cfg$conception.alpha_base <- lhs.df$conception.alpha_base[sim.id]
  cfg$person.eagerness.dist.gamma.a <- lhs.df$person.eagerness.dist.gamma.a[sim.id]
  cfg$person.eagerness.dist.gamma.b <- lhs.df$person.eagerness.dist.gamma.b[sim.id]
  cfg$formation.hazard.agegapry.gap_factor_man_exp <- lhs.df$formation.hazard.agegapry.gap_factor_man_exp[sim.id]
  cfg$formation.hazard.agegapry.gap_factor_woman_exp <- lhs.df$formation.hazard.agegapry.gap_factor_woman_exp[sim.id]
  cfg$person.agegap.man.dist.normal.mu <- lhs.df$person.agegap.man.dist.normal.mu[sim.id]
  cfg$person.agegap.woman.dist.normal.mu <- lhs.df$person.agegap.woman.dist.normal.mu[sim.id]
  cfg$person.agegap.man.dist.normal.sigma <- lhs.df$person.agegap.man.dist.normal.sigma[sim.id]
  cfg$person.agegap.woman.dist.normal.sigma <- lhs.df$person.agegap.woman.dist.normal.sigma[sim.id]
  simpact.seed.id <- sim.id

  testoutput <- simpact.run(configParams = cfg, destDir = "temp", agedist = agedist.data.frame, intervention = iv,
                            seed = simpact.seed.id)


  if (testoutput$simulationtime < cfg$population.simtime)
  {
    # Ik kan op dit moment twee redenen bedenken waarom de simulatie te vroeg zou stoppen
    #  - maximaal aantal events is bereikt, kan op gecheckt worden dmv ret["eventsexecuted"]
    #  - geen events meer, gebeurt bvb als populatie uitsterft.
    if (testoutput$eventsexecuted >= cfg$population.maxevents-1)
    {
      # Ik doe hier een -1 omdat in R getallen standaard voorgesteld worden als floating
      # point getallen, en echte gelijkheden daarmee nogal gevaarlijk zijn.
      stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
    }
    else
    {
      # Misschien moet er in dit geval niet gestopt worden
      stop("Simulation stopped prematurely, probably ran out of events")
    }
  }

  #Compute the target statistics from each run of the simulation.
  datalist.test <- readthedata(testoutput)

}


succInANDOut.df<- function(design.points=10){

  # Running the simpact4emulation function with error catcher in a loop
  for (sim.id in inANDout.df$sim.id){

    if(sim.id > 832){
      #create a df to collect the repetead runs to be averaged
      outStats <- as.data.frame(cbind(growth.rate = rep(NA, 10), median.AD = rep(NA, 10), Q1.AD = rep(NA, 10), Q3.AD = rep(NA, 10),
                        prev.15.25 = rep(NA, 10), prev.25.50 = rep(NA, 10), ART.cov.15.50 = rep(NA, 10),
                        incid.wom.15.30 = rep(NA, 10), frac.degreeGT1.wom.15.30 = rep(NA, 10)))

      for (j in 1:10){
        simulation.number.count <- simulation.number.count + 1
        print(paste("Simulation Count:",simulation.number.count,"Design number:", sim.id, sep = " "))

        datalist.test <- tryCatch(simpact4emulation(sim.id, lhs.df, cfg, agedist.data.frame, iv), error = errFunction)
        if(length(datalist.test)>1){out.test <- summary.statistics(datalist.test)}else{out.test <- rep(NA,9)}

        outStats[j,] <- out.test

      }


      out.test <- c(mean(as.numeric(outStats$growth.rate), na.rm=TRUE), mean(as.numeric(outStats$median.AD), na.rm=TRUE),
                    mean(as.numeric(outStats$Q1.AD), na.rm=TRUE), mean(as.numeric(outStats$Q3.AD), na.rm=TRUE),
                    mean(as.numeric(outStats$prev.15.25), na.rm=TRUE), mean(as.numeric(outStats$prev.25.50), na.rm=TRUE),
                    mean(as.numeric(outStats$ART.cov.15.50), na.rm=TRUE), mean(as.numeric(outStats$incid.wom.15.30), na.rm=TRUE),
                    mean(as.numeric(outStats$frac.degreeGT1.wom.15.30), na.rm=TRUE))



      # Inserting the output to the inANDout dataframe


      big.insert <- length(inANDout.df) - 8

      inANDout.df[sim.id, big.insert:length(inANDout.df)] <- out.test

      write.csv(inANDout.df, file =paste("RowUpdate","-",design.points,"Points",variables,"Par_Partial",Sys.Date(), ".csv", sep=""), row.names = FALSE)
      unlink(paste(getwd(),"/temp/*", sep=""))
    }
    #write the data stating the number of design points and the number of var parameters
  }

  write.csv(inANDout.df, file =paste("inANDout.df","-",design.points,"Points",variables,"Par",Sys.Date(), ".csv", sep=""), row.names = FALSE)

  return(inANDout.df)

}

#succRows.df <- success
start.time = proc.time()
inANDout.df <- succInANDOut.df(design.points)
end.time = proc.time() - start.time


