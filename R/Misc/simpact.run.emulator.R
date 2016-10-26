##neatly load all the packages listed in p_load()
rm(list=ls())
pacman::p_load(RSimpactCyan, RSimpactHelper, data.table, magrittr, dplyr, exactci,nlme,
               ggplot2, readcsvcolumns, survival, KMsurv, tidyr, lhs)

## First set up the runs from Simpact with the parameters that you want to vary.
simpact.set.simulation("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)

#### Set input params
##Specifying the initially chosen values for the simulation.
cfg <- input.params.creator(population.simtime = 40, population.numwomen = 500, population.nummen = 500)

#number of simulations repeats.
design.points <- 50
simulation.number.count <- 0
#intervention introduced
# Simulation starts in 1977. After 27 years (in 2004), ART is introduced.
iv <- intervention.introduced(list(27,0,100,2),list(30,200,1.5), list(33,350,1),list(36,500,0.5))

input.varied.params <- c("person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b", "conception.alpha_base",
                 "formation.hazard.agegapry.numrel_man", "formation.hazard.agegapry.eagerness_diff",
                 "formation.hazard.agegapry.gap_factor_man_exp", "person.agegap.man.dist.normal.mu",
                 "person.agegap.woman.dist.normal.mu","person.agegap.man.dist.normal.sigma",
                 "person.agegap.woman.dist.normal.sigma")

input.varied.params.boundaries <- list(person.eagerness.man.dist.gamma.a.min =0.1, person.eagerness.man.dist.gamma.a.max =2,
                               person.eagerness.man.dist.gamma.b.min = 5, person.eagerness.man.dist.gamma.b.max = 60,
                               conception.alpha_base.min = -3.6, conception.alpha_base.max = -1.2,
                               formation.hazard.agegapry.numrel_man.min = -1.5, formation.hazard.agegapry.numrel_man.max = -0.1,
                               formation.hazard.agegapry.eagerness_diff.min = -0.1, formation.hazard.agegapry.eagerness_diff.max = 0,
                               formation.hazard.agegapry.gap_factor_woman_exp.min = -1.5, formation.hazard.agegapry.gap_factor_woman_exp.max =-0.4,
                               person.agegap.man.dist.normal.mu.min = 0, person.agegap.man.dist.normal.mu.max = 4,
                               person.agegap.woman.dist.normal.mu.min =0, person.agegap.woman.dist.normal.mu.max = 4,
                               person.agegap.man.dist.normal.sigma.min = 0.5, person.agegap.man.dist.normal.sigma.max =2,
                               person.agegap.woman.dist.normal.sigma.min =0.5, person.agegap.woman.dist.normal.sigma.max =2)

target.variables <-c("growth.rate", "median.AD", "Q1.AD", "Q3.AD", "prev.men.15.25", "prev.men.25.50",
                       "ART.cov.15.50", "incid.wom.15.30", "frac.degreeGT1.wom.15.30", "mean.degree",
                     "median.degree", "Q1.degree", "Q3.degree")

##Rows repeat statistics for each run (saving for house keeping)
repeat.sum.stats.df <- data.frame(matrix(NA, nrow = 1, ncol = length(target.variables)))
names(repeat.sum.stats.df) <- target.variables
repeat.sum.stats.df$sim.id <- NA

## Create the list of parameters to be stimated with their min and max

# Creating the LHS over the 0-1 uniform parameter space for the parameters to be estimated
variables <- length(input.varied.params)
set.seed(1)
rlhs <- randomLHS(design.points, variables)

## run one simulation to collect output parameters
testoutput.headers <- simpact.run(configParams = cfg, destDir = "temp", agedist = agedist.data.frame,
                                  intervention = iv, seed = 1)
datalist.headers <- readthedata(testoutput.headers)

lhs.df <- as.data.frame(head(datalist.headers$itable,1))
unlink(paste(getwd(),"/temp/*", sep=""))

#rm all the columns that have text + t
lhs.df <- within(lhs.df, rm(t, population.agedistfile, periodiclogging.outfile.logperiodic, logsystem.outfile.logevents,
                            logsystem.outfile.loglocation, logsystem.outfile.logpersons, logsystem.outfile.logrelations,
                            logsystem.outfile.logsettings, logsystem.outfile.logtreatments))

lhs.df <- as.data.frame(lapply(lhs.df, rep, design.points))

#Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)

x.index <- 0
for (j in input.varied.params){
  x.index <- x.index + 1

  min.var <- input.varied.params.boundaries[paste(j,".min",sep = "")][[1]]
  max.var <- input.varied.params.boundaries[paste(j,".max",sep = "")][[1]]
  col.index <- which(colnames(lhs.df)==j)
  lhs.df[col.index] <- qunif(rlhs[ , x.index], min = as.numeric(min.var), max = as.numeric(max.var))
}

#Set those that need to use the same values
lhs.df$person.eagerness.woman.dist.gamma.a <- lhs.df$person.eagerness.man.dist.gamma.a
lhs.df$person.eagerness.woman.dist.gamma.b <- lhs.df$person.eagerness.man.dist.gamma.b
lhs.df$formation.hazard.agegapry.numrel_woman <- lhs.df$formation.hazard.agegapry.numrel_man
lhs.df$formation.hazard.agegapry.gap_factor_man_exp <- lhs.df$formation.hazard.agegapry.gap_factor_woman_exp
lhs.df$success.rows <- NA

#Create a dataframe with NA for the summary statistics
summary.stats.df <- data.frame(matrix(NA, nrow = design.points, ncol = length(target.variables)))
names(summary.stats.df) <- target.variables
# Creating a dataframe for input AND output
inANDout.df <- cbind.data.frame(sim.id = 1:design.points, rlhs, lhs.df, summary.stats.df)

# Creating a new Simpact4emulation function
# Generating New function from here
simpact4emulation <- function(sim.id){

  #changing the parameters to indicate the row in the simulation run
  for (cfg.par in input.varied.params){
    assign.cfg.value <- lhs.df[sim.id,cfg.par]
    cfg[cfg.par][[1]] <- assign.cfg.value
  }

  ## Set those parameters that need to be the same
  cfg$person.eagerness.woman.dist.gamma.a <- cfg$person.eagerness.man.dist.gamma.a
  cfg$person.eagerness.woman.dist.gamma.b <- cfg$person.eagerness.man.dist.gamma.b
  cfg$formation.hazard.agegapry.numrel_woman <- cfg$formation.hazard.agegapry.numrel_man
  cfg$formation.hazard.agegapry.gap_factor_man_exp <- cfg$formation.hazard.agegapry.gap_factor_woman_exp

  simpact.seed.id <- sim.id

  #testoutput <- simpact.run(configParams = cfg, destDir = "temp", agedist = agedist.data.frame, intervention = iv,
  #                          seed = simpact.seed.id)
  #Allow for different summarise statistics per/run
  testoutput <- simpact.run(configParams = cfg, destDir = "temp", agedist = agedist.data.frame, intervention = iv)


  if (testoutput$simulationtime < cfg$population.simtime)
  {
    if (testoutput$eventsexecuted >= cfg$population.maxevents-1)  #use ifelse
    {
      stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
    }
    else
    {
      stop("Simulation stopped prematurely, probably ran out of events")
    }
  }
  datalist.test <- readthedata(testoutput)
}



succInANDOut.df<- function(design.points=10){

  # Running the simpact4emulation function with error catcher in a loop
  for (sim.id in inANDout.df$sim.id){
    set.average.number <- 10
    #create a df to collect the repetead runs to be averaged
    outStats.df <- data.frame(matrix(NA, nrow = 1, ncol = length(target.variables)))
    names(outStats.df) <- target.variables


    for (rep.sim.id in 1:set.average.number){
      simulation.number.count <- simulation.number.count + 1
      print(paste("Simulation Count:",simulation.number.count,"Design number:", sim.id,"/",design.points, sep = " "))

      datalist.test <- tryCatch(simpact4emulation(sim.id), error = errFunction)

      if(length(datalist.test)>1){
        #get the summary statistics for each run
        out.test <- output.summary.maker(datalist = datalist.test, growth.rate=list(timewindow.min = 0, timewindow.max = unique(datalist.test$itable$population.simtime)),
                                           agemix.maker=list(agegroup.min = 15, agegroup.max=30, timepoint =30,
                                                             timewindow = 1, start=FALSE, gender = "female"),
                                           prev.15.25 = list(age.group.min=15, age.group.max=25, timepoint = 35, gender = "men"),
                                           prev.25.50 = list(age.group.min=25, age.group.max=50, timepoint = 35, gender = "men"),
                                           art.coverage = list(age.group.min=15, age.group.max=50, timepoint = 34, gender = "men"),
                                           inc.15.30 = list(age.group.min=15, age.group.max=30, timewindow.min = 30,
                                                            timewindow.max = 40, gender = "women", only.active = "Harling"),
                                           partner.degree = list(age.group.min=15, age.group.max=30, hivstatus = 0, survey.time = 30,
                                                                 window.width = 1, gender="female", only.new = FALSE))
        ##get the summary statistics as indicated by target.variables
        out.test <- out.test[,target.variables]
      }else{out.test <- rep(NA,length(target.variables))}

      #print(out.test)
      outStats.df[rep.sim.id,] <- out.test

    }
    #average the ten runs to get the averal summary statistic
    success.complete.df <- dplyr::filter(outStats.df, complete.cases(outStats.df))
    success.rows <- nrow(success.complete.df)
    out.test <- colMeans(outStats.df, na.rm = TRUE)

    outStats.df$sim.id <- sim.id
    repeat.sum.stats.df <- rbind(repeat.sum.stats.df,outStats.df)

    # Inserting the output to the inANDout dataframe
    big.insert <- length(inANDout.df) - length(target.variables) + 1

    inANDout.df[sim.id, big.insert:length(inANDout.df)] <- out.test
    inANDout.df$success.rows[sim.id] <- success.rows

    write.csv(inANDout.df, file =paste("RowUpdate","-",design.points,"Points",variables,"Par_Partial",Sys.Date(), ".csv", sep=""), row.names = FALSE)
    write.csv(repeat.sum.stats.df, file =paste("RepeatAverage.df","-",design.points,"Points",variables,"Par_Partial",Sys.Date(), ".csv", sep=""), row.names = FALSE)
    #unlink(paste(getwd(),"/temp/*", sep="")) #Delete all the files produced by simpact.run()

    #write the data stating the number of design points and the number of var parameters
  }

  write.csv(inANDout.df, file =paste("inANDout.df","-",design.points,"Points",variables,"Par",Sys.Date(), ".csv", sep=""), row.names = FALSE)
  write.csv(repeat.sum.stats.df, file =paste("RepeatAverage.df","-",design.points,"Points",variables,"Par",Sys.Date(), ".csv", sep=""), row.names = FALSE)

  return(inANDout.df)
}

#succRows.df <- success (Run the whole file).
start.time = proc.time()
inANDout.df <- succInANDOut.df(design.points)
end.time = proc.time() - start.time


##### either read the data from before or run some data here and then Do the following










