#!/usr/bin/env/ Rscript
#get the necessary libraries
pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper)
#data file to read

comp <- "win" #lin #mac

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
  dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else{dirname <- "~/Documents/RSimpactHelp"  #mac directory here
}

# There may be ways to make main.filename dynamic: so the file name does not need to be manually updated
main.filename <- "simpactInputParams.df-200Points11Par2017-01-13.csv" #Read the file produced by varying parameters *design.points
file.chunk.name.csv <-paste0(dirname, "/", main.filename) #### Input file name is produced from the .sh script
inPUT.df.complete <- read.csv(file = file.chunk.name.csv, header = TRUE, sep = ",")

################################# YOU CAN EITHER RUNN THIS LINE BELOW or READ THE FILE SAVED ALREADY ############################
# inPUT.df.complete <- simpact.config.inputs(design.points = 200, conception.alpha_base = c(-3.6, -1.2),
#                                            formation.hazard.agegapry.baseline = c(1.5, 3),
#                                            person.eagerness.man.dist.gamma.a = c(0.1, 2),
#                                            person.eagerness.man.dist.gamma.b = c(5, 60),
#                                            formation.hazard.agegapry.numrel_man = c(-1.5, -0.1),
#                                            formation.hazard.agegapry.eagerness_diff = c(-0.1, 0),
#                                            formation.hazard.agegapry.gap_factor_man_exp = c(-1.5, -0.4),
#                                            person.agegap.man.dist.normal.mu = c(0, 4),
#                                            person.agegap.woman.dist.normal.mu = c(0, 4),
#                                            person.agegap.man.dist.normal.sigma = c(0.5, 2),
#                                            person.agegap.woman.dist.normal.sigma = c(0.5, 2) )

##################################################################################################################################

#Select a chunk to send to process
min.chunk <- 1
max.chunk <- 200
inANDout.df.chunk <- inPUT.df.complete[min.chunk:max.chunk,]

#make sure there are no empty rows
inANDout.df.chunk <- inANDout.df.chunk[!is.na(inANDout.df.chunk$sim.id),]

sim_repeat <- 5
ncluster.use <- 4 # number of cores per node

#indicate the target statitics that you want to hit
target.variables <- c("growth.rate", "inc.men.15.30", "inc.wom.15.30", "prev.men.15.25",
                      "prev.wom.15.25","prev.men.25.50", "prev.wom.25.50", "ART.cov.men.15.50", "ART.cov.wom.15.50")

##Each of these should be calculated after each run, else we give an NA

#set the prior names - varied parameters
preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

#rbind all the results for this chunk to be merged after
#Create a dataframe with NA for the summary statistics Will collect all the chunks with the sim.id to link back
chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 0, ncol = length(target.variables)+1))
names(chunk.summary.stats.df) <- c(target.variables, "sim.id")


simpact4ABC.chunk.wrapper <- function(simpact.chunk.prior){

  #This needs to be read by each processor
  pacman::p_load(RSimpactHelper)
  target.variables <- c("growth.rate", "inc.men.15.30", "inc.wom.15.30", "prev.men.15.25",
                        "prev.wom.15.25","prev.men.25.50", "prev.wom.25.50", "ART.cov.men.15.50", "ART.cov.wom.15.50")
  err.functionGEN <- function(e){
    if (length(grep("MAXEVENTS",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    if (length(grep("internal event time",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    stop(e)
  }

  simpact.chunk.run <- function(input.chunk.params){

    pacman::p_load(RSimpactCyan, RSimpactHelper, dplyr,lhs,data.table, dplyr, magrittr, exactci,
                   nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
                   igraph,lhs, GGally, tidyr)

    ## Run preprior.names.chunk and copy the results here.
    input.varied.params.plus <- c("conception.alpha_base", "formation.hazard.agegapry.baseline",
                                  "person.eagerness.man.dist.gamma.a",  "person.eagerness.man.dist.gamma.b")

    target.variables <- c("growth.rate", "inc.men.15.30", "inc.wom.15.30", "prev.men.15.25", "prev.wom.15.25","prev.men.25.50",
                          "prev.wom.25.50", "ART.cov.men.15.50", "ART.cov.wom.15.50")

    simulation.type <- ("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
    simpact.set.simulation(simulation.type)
    agedist.chunk.data.frame <- agedistr.creator(shape = 5, scale = 65)

    #### Set input params
    ##Specifying the initially chosen values for the simulation.
    cfg.chunk <- input.params.creator(population.simtime = 40, population.numwomen = 500, population.nummen = 500,
                                      simulation.type = simulation.type)

    #intervention introduced See the intervention.introduced
    # Simulation starts in 1977. After 27 years (in 2004), ART is introduced.
    iv.chunk <- intervention.introduced(simulation.type = simulation.type)

    #The first parameter is set to be the seed value
    seed.chunk.id <- input.chunk.params[1]

    #set up the parameters to be varied in the model starting from 2:length of the varied params.
    j <- 1
    for (cfg.chunk.par in input.varied.params.plus){
      j <- j + 1
      assign.chunk.cfg.value <- input.chunk.params[j]
      cfg.chunk[cfg.chunk.par][[1]] <- assign.chunk.cfg.value
    }

    ## Keep the files produced in subfolders
    generate.filename <- function(how.long){
      chars <- c(letters, LETTERS)
      paste0(sample(chars,how.long), collapse = "")
    }

    print(cfg.chunk)

    sub.dir.sim.id <- generate.filename(8)
    sub.dir.rename <- paste0("temp/",sub.dir.sim.id,"/")

    testoutput <- simpact.run(configParams = cfg.chunk,
                              destDir = sub.dir.rename,
                              agedist = agedist.chunk.data.frame,
                              intervention = iv.chunk,
                              identifierFormat = paste0("%T-%y-%m-%d-%H-%M-%S_%p_%r%r%r%r%r%r%r%r_",sub.dir.sim.id,"-"),
                              seed = seed.chunk.id)

    if(testoutput$simulationtime < cfg.chunk$population.simtime)#{chunk.datalist.test <- "Error"}else{chunk.datalist.test <- readthedata(testoutput)}
    {
      if (testoutput$eventsexecuted >= cfg.chunk$population.maxevents-1)  #use ifelse
      {
        stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
      }
      else
      {
        stop("Simulation stopped prematurely, probably ran out of events")
      }
    }
    chunk.datalist.test <- readthedata(testoutput)

    if(length(chunk.datalist.test)>1){
      #get the summary statistics for each run
      growth.rate = pop.growth.calculator(datalist = chunk.datalist.test,
                                          timewindow = c(0, timewindow.max=unique(chunk.datalist.test$itable$population.simtime)))

      inc.15.30 <- incidence.calculator(datalist = chunk.datalist.test, agegroup = c(15, 30), timewindow = c(20, 30), only.active = "No")
      inc.men.15.30 <- inc.15.30$incidence[1]
      inc.wom.15.30 <- inc.15.30$incidence[2]
      prev.15.25 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(15, 25), timepoint = 35)
      prev.men.15.25 = prev.15.25$pointprevalence[1]
      prev.wom.15.25 = prev.15.25$pointprevalence[2]
      prev.25.50 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(25, 50), timepoint = 35)
      prev.men.25.50 = prev.25.50$pointprevalence[1]
      prev.wom.25.50 = prev.25.50$pointprevalence[2]
      ARTcov <- ART.coverage.calculator(datalist = chunk.datalist.test, agegroup = c(15, 50), timepoint =35, site="All")
      ART.cov.men.15.50 <- ARTcov$ART.coverage[1]
      ART.cov.wom.15.50 <- ARTcov$ART.coverage[2]

      ##get the summary statistics as indicated by target.variables
      out.statistic <- c(growth.rate, inc.men.15.30, inc.wom.15.30, prev.men.15.25, prev.wom.15.25,
                         prev.men.25.50, prev.wom.25.50, ART.cov.men.15.50, ART.cov.wom.15.50)
      ##out.test.degree <- out.statistic[[2]]
    }else{
      out.statistic <- rep(NA,length(target.variables))
      ##out.statistic.degree <- NA
    }

    chunk.summary.stats <- out.statistic

    return(chunk.summary.stats)
  }

  chunk.summary.stats <- tryCatch(simpact.chunk.run(simpact.chunk.prior),
                                  error = err.functionGEN)
}


start.chunk.time <- proc.time()
for (chunk.sim.id in inANDout.df.chunk$sim.id){

  simpact.chunk.prior = list()

  for (i in preprior.names.chunk){

    #col.index <- which(colnames(preprior.names.chunk)==i)

    prior.chunk.val <- list(c("runif",1,as.numeric(inANDout.df.chunk[chunk.sim.id,i]),
                              as.numeric(inANDout.df.chunk[chunk.sim.id,i])), c("dunif",0,1))
    simpact.chunk.prior[[length(simpact.chunk.prior)+1]] <- prior.chunk.val
  }

  print(paste("Working on simulation number: ", chunk.sim.id, sep=" "))
  #invoke the ABC_rejection method repeating the number of simulation X* for each chunk row.
  ABC.chunk.result <- ABC_rejection(model = simpact4ABC.chunk.wrapper,
                                        prior = simpact.chunk.prior,
                                        nb_simul= sim_repeat,
                                        use_seed = TRUE,
                                        seed_count = 0,
                                        n_cluster = ncluster.use)


  #Save the statistics results with the chunk row sim.id repeated X* from the ABC_rejection method
  ABC.results.chunk.statistics <- data.frame(ABC.chunk.result$stats)
  names(ABC.results.chunk.statistics) <- target.variables
  ABC.results.chunk.statistics$sim.id <- chunk.sim.id

  chunk.summary.stats.df <- rbind(chunk.summary.stats.df, ABC.results.chunk.statistics)

}

inputANDoutput.chunk.df  <- left_join(chunk.summary.stats.df, inANDout.df.chunk, by = "sim.id")


write.csv(inputANDoutput.chunk.df, file =paste0("SummaryOutPut-inANDout.df.chunk-",min.chunk,"-",max.chunk,"-",Sys.Date(),
                                               ".csv"), row.names = FALSE)
end.chunk.time <- proc.time() - start.chunk.time


