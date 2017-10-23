#!/usr/bin/env/ Rscript
#get the necessary libraries
pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper, lhs)

comp <- "lin" #lin #mac #chpc #gent

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
    dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"}else if(comp=="gent"){
    dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"}else{
    dirname <- "~/Documents/RSimpactHelp"  #mac directory here
}

all.sim.start <- as.numeric(proc.time()[3])

#source simpact set parameters
source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.simpact.parameters.R")
inPUT.df.complete <- source.simpact.parameters(init.design.points = 10, resample = 1)

#select a chunk to process
sel.id.list <- c(1:10)  # or c(1,3,5)

#set min and max to for result save
min.chunk <- min(sel.id.list)
max.chunk <- max(sel.id.list)

#set how many time the single row will be repeated
sim_repeat <- 2
#number of cores per node
ncluster.use <- 4

#if you are doing many simulation you can also use a prepared file amd read as csv file.
#inPUT.df.complete <- data.frame(read.csv(file = paste0(dirname,"PARAMETER_FILE.csv"),
#                                         header = TRUE, stringsAsFactors = FALSE) )

#select a subset of the parameter set
inANDout.df.chunk <- as.data.frame(subset(inPUT.df.complete, sim.id %in% sel.id.list))
inANDout.df.chunk <- inANDout.df.chunk[!is.na(inANDout.df.chunk$sim.id),]

#indicate the target statitics that you want to hit
source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.summary.statistics.creator.R")
##Each of these should be calculated after each run, else we give an NA

#set the prior names - varied parameters
preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

#rbind all the results for this chunk to be merged after
#Create a dataframe with NA for the summary statistics Will collect all the chunks with the sim.id to link back
chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 0, ncol = length(target.variables)+1))
names(chunk.summary.stats.df) <- c(target.variables, "sim.id")

############   MAIN Simulation is here #######################
simpact4ABC.chunk.wrapper <- function(simpact.chunk.prior){

  #This needs to be read by each processor
  pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper, lhs)

  source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.summary.statistics.creator.R")

  err.functionGEN <- function(e){
    if (length(grep("MAXEVENTS",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    if (length(grep("internal event time",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    stop(e)
  }

  simpact.chunk.run <- function(input.chunk.params){

    pacman::p_load(RSimpactCyan, RSimpactHelper, dplyr, data.table, magrittr, exactci, tidyr, lhs)

    source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.simpact.parameters.R")
    getparam.names <- head(source.simpact.parameters(),1)

    #if you are doing many simulation you can also use a pre-prepared file
    #and read in a csv file.
    #getparam.names <- data.frame(read.csv(file = paste0(dirname,"PARAMETER_FILE.csv"),
    #                                         header = TRUE, stringsAsFactors = FALSE) )

    getparam.names <- names(dplyr::select(getparam.names, contains(".")))
    input.varied.params.plus <- getparam.names[2:length(getparam.names)]

    source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.summary.statistics.creator.R")

    #Set MaxART Simpact simulation?
    simulation.type <- "maxart"
    simpact.set.simulation(simulation.type)

    #This is largerly motivated by the 1970 Swaziland UN population
    #agedistr.creator(shape = 2, scale = 25)
    #sticking with our default.
    agedist.chunk.data.frame <- agedistr.creator(shape = 5, scale = 65)

    #When is HIV seeded and for how long the simulation will be.
    seed.hiv.time <- round(as.numeric(difftime(seed.hiv.date, sim.start.full, units = "days")/365.242),0)
    sim.duration <- round(as.numeric(difftime(sim.end.full,sim.start.full, units = "days")/365.242),0)

    #initial population
    num.women.prop <- round(women.frac * init.population.total, 0)
    num.men.prop <- init.population.total - num.women.prop

    ##Specifying the initially chosen values for the simulation.
    cfg.chunk <- input.params.creator(population.simtime = sim.duration,
                                      population.numwomen = num.women.prop,
                                      population.nummen = num.men.prop,
                                      hivseed.time = seed.hiv.time,
                                      debut.debutage = 12,  #maxart
                                      simulation.type = simulation.type)

    #add some of the MaxART specific parameters
    if(simulation.type == "maxart"){
      cfg.chunk$facilities.randomization <- "${SIMPACT_DATA_DIR}maxart-randomization-fake_5.csv"
      cfg.chunk$maxart.starttime <- round(as.numeric(difftime(maxart.starttime ,sim.start.full, units = "days")/365),0)
      #cfg.chunk$person.geo.dist2d.discrete.maskfile <-  "" #Defaults to Hhohho  region
    }

    #intervention introduced see the intervention.introduced
    source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.intervention.introduced.R")
    iv.chunk <- intervention.introduced(simulation.type = "pre.hhohho", simulation.start = sim.start.full)

    #The first parameter is set to be the seed value
    seed.chunk.id <- input.chunk.params[1]

    #set up the parameters to be varied in the model starting from 2:length of the varied params.
    j <- 1
    for (cfg.chunk.par in input.varied.params.plus){
      j <- j + 1
      assign.chunk.cfg.value <- input.chunk.params[j]
      cfg.chunk[cfg.chunk.par][[1]] <- assign.chunk.cfg.value
      #setting up a value that is depended on the other input (we can do this for many other as needed)
      if(cfg.chunk.par == "hivtransmission.param.f1"){
        f2.num <- log((1+assign.chunk.cfg.value)/2)
        f2.den <- log(assign.chunk.cfg.value)
        cfg.chunk["hivtransmission.param.f2"][[1]] <- log(f2.num / f2.den)/5
      }
    }

    ## Keep the files produced in subfolders
    generate.filename <- function(how.long){

      rn <- sample(1:100,1)
      t <- as.numeric(Sys.time())
      set.seed((t - floor(t)) * 1e8)
      chars <- c(letters, LETTERS)
      sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")

      noise.sample1 <- sample(8:15,1, replace = TRUE)
      sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
      noise.sample <- sample(1:1000,1)
      noise.sample2 <- sample(8:17,1, replace = TRUE)
      sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                               paste0(sample(chars,noise.sample2), collapse = ""),
                               noise.sample, rn, Sys.getpid())

      return(sub.dir.sim.id)

    }

    sub.dir.rename <- paste0("temp/",generate.filename(10))


    testoutput <- simpact.run(configParams = cfg.chunk,
                              destDir = sub.dir.rename,
                              agedist = agedist.chunk.data.frame,
                              intervention = iv.chunk,
                              identifierFormat = paste0("%T-%y-%m-%d-%H-%M-%S_%p_%r%r%r%r%r%r%r%r_"),
                              seed = seed.chunk.id)

    if(testoutput$simulationtime < cfg.chunk$population.simtime)
    {
      if (testoutput$eventsexecuted >= cfg.chunk$population.maxevents-1)
      {
        stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
      }
      else
      {
        stop("Simulation stopped prematurely, probably ran out of events")
      }
    }
    chunk.datalist.test <- readthedata(testoutput)

    #save each of the run output.
    #save(chunk.datalist.test, file = paste0("temp/chunk.datalist.",substring(sub.dir.rename, 6, 15),".rda"))

    #delete all the file created during the current simulation
    unlink(paste0(sub.dir.rename,"/"), recursive = TRUE)


    if(length(chunk.datalist.test)>1){
      #get the summary statistics for each run

      source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/pre.hho.sim.summary.targets.creator.R")
      source("R/Misc/Pre.hhohhoMaxARTFinal.Sim/source.ART.retention.R") #dummy cluster not updated
      ##get the summary statistics as indicated by target.variables
      out.statistic <- pre.hhohho.sim.summary.creator(chunk.datalist.test)
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


start.chunk.time <- as.numeric(proc.time()[3])
for (chunk.sim.id in inANDout.df.chunk$sim.id){

  simpact.chunk.prior <- list()

  for (i in preprior.names.chunk){

    #col.index <- which(colnames(preprior.names.chunk)==i)

    prior.chunk.val <- list(c("runif",1,as.numeric(inANDout.df.chunk[inANDout.df.chunk$sim.id==chunk.sim.id,i]),
                              as.numeric(inANDout.df.chunk[inANDout.df.chunk$sim.id==chunk.sim.id,i])), c("dunif",0,1))
    simpact.chunk.prior[[length(simpact.chunk.prior)+1]] <- prior.chunk.val
  }

  print(paste("Working on simulation number: ", chunk.sim.id, sep=" "))
  #invoke the ABC_rejection method repeating the number of simulation X* for each chunk row.
  ABC.chunk.result <- ABC_rejection(model = simpact4ABC.chunk.wrapper,
                                        prior = simpact.chunk.prior,
                                        nb_simul= sim_repeat,
                                        use_seed = TRUE,
                                        seed_count = 1,
                                        n_cluster = ncluster.use)


  #Save the statistics results with the chunk row sim.id repeated X* from the ABC_rejection method
  ABC.results.chunk.statistics <- data.frame(ABC.chunk.result$stats)
  names(ABC.results.chunk.statistics) <- target.variables
  ABC.results.chunk.statistics$sim.id <- chunk.sim.id

  chunk.summary.stats.df <- rbind(chunk.summary.stats.df, ABC.results.chunk.statistics)

}

#collect the parameters and the target statistics
inputANDoutput.chunk.df  <- left_join(chunk.summary.stats.df, inANDout.df.chunk, by = "sim.id")

rand.string <- paste0(sample(c(LETTERS,letters), 10), collapse = "")

filename.run <- paste0(dirname,"/","SummaryOutPut-df-",min.chunk,"-",max.chunk,".csv")

write.csv(inputANDoutput.chunk.df, file = filename.run, row.names = FALSE)

end.chunk.time <- as.numeric(proc.time()[3]) - start.chunk.time
print(paste0("Total time to compute simulation: ", end.chunk.time))

all.sim.end <- as.numeric(proc.time()[3]) - all.sim.start
print(paste0("Total time to finish simulation: ", all.sim.end))





