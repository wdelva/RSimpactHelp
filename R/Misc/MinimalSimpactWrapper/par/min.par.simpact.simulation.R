#!/usr/bin/env/ Rscript
min.par.simulation.wrapper <- function(set.seed.id = 1,   #what seed to use
                                   design.points = 16,
                                   par.repeat = 1,    #each row is repeated once.
                                   ncluster.use = 4,  #Not useful in this
                                   min.sim = 1, max.sim = 4,
                                   datalist = datalist,
                                   cal.simulation = FALSE){

  start.sim.time <- as.numeric(proc.time()[3])
  #Source all latest RSimpactHelper functions
  source("R/Misc/MinimalSimpactWrapper/min.rsimpacthelper.R")

  #source simpact set parameters
  source("R/Misc/MinimalSimpactWrapper/min.simpact.parameters.R")

  #what are the varied parameters
  #datalist obtained from the calibration for simulation
  inPUT.df.complete  <-  head(source.simpact.parameters(),1)

  #set varied parameters names
  preprior.chunk <- names(dplyr::select(inPUT.df.complete, contains(".")))
  preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

  if(cal.simulation == FALSE ){

    #This will get the initial sampled parameters
    inPUT.df.complete <- source.simpact.parameters(init.design.points = design.points, resample = 1)

    #select only simpact parameters
    inANDout.df.chunk <- subset(inPUT.df.complete, select = c("sim.id", preprior.names.chunk))

    #check if min and max makes sense.
    if(max.sim > nrow(inANDout.df.chunk)){max.sim <- nrow(inANDout.df.chunk)}
    if(min.sim < 1){min.sim <- 1}

    #indicate the simulation ids that you require.
    sel.id.list <- c(min.sim:max.sim) #or c(1, 4, 19)

    #select a subset of the parameter set
    inANDout.df.chunk <- as.data.frame(subset(inANDout.df.chunk, sim.id %in% sel.id.list))

    #preserve the original inputfile
    inANDout.df.chunk.uni <- inANDout.df.chunk

  }else{

    inPUT.df.complete <- datalist
    #select only simpact parameters
    inANDout.df.chunk <- subset(datalist, select = c("sim.id", preprior.names.chunk))

    #incase there are dublicate sim.id return only the unique ones.
    inANDout.df.chunk <- inANDout.df.chunk[!duplicated(inANDout.df.chunk[,c("sim.id")]),]

    #check if method exit, then then include it.
    if(!("method" %in% colnames(inPUT.df.complete))) { inPUT.df.complete$method <- "Init" }

    #get the method and the id number.
    method.id.name <- subset(inPUT.df.complete, select = c("sim.id","method"))

    #preserve the original inputfile
    inANDout.df.chunk.uni <- left_join(inANDout.df.chunk, method.id.name, by = "sim.id")

    min.sim <- 1
    max.sim <- nrow(inANDout.df.chunk)

  }

  #If you want to repeat the same parameter set many times default is just one.
  if(par.repeat != 1){

    inANDout.df.chunk <- bind_rows(replicate(par.repeat, inANDout.df.chunk, simplify = FALSE))

    inANDout.df.chunk <- inANDout.df.chunk %>%
      group_by(sim.id) %>%
      mutate(sim.seed = set.seed.id:(set.seed.id+n()-1),
             sim.sim.id = sim.id) %>%
      as.data.frame

    inANDout.df.chunk$sim.id <- 1:nrow(inANDout.df.chunk)
    inANDout.df.chunk$sim.id <- as.numeric(inANDout.df.chunk$sim.id)

  }

  #check if your df is being created correctly
  file.name.test.df <- paste0(dirname,"/", min.sim, "-",
                              max.sim,"-df.repeat-par.csv")

  write.csv(inANDout.df.chunk, file = file.name.test.df, row.names = FALSE)

  #indicate the target statitics that you want to hit
  source("R/Misc/MinimalSimpactWrapper/min.summary.statistics.creator.R")

  ##Each of these should be calculated after each run, else we give an NA

  ############   MAIN Simulation is here #######################
  simpact4ABC.chunk.wrapper <- function(simpact.chunk.prior){

    #This needs to be read by each processor
    pacman::p_load(data.table, tidyr, dplyr, magrittr, exactci,
                   fBasics, exactci, RSimpactCyan, readcsvcolumns, lhs)

    source("R/Misc/MinimalSimpactWrapper/min.summary.statistics.creator.R")

    err.functionGEN <- function(e){
      if (length(grep("MAXEVENTS",e$message)) != 0)
        return(chunk.summary.stats = rep(NA,length(target.variables)))
      if (length(grep("internal event time",e$message)) != 0)
        return(chunk.summary.stats = rep(NA,length(target.variables)))
      stop(e)
    }

    simpact.chunk.run <- function(input.chunk.params){

      pacman::p_load(data.table, tidyr, dplyr, magrittr, exactci,
                     fBasics, exactci, RSimpactCyan, readcsvcolumns, lhs)

      #Source all latest RSimpactHelper functions
      source("R/Misc/MinimalSimpactWrapper/min.rsimpacthelper.R")

      source("R/Misc/MinimalSimpactWrapper/min.simpact.parameters.R")
      getparam.names <- head(source.simpact.parameters(),1)
      getparam.names <- names(dplyr::select(getparam.names, contains(".")))
      input.varied.params.plus <- getparam.names[2:length(getparam.names)]

      source("R/Misc/MinimalSimpactWrapper/min.summary.statistics.creator.R")

      #Set Simpact simulation type
      simpact.set.simulation(simulation.type)

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
        cfg.chunk$facilities.randomization <- "${SIMPACT_DATA_DIR}maxart-randomization.csv"
        cfg.chunk$maxart.starttime <- round(as.numeric(difftime(maxart.starttime, sim.start.full, units = "days")/365),0)
        #cfg.chunk$person.geo.dist2d.discrete.maskfile <-  "" #Defaults to Hhohho  region
      }

      #intervention introduced see the intervention.introduced
      iv.chunk <- intervention.introduced(simulation.type = "maxart", simulation.start = sim.start.full)

      #The first parameter is set to be the seed value
      seed.chunk.id <- input.chunk.params$seed.value[1]

      #set up the parameters to be varied in the model starting from 2:length of the varied params.
      #j <- 1
      for (cfg.chunk.par in input.varied.params.plus){
        #j <- j + 1
        assign.chunk.cfg.value <- input.chunk.params[,cfg.chunk.par][1]
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

        noise.sample1 <- sample(4:10,1, replace = TRUE)
        sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
        noise.sample <- sample(1:50,1)
        noise.sample2 <- sample(8:17,1, replace = TRUE)
        sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                                 paste0(sample(chars,noise.sample2), collapse = ""),
                                 noise.sample, rn, Sys.getpid())

        return(sub.dir.sim.id)

      }

      sub.dir.rename <- paste0("temp/", generate.filename(10))

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

      #save each of the run output. (NICE to have)
      #save(chunk.datalist.test, file = paste0("temp/TESTPARAM.chunk.datalist.",
      #                               substring(sub.dir.rename, 6, 15),".rda"))

      #delete all the file created during the current simulation
      #unlink(paste0(sub.dir.rename,"/"), recursive = TRUE)

      if(length(chunk.datalist.test)>1){

        #Source all latest RSimpactHelper functions
        source("R/Misc/MinimalSimpactWrapper/min.rsimpacthelper.R")

        #get the summary statistics for each run
        source("R/Misc/MinimalSimpactWrapper/min.summary.targets.creator.R")

        ##get the summary statistics as indicated by target.variables
        out.statistic <- sim.summary.creator.all(chunk.datalist.test)

      }else{
        #if simulation producing data list create NA results
        out.statistic <- rep(NA,length(target.variables))

      }

      chunk.summary.stats <- out.statistic

      return(chunk.summary.stats)
    }

    chunk.summary.stats <- tryCatch(simpact.chunk.run(simpact.chunk.prior),
                                    error = err.functionGEN)
  }

  df.names.chunk <- names(inANDout.df.chunk)

  sim.run <- function(df.chunk.split){

    #Set each row as.data.frame with col names
    df.chunk.split <- as.data.frame(df.chunk.split)
    names(df.chunk.split) <- df.names.chunk  #it has names already??

    #select only the varied parameters.
    simpact.chunk.par <- subset(df.chunk.split, select = preprior.names.chunk)

    #set seed
    if(par.repeat != 1){
      seed.value <- as.numeric(df.chunk.split$sim.seed)
      chunk.sim.id <- df.chunk.split$sim.sim.id
    }else{
      seed.value <- set.seed.id
      chunk.sim.id <- df.chunk.split$sim.id
    }

    simpact.chunk.par$seed.value <- seed.value

    print(paste("Working on simulation number: ", chunk.sim.id, " SeeId: ", seed.value, sep=" "))

    ABC.results.chunk.statistics <- simpact4ABC.chunk.wrapper(simpact.chunk.par)

    #Save the statistics results with the chunk row sim.id
    ABC.results.chunk.statistics <- c(ABC.results.chunk.statistics, chunk.sim.id)

    return(ABC.results.chunk.statistics)

  }

  cl <- makeCluster(getOption("cl.cores", ncluster.use))
  chunk.summary.stats.df <-  NULL
  list_param <- list(NULL)

  nb_simul <- nrow(inANDout.df.chunk)

  for (i in 1:nb_simul) {
    list_param[[i]] <- inANDout.df.chunk[i, ]
  }

  cl.summary.stats.df <- parLapplyLB(cl, list_param, sim.run)

  chunk.summary.stats.df <- do.call(rbind, cl.summary.stats.df)

  stopCluster(cl)

  chunk.summary.stats.df <- as.data.frame(chunk.summary.stats.df)

  #give the calibrated values names
  names(chunk.summary.stats.df) <- c(target.variables,"sim.id")

  #collect the parameters and the target statistics
  inputANDoutput.chunk.df  <- left_join(chunk.summary.stats.df, inANDout.df.chunk.uni, by = "sim.id")

  if(par.repeat != 1){
    inputANDoutput.chunk.df <- inputANDoutput.chunk.df %>% arrange(sim.id) #order by sim.id
  }

  #Save the results to file.
  rand.string <- paste0(sample(c(LETTERS,letters), 10), collapse = "")
  filename.run <- paste0(dirname,"/", save.name, min.sim,"-", max.sim,
                         "-SummaryOutPut-",rand.string,"-par.csv")

  write.csv(inputANDoutput.chunk.df, file = filename.run, row.names = FALSE)

  #check how long the simulation took in seconds.
  end.sim.time <- as.numeric(proc.time()[3]) - start.sim.time
  print(end.sim.time)

  #Returns the df that contains the simulated statistics with
  #respective simpact parameters
  return(inputANDoutput.chunk.df)

}

