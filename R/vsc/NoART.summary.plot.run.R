#get the necessary libraries
#pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper, lhs, gtools, ggplot2, gridExtra)


#============================================================================================
#set the prior names - varied parameters
preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

############   MAIN Simulation is here #######################

simpact4ABC.chunk.wrapper <- function(simpact.chunk.prior){

  #This needs to be read by each processor
  pacman::p_load(RSimpactHelper)
  target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
                        "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
                        "ART.cov.wom.18.50", "median.wom.18.50.AD")
  err.functionGEN <- function(e){
    if (length(grep("MAXEVENTS",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    if (length(grep("internal event time",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    stop(e)
  }

  simpact.chunk.run <- function(input.chunk.params){

    pacman::p_load(RSimpactCyan, RSimpactHelper, dplyr, data.table, magrittr, exactci, tidyr)

    ## Run preprior.names.chunk and copy the results here.
    input.varied.params.plus <- c("conception.alpha_base",
                                  #"person.art.accept.threshold.dist.fixed.value",
                                  "person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b",
                                  "person.eagerness.woman.dist.gamma.a", "formation.hazard.agegapry.eagerness_diff",
                                  "person.eagerness.woman.dist.gamma.b", #"formation.hazard.agegapry.numrel_man",
                                  "formation.hazard.agegapry.numrel_woman", "formation.hazard.agegapry.gap_factor_man_exp",
                                  "formation.hazard.agegapry.gap_factor_woman_exp", #"person.agegap.man.dist.normal.mu",
                                  "person.agegap.woman.dist.normal.mu", "person.agegap.man.dist.normal.sigma",
                                  #"person.agegap.woman.dist.normal.sigma",
                                  "hivtransmission.param.f1")

    target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
                          "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
                          "ART.cov.wom.18.50", "median.wom.18.50.AD")

    simulation.type <- ("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
    simpact.set.simulation(simulation.type)
    agedist.chunk.data.frame <- agedistr.creator(shape = 5, scale = 65)

    #### Set input params
    ##Specifying the initially chosen values for the simulation.
    cfg.chunk <- input.params.creator(population.simtime = 40, population.numwomen = 1000,
                                      population.nummen = 1000, simulation.type = simulation.type)

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
      #setting up a value that is depended on the other input
      if(cfg.chunk.par == "hivtransmission.param.f1"){
        f2.num <- log((1+assign.chunk.cfg.value)/2)
        f2.den <- log(assign.chunk.cfg.value)
        cfg.chunk["hivtransmission.param.f2"][[1]] <- log(f2.num / f2.den)/5
      }
    }

    cfg.chunk["person.art.accept.threshold.dist.fixed.value"][[1]] <- 0
    cfg.chunk["mortality.aids.survtime.C"][[1]] <- 62
    cfg.chunk["mortality.aids.survtime.k"][[1]] <- -0.2

    ## Keep the files produced in subfolders
    generate.filename <- function(how.long){
      chars <- c(letters, LETTERS)
      paste0(sample(chars,how.long), collapse = "")
    }

    sub.dir.sim.id <- generate.filename(8)
    sub.dir.rename <- paste0("temp/",sub.dir.sim.id,"/")

    testoutput <- simpact.run(configParams = cfg.chunk,
                              destDir = sub.dir.rename,
                              agedist = agedist.chunk.data.frame,
                              intervention = iv.chunk,
                              identifierFormat = paste0("%T-%y-%m-%d-%H-%M-%S_%p_%r%r%r%r%r%r%r%r_",
                                                        sub.dir.sim.id,"-"),
                              seed = seed.chunk.id)

    if(testoutput$simulationtime < cfg.chunk$population.simtime)
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

    #save each of the run output.
    save(chunk.datalist.test, file = paste0("temp/","chunk.datalist.",sub.dir.sim.id,".rda"))

    #delete all the file created during the current simulation
    unlink(paste0("temp/",sub.dir.sim.id), recursive = TRUE, force = TRUE)

    if(length(chunk.datalist.test)>1){
      #get the summary statistics for each run
      growth.rate <- pop.growth.calculator(datalist = chunk.datalist.test,
                                          timewindow = c(0,
                                          timewindow.max=unique(chunk.datalist.test$itable$population.simtime)))

      inc.20.25 <- incidence.calculator(datalist = chunk.datalist.test, agegroup = c(20, 25),
                                        timewindow = c(32, 34), only.active = "No")
      inc.men.20.25 <- inc.20.25$incidence[1]
      inc.wom.20.25 <- inc.20.25$incidence[2]
      prev.25.30 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(25, 30),
                                         timepoint = 34)
      prev.men.25.30 = prev.25.30$pointprevalence[1]
      prev.wom.25.30 = prev.25.30$pointprevalence[2]
      prev.30.35 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(30, 35),
                                         timepoint = 34)
      prev.men.30.35 = prev.30.35$pointprevalence[1]
      prev.wom.30.35 = prev.30.35$pointprevalence[2]
      ARTcov <- ART.coverage.calculator(datalist = chunk.datalist.test, agegroup = c(18, 50),
                                        timepoint = 34, site="All")
      ART.cov.men.18.50 <- ARTcov$ART.coverage[1]
      ART.cov.wom.18.50 <- ARTcov$ART.coverage[2]

      agemix.df <- agemix.df.maker(chunk.datalist.test)
      pattern <- pattern.modeller(dataframe = agemix.df, agegroup = c(18, 50),
                                  timepoint = 34, timewindow = 1, start = FALSE)
      median.wom.18.50.AD <- as.numeric(median(pattern[[1]]$AgeGap[pattern[[1]]$Gender == "female"]))

      ##get the summary statistics as indicated by target.variables
      out.statistic <- c(growth.rate, inc.men.20.25, inc.wom.20.25, prev.men.25.30,
                         prev.wom.25.30, prev.men.30.35, prev.wom.30.35, ART.cov.men.18.50,
                         ART.cov.wom.18.50, median.wom.18.50.AD)
      ##out.test.degree <- out.statistic[[2]]

      #Process the summary stats for different year points for ploting.
      start.year <- 0 #chunk.datalist.test$itable$hivseed.time[1]
      end.year <- unique(chunk.datalist.test$itable$population.simtime)
      year.time.int <- start.year:end.year

      plot.sum.res <- c() #collect all the results from the wide summary calc
      plot.age.group <- c(15,50)

      for(time.point in year.time.int){
        prev.15.50.sum <- prevalence.calculator(datalist = chunk.datalist.test, agegroup = plot.age.group,
                                                timepoint = time.point)

        inc.row.index <- time.point + 3 #calc within three years
        inc.15.50.sum <- incidence.calculator(datalist = chunk.datalist.test, agegroup = plot.age.group,
                                              timewindow = c(time.point, inc.row.index),
                                              only.active = "No")

        prev.loop.result <- c()
        inci.loop.result <- c()
        #3=overal, 2=woman, 1=man
        for(k in 3:1){
          prev.loop.result <- c(prev.loop.result, as.numeric(prev.15.50.sum$pointprevalence[k]),
                                as.numeric(prev.15.50.sum$pointprevalence.95.ll[k]),
                                as.numeric(prev.15.50.sum$pointprevalence.95.ul[k]))

          inci.loop.result <- c(inci.loop.result, as.numeric(inc.15.50.sum$incidence[k]),
                                as.numeric(inc.15.50.sum$incidence.95.ll[k]),
                                as.numeric(inc.15.50.sum$incidence.95.ul[k]))
        }
        #join the results from start to end time
        plot.sum.res <- c(plot.sum.res, prev.loop.result, inci.loop.result)

      }

      out.statistic <- c(out.statistic, plot.sum.res)

    }else{
      plot.sum.res <- rep(NA,18 * length(year.time.int)) #point + CI each time step.
      out.statistic <- rep(NA,length(target.variables))
      out.statistic <- c(out.statistic, plot.sum.res)
      ##out.statistic.degree <- NA
    }

    chunk.summary.stats <- out.statistic

    return(chunk.summary.stats)
  }

  chunk.summary.stats <- tryCatch(simpact.chunk.run(simpact.chunk.prior),
                                  error = err.functionGEN)
}


start.chunk.time <- proc.time()
set.init.res <- 0
for (chunk.sim.id in inANDout.df.chunk$sim.id){

  simpact.chunk.prior = list()

  for (i in preprior.names.chunk){

    prior.chunk.val <- list(c("runif",1,
                            as.numeric(inANDout.df.chunk[inANDout.df.chunk$sim.id==chunk.sim.id,i]),
                            as.numeric(inANDout.df.chunk[inANDout.df.chunk$sim.id==chunk.sim.id,i])),
                            c("dunif",0,1))
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

  length(ABC.results.chunk.statistics)

  #get the year start = 1 to end to be ploted against
  plot.res.size <- length(ABC.results.chunk.statistics) -  length(target.variables)
  plot.sum.rep <- plot.res.size/length(plot.summary.names)

  year.time.int <- 1:plot.sum.rep

  plot.index.names <- rep(sprintf("%02d",year.time.int), each = length(plot.summary.names))
  plot.colnames.summary <- paste0(plot.summary.names, plot.index.names) #names for plot

  target.plot.target.names <- c(target.variables, plot.colnames.summary )

  names(ABC.results.chunk.statistics) <- target.plot.target.names
  ABC.results.chunk.statistics$sim.id <- chunk.sim.id

  if(set.init.res == 0){
    #rbind all the results for this chunk to be merged after
    #Create a dataframe with NA for the summary statistics for all chunks with similar sim.id
    chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 0,
                                                ncol = length(target.plot.target.names)+2))
    names(chunk.summary.stats.df) <- c(target.plot.target.names, "sim.id")
    set.init.res <- set.init.res + 1
  }

  chunk.summary.stats.df <- rbind(chunk.summary.stats.df, ABC.results.chunk.statistics)

}
chunk.summary.stats.df$ArtYesNo <- "NoART"


