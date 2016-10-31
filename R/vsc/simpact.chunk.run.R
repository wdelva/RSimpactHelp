simpact.chunk.run <- function(input.chunk.params){
  pacman::p_load(RSimpactCyan, RSimpactHelper, dplyr,lhs,data.table, dplyr, magrittr, exactci,
                 nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
                 igraph,lhs, GGally, emulator, multivator, tidyr)
  ### If more varaiables are added to the lhs.df (Latin Hypercude Sample), This is input variable needs is needed
  input.varied.params.plus <- c("person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b", "conception.alpha_base",
                           "formation.hazard.agegapry.numrel_man", "formation.hazard.agegapry.eagerness_diff",
                           "formation.hazard.agegapry.gap_factor_man_exp", "person.agegap.man.dist.normal.mu",
                           "person.agegap.woman.dist.normal.mu","person.agegap.man.dist.normal.sigma",
                           "person.agegap.woman.dist.normal.sigma", "person.eagerness.woman.dist.gamma.a",
                           "person.eagerness.woman.dist.gamma.b", "formation.hazard.agegapry.numrel_woman",
                           "formation.hazard.agegapry.gap_factor_woman_exp")
  
  ## Same here not needed, in case you do not need all the target statistics
  target.variables <-c("growth.rate", "median.AD", "Q1.AD", "Q3.AD", "prev.men.15.25", "prev.men.25.50",
                       "ART.cov.15.50", "incid.wom.15.30", "frac.degreeGT1.wom.15.30", "mean.degree",
                       "median.degree", "Q1.degree", "Q3.degree")
  
  simpact.set.simulation("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
  agedist.chunk.data.frame <- agedistr.creator(shape = 5, scale = 65)
  
  #### Set input params
  ##Specifying the initially chosen values for the simulation.
  cfg.chunk <- input.params.creator(population.simtime = 40, population.numwomen = 500, population.nummen = 500)
  
  #intervention introduced See the intervention.introduced
  # Simulation starts in 1977. After 27 years (in 2004), ART is introduced.
  iv.chunk <- intervention.introduced(list(27,0,100,2),list(30,200,1.5), list(33,350,1),list(36,500,0.5))
  
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
  
    #if(testoutput$simulationtime < cfg.chunk$population.simtime){chunk.datalist.test <- "Error"}else{chunk.datalist.test <- readthedata(testoutput)}
    # #{
    #   if (testoutput$eventsexecuted >= cfg.chunk$population.maxevents-1)  #use ifelse
    #   {
    #     stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
    #   }
    #   else
    #   {
    #     stop("Simulation stopped prematurely, probably ran out of events")
    #   }
    # }
  chunk.datalist.test <- readthedata(testoutput)
  
  
  
  if(length(chunk.datalist.test)>1){
    #get the summary statistics for each run
    out.statistics <- output.summary.maker(datalist = chunk.datalist.test, 
                                           growth.rate=list(timewindow.min = 0, timewindow.max = unique(chunk.datalist.test$itable$population.simtime)),
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
    out.statistic.no.degree <- out.statistics[,target.variables]
    ##out.test.degree <- out.statistic[[2]]
  }else{
    out.statistic.no.degree <- rep(NA,length(target.variables))
    ##out.statistic.degree <- NA
  }
  
  chunk.summary.stats <- out.statistic.no.degree
  return(chunk.summary.stats)
}