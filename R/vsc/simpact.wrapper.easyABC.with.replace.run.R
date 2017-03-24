#!/usr/bin/env/ Rscript

#get the necessary libraries
pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper, lhs)
#data file to read

comp <- "lin" #lin #mac #chpc #gent

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
    dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"}else if(comp=="gent"){
    dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"}else{
    dirname <- "~/Documents/RSimpactHelp"  #mac directory here
}

all.sim.start <- proc.time()

###### Generate the input parameters for the simulation ###############################################
inPUT.df.complete <- simpact.config.inputs(design.points = 10, resample.count = 1,
                                           conception.alpha_base = c(-5, -1), #c(-4, -1.5)
                                           person.art.accept.threshold.dist.fixed.value = c(0.65, 0.99), #good
                                           person.eagerness.man.dist.gamma.a = c(0.3, 1.5), #good
                                           person.eagerness.man.dist.gamma.b = c(10, 50), #c(10,50)
                                           person.eagerness.woman.dist.gamma.a = c(0.2, 1.2), #c(0.3, 1.5)
                                           formation.hazard.agegapry.eagerness_diff =  c(-0.02, 0), #c(-0.1, 0),
                                           person.eagerness.woman.dist.gamma.b = c(10, 50), #good
                                           formation.hazard.agegapry.numrel_man = c(-1.3, -0.6), #c(-2, -0.1)
                                           formation.hazard.agegapry.numrel_woman = c(-1.3, -0.6), #c(-2, -0.1),
                                           formation.hazard.agegapry.gap_factor_man_exp = c(-0.4, -0.1),#c(-2, -0.1),
                                           formation.hazard.agegapry.gap_factor_woman_exp = c(-1.5, -0.1),#c(-2, -0.1),
                                           person.agegap.man.dist.normal.mu = c(1, 4), #c(1, 5),
                                           person.agegap.woman.dist.normal.mu = c(2, 3.5), #c(1, 5),
                                           person.agegap.man.dist.normal.sigma = c(2.5, 3.5), #c(0.5, 3.5),
                                           person.agegap.woman.dist.normal.sigma = c(1.5, 2.6), #c(0.5, 3.5),
                                           hivtransmission.param.f1 = c(log(2), log(3.5))
                                           )
#
#######################################################################################################

#initial set up with no replacement
count.match <- 0

rep.sample.number <- nrow(inPUT.df.complete)

set.new.seed <- 1

while(count.match < 150){
  #Select a chunk to send to process
  min.chunk <- 1
  max.chunk <- 10

  if(count.match != 0){

    inPUT.df.complete <- simpact.config.inputs.add.sample(datalist = inputANDoutput.chunk.df,
                                        resample.points = rep.sample.number, set.seed.new = set.new.seed,
                                        conception.alpha_base = c(-5, -0.1), #c(-4, -1.5)
                                        person.art.accept.threshold.dist.fixed.value = c(0.65, 0.99), #good
                                        person.eagerness.man.dist.gamma.a = c(0.3, 1.5), #good
                                        person.eagerness.man.dist.gamma.b = c(10, 60), #c(10,50)
                                        person.eagerness.woman.dist.gamma.a = c(0.8, 1.5), #c(0.3, 1.5)
                                        formation.hazard.agegapry.eagerness_diff = c(-0.02, -0.013), #c(-0.1, 0),
                                        person.eagerness.woman.dist.gamma.b = c(10, 60), #good
                                        formation.hazard.agegapry.numrel_man = c(-1.3, -0.6), #c(-2, -0.1)
                                        formation.hazard.agegapry.numrel_woman = c(-1.3, -1.0), #c(-2, -0.1),
                                        formation.hazard.agegapry.gap_factor_man_exp = c(-0.4, -0.1),#c(-2, -0.1),
                                        formation.hazard.agegapry.gap_factor_woman_exp = c(-0.4, -0.01),#c(-2, -0.1),
                                        person.agegap.man.dist.normal.mu = c(2.7, 5), #c(1, 5),
                                        person.agegap.woman.dist.normal.mu = c(3.2, 4), #c(1, 5),
                                        person.agegap.man.dist.normal.sigma = c(1.5, 2.9), #c(0.5, 3.5),
                                        person.agegap.woman.dist.normal.sigma = c(1.5, 2.6), #c(0.5, 3.5),
                                        hivtransmission.param.f1 = c(log(2), log(3.5))
                                        )

  }


  if(max.chunk > nrow(inPUT.df.complete)){max.chunk <- nrow(inPUT.df.complete)}
  if(min.chunk > nrow(inPUT.df.complete) || min.chunk < 1){min.chunk <- max.chunk}

  inANDout.df.chunk <- inPUT.df.complete[min.chunk:max.chunk,]

  #make sure there are no empty rows
  inANDout.df.chunk <- inANDout.df.chunk[!is.na(inANDout.df.chunk$sim.id),]

  #set how many time the single row will be repeated
  sim_repeat <- 4

  # number of cores per node
  ncluster.use <- 4

  #indicate the target statitics that you want to hit
  target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
                        "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
                        "ART.cov.wom.18.50", "median.wom.18.50.AD")

  ##Each of these should be calculated after each run, else we give an NA

  #set the prior names - varied parameters
  preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
  preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

  #rbind all the results for this chunk to be merged after
  #Create a dataframe with NA for the summary statistics Will collect all the chunks with the sim.id to link back
  chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 0, ncol = length(target.variables)+2))
  names(chunk.summary.stats.df) <- c(target.variables, "sim.id")

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
      input.varied.params.plus <- c("conception.alpha_base", "person.art.accept.threshold.dist.fixed.value",
                                    "person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b",
                                    "person.eagerness.woman.dist.gamma.a", "formation.hazard.agegapry.eagerness_diff",
                                    "person.eagerness.woman.dist.gamma.b", "formation.hazard.agegapry.numrel_man",
                                    "formation.hazard.agegapry.numrel_woman", "formation.hazard.agegapry.gap_factor_man_exp",
                                    "formation.hazard.agegapry.gap_factor_woman_exp", "person.agegap.man.dist.normal.mu",
                                    "person.agegap.woman.dist.normal.mu", "person.agegap.man.dist.normal.sigma",
                                    "person.agegap.woman.dist.normal.sigma", "hivtransmission.param.f1")

      target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
                            "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
                            "ART.cov.wom.18.50", "median.wom.18.50.AD")

      simulation.type <- ("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
      simpact.set.simulation(simulation.type)
      agedist.chunk.data.frame <- agedistr.creator(shape = 5, scale = 65)

      #### Set input params
      ##Specifying the initially chosen values for the simulation.
      cfg.chunk <- input.params.creator(population.simtime = 40, population.numwomen = 1000, population.nummen = 1000,
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
        #setting up a value that is depended on the other input (we can do this for many other as needed)
        if(cfg.chunk.par == "hivtransmission.param.f1"){
          f2.num <- log((1+assign.chunk.cfg.value)/2)
          f2.den <- log(assign.chunk.cfg.value)
          cfg.chunk["hivtransmission.param.f2"][[1]] <- log(f2.num / f2.den)/5
        }
      }

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
      #save(chunk.datalist.test, file = paste0("temp/","chunk.datalist.",sub.dir.sim.id,".rda"))

      #delete all the file created during the current simulation
      unlink(paste0("temp/",sub.dir.sim.id), recursive = TRUE, force = TRUE)

      if(length(chunk.datalist.test)>1){
        #get the summary statistics for each run
        end.time.wind <- unique(chunk.datalist.test$itable$population.simtime)
        growth.rate <- pop.growth.calculator(datalist = chunk.datalist.test,
                                            timewindow = c(0, timewindow.max = end.time.wind))

        inc.20.25 <- incidence.calculator(datalist = chunk.datalist.test, agegroup = c(20, 25),
                                          timewindow = c(32, 35), only.active = "No")
        inc.men.20.25 <- inc.20.25$incidence[1]
        inc.wom.20.25 <- inc.20.25$incidence[2]
        prev.25.30 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(25, 30),
                                           timepoint = 35)
        prev.men.25.30 = prev.25.30$pointprevalence[1]
        prev.wom.25.30 = prev.25.30$pointprevalence[2]
        prev.30.35 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(30, 35),
                                           timepoint = 35)
        prev.men.30.35 = prev.30.35$pointprevalence[1]
        prev.wom.30.35 = prev.30.35$pointprevalence[2]
        ARTcov <- ART.coverage.calculator(datalist = chunk.datalist.test, agegroup = c(18, 50),
                                          timepoint = 35, site="All")
        ART.cov.men.18.50 <- ARTcov$ART.coverage[1]
        ART.cov.wom.18.50 <- ARTcov$ART.coverage[2]

        agemix.df <- agemix.df.maker(chunk.datalist.test)
        pattern <- pattern.modeller(dataframe = agemix.df, agegroup = c(18, 50),
                                    timepoint = 35, timewindow = 1, start = FALSE)
        median.wom.18.50.AD <- as.numeric(median(pattern[[1]]$AgeGap[pattern[[1]]$Gender == "female"]))

        ##get the summary statistics as indicated by target.variables
        out.statistic <- c(growth.rate, inc.men.20.25, inc.wom.20.25, prev.men.25.30, prev.wom.25.30,
                           prev.men.30.35, prev.wom.30.35, ART.cov.men.18.50,
                           ART.cov.wom.18.50, median.wom.18.50.AD)
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
                                          seed_count = 0,
                                          n_cluster = ncluster.use)


    #Save the statistics results with the chunk row sim.id repeated X* from the ABC_rejection method
    ABC.results.chunk.statistics <- data.frame(ABC.chunk.result$stats)
    names(ABC.results.chunk.statistics) <- target.variables
    ABC.results.chunk.statistics$sim.id <- chunk.sim.id

    chunk.summary.stats.df <- rbind(chunk.summary.stats.df, ABC.results.chunk.statistics)

  }

  #target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
  #                      "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
  #                      "ART.cov.wom.18.50", "median.wom.18.50.AD")
  target.stats <- c(0.015, 0.016, 0.043, 0.21, 0.47, 0.37, 0.54, 0.33, 0.34, 5)


  #chunk.summary.stats.df <- aggregate(chunk.summary.stats.df,
  #                               by = list(chunk.summary.stats.df$sim.id), FUN = "mean", na.rm = TRUE)

  chunk.summary.stats.df <- subset(chunk.summary.stats.df, select = -c(Group.1))

  sum.square.df <- subset(chunk.summary.stats.df, select=-c(sim.id))

  #compute sum differnce of squares
  sum.square.df <- rowSums(as.data.frame(t(apply(sum.square.df, 1, function(x) (x - t(target.stats))^2))))
  sum.square.df <- as.data.frame(sum.square.df)

  #Discard Low sum difference of squares due to zeros in summary of results
  less.zeros.df <- rowSums(chunk.summary.stats.df==0)  < 2

  chunk.summary.stats.df <- cbind(chunk.summary.stats.df, sum.square.df)

  chunk.summary.stats.df$match <- chunk.summary.stats.df$sum.square.df <= 0.25 & less.zeros.df


  inputANDoutput.chunk.df  <- left_join(chunk.summary.stats.df, inANDout.df.chunk, by = "sim.id")


  if(count.match ==0) {
    rand.string <- paste0(sample(c(LETTERS,letters), 6), collapse = "")

    filename.run <- paste0(dirname,"/","SummaryOutPut.df.ReSample-",
                           rand.string,"-",set.new.seed,".csv")

    write.csv(inputANDoutput.chunk.df, file = filename.run, row.names = FALSE)

  }else{

    append.results <- read.csv(file = filename.run, header = TRUE, sep = ",")

    inputANDoutput.chunk.df <- rbind(append.results, inputANDoutput.chunk.df)

    write.csv(inputANDoutput.chunk.df, file = filename.run, row.names = FALSE)
  }


  if(sum(chunk.summary.stats.df$match, na.rm = TRUE) == 0 & count.match == 0){
    count.match <- count.match + 1
  }else{
    count.match <- count.match + sum(chunk.summary.stats.df$match, na.rm = TRUE)
  }

}

end.chunk.time <- proc.time() - start.chunk.time

all.sim.end <- proc.time() - all.sim.start

#use when you have a file that was generated from wrapper with.replace.run
file.name.csv <- paste0(dirname,"/","SummaryOutPut.df.ReSample-iPUwps-pairs.csv")
complete.results.res <- data.frame(read.csv(file = file.name.csv, header = TRUE))


#Analysis of results
#complete.results.res <- subset(complete.results.res, complete.results.res$match == TRUE)

start.pair <- length(complete.results.res)-29
end.pair <- length(complete.results.res)-31

#use this to check the range of the varied parameters
pairs(complete.results.res[start.pair:end.pair],  col = 1+complete.results.res$match, pch = 16, cex = 2)












