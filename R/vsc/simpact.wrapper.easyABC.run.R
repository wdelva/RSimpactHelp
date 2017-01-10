#!/usr/bin/env/ Rscript
#get the necessary libraries
pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper)
#data file to read
dirname <- getwd()
# There may be ways to make line 6 dynamic: so the file name does not need to be manually updated
main.filename <- "INPUT.df-1000Points10Par2016-11-04.csv" #Read the file produced by varying parameters *design.points
file.chunk.name.csv <-paste0(dirname, "/", main.filename) #### Input file name is produced from the .sh script
inPUT.df.complete <- read.csv(file = file.chunk.name.csv, header = TRUE, sep = ",")


#Select a chunk to send to process
min.chunk <- 1
max.chunk <- 1000
inANDout.df.chunk <- inPUT.df.complete[min.chunk:max.chunk,]

#make sure there are no empty rows
inANDout.df.chunk <- inANDout.df.chunk[!is.na(inANDout.df.chunk$sim.id),]

sim_repeat <- 2
ncluster.use <- 8 # number of cores per node

## In case you need more target statistics you can add here. The default are 13. Remember to change
## in the target.variables in the main file as well.
target.variables <- select.summary.params()[[1]]


#set the prior names - varied parameters
preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

#rbind all the results for this chunk to be merged after
#Create a dataframe with NA for the summary statistics Will collect all the chunks with the sim.id to link back
chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 0, ncol = length(target.variables)+1))
names(chunk.summary.stats.df) <- c(target.variables, "sim.id")


simpact4ABC.chunk.wrapper <- function(simpact.chunk.prior){
  pacman::p_load(RSimpactHelper)
  simpact.chunk.run <- function(input.chunk.params){

    pacman::p_load(RSimpactCyan, RSimpactHelper, dplyr,lhs,data.table, dplyr, magrittr, exactci,
                   nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
                   igraph,lhs, GGally, emulator, multivator, tidyr)


    input.varied.params.plus <- varied.simpact.params()[[1]]

    target.variables <- select.summary.params()[[1]]

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
      out.statistics <- output.summary.maker(datalist = chunk.datalist.test,
                                             growth.rate=list(timewindow=c(0, timewindow.max=unique(chunk.datalist.test$itable$population.simtime))),
                                             agemix.maker=list(agegroup=c(15,30), timepoint=30, timewindow=1, start=FALSE, gender="female"),
                                             prev.15.25=list(age.group=c(15,25), timepoint=35, gender="men"),
                                             prev.25.50=list(age.group=c(25,50), timepoint=35, gender="men"),
                                             art.coverage=list(age.group=c(15,50), timepoint=34, gender="men"),
                                             inc.15.30=list(age.group=c(15,30), timewindow=c(30,40), gender="women", only.active="Harling"),
                                             partner.degree=list(age.group=c(15,30), hivstatus=0, survey.time=30,
                                                                   window.width=1, gender="female", only.new=FALSE))
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

  chunk.summary.stats <- tryCatch(simpact.chunk.run(simpact.chunk.prior),
                                  error = err.function)
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


