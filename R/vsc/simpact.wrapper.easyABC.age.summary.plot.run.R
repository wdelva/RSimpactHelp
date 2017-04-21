#get the necessary libraries
pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper, lhs, gtools, ggplot2, gridExtra)
#data file to read

comp <- "lin" #lin #mac #cluster

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
    dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"}else if(comp=="gent"){
    dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"
    }else{dirname <- "~/Documents/RSimpactHelp"  #mac directory here
}

####  This will create the varied parameter space ###################################################
inPUT.df.complete <- simpact.config.inputs(design.points = 20000,
                                           conception.alpha_base = c(-5, -0), #c(-4, -1.5)
                                           person.art.accept.threshold.dist.fixed.value = c(0.65, 0.85), #good
                                           person.eagerness.man.dist.gamma.a = c(0.4, 1.3), #good
                                           person.eagerness.man.dist.gamma.b = c(30, 60), #c(10,50)
                                           person.eagerness.woman.dist.gamma.a = c(0.9, 1.4), #c(0.3, 1.5)
                                           formation.hazard.agegapry.eagerness_diff =  c(-0.02, -0.013), #c(-0.1, 0),
                                           person.eagerness.woman.dist.gamma.b = c(10, 45), #good
                                           formation.hazard.agegapry.numrel_man = c(-1.2, -0.6), #c(-2, -0.1)
                                           formation.hazard.agegapry.numrel_woman = c(-1.3, -1.0), #c(-2, -0.1),
                                           formation.hazard.agegapry.gap_factor_man_exp = c(-0.41, -0.2),#c(-2, -0.1),
                                           formation.hazard.agegapry.gap_factor_woman_exp = c(-0.4, -0.1),#c(-2, -0.1),
                                           person.agegap.man.dist.normal.mu = c(3, 5), #c(1, 5),
                                           person.agegap.woman.dist.normal.mu = c(3.2, 3.8), #c(1, 5),
                                           person.agegap.man.dist.normal.sigma = c(1.9, 2.9), #c(0.5, 3.5),
                                           person.agegap.woman.dist.normal.sigma = c(1.7, 2.5), #c(0.5, 3.5),
                                           hivtransmission.param.f1 = c(log(2), log(3.5))
                                           )

################################################################################################
#use when you have a file that was generated from wrapper with.replace.run
file.name.csv <- paste0(dirname,"/","inputANDoutput.SSE.df.csv")
complete.results <- data.frame(read.csv(file = file.name.csv, header = TRUE))

#### WIM Course PLACE holder
complete.results <- complete.results[complete.results$sim.id==6994,]
complete.results <- aggregate(complete.results,
                                    by = list(complete.results$sim.id), FUN = "mean", na.rm = TRUE)
#remove the name produced by mean grouping
complete.results <- subset(complete.results, select=-c(Group.1))

complete.results$match <- TRUE
complete.results$sum.square.df <- 0.006369786
##################

################ MICE METHOD  #######
complete.results <- mice.df.out[,15:26]
xdesign <- data.frame(matrix(NA, nrow = nrow(complete.results), ncol = ncol(complete.results) + 1))
names(xdesign) <- c("sim.id", paste0("xdesign",1:(length(xdesign)-1)))
xdesign$sim.id <- 1:nrow(complete.results)
complete.results <- cbind(xdesign, complete.results)

################
complete.results$match <- TRUE
complete.results$sum.square.df <- 0.006369786
#complete.results <- complete.results[order(complete.results$sum.square.df),]
match.true <- subset(complete.results, complete.results$match == TRUE)
inPUT.df.complete <- simpact.config.inputs.par.select(datalist = match.true)
################################################################################################

#Select a chunk to send to process
min.chunk <- 1
max.chunk <- 100

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
##Each of these should be calculated after each run, else we give an NA
target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
                      "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
                      "ART.cov.wom.18.50", "median.wom.18.50.AD")

#these are the summary we want to plot
#IF ANY change here, you need to change the calculation step as well
plot.summary.names <- c("prev.overalm.pt","prev.overalm.ll","prev.overalm.ul",
                        "inci.overalm.pt","inci.overalm.ll","inci.overalm.ul",
                        "prev.overalw.pt","prev.overalw.ll","prev.overalw.ul",
                        "inci.overalw.pt","inci.overalw.ll","inci.overalw.ul",
                        "prev.overalt.pt","prev.overalt.ll","prev.overalt.ul",
                        "inci.overalt.pt","inci.overalt.ll","inci.overalt.ul"
                        )
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
    unlink(paste0("temp/",sub.dir.sim.id), recursive = TRUE)

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

      #Process the summary stats for different year points for ploting.
      start.year <- chunk.datalist.test$itable$hivseed.time[1]
      end.year <- chunk.datalist.test$itable$population.simtime[1]
      year.time.int <- start.year:end.year

      plot.sum.res <- c() #collect all the results from the wide summary calc
      plot.age.group <- list(c(18,20), c(20,25), c(25,30), c(30,35), c(35,40), c(40,45),c(45,50))

      for(time.point in year.time.int){
        for(age in plot.age.group){
		     prev.sum.stat <- prevalence.calculator(datalist = chunk.datalist.test, agegroup = age,
		                                             timepoint = time.point)

		     inc.row.index <- time.point + 3#calc within three years
		     inc.sum.stat <- incidence.calculator(datalist = chunk.datalist.test, agegroup = age,
		                                           timewindow = c(time.point, inc.row.index),
		                                           only.active = "No")

		     #prev.loop.result <- c()
		     #inci.loop.result <- c()
		     #3=overal, 2=woman, 1=man (0 <- Men, 1 <- Women SIMPACT)
		     for(k in 1:3){
		       prev.loop.result <- c(as.numeric(prev.sum.stat$pointprevalence[k]),
		                             as.numeric(prev.sum.stat$pointprevalence.95.ll[k]),
		                             as.numeric(prev.sum.stat$pointprevalence.95.ul[k]))

		       inci.loop.result <- c(as.numeric(inc.sum.stat$incidence[k]),
		                             as.numeric(inc.sum.stat$incidence.95.ll[k]),
		                             as.numeric(inc.sum.stat$incidence.95.ul[k]))
		       #join the results from start to end time
		       plot.sum.res <- c(plot.sum.res, prev.loop.result, inci.loop.result)
		       }
		    }

      }

      out.statistic <- c(out.statistic, plot.sum.res)

    }else{
      plot.sum.res <- rep(NA,18*7*length(year.time.int)) #point+CI @time step & age grp
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

  #get the year start = 1 to end to be ploted against
  plot.age.group <- list(c(18,201), c(20,252), c(25,303), c(30,354), c(35,405), c(40,456),c(45,507))
  plot.res.size <- length(ABC.results.chunk.statistics) -  length(target.variables)
  plot.sum.rep <- plot.res.size/(7 * 18)

  year.time.int <- 1:plot.sum.rep

  age.plot.names <- c()
  for(i.time in year.time.int){

    for(z in plot.age.group){

      for(sum.nam in plot.summary.names){
        #print(z)
        min.age <- z[1]
        max.age <- z[2]
        age.plot.names <- c(age.plot.names, paste0(sum.nam, min.age,max.age, sprintf("%02d",i.time)))
      }

    }
  }
  #plot.index.names <- rep(sprintf("%02d",year.time.int), each = 7*18, times = length(plot.sum.rep))
  #name of the plot with an indication of the run
  #plot.colnames.summary <- paste0(age.plot.names, plot.index.names)

  target.plot.target.names <- c(target.variables, age.plot.names )

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

#Add the summarries to the simpact parameters
inputANDoutput.chunk.df  <- left_join(chunk.summary.stats.df, inANDout.df.chunk, by = "sim.id")

#Calculate how much time the simulation took
end.chunk.time <- proc.time() - start.chunk.time

rand.string <- paste0(sample(c(LETTERS,letters), 5), collapse="")

write.csv(inputANDoutput.chunk.df, file = paste0(dirname,"/","SummaryOutPut-inANDout.df.chunk-Age-WIMCourse",
                                                 rand.string,".csv"),
                                                row.names = FALSE)

#Read the file for summary plot UNCOMMENT ONCE TESTING IS DONE
#file.name.csv <- paste0(dirname,"/","SummaryOutPut.df.ReSample-iPUwps-1.csv")
#inputANDoutput.chunk.plot.df <- data.frame(read.csv(file = file.name.csv, header = TRUE))


inputANDoutput.chunk.plot.df <- inputANDoutput.chunk.df
inputANDoutput.chunk.plot.df <- subset(inputANDoutput.chunk.df, sim.id == 59)


############ Doing the mean comment below if needed  ######################################
inputANDoutput.chunk.plot.df <- aggregate(inputANDoutput.chunk.plot.df,
                               by = list(inputANDoutput.chunk.plot.df$sim.id), FUN = "mean")

inputANDoutput.chunk.plot.df <- subset(inputANDoutput.chunk.plot.df, select=-c(Group.1))
##################################################################################################

#Set the targets for the summary statistics.
target.stats <- c(0.015, 0.016, 0.043, 0.21, 0.47, 0.37, 0.54, 0.33, 0.34, 5)

mean.sum.square.df <- inputANDoutput.chunk.plot.df[,1:10]
#compute mean sum of difference squared
mean.sse <- rowSums(as.data.frame(t(apply(mean.sum.square.df[,4:7], 1,
                                          function(x) (((x - t(target.stats[4:7]) )^2)  ) ) ) ) )
mean.sse <- as.data.frame(mean.sse)
##################################################################################################


# Prepare and generate inc and prev plots
prev.inci.select.sum <- dplyr::select(inputANDoutput.chunk.plot.df, contains(".overal"))
sim.id.select <- dplyr::select(inputANDoutput.chunk.plot.df, contains("sim"))
prev.inci.select.sum <- cbind(prev.inci.select.sum, sim.id.select)

time.sim.id <- length(unique(prev.inci.select.sum$sim.id))
rep.sim.id <- sum(prev.inci.select.sum$sim.id==unique(prev.inci.select.sum$sim.id)[1])

## create a unique sim.id from repeated sim id's
prev.inci.select.sum <- prev.inci.select.sum %>%
  mutate(sim.id.unique = paste0(sim.id, rep(letters[1:rep.sim.id], times=time.sim.id)))

## wide to long format
prev.inci.select.sum <- prev.inci.select.sum %>%
  tidyr::gather(point.est, prev.inc, 1:(length(prev.inci.select.sum)-2))

#age.group vector
prev.inci.select.sum <- prev.inci.select.sum %>%
  mutate(age.group.point = as.numeric(substr(point.est, nchar(point.est)-2, nchar(point.est)-2)))

#age group vector label
prev.inci.select.sum <- prev.inci.select.sum %>%
  mutate(age.group.label = as.numeric(substr(point.est, nchar(point.est)-6, nchar(point.est)-3)))

#time vector
prev.inci.select.sum <- prev.inci.select.sum %>%
  mutate(year.time = as.numeric(substr(point.est, nchar(point.est)-1, nchar(point.est))))

#indicate the plot type
prev.inci.select.sum <- prev.inci.select.sum %>%
  mutate(plot.type = substr(point.est, 1, 4))

#change observation to real %
prev.inci.select.sum$prev.inc <- round(100* prev.inci.select.sum$prev.inc, 2)

#gender
prev.inci.select.sum$Gender <- NA
prev.inci.select.sum$Gender[grep(".overalt",prev.inci.select.sum$point.est) ] <- "Total"
prev.inci.select.sum$Gender[grep(".overalm",prev.inci.select.sum$point.est) ] <- "Men"
prev.inci.select.sum$Gender[grep(".overalw",prev.inci.select.sum$point.est) ] <- "Women"

#point, lci, uci
prev.inci.select.sum$sum.type <- NA
prev.inci.select.sum$sum.type[grep(".pt",prev.inci.select.sum$point.est)] <- "Point"
prev.inci.select.sum$sum.type[grep(".ll",prev.inci.select.sum$point.est)] <- "Lower"
prev.inci.select.sum$sum.type[grep(".ul",prev.inci.select.sum$point.est)] <- "Upper"

#Select year to plot the age distribution start year 1977 <- 10 - 40 yrs
# then 2011 is year 33
selected.year <- 25
target.plot.list <- c("Women")

#prev select for plot # Gender != "Total"  ## sum.type == "Point",
prev.select.sum.plot <- dplyr::filter(prev.inci.select.sum, plot.type == "prev",
                                      sum.type == "Point", year.time == selected.year,
                                      Gender  %in% target.plot.list)
#inci select for plot
inci.select.sum.plot <- dplyr::filter(prev.inci.select.sum, plot.type=="inci",
                                      sum.type == "Point", year.time == selected.year,
                                      Gender  %in% target.plot.list)

if(selected.year == 25){
  shims.prev.men <- list(c(rep(prev.select.sum.plot$sim.id[1],7)),
                     c(rep(prev.select.sum.plot$sim.id.unique[1],7)),
                     c(rep("prev.age.n",7)),
                     c(0.8, 6.6, 21.3, 36.6,47,45.5,42.5), c(1:7),
                     c(1820, 2025, 2530, 3035, 3540, 4045, 4550),
                     c(rep(25, 7)), c(rep("prev",7)), c(rep("SHIMS-Men", 7)),
                     c(rep("Point",7)))

  shims.prev.men <- as.data.frame(shims.prev.men)
  names(shims.prev.men) <- names(prev.select.sum.plot)

  shims.prev.women <- list(c(rep(prev.select.sum.plot$sim.id[1],7)),
                         c(rep(prev.select.sum.plot$sim.id.unique[1],7)),
                         c(rep("prev.age.n",7)),
                         c(14.3, 31.5, 46.7, 53.8, 49.1,39.7, 31.6), c(1:7),
                         c(1820, 2025, 2530, 3035, 3540, 4045, 4550),
                         c(rep(25, 7)), c(rep("prev",7)), c(rep("SHIMS-Wom", 7)),
                         c(rep("Point",7)))

  shims.prev.women <- as.data.frame(shims.prev.women)
  names(shims.prev.women) <- names(prev.select.sum.plot)

  prev.select.sum.plot <- rbind(prev.select.sum.plot, shims.prev.women, shims.prev.men)

  shims.inc.men <- list(c(rep(inci.select.sum.plot$sim.id[1],7)),
                        c(rep(inci.select.sum.plot$sim.id.unique[1],7)),
                        c(rep("inci.age.n",7)),
                        c(0.8, 1.6, 2.6, 3.1, 0.4, 1.2, NA), c(1:7),
                        c(1820, 2025, 2530, 3035, 3540, 4045, 4550),
                        c(rep(25, 7)), c(rep("inci",7)), c(rep("SHIMS-Men", 7)),
                        c(rep("Point",7)))

  shims.inc.men <- as.data.frame(shims.inc.men)
  names(shims.inc.men) <- names(inci.select.sum.plot)


  shims.inc.women <- list(c(rep(inci.select.sum.plot$sim.id[1],7)),
                        c(rep(inci.select.sum.plot$sim.id.unique[1],7)),
                        c(rep("inci.age.n",7)),
                        c(3.8, 4.3, 2.0, 2.7, 4.0, 2.1, 1.2), c(1:7),
                        c(1820, 2025, 2530, 3035, 3540, 4045, 4550),
                        c(rep(25, 7)), c(rep("inci",7)), c(rep("SHIMS-Wom", 7)),
                        c(rep("Point",7)))

  shims.inc.women <- as.data.frame(shims.inc.women)
  names(shims.inc.women) <- names(inci.select.sum.plot)

  inci.select.sum.plot <- rbind(inci.select.sum.plot, shims.inc.women, shims.inc.men)
}



year.plot.selection <- 1976 + selected.year + 10 #start of the simulation plus yrs after start

inc.plot <- ggplot(inci.select.sum.plot,
                   aes(x=age.group.point, y=prev.inc, group=interaction(sim.id.unique, Gender))) +
  geom_point() + geom_line() + aes(colour = Gender)  +
  xlim("18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49") +
  #facet_grid(sim.id~plot.type) + xlab(paste0(year.plot.selection," Age Group")) +
  #ylab(" HIV Incidence") +
  labs(title = paste0(year.plot.selection, " Model HIV Incidence in Swaziland"),
       x = "Age Group", y = " HIV Incidence") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x  = element_text(vjust=0.5, size=10),
        axis.title.x = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=10),
        axis.title.y = element_text(size=16))


prev.plot <- ggplot(prev.select.sum.plot,
                    aes(x=age.group.point, y=prev.inc, group=interaction(sim.id.unique, Gender))) +
  geom_point() + geom_line() + aes(colour = Gender) +
  geom_text(aes(label=prev.inc), size=4, vjust=-.8, show.legend = FALSE) +
  xlim("18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49") +
  ylim(0,65) +
  #facet_grid(sim.id~plot.type) +
  #xlab(paste0(year.plot.selection," Age Group")) +
  #ylab(" HIV Prevalence") +
  labs(title = paste0(year.plot.selection, " Model HIV Prevalence in Swaziland"),
       x = "Age Group", y = " HIV Prevalence") +
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x  = element_text(vjust=0.5, size=10),
        axis.title.x = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=10),
        axis.title.y = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5, size=16))
  #+ theme(legend.position="none")


#have prev and inc on the same plot.
grid.arrange(inc.plot, prev.plot, ncol = 1 )

# Use semi-transparent fill
inci.bar.plot <- ggplot(inci.select.sum.plot, aes(x=age.group.point, y=prev.inc, fill=Gender, color=Gender))

inci.bar.plot <- inci.bar.plot + geom_bar(stat = "identity", position = position_dodge(0.9))

#geom_errorbar(aes(ymax = ), position = position_dodge(0.9),
#              width = 0.25)

inci.bar.plot <- inci.bar.plot + xlim("18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49") +
  ylim(0,6) + labs(title = paste0(year.plot.selection, " Model HIV Incidence in Swaziland"),
       x = "Age Group", y = " HIV Incidence") +
  theme_bw()+
  theme( panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x  = element_text(vjust=0.5, size=10),
        axis.title.x = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=10),
        axis.title.y = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5, size=16)) +
  theme(
    legend.position = c(0.99, 0.9),
    legend.justification = c("right", "top")
  )


#PrevHist
target.plot.hist <- c("Men", "SHIMS-Men")

#prev select for plot # Gender != "Total"  ## sum.type == "Point",
prev.select.hist.plot <- dplyr::filter(prev.select.sum.plot,
                                      sum.type == "Point",
                                      Gender  %in% target.plot.hist)

# Use semi-transparent fill
prev.bar.plot <- ggplot(prev.select.hist.plot, aes(x=age.group.point, y=prev.inc, fill=Gender, color=Gender))

prev.bar.plot <- prev.bar.plot + geom_bar(stat = "identity", position = position_dodge(0.9))

prev.bar.plot <- prev.bar.plot + xlim("18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49") +
  labs(title = paste0(year.plot.selection, " Model HIV Prevalence in Swaziland"),
                   x = "Age Group", y = " HIV Prevalence") +
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x  = element_text(vjust=0.5, size=10),
        axis.title.x = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=10),
        axis.title.y = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5, size=16)) +
  theme(
    legend.position = c(0.9, 0.3),
    legend.justification = c("right", "top")
  )





#one with confident interval
#geom_ribbon(data = prev.inci.select.sum, aes(ymin=Lower, ymax=upper), alpha=0.3)
#


