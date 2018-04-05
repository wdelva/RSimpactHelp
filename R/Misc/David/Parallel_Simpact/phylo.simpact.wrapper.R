
phylo.simpact.wrapper <- function(inputvector = inputvector){
  
  
  # ## For sub-optimal sequence coverage
  
  # Kick-start the missing sequence data simulation study:
  
  # Let's start with verifying that we have parameter combinations that produce "useful" output.
  # i.e. output that is more or less in line with the demography, behaviour and epidemic that we want to consider.
  library(RSimpactCyan)
  library(RSimpactHelper)
  
  library(devtools)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(readr)
  library(phangorn)
  library(dplyr)
  library(adephylo)
  library(treedater)
  library(geiger)
  library(picante)
  library(igraph)
  library(ggplot2)
  
  library(phyloTop) # these may be new to install
  library(phytools)
  library(lme4)
  library(data.table)
  
  destDir <- "/home/david/Dropbox/1.Test.Master.Model.2018/Parallel_Simpact/temp" # laptop
  
  
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                   population.msm = "no",
                                   population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
                                   population.nummen = 600, #3000, #600, # 3800, #2500,
                                   population.numwomen = 600, # 3000, #600, #4200, #2500,
                                   hivseed.time = 20, # 10,
                                   hivseed.type = "amount",
                                   hivseed.amount = 20, #30,
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   hivtransmission.param.a = -1, # -1,
                                   hivtransmission.param.b = -90,
                                   hivtransmission.param.c = 0.5,
                                   hivtransmission.param.f1 = log(2), #log(inputvector[2]) , #log(2),
                                   hivtransmission.param.f2 = log(log(1.4) / log(2)) / 5, #log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5, #log(log(1.4) / log(2)) / 5,
                                   formation.hazard.agegapry.gap_factor_man_age = -0.01, #-0.01472653928518528523251061,
                                   formation.hazard.agegapry.gap_factor_woman_age = -0.01, #-0.0726539285185285232510561,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   formation.hazard.agegapry.gap_factor_man_const = 0,
                                   formation.hazard.agegapry.gap_factor_woman_const = 0,
                                   formation.hazard.agegapry.gap_factor_man_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_factor_woman_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_agescale_man = 0.25, #inputvector[3], # 0.25,
                                   formation.hazard.agegapry.gap_agescale_woman = 0.25, #inputvector[3], # 0.25,#-0.30000007,#-0.03,
                                   debut.debutage = 15,
                                   conception.alpha_base = -2.5#inputvector[14]#-2.5#,
                                   #person.art.accept.threshold.dist.fixed.value = 0
  )
  
  
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3 
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  
  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.4
  cfg.list["diagnosis.baseline"] <- -2
  
  
  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- -2 # 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  
  ### add something about diagnosis
  art.intro["diagnosis.agefactor"] <- 0
  art.intro["diagnosis.genderfactor"] <- 0
  art.intro["diagnosis.diagpartnersfactor"] <- 0
  art.intro["diagnosis.isdiagnosedfactor"] <- 0
  ### end of add-on about diagnosis
  
  
  
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
  
  
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- -2 # 0#100
  art.intro1["monitoring.cd4.threshold"] <- 150 # 1200
  
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 36
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
  
  interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist # scenario(interventionlist, tasp.indicator)
  
  
  
  cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
  cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[3]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[3]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]
  
  
  
  
  cfg <- cfg.list
  
  cfg["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  # cfg["monitoring.fraction.log_viralload"] <- 0.3
  cfg["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine
  
  seedid <- inputvector[1]
  #cfg["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
  #cfg["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
  cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10] ######### -0.5
  cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10] ######### -0.5
  cfg["formation.hazard.agegapry.baseline"] <- inputvector[11]
  
  cfg["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
  cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
  cfg["conception.alpha_base"] <- inputvector[14] #is conception.alpha.base (higher up)
  cfg["dissolution.alpha_0"] <- inputvector[15]
  cfg["dissolution.alpha_4"] <- inputvector[16]
  
  
  #
  # # # # Run Simpact
  results <- simpact.run(configParams = cfg,
                         destDir = "temp",
                         agedist = age.distr,
                         seed = seedid,
                         intervention = intervention)
  
  datalist <- readthedata(results)
  
  
  # Let's examine the output:
  # source("/home/david/Dropbox/Stress_Testing/00-Functions.R") # Laptop
  
  source("/home/david/Dropbox/Stress_Testing/00-Functions.R")
  
  # source("/user/data/gent/vsc400/vsc40070/phylo/00-Functions.R")
  
  ###
  # Firstly, let's look at the age mixing pattern
  ###
  
  datalist.agemix <- readthedata(results)
  
  agemix.df <- agemix.df.maker(datalist.agemix)
  
  agemix.model <- pattern.modeller(dataframe = agemix.df,
                                   agegroup = c(15, 50),
                                   timepoint = datalist.agemix$itable$population.simtime[1],
                                   timewindow = 3)#1)#3)
  
  # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
  #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  
  men.lmer <- tryCatch(ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
                       error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  
  bignumber <- NA # let's try if NA works (instead of 9999 for example)
  
  AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
  slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
  WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
  
  BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
  
  intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
  
  
  age.scatter.df <- agemix.model[[1]]
  
  # # Scatter plot of age-mixing pattern
  # ggplot(data = dplyr::filter(age.scatter.df, Gender == "male", agerelform >=18, agerelform < 50),
  #        aes(x = agerelform, y = pagerelform)) +
  #   geom_point(alpha = 0.25) +
  #   geom_abline(size = 1,
  #               aes(intercept = 0, slope = 1, linetype = "Same age"),
  #               show.legend = FALSE) +
  #   #facet_grid(. ~ Gender) +
  #   scale_y_continuous(name = "Age of woman") +
  #   # scale_linetype_manual('Lines',
  #   #                        values = c("Population mean" = 1, "Same age" = 2)) +
  #   xlab("Age of man") +
  #   guides(linetype = guide_legend(keywidth = 2, keyheight = 1)) +
  #   coord_fixed(ratio = 1) #+ theme
  # 
  # # Histogram of age differences
  # ggplot(data = dplyr::filter(age.scatter.df, Gender == "male", agerelform >=18, agerelform < 50),
  #        aes(AgeGap, ..count..)) +
  #   geom_histogram(bins = 20) +
  #   xlab("Age difference") +
  #   ylab("Count")
  
  
  
  ###
  # Now we fit the negative binomial distribution
  ###
  degrees.df <- degree.df.maker(dataframe.df = agemix.df,
                                agegroup = c(18, 50),
                                hivstatus = 2,
                                survey.time = datalist.agemix$itable$population.simtime[1],
                                window.width = 1,
                                gender.degree = "male",
                                only.new = TRUE)
  
  
  # If we want to know % of men who had 0 partners in the past 12 months,
  # We need to compare nrow(degree.df) with the number of 18-50 year old men #HIV negative women
  # that were alive at the time of the survey
  allsurveymen <- dplyr::filter(datalist.agemix$ptable,
                                Gender == 1, # Male
                                TOD > datalist.agemix$itable$population.simtime[1], # Still alive at the time of the survey
                                TOB <= datalist.agemix$itable$population.simtime[1] - 18, # Not younger than 18 at the time of the survey
                                TOB > datalist.agemix$itable$population.simtime[1] - 50)#, # Not older than 50 at the time of the survey
  #InfectTime > datalist.agemix$itable$population.simtime[1]) # HIV negative at the time of the survey
  
  # Creating the vector of degrees
  degree.vector <- c(rep(0, times = (nrow(allsurveymen) - nrow(degrees.df))), degrees.df$Degree)
  meandegree.male <- mean(degree.vector)
  # hist(degree.vector, 10)
  # table(degree.vector)
  # degree.vector <- rnbinom(n = 100, size =1.28, mu = 0.66) # PLACEHOLDER FOR NOW
  # Fitting the negative binomial distribution to this vector
  
  fit.negbin <- tryCatch(fitdist(degree.vector, "nbinom"), error = agemixing.lme.errFunction)
  shape.nb.male <- ifelse(length(fit.negbin) > 0, as.numeric(fit.negbin$estimate[2]), bignumber)
  scale.nb.male <- ifelse(length(fit.negbin) > 0, as.numeric(fit.negbin$estimate[1]), bignumber) #(theta = p/(1-p))
  
  
  # Concurrency point prevalence 6 months before a survey, among men
  pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                     timepoint = datalist$itable$population.simtime[1] - 0.5)
  
  
  hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25),
                                               timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
  hiv.prev.25.34.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(25, 35),
                                                timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  hiv.prev.25.34.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(25, 35),
                                              timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
  hiv.prev.35.44.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(35, 45),
                                                timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  hiv.prev.35.44.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(35, 45),
                                              timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
  
  # 
  # prev.plot <- prevalence.plotter(datalist = datalist.agemix,
  #                                 agegroup = c(15, 50))
  
  
  
  growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                                      timewindow = c(0, datalist.agemix$itable$population.simtime[1]))
  
  
  
  
  
  ART.coverage.vector.creator <- function(datalist = datalist,
                                          agegroup = c(15, 50)){
    ART.cov.eval.timepoints <- seq(from = datalist$itable$t[2],
                                   to = datalist$itable$population.simtime[1])
    ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
    for (art.cov.index in 1:length(ART.cov.vector)){
      ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist,
                                                                   agegroup = agegroup,
                                                                   timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
        sum(ART.coverage.calculator(datalist = datalist,
                                    agegroup = agegroup,
                                    timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
    }
    return(ART.cov.vector)
  }
  
  # 
  # cov.vector <- ART.coverage.vector.creator(datalist = datalist.agemix,
  #                                           agegroup = c(15, 50))
  # plot(cov.vector)
  
  #######
  # And now we add HIV incidence at two selected time points
  #######
  # inc.times <- seq(from = 11, #datalist.agemix$itable$population.simtime[1] - 5,
  #                  by = 1, #5,
  #                  to = datalist.agemix$itable$population.simtime[1]) # 25, which is a 5-year forward projection
  # inc.vector <- rep(NA, length(inc.times))
  # for (inc.time in (1:length(inc.times))){
  #
  #   inc.tibble <- incidence.calculator(datalist = datalist.agemix,
  #                                                agegroup = c(0, 500), # Essentially incidence in the entire population
  #                                                timewindow = c(inc.times[inc.time] - 1, inc.times[inc.time]),
  #                                                only.active = "No")
  #   inc.vector[inc.time] <- inc.tibble$incidence[3]
  # }
  
  
  #      inc.time <- datalist.agemix$itable$population.simtime[1]
  #      inc.tibble <- incidence.calculator(datalist = datalist.agemix,
  #                                         agegroup = c(0, 500), # Essentially incidence in the entire population
  #                                         timewindow = c(inc.time - 1, inc.time),
  #                                         only.active = "No")
  #      inc.end <- inc.tibble$incidence[3]
  
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  # Resource required RSimpactHelp function in my branch
  # source("/user/data/gent/vsc400/vsc40070/phylo/transmNetworkBuilder.diff2.R")
  # source("/user/data/gent/vsc400/vsc40070/phylo/trans.network2tree.R")
  # source("/user/data/gent/vsc400/vsc40070/phylo/epi2tree2.R")
  
  # Resource required RSimpactHelp function in my branch
  source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff3.R")
  source("/home/david/RSimpactHelp/R/trans.network2tree.R")
  source("/home/david/RSimpactHelp/R/epi2tree2.R")
  
  
  simpact.trans.net <- transmNetworkBuilder.diff3(datalist = datalist.agemix, endpoint = 40)
  
  
  # save(simpact.trans.net, file = "simpact.trans.net.RData")
  
  # Read saved output data set
  # simpact.trans.net <- get(load("simpact.trans.net.RData"))
  
  # Check trees with negative branch lengths and at least with 3 leaves
  
  # smallest.branches <- rep(NA, times = length(simpact.trans.net))
  # for (list.element in 1:length(simpact.trans.net)){
  #   net.list <- simpact.trans.net[[list.element]]
  #   if(length(net.list$id) > 2){
  #     tree.tr <- epi2tree2(net.list)
  #     smallest.branch <- min(tree.tr$edge.length)
  #     smallest.branches[list.element] <- smallest.branch
  #   }
  # }
  # min(smallest.branches, na.rm = TRUE) #
  # ## seeds and transmission network sizes:
  # # 2>>10,  3>>212,  8>>20, 10>>10, 12>>3, 16>>13, 19>>1741
  # # which(smallest.branches!="NA")
  
  
  
  ###############################
  # Step 3: Sequence simulation #
  ###############################
  
  
  # Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp
  
  # Note: transmission network with less than 3 individuals will not be considered
  
  
  seed <- inputvector[1]
  
  trans.net <- simpact.trans.net # all transmission networks
  num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
  # constrained to rename IDs to -1, 0, 1, 2, ...
  num.i <- vector() # i_th seed in the list of seeds
  
  for(i in 1:length(trans.net)){
    
    tree.n <- trans.net[[i]] # transmission network for i^th seed
    
    if(nrow(as.data.frame(tree.n)) >= 3){
      
      # Construct transmission trees
      tree.i <- trans.network2tree(transnetwork = tree.n)
      num.trees <- c(num.trees,tree.n$id[1])
      num.i <- c(num.i,i)
      
      tree.j <- tree.i
      tree.j$tip.label <- paste(i,".", tree.j$tip.label, ".A", sep = "")
      
      # Save the transmission tree
      write.tree(tree.j, file = paste("tree.model1.seed",i,".nwk", sep = ""))
      
      # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = "")
      
      
      # Keep sampling dates, add "A" on ID in order to handle all ID's as characters
      id.samplingtime <- as.data.frame(cbind(paste(i,".",tree.n$id, ".A", sep = ""), tree.n$dtimes)) # IDs and their sampling times in the transmission network
      
      write.csv(id.samplingtime,file=paste("samplingtimes_seed_",i,".csv", sep = ""))
      
    }
  }
  
  IDs.transm <- num.i # vector of seeds with at least 2 transmissions
  
  
  
  ### Binding all transmission trees together ###
  ###############################################
  
  # Make a list of transmission trees
  
  trees <- list() # list of all transmission trees
  
  for(j in 1:length(IDs.transm)){
    p <- IDs.transm[j]
    tr <- read.tree(file = paste("tree.model1.seed",p,".nwk", sep = ""))
    trees[[j]] <- tr
  }
  
  class(trees)<-"multiPhylo"
  
  # print(trees,details=TRUE)
  
  
  
  
  # Function to bind all transmission trees in the list above
  
  bind.trees<-function(trees){
    if(length(trees)==2) return(bind.tree(trees[[1]],trees[[2]], where = "root", position = 0))
    else {
      trees<-c(bind.tree(trees[[1]],trees[[2]]),
               if(length(trees)>2) trees[3:length(trees)])
      trees<-bind.trees(trees) ## this is the recursive part
      return(trees)
    }
  }
  
  combined.tree<-bind.trees(trees) # This is a polytomy tree
  
  # is.binary.tree(combined.tree)
  
  # resolve polytomies of the combined transmission tree
  
  resolved.combined.tree <- multi2di(combined.tree) # resolve polytomies
  
  # is.binary.tree(resolved.combined.tree)
  
  
  ### Bind sampling dates ###
  ###########################
  
  for (i in IDs.transm){
    if(i==IDs.transm[1]){
      colitem <- read.csv(paste("samplingtimes_seed_",i,".csv", sep = ""))
    }
    else{
      
      ritem <- read.csv(paste("samplingtimes_seed_",i,".csv", sep = ""))
      colitem <- rbind(colitem, ritem)
    }
  }
  
  
  # Prepare to run the sequence simulation
  
  seq.rand <- 1 # first sequence
  
  n.tr <- 1 # number of transmission tree
  
  seed.id <- inputvector[1] # seed for simulation
  
  # # call the seed sequences - pool of viruses and rename the file
  file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("seed.seq.bis.nwk", sep = ""))
  # add the number of tree in the file and
  write(n.tr,file = paste("seed.seq.bis.nwk",sep = ""), append = TRUE)  # n.tr
  # the tree, to prepare the file to simulate the evolution of the virus across the tree
  
  write.tree(resolved.combined.tree,file = paste("seed.seq.bis.nwk", sep = ""), append = TRUE)
  
  file.rename(from = paste("seed.seq.bis.nwk", sep = ""), to = paste("seed.seq.bis.sim.nwk", sep = ""))
  
  system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230  -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< seed.seq.bis.sim.nwk -z",seed.id," > B.Epidemic_seed.seq.bis.sim.nwk.fasta",sep = ""))
  
  # a: shape parameter of Gamma > Gamma Rate Heterogeneity
  # g: category of Gamma > Discrete Gamma Rate Heterogeneity
  # r: rate matrix
  # s: scale which is the substitution rate of pol gene
  # z: seed
  
  
  ##################################################################################################
  # Step 4: Construct time stamped phylogenetic tree of the epidemic with simulated sequences data #
  ##################################################################################################
  
  
  # Function to tranform dates in named vector to be handled by treedater
  
  dates.Transform.NamedVector  <- function(dates=dates){
    
    dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # dates
    names(dates.val) <- as.character(dates$V1) # names are the names of the tips
    
    return(dates.val)
  }
  
  # 4.1. Construct phylogenetic tree
  
  # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree
  
  system(paste("./FastTree  -gtr -nt <", paste("B.Epidemic_seed.seq.bis.sim.nwk.fasta"), paste(">B.Epidemic_seed.seq.bis.sim.nwk.fasta.tree", sep = "")))
  
  
  samp.dates <- colitem # call the data table of dates
  
  tree.const <- read.tree(paste("B.Epidemic_seed.seq.bis.sim.nwk.fasta.tree", sep = ""))
  
  time.samp <- dates.Transform.NamedVector(dates=samp.dates) # name the dates
  
  
  ### Match dates and phylogenetic tree leaves ###
  ################################################
  
  tree.tips <- as.character(tree.const$tip.label)
  
  Ord.tree.dates <- vector()
  
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }
  
  # 4.2. Calibrate the phylogenetic tree
  
  # Use of library(treedater) to calibrate internal nodes
  
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 3012) # s is the length of sequence
  
  # save(dater.tree, file = paste("dated.tree.object.Rdata", sep = ""))
  
  # test.t.dated <- dater.tree
  
  # dater.tree <- get(load("dated.tree.object.Rdata"))
  
  N <- node.age(dater.tree)
  
  int.node.age <- N$Ti # internal nodes ages
  
  latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
  
  
  write.tree(dater.tree, file = "calibrated.tree.nwk") # for further manipulation in an easy format
  
  ### Lineage Through Time ###
  ############################
  
  # Estimating confidence intervals for rates and dates using a parametric bootstrap
  
  # pb <- parboot.treedater(dater.tree,
  #                         nreps = 100,
  #                         level = .95 ) # Lineage Through Time
  
  # plot.parboot.ltt( pb ) # This function computes the lineages through time given bootstrap replicate trees.
  # The pseudo-maximum likelihood estimate is plotted alongside CIs based on bootstrap trees.
  
  # Function to compute lineages through time with confidence intervals
  # 
  # plot.parboot.ltt.dat <- function (pbtd, t0 = NA, res = 100, ...)
  # {
  #   t1 <- max(pbtd$td$sts, na.rm = T)
  #   if (is.na(t0))
  #     t0 <- min(sapply(pbtd$trees, function(tr) tr$timeOf))
  #   times <- seq(t0, t1, l = res)
  #   ltt <- cbind(times = times, t(sapply(times, function(t) {
  #     c(pml = sum(pbtd$td$sts > t) - sum(pbtd$td$Ti > t), setNames(quantile(sapply(pbtd$trees,
  #                                                                                  function(tre) sum(tre$sts > t) - sum(tre$Ti > t)),
  #                                                                           probs = c(0.025, 0.5, 0.975)), c("lb", "median",
  #                                                                                                            "ub")))
  #   })))
  #   pl.df <- as.data.frame(ltt)
  #   return(pl.df)
  #   # p <- ggplot(pl.df) + geom_ribbon(aes(x = times, ymin = lb,
  #   #                                      ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
  #   # p <- p + geom_path(aes(x = times, y = pml))
  #   # (p <- p + ylab("Lineages through time") + xlab("Time"))
  # }
  # 
  # LTT <- plot.parboot.ltt.dat(pb)
  # 
  
  
  ###################################################
  # Step 5: Visualise and Use of simulation outputs #
  ###################################################
  
  
  
  ## Count transmission events per each transmission network ##
  #############################################################
  
  
  trans.net <- simpact.trans.net
  
  transmissions.vec.i <- list()
  time.calendar.vec.i <- list()
  
  for (i in 1:length(IDs.transm)){
    
    p <- IDs.transm[i]
    
    trans.net.i <- simpact.trans.net[[p]]
    
    Infec.time.i <- trans.net.i$InfecTime + 1987
    
    min.val = 1987
    max.val = round(max(trans.net.i$itimes)) + 1987
    
    step.int=1
    d <- (max.val-min.val)/step.int
    
    
    dat.f.trans.i <- as.data.frame(trans.net.i)
    dat.f.trans.i$itimes <- abs(dat.f.trans.i$itimes-40)+1977
    
    numb.tra <- vector()
    i.vec <- vector()
    for (j in 1:d){
      inf <- 1986+j
      sup <- 1987+j
      dat.f.trans.count <- dat.f.trans.i[which(dat.f.trans.i$itimes <= sup & dat.f.trans.i$itimes  >= inf),]
      numb.i <- nrow(dat.f.trans.count)
      numb.tra <- c(numb.tra, numb.i)
      i.vec <- c(i.vec, sup)
    }
    
    transmissions.vec.i[[i]] <- numb.tra
    time.calendar.vec.i[[i]] <- i.vec
    
  }
  
  ## Summing transmission events in all transmission networks ##
  ##############################################################
  
  for (i in 1:length(IDs.transm)){
    if(i==1){
      trans.sum <- transmissions.vec.i[[i]]
    }
    else{
      
      read.trans.sum <- transmissions.vec.i[[i]]
      trans.sum <- trans.sum + read.trans.sum
    }
  }
  
  
  
  ## Summing internal nodes in same age interval ##
  #################################################
  
  dt.node.age.dt <- int.node.age
  
  min.val = 1987
  max.val = 2017
  
  step.int=1
  d <- (max.val-min.val)/step.int
  
  int.node.vec <- vector()
  for (i in 1:d) {
    inf <- 1986+i
    sup <- 1987+i
    int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt >= inf]
    int.node.vec <- c(int.node.vec,length(int.node.age.i))
  }
  
  
  # Transmissison events and internal nodes within same time intervals
  
  numb.tra <- trans.sum
  int.node.vec <- int.node.vec
  
  
  ## Transmission network plotting: Union of all transmission networks ##
  #######################################################################
  
  network.list <- list()
  
  for (i in 1:length(IDs.transm)){
    
    p <- IDs.transm[i]
    
    trans.net.i <- as.data.frame(simpact.trans.net[[p]])
    
    trans.net.i <- trans.net.i[-1,]
    
    trans.net.i$id <- paste("A.",p,".",trans.net.i$id, sep = "")
    trans.net.i$parent <- paste("A.",p,".",trans.net.i$parent, sep = "")
    
    graph.build <- trans.net.i
    
    graph.build[,4] <- as.character(graph.build$parent) # donors
    graph.build[,3] <- as.character(graph.build$id) # recipients
    gag = as.matrix(graph.build)
    ga.graph = graph.edgelist(gag[,4:3])
    
    V(ga.graph)$color <- "red"
    
    # transNet.yrs.Old <- delete.vertices(ga.graph, "-1")
    
    network.list[[i]] <- ga.graph
    
  }
  
  
  for (i in 1:length(IDs.transm)){
    if(i==1){
      trans.network.union <- network.list[[i]]
    }
    else{
      
      read.trans.network.union <- network.list[[i]]
      trans.network.union <- trans.network.union + read.trans.network.union
    }
  }
  
  
  # 1. Transmission network from simpact                           # 1 #
  # plot.igraph(trans.network.union, edge.arrow.size=0.1, vertex.size=5,
  #             edge.color="black",
  #             asp = 1,
  #             xlim = c(-1, 2),
  #             ylim = c(-0.5,0.5),
  #             vertex.frame.color="black",
  #             vertex.label.color="black",
  #             #vertex.label = NULL,
  #             layout = layout_with_kk,
  #             edge.width = 1,
  #             vertex.label.cex=0.1,
  #             vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0
  #             #main = "True transmission network"
  # )
  
  # 2. Phylogenetic tree
  
  # plot(dater.tree, show.tip.label=FALSE,
  #      edge.width=1,
  #      edge.color="blue") # Try a few different settings!
  # axisPhylo(backward = FALSE)
  
  
  # 3. Transmission events and internal nodes
  # x <- i.vec
  # plot(x, int.node.vec, type="b", col="red", lwd=2,
  #      xlab = "Calendar time",
  #      ylab = "Count") # 1 > 1
  # lines(x, numb.tra, col='green3', type='b', lwd=2)
  # legend("topleft", legend = c("Internal nodes", "Transmission events"),
  #        col=c("red","green3"), pch=1)
  
  # 4. Plot of lineages through time with confidence intervals
  
  # plot.parboot.ltt( pb )
  
  # Or
  # 
  # pl.df <- LTT
  # 
  # p <- ggplot(pl.df) + geom_ribbon(aes(x = times, ymin = lb,
  #                                      ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
  # p <- p + geom_path(aes(x = times, y = pml))
  # (p <- p + ylab("Lineages through time") + xlab("Calendar time"))
  # 
  
  # # Object for plotting by ggplot
  # SimpactPaperPhyloExample <- list()
  # SimpactPaperPhyloExample$transm.network <- trans.network.union
  # SimpactPaperPhyloExample$dater.tree <- dater.tree
  # SimpactPaperPhyloExample$years.vec <- i.vec
  # SimpactPaperPhyloExample$int.node.vec <- int.node.vec
  # SimpactPaperPhyloExample$numb.trans.vec <- numb.tra
  # save(SimpactPaperPhyloExample, file = "SimpactPaperPhyloExample.RData")
  
  ########################################
  ### FEATURES FROM PHYLOGENETIC TREE ####
  ########################################
  
  # 1. Mean height of internal nodes & Maximum height
  
  # The height of a node is the number of edges on the longest path from the node to a leaf.
  # A leaf node will have a height of 0.
  
  # library(phytools)
  
  tree.cal <- read.tree("calibrated.tree.nwk") 
  
  H <- nodeHeights(dater.tree) # similar to node.depth.edgelength(tree) 
  
  
  # It's clear from a casual inspection of the matrix that each parent node height (in the right column) 
  # is represented twice and only twice. Thus, if we exclude the root node (zero height), 
  # we can just take the mean of H[,1].
  
  feature3 <- mean(sort(H[,1])[3:nrow(H)]) # Statistics
  
  # library(phyloTop)
  
  feature4 <- maxHeight(tree.cal, normalise = FALSE) # Statistics
  
  
  
  # Phylogenetic tree balace  features
  
  # library(phyloTop)
  
  feature5 <- colless.phylo(tree.cal, normalise = TRUE) # Statistics
  
  feature6 <- sackin.phylo(tree.cal, normalise = FALSE) # Statistics
  
  ###########################################################
  ######    TRANSMISSION NETWORK CHARACTERISTICS    #########
  ###########################################################
  
  
  ###########################
  # 1. Time trend incidence #
  ###########################
  
  incidence.df.15.24 <- as.data.frame(incidence.calculator(datalist = datalist.agemix,
                                                           agegroup = c(15, 25), timewindow = c(10, 40)))
  
  incidence.df.15.24.men <- incidence.df.15.24$incidence[1] # Statistics to save
  incidence.df.15.24.women <- incidence.df.15.24$incidence[2] # Statistics to save
  
  incidence.df.15.24.men.95.ll <- incidence.df.15.24$incidence.95.ll[1] # Statistics to save
  incidence.df.15.24.women.95.ll <- incidence.df.15.24$incidence.95.ll[2] # Statistics to save
  
  incidence.df.15.24.men.95.ul <- incidence.df.15.24$incidence.95.ul[1] # Statistics to save
  incidence.df.15.24.women.95.ul <- incidence.df.15.24$incidence.95.ul[2] # Statistics to save
  
  
  incidence.df.25.34 <- as.data.frame(incidence.calculator(datalist = datalist.agemix,
                                                           agegroup = c(25, 35), timewindow = c(10, 40)))
  
  incidence.df.25.34.men <- incidence.df.25.34$incidence[1] # Statistics to save
  incidence.df.25.34.women <- incidence.df.15.24$incidence[2] # Statistics to save
  
  incidence.df.25.34.men.95.ll <- incidence.df.25.34$incidence.95.ll[1] # Statistics to save
  incidence.df.25.34.women.95.ll <- incidence.df.25.34$incidence.95.ll[2] # Statistics to save
  
  incidence.df.25.34.men.95.ul <- incidence.df.25.34$incidence.95.ul[1] # Statistics to save
  incidence.df.25.34.women.95.ul <- incidence.df.25.34$incidence.95.ul[2] # Statistics to save
  
  
  incidence.df.35.44 <- as.data.frame(incidence.calculator(datalist = datalist.agemix,
                                                           agegroup = c(35, 45), timewindow = c(10, 40)))
  
  incidence.df.35.44.men <- incidence.df.35.44$incidence[1] # Statistics to save
  incidence.df.35.44.women <- incidence.df.35.44$incidence[2] # Statistics to save
  
  incidence.df.35.44.men.95.ll <- incidence.df.35.44$incidence.95.ll[1] # Statistics to save
  incidence.df.35.44.women.95.ll <- incidence.df.35.44$incidence.95.ll[2] # Statistics to save
  
  incidence.df.35.44.men.95.ul <- incidence.df.35.44$incidence.95.ul[1] # Statistics to save
  incidence.df.35.44.women.95.ul <- incidence.df.35.44$incidence.95.ul[2] # Statistics to save
  
  
  ##########################################
  # 2. Age mixing in transmission networks #
  ##########################################
  
  # Source transmission networks builder function
  
  # source("/home/david/Dropbox/Niyukuri/Stress_Testing_Cluster_Code/transmNetworkBuilder.diff2.R")
  # 
  
  # Function of age mixing data: agemixing.trans.df(datalist = datalist, 
  #                                                 trans.network = trans.network)
  
  # require library(data.table)
  
  #' return a data table of all transmissions 
  #' ID1 : ID of men  
  #' ID2 : ID of women    
  #' TOBID1 : Time of Birth of men     
  #' TOBID2 : Time of Birth of women   
  #' AgeID1 : Age of men at transmission time 
  #' AgeID2 : Age of women at transmission time
  #' AgeGap : Age gap between man & woman
  #' infecttime : Infection time
  #' samptime : Removal time
  
  # Function to compute the data table of age mixing in transmission
  ####################################################################
  
  # Due to much computation time of transmNetworkBuilder.diff2() we put it external to age-mixing function
  
  
  # 
  # trans.network <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = datalist$itable$population.simtime[1])
  # 
  
  trans.network <- simpact.trans.net
  
  
  # Define function for age-mixing in transmission networks
  
  agemixing.trans.df <- function(datalist = datalist, 
                                 trans.network = trans.network){
    
    pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
    
    pers.infec <- pers.infec.raw[which(pers.infec.raw$InfectTime <= datalist$itable$population.simtime[1]),]
    
    # person table of infected individuals by seed event
    pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)
    
    # id of people who got infection by seed event: seeds.id
    seeds.id <- pers.table.seed$ID # do
    
    infectionTable <- vector("list", length(seeds.id))
    
    for (i in 1: length(seeds.id)) {
      
      trans.network.i <- as.data.frame(trans.network[[i]])
      
      if(nrow(trans.network.i) >= 3){
        
        trans.network.i <- trans.network.i[-1,]
        
        rrtable <- as.data.frame(cbind(trans.network.i$RecId, trans.network.i$DonId, 
                                       trans.network.i$GenderRec, trans.network.i$GenderDon,
                                       trans.network.i$TOBRec, trans.network.i$TOBDon,
                                       trans.network.i$InfecTime, trans.network.i$SampTime))
        
        names(rrtable) <- c("RecId", "DonId", "GenderRec",
                            "GenderDon", "TOBRec", "TOBDon", "InfecTime", "SampTime")
        
        rrtable.men <- subset(rrtable, rrtable$GenderDon=="0")
        rrtable.women <- subset(rrtable, rrtable$GenderDon=="1")
        
        ids.men <- rrtable.men$DonId
        ids.men.part.w <- rrtable.men$RecId
        age.gap.ID1 <- abs(rrtable.men$TOBDon) - abs(rrtable.men$TOBRec) # men are donors
        tob.men.ID1 <- rrtable.men$TOBDon
        tob.women.ID1 <- rrtable.men$TOBRec
        
        age.men.ID1 <- abs(rrtable.men$TOBDon) + rrtable.men$InfecTime
        age.women.ID1 <- abs(rrtable.men$TOBRec) + rrtable.men$InfecTime
        
        infectime.m <- rrtable.men$InfecTime
        samptime.m <- rrtable.men$SampTime
        
        ID1.m <- ids.men
        ID2.m <- ids.men.part.w
        age.gap.m <- age.gap.ID1
        infectime.m <- infectime.m 
        
        infectable.m <- cbind(ID1.m, ID2.m, tob.men.ID1, tob.women.ID1, age.men.ID1, age.women.ID1, age.gap.m, infectime.m, samptime.m)
        
        ids.women <- rrtable.women$DonId
        ids.women.part.m <- rrtable.women$RecId
        age.gap.ID2 <- abs(rrtable.women$TOBRec) - abs(rrtable.women$TOBDon) # men are receptors
        tob.men.ID2 <- rrtable.women$TOBRec
        tob.women.ID2 <- rrtable.women$TOBDon
        
        age.men.ID2 <- abs(rrtable.women$TOBDon) + rrtable.women$InfecTime
        age.women.ID2 <- abs(rrtable.women$TOBRec) + rrtable.women$InfecTime
        
        infectime.w <- rrtable.women$InfecTime
        samptime.w <- rrtable.women$SampTime
        
        ID1.w <- ids.women.part.m
        ID2.w <- ids.women
        age.gap.w <- age.gap.ID2
        infectime.w <- infectime.w
        
        infectable.w <- cbind(ID1.w, ID2.w, tob.men.ID2, tob.women.ID2, age.men.ID2, age.women.ID2, age.gap.w, infectime.w, samptime.w)
        
        infectable.i <- as.data.frame(rbind(infectable.m, infectable.w))
        
        names(infectable.i) <- c("ID1", "ID2", "TOBID1", "TOBID2", "AgeID1", "AgeID2", "AgeGap", "infecttime", "samptime")
        infectionTable[[i]] <- infectable.i
      }
      
      
    }
    
    
    infecttable <- rbindlist(infectionTable) 
    
    return(infecttable)
    
  }
  
  
  # Call the function of age-mixing
  
  agemix.df <- agemixing.trans.df(datalist = datalist.agemix, 
                                  trans.network = trans.network)
  
  
  # Fit mixed effects model to age mixing data
  
  # require  library(lme4)
  
  fit.agemix.trans <- function(datatable = agemix.df){
    
    datatable <- agemix.df
    
    agemix.inter <- lmer(AgeID2 ~ AgeID1 + (1|ID1), data = datatable) # choice of age_woman by a man
    
    return(agemix.inter)
    
  }
  
  
  agemix.fit <- fit.agemix.trans(datatable = agemix.df)
  
  
  
  # Visualisation of age mixing in transmission
  
  coef.inter <- fixef(agemix.fit)
  
  age.mix.intercept <- coef.inter[1] # Statistics to save
  
  age.mix.slope <- coef.inter[2] # Statistics to save
  
  # 
  # x=agemix.df$AgeID1
  # y=agemix.df$AgeID2
  # plot(x, y, lwd=1, col = "blue",
  #      xlab = "Man Age",
  #      ylab = "Woman Age")
  # abline(coef = c(coef.inter[1], coef.inter[2]))
  
  
  ###########################
  # 3. Onward transmissions # considering people who already died > they had full time to be able to transmit the infection
  ###########################
  
  # Function to return a list of how many onward infections produce by donors (after cquiring the infection)
  
  # requires transNetBuilder3() function
  
  onwardtransmissions.dat <- function(datalist = datalist, 
                                      trans.network = trans.network){
    
    
    pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
    
    pers.infec.raw.died <- pers.infec.raw[pers.infec.raw$TOD != "Inf", ] # consider only these who had full time to transmit before they die
    
    pers.infec <- pers.infec.raw.died[which(pers.infec.raw.died$InfectTime <= datalist$itable$population.simtime[1]),]
    
    # person table of infected individuals by seed event
    pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)
    
    # id of people who got infection by seed event: seeds.id
    seeds.id <- pers.table.seed$ID # do
    
    
    # Onward transmissions in each transmission network
    
    onwardtransm <- vector("list", length(seeds.id))
    
    for (j in 1: length(seeds.id)) {
      
      trans.network.j <- as.data.frame(trans.network[[j]])
      
      trans.network.j <- trans.network.j[-1,] # remove the universal infector
      
      if(nrow(trans.network.j) > 1){ # consider transmission networks with at least one onward transmission
        
        d.j <- table(trans.network.j$DonId) # in the transmission table, the number of times DonId appears is the number of Onward transmissions after acuiring the infection 
        num.j <- as.data.frame(as.numeric(d.j))
        names(num.j) <- c("TransCount")
        onwardtransm[[j]] <- num.j
        
      }
      
    }
    
    onwardtransmissions <- rbindlist(onwardtransm) 
    
    count.dat <- onwardtransmissions$TransCount
    
    return(count.dat) # count.dat = all infections - seeds which didn;t produce at least one transmission
    
  }
  
  
  transm.count <- onwardtransmissions.dat(datalist = datalist.agemix, 
                                          trans.network = trans.network)
  
  
  mean.transm.count <- mean(transm.count) # Statistics to save
  sd.transm.count <- sd(transm.count) # Statistics to save
  
  # Advise for which distribution: power-law???
  
  # Visualisation of onward distribution
  
  # hist(transm.count)
  
  # Power alw fitting of onward transmissions
  
  fit.pow <- power.law.fit(transm.count)
  
  alpha.stat <- fit.pow$alpha # Statistics to save
  min.stat <- fit.pow$xmin # Statistics to save
  logLik.stat <- fit.pow$logLik # Statistics to save
  ks.test.stat <- fit.pow$KS.stat # Statistics to save
  ks.p.stat <- fit.pow$KS.p # Statistics to save
  
  
  
  ## Relationship rate calculator function
  
  relationship.rate.calculator <- function(datalist = datalist,
                                           timewindow = c(20, 40), int = FALSE, by=1){
    
    if(int==FALSE){
      Rels.table <- datalist$rtable
      
      Rels.table.window <- Rels.table %>%
        subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])
      
      numb.rels <- nrow(Rels.table.window)
      
      numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women
      
      numb.rels.women <- length(unique(Rels.table.window$ID2))
      
      rels.rate <- (numb.rels)/ diff(timewindow)
      
      return(rels.rate)
    }
    
    if(int==TRUE){
      
      Rels.table <- datalist$rtable
      
      upper.limit <- ceiling(diff(timewindow)/by)
      
      interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)
      
      rels.int <- vector()
      rels.rate.int <- vector()
      
      
      for(i in 0:(upper.limit-2)){
        
        timewindow.int <- c(interv.time[1+i], interv.time[2+i])
        
        rels.numb <- nrow(Rels.table %>%
                            subset(FormTime <=timewindow.int[2] & FormTime >= timewindow.int[1]))
        
        
        rels.int <- c(rels.int,  rels.numb)
        
        
        rels.rate.int <- c(rels.rate.int, (rels.numb / diff(timewindow.int)))
        
      }
      
      return(rels.rate.int)
    }
    
    
  }
  
  
  # call the relationship rate function 
  rels.rate <- relationship.rate.calculator(datalist = datalist.agemix,
                                            timewindow = c(datalist$itable$hivseed.time[[1]], 
                                                           datalist$itable$population.simtime[1]), int = FALSE, by=1)
  
  
  ## Transmission rate function
  
  transmission.rate.calculator <- function(datalist = datalist,
                                           timewindow = c(20, 40), int = FALSE, by=1){
    
    if(int==FALSE){
      
      Infec.pop.table <- datalist$ptable[InfectType==1]
      
      numb.infec.pop <- nrow(Infec.pop.table %>%
                               subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1]))
      
      
      transm.rate <- numb.infec.pop / diff(timewindow)
      
      return(transm.rate)
    }
    
    if(int==TRUE){
      
      Infec.pop.table <- datalist$ptable[InfectType==1]
      
      upper.limit <- ceiling(diff(timewindow)/by)
      
      interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)
      
      infec.pop.int <- vector()
      trans.rate.int <- vector()
      
      
      for(i in 0:(upper.limit-2)){
        
        timewindow.int <- c(interv.time[1+i], interv.time[2+i])
        
        infec.pop.numb <- nrow(Infec.pop.table %>%
                                 subset(InfectTime <=timewindow.int[2] & InfectTime >= timewindow.int[1]))
        
        
        infec.pop.int <- c(infec.pop.int,  infec.pop.numb)
        
        
        trans.rate.int <- c(trans.rate.int, (infec.pop.numb / diff(timewindow.int)))
        
      }
      
      return(trans.rate.int)
    }
    
  }
  
  
  # call the transmission rate function 
  trans.rate <- transmission.rate.calculator(datalist = datalist.agemix,
                                             timewindow = c(datalist$itable$hivseed.time[[1]], 
                                                            datalist$itable$population.simtime[1]), int = FALSE, by=1)
  
  
  
  #### Combined all together  in a vector ###
  
  summary.df <- c(
    
    pp.cp.6months.male, # <- concurr.pointprev.calculator(datalist = datalist.agemix,
    #                              timepoint = datalist$itable$population.simtime[1] - 0.5)
    hiv.prev.lt25.women, # <- prevalence.calculator(datalist = datalist.agemix,
    #                       agegroup = c(15, 25),
    #                      timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
    hiv.prev.lt25.men, # <- prevalence.calculator(datalist = datalist.agemix,
    #                       agegroup = c(15, 25),
    #                      timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
    hiv.prev.25.34.women, # <- prevalence.calculator(datalist = datalist.agemix,
    #                       agegroup = c(25, 35),
    #                      timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
    hiv.prev.25.34.men, # <- prevalence.calculator(datalist = datalist.agemix,
    #                       agegroup = c(25, 35),
    #                      timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
    hiv.prev.35.44.women, # <- prevalence.calculator(datalist = datalist.agemix,
    #                       agegroup = c(35, 45),
    #                      timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
    hiv.prev.35.44.men, # <- prevalence.calculator(datalist = datalist.agemix,
    #                       agegroup = c(35, 45),
    #                      timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
    growthrate, # <- pop.growth.calculator(datalist = datalist.agemix,
    #                       timewindow = c(0, datalist.agemix$itable$population.simtime[1]))
    
    feature3, # <- mean(sort(H[,1])[3:nrow(H)]) # Statistics
    feature4, # <- maxHeight(tree.cal, normalise = FALSE) # Statistics
    feature5, # <- colless.phylo(tree.cal, normalise = TRUE) # Statistics
    feature6, # <- sackin.phylo(tree.cal, normalise = FALSE) # Statistics
    
    incidence.df.15.24.men, # <- incidence.df.15.24$incidence[1] # Statistics to save
    incidence.df.15.24.women, # <- incidence.df.15.24$incidence[2] # Statistics to save
    # incidence.df.15.24.men.95.ll <- incidence.df.15.24$incidence.95.ll[1] # Statistics to save
    # incidence.df.15.24.women.95.ll <- incidence.df.15.24$incidence.95.ll[2] # Statistics to save
    # incidence.df.15.24.men.95.ul <- incidence.df.15.24$incidence.95.ul[1] # Statistics to save
    # incidence.df.15.24.women.95.ul <- incidence.df.15.24$incidence.95.ul[2] # Statistics to save
    incidence.df.25.34.men, # <- incidence.df.25.34$incidence[1] # Statistics to save
    incidence.df.25.34.women, # <- incidence.df.15.24$incidence[2] # Statistics to save
    # incidence.df.25.34.men.95.ll <- incidence.df.25.34$incidence.95.ll[1] # Statistics to save
    # incidence.df.25.34.women.95.ll <- incidence.df.25.34$incidence.95.ll[2] # Statistics to save
    # incidence.df.25.34.men.95.ul <- incidence.df.25.34$incidence.95.ul[1] # Statistics to save
    # incidence.df.25.34.women.95.ul <- incidence.df.25.34$incidence.95.ul[2] # Statistics to save
    incidence.df.35.44.men, # <- incidence.df.35.44$incidence[1] # Statistics to save
    incidence.df.35.44.women, # <- incidence.df.35.44$incidence[2] # Statistics to save
    # incidence.df.35.44.men.95.ll <- incidence.df.35.44$incidence.95.ll[1] # Statistics to save
    # incidence.df.35.44.women.95.ll <- incidence.df.35.44$incidence.95.ll[2] # Statistics to save
    # incidence.df.35.44.men.95.ul <- incidence.df.35.44$incidence.95.ul[1] # Statistics to save
    # incidence.df.35.44.women.95.ul <- incidence.df.35.44$incidence.95.ul[2] # Statistics to save
    # 
    
    
    age.mix.intercept[[1]], # <- coef.inter[1] # Statistics to save
    age.mix.slope[[1]], # <- coef.inter[2] # Statistics to save
    
    mean.transm.count, # <- mean(transm.count) # Statistics to save
    sd.transm.count, # <- sd(transm.count) # Statistics to save
    
    alpha.stat, # <- fit.pow$alpha # Statistics to save
    min.stat, # <- fit.pow$xmin # Statistics to save
    logLik.stat, # <- fit.pow$logLik # Statistics to save
    ks.test.stat, # <- fit.pow$KS.stat # Statistics to save
    ks.p.stat, # <- fit.pow$KS.p # Statistics to save
    
    rels.rate,
    trans.rate
    
  )
  
  
  
  
  summary.features.df <- summary.df
  
  
  
  # Name the columns
  
  features.names <- c("concur.prev.6months.male", 
                      
                      "hiv.prev.lt25.women",
                      "hiv.prev.lt25.men",
                      "hiv.prev.25.34.women",
                      "hiv.prev.25.34.men",
                      "hiv.prev.35.44.women",
                      "hiv.prev.35.44.men", 
                      
                      "Pop.growthrate", 
                      
                      "meanHeightTree", 
                      "maxHeightTree", 
                      "CollesIndex", 
                      "SackinIndex", 
                      
                      "incidence.df.15.24.men", 
                      "incidence.df.15.24.women", 
                      "incidence.df.25.34.men", 
                      "incidence.df.25.34.women", 
                      "incidence.df.35.44.men",
                      "incidence.df.35.44.women", 
                      
                      "trans.age.mix.intercept",
                      "trans.age.mix.slope",
                      
                      "mean.transm.count", 
                      "sd.transm.count",
                      
                      "alpha.stat", 
                      "min.stat", 
                      "logLik.stat",
                      "ks.test.stat", 
                      "ks.p.stat",
                      "relationship.rate",
                      "transmission.rate")
  
  names(summary.features.df) <- c(features.names)
  
  outputvector <- summary.features.df
  
  return(outputvector)
  
}

