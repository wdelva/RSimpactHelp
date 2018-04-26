
#' Wrapper function for simpact and phylodynamics
#'
#' This function simulates an HIV epidemic and molecular evolution of HIV across
#' transmission network
#'
#' @param inputvector Vector of model input parameter values
#' @return A vector of model features
#' @import RSimpactCyan
#' @importFrom phytools nodeHeights
#' @importFrom picante node.age
#' @export
#'
wrapper.phylo.simpact.parallel <- function(inputvector){


  # Working directory where we find seq-gen, FastTree, and HIV seed file hiv.seq.C.pol.j.fasta

  # work.dir <- "/home/david/Desktop/TEST_20_4_2018" # on laptop

  work.dir <- "/user/data/gent/vsc400/vsc40070/phylo" # on cluster



  # If we do not include these packages we still get errors

  # pacman::p_load(RSimpactCyan, RSimpactHelper, devtools, Rcpp, ape, expoTree,
  #                data.table, tidyr, phylosim, readr, dplyr, adephylo, geiger,
  #                picante, igraph, ggplot2, magrittr, lmtest, gsubfn, utils, pcaPP,
  #                phyloTop, phytools, lme4, data.table, treedater, phangorn, nlme,
  #                fitdistrplus, apTreeshape)

  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  library(tidyr)
  library(readr)
  library(dplyr)
  library(adephylo)
  library(geiger)
  library(picante)
  library(igraph)
  library(magrittr)
  library(lmtest)
  library(gsubfn)
  library(utils)
  library(pcaPP)
  library(phyloTop)
  library(phytools)
  library(lme4)
  library(treedater)
  library(nlme)
  library(fitdistrplus)
  library(apTreeshape)

  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################


  age.distr <- agedistr.creator(shape = 5, scale = 65)

  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                   # population.msm = "no",
                                   population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
                                   population.nummen = 100, #3000, #600, # 3800, #2500,
                                   population.numwomen = 100, # 3000, #600, #4200, #2500,
                                   hivseed.time = 10, # 20,
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

  cfg.list["hivtransmission.param.f1"] <- log(inputvector[2])
  cfg.list["hivtransmission.param.f2"] <- log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[3]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[3]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]

  #cfg <- cfg.list

  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  # cfg["monitoring.fraction.log_viralload"] <- 0.3
  cfg.list["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

  seedid <- inputvector[1]
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
  cfg.list["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10] ######### -0.5
  cfg.list["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10] ######### -0.5
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[11]

  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
  cfg.list["conception.alpha_base"] <- inputvector[14] #is conception.alpha.base (higher up)
  cfg.list["dissolution.alpha_0"] <- inputvector[15]
  cfg.list["dissolution.alpha_4"] <- inputvector[16]

  #creating subfolder with unique name for each simulation
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
                             paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)

    return(sub.dir.sim.id)
  }

  sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))

  # # # # Run Simpact

  seedid <- inputvector[1]

  results <- simpact.run(configParams = cfg.list,
                         destDir = sub.dir.rename,
                         agedist = age.distr,
                         seed = seedid,
                         intervention = intervention)

  datalist <- readthedata(results)



  growthrate <- pop.growth.calculator(datalist = datalist,
                                      timewindow = c(0, datalist$itable$population.simtime[1]))



  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################



  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)




  ###############################
  # Step 3: Sequence simulation #
  ###############################


  trans.net <- simpact.trans.net # all transmission networks


  dirseqgen <- work.dir

  seeds.num <- inputvector[1]

  # Sequence simulation is done for at least a transmission network with 6 individuals
  # This means that limitTransmEvents equal at least 7

  sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                 sub.dir.rename = sub.dir.rename,
                                 seq.gen.tool = "seq-gen",
                                 datalist = datalist,
                                 seeds.num = seeds.num,
                                 endpoint = 40,
                                 limitTransmEvents = 7, # no less than 7
                                 hiv.seq.file = "hiv.seq.C.pol.j.fasta") # hiv.seq.file lodged in work.dir



  #####################################################
  # Step 4: Construct time stamped phylogenetic trees #
  #####################################################


  dirfasttree <- work.dir


  if(file.exists(paste0(sub.dir.rename, "/C.Epidemic_seed.seq.bis.sim.nwk.fasta"))==TRUE){

    # check if the run has simulated sequence,
    # in other words: e have a transmission network with at least 6 individuals


    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = sub.dir.rename,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40)


    ape::write.tree(tree.calib, file = paste0(sub.dir.rename,"/calibrated.tree.nwk"))


    N <- picante::node.age(tree.calib)

    int.node.age <- N$Ti # internal nodes ages

    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


    ########################################
    ### FEATURES FROM PHYLOGENETIC TREE ####
    ########################################

    # 1. Mean height of internal nodes & Maximum height

    # The height of a node is the number of edges on the longest path from the node to a leaf.
    # A leaf node will have a height of 0.

    # library(phytools)

    tree.cal <- ape::read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))


    H <- phytools::nodeHeights(tree.cal) # similar to node.depth.edgelength(tree)


    # It's clear from a casual inspection of the matrix that each parent node height (in the right column)
    # is represented twice and only twice. Thus, if we exclude the root node (zero height),
    # we can just take the mean of H[,1].

    feature3 <- mean(sort(H[,1])[3:nrow(H)]) # Statistics

    # # library(phyloTop)
    #
    # feature4 <- maxHeight(tree.cal, normalise = FALSE) # Statistics
    #
    #
    #
    # # Phylogenetic tree balace  features
    #
    # # library(phyloTop)
    #
    # feature5 <- colless.phylo(tree.cal, normalise = TRUE) # Statistics
    #
    # feature6 <- sackin.phylo(tree.cal, normalise = FALSE) # Statistics

    summary.df <- c(growthrate, feature3)

    features.names <- c("growth.rate", "mean.height.tree")

    names(summary.df) <- features.names

    return(summary.df)


  }else{

    features.names <- c("growth.rate", "mean.height.tree")

    summary.df <- rep(NA,length(features.names))

    names(summary.df) <- features.names

    return(summary.df)

  }


  #  unlink(paste0(sub.dir.rename, "/"), recursive = TRUE)

}

