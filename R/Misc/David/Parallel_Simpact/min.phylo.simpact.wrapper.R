
phylo.simpact.wrapper <- function(inputvector){

  pacman::p_load(RSimpactCyan, RSimpactHelper, devtools, Rcpp, ape, expoTree,
                 data.table, tidyr, phylosim, readr, dplyr, adephylo, geiger,
                 picante, igraph, ggplot2, magrittr, lmtest, gsubfn, utils, pcaPP,
                 phyloTop, phytools, lme4, data.table, treedater, phangorn, nlme,
                 fitdistrplus, apTreeshape)

  # ## For sub-optimal sequence coverage

  # Kick-start the missing sequence data simulation study:

  # Let's start with verifying that we have parameter combinations that produce "useful" output.
  # i.e. output that is more or less in line with the demography, behaviour and epidemic that we want to consider.


  age.distr <- agedistr.creator(shape = 5, scale = 65)

  #source the input files var name <- cfg.list <- cfg
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

  sub.dir.rename <- paste0("temp/",generate.filename(10))

  # # # # Run Simpact
  results <- simpact.run(configParams = cfg.list,
                         destDir = sub.dir.rename,
                         agedist = age.distr,
                         seed = seedid,
                         intervention = intervention)

  datalist <- readthedata(results)

  # Let's examine the output:
  #source the summary functions
  source("R/transmNetworkBuilder.diff3.R")
  source("R/Misc/David/Parallel_Simpact/00-Functions.R")

  growthrate <- pop.growth.calculator(datalist = datalist,
                                      timewindow = c(0, datalist$itable$population.simtime[1]))


  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40) # transmNetworkBuilder.diff3

  trans.net <- simpact.trans.net # all transmission networks

  limitTransmEvents <- 7

  # Vector of number of individuals in existing transmission networks
  # with at least limitTransmEvents-1 individuals

  TransmEventsCountVector <- vector()

  for(i in 1:length(trans.net)){
    trans.net.i.check <- as.data.frame(trans.net[[i]])
    if(nrow(trans.net.i.check)>=limitTransmEvents){
      TransmEventsCountVector <- c(TransmEventsCountVector, nrow(trans.net.i.check))
    }
  }


  ##### ######################
  #### Only One transmission network with at least limitTransmEvents each of them
  ############################
  if(length(TransmEventsCountVector)==1){



    num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
    # constrained to rename IDs to -1, 0, 1, 2, ...
    num.i <- vector() # i_th seed in the list of seeds

    for(i in 1:length(trans.net)){

      tree.n <- trans.net[[i]] # transmission network for i^th seed

      if(nrow(as.data.frame(tree.n)) >= limitTransmEvents){

        # Construct transmission trees
        tree.i <- trans.network2tree(transnetwork = tree.n)
        num.trees <- c(num.trees,tree.n$id[1])
        num.i <- c(num.i,i)

        tree.j <- tree.i
        tree.j$tip.label <- paste(i,".", tree.j$tip.label, ".A", sep = "")

        # Save the transmission tree
        write.tree(tree.j, file = paste0(sub.dir.rename,"/tree.model1.seed_",i,".nwk"))

        # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = "")


        # Keep sampling dates, add "A" on ID in order to handle all ID's as characters
        id.samplingtime <- as.data.frame(cbind(paste0(i,".",tree.n$id, ".A"), tree.n$dtimes)) # IDs and their sampling times in the transmission network

        write.csv(id.samplingtime,file=paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))

      }
    }

    IDs.transm <- num.i # vector of seeds with at least 2 transmissions



    ### Binding all transmission trees together ### DESACTIVATE THIS >>>
    ############################################### it's applicale whenw e have at least 2 transmission networks


    # Make a list of transmission trees
#
#     trees <- list() # list of all transmission trees
#
#     for(j in 1:length(IDs.transm)){
#       p <- IDs.transm[j]
#       tr <- read.tree(file = paste0(sub.dir.rename,"/tree.model1.seed_",p,".nwk"))
#       trees[[j]] <- tr
#     }
#
#     class(trees)<-"multiPhylo"
#
#     # print(trees,details=TRUE)
#
#     # Function to bind all transmission trees in the list above
#
#     bind.trees<-function(trees){
#       if(length(trees)==2) return(bind.tree(trees[[1]],trees[[2]], where = "root", position = 0))
#       else {
#         trees<-c(bind.tree(trees[[1]],trees[[2]]),
#                  if(length(trees)>2) trees[3:length(trees)])
#         trees<-bind.trees(trees) ## this is the recursive part
#         return(trees)
#       }
#     }
#
#     combined.tree<-bind.trees(trees) # This is a polytomy tree

    # is.binary.tree(combined.tree)

    # resolve polytomies of the combined transmission tree

    # resolved.combined.tree <- multi2di(combined.tree) # resolve polytomies

    p <- IDs.transm[1]

    resolved.combined.tree <- read.tree(file = paste0(sub.dir.rename,"/tree.model1.seed_",p,".nwk"))

    # is.binary.tree(resolved.combined.tree)


    ### Bind sampling dates ###  DESACTIVATE THIS >>> it's applicale whenw e have at least 2 transmission networks
    ###########################

    # for (i in IDs.transm){
    #   if(i==IDs.transm[1]){
    #     colitem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))
    #   }
    #   else{
    #
    #     ritem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))
    #     colitem <- rbind(colitem, ritem)
    #   }
    # }

    colitem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",p,".csv"))


    # Prepare to run the sequence simulation

    seq.rand <- 1 # first sequence

    n.tr <- 1 # number of transmission tree

    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste0("R/Misc/David/Parallel_Simpact/hiv.seq.C.pol.j.fasta"),paste0(sub.dir.rename,"/seed.seq.bis.nwk"))
    # add the number of tree in the file and
    write(n.tr,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree

    write.tree(resolved.combined.tree,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)

    file.rename(from = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), to = paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk"))

    out.seq.gen.file <- paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta")
    in.seq.gen.file <- paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk")

    system(paste0("R/Misc/David/Parallel_Simpact/seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230  -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< ",in.seq.gen.file," -z",seedid," > ", out.seq.gen.file))

    # Function to tranform dates in named vector to be handled by treedater

    dates.Transform.NamedVector  <- function(dates=dates){

      dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # dates
      names(dates.val) <- as.character(dates$V1) # names are the names of the tips

      return(dates.val)
    }

    # 4.1. Construct phylogenetic tree

    # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree


    out.fast.tree.file <- paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta.tree")

    system(paste("R/Misc/David/Parallel_Simpact/FastTree -gtr -nt < ", out.seq.gen.file, paste0("> ", out.fast.tree.file)))


    samp.dates <- colitem # call the data table of dates

    tree.const <- read.tree(out.fast.tree.file)

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


    write.tree(dater.tree, file = paste0(sub.dir.rename,"/calibrated.tree.nwk")) # for further manipulation in an easy format


    ## Count transmission events per each transmission network ##
    #############################################################


    trans.net <- simpact.trans.net

    transmissions.vec.i <- list()
    time.calendar.vec.i <- list()

    for (i in 1:length(IDs.transm)){

      p <- IDs.transm[i]

      trans.net.i <- simpact.trans.net[[p]]

      Infec.time.i <- trans.net.i$InfecTime + 1977

      min.val = 1987
      max.val = datalist$itable$population.simtime[1] + 1977 # round(max(trans.net.i$InfecTime)) + 1977

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
        if(is.null(numb.i)==TRUE){
          numb.tra <- c(numb.tra, 0)
        }else{
          numb.tra <- c(numb.tra, numb.i)
        }
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
    max.val = datalist$itable$population.simtime[1] + 1977 # round(max(trans.net.i$InfecTime)) + 1977

    step.int=1

    d <- (max.val-min.val)/step.int

    int.node.vec <- vector()
    for (i in 1:d) {
      inf <- 1986+i
      sup <- 1987+i
      int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt >= inf]
      if(is.null(int.node.age.i)==TRUE){
        int.node.vec <- c(int.node.vec,0)
      }else{
        int.node.vec <- c(int.node.vec,length(int.node.age.i))
      }

    }


    # Transmissison events and internal nodes within same time intervals

    numb.tra <- trans.sum
    int.node.vec <- int.node.vec


    # x <- i.vec
    # plot(x, int.node.vec, type="b", col="red", lwd=2,
    #      xlab = "Calendar time",
    #      ylab = "Count") # 1 > 1
    # lines(x, numb.tra, col='green3', type='b', lwd=2)
    # legend("topleft", legend = c("Internal nodes", "Transmission events"),
    #        col=c("red","green3"), pch=1)

    ## Transmission network plotting: Union of all transmission networks ##
    #######################################################################

    network.list <- list()

    for (i in 1:length(IDs.transm)){

      p <- IDs.transm[i]

      trans.net.i <- as.data.frame(simpact.trans.net[[p]])

      trans.net.i <- trans.net.i[-1,]

      trans.net.i$id <- paste("C.",p,".",trans.net.i$id, sep = "")
      trans.net.i$parent <- paste("C.",p,".",trans.net.i$parent, sep = "")

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


    # library(phytools)

    tree.cal <- read.tree(sub.dir.rename,"/calibrated.tree.nwk")

    H <- nodeHeights(dater.tree) # similar to node.depth.edgelength(tree)


    # It's clear from a casual inspection of the matrix that each parent node height (in the right column)
    # is represented twice and only twice. Thus, if we exclude the root node (zero height),
    # we can just take the mean of H[,1].

    feature3 <- mean(sort(H[,1])[3:nrow(H)]) # Statistics


    summary.df <- c(growthrate,
                    feature3) # <- mean(sort(H[,1])[3:nrow(H)]) # Statistics  )


    # Name the columns

    features.names <- c("Pop.growthrate","meanHeightTree")

    names(summary.df) <- c(features.names)

    return(summary.df)

  }



  ##### ######################
  #### Two or more transmissions networks with at least limitTransmEvents each of them
  ############################

  if(length(TransmEventsCountVector)>=2){ # at least two transmissions networks ith at least 6
                                          # transmissions each


    num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
    # constrained to rename IDs to -1, 0, 1, 2, ...
    num.i <- vector() # i_th seed in the list of seeds

    for(i in 1:length(trans.net)){

      tree.n <- trans.net[[i]] # transmission network for i^th seed

      if(nrow(as.data.frame(tree.n)) >= limitTransmEvents){

        # Construct transmission trees
        tree.i <- trans.network2tree(transnetwork = tree.n)
        num.trees <- c(num.trees,tree.n$id[1])
        num.i <- c(num.i,i)

        tree.j <- tree.i
        tree.j$tip.label <- paste(i,".", tree.j$tip.label, ".A", sep = "")

        # Save the transmission tree
        write.tree(tree.j, file = paste0(sub.dir.rename,"/tree.model1.seed_",i,".nwk"))

        # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = "")


        # Keep sampling dates, add "A" on ID in order to handle all ID's as characters
        id.samplingtime <- as.data.frame(cbind(paste0(i,".",tree.n$id, ".A"), tree.n$dtimes)) # IDs and their sampling times in the transmission network

        write.csv(id.samplingtime,file=paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))

      }
    }

    IDs.transm <- num.i # vector of seeds with at least 2 transmissions



    ### Binding all transmission trees together ###
    ###############################################

    # Make a list of transmission trees

    trees <- list() # list of all transmission trees

    for(j in 1:length(IDs.transm)){
      p <- IDs.transm[j]
      tr <- read.tree(file = paste0(sub.dir.rename,"/tree.model1.seed_",p,".nwk"))
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
        colitem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))
      }
      else{

        ritem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))
        colitem <- rbind(colitem, ritem)
      }
    }


    # Prepare to run the sequence simulation

    seq.rand <- 1 # first sequence

    n.tr <- 1 # number of transmission tree

    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste0("R/Misc/David/Parallel_Simpact/hiv.seq.C.pol.j.fasta"),paste0(sub.dir.rename,"/seed.seq.bis.nwk"))
    # add the number of tree in the file and
    write(n.tr,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree

    write.tree(resolved.combined.tree,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)

    file.rename(from = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), to = paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk"))

    out.seq.gen.file <- paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta")
    in.seq.gen.file <- paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk")

    system(paste0("R/Misc/David/Parallel_Simpact/seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230  -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< ",in.seq.gen.file," -z",seedid," > ", out.seq.gen.file))

    # Function to tranform dates in named vector to be handled by treedater

    dates.Transform.NamedVector  <- function(dates=dates){

      dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # dates
      names(dates.val) <- as.character(dates$V1) # names are the names of the tips

      return(dates.val)
    }

    # 4.1. Construct phylogenetic tree

    # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree


    out.fast.tree.file <- paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta.tree")

    system(paste("R/Misc/David/Parallel_Simpact/FastTree -gtr -nt < ", out.seq.gen.file, paste0("> ", out.fast.tree.file)))


    samp.dates <- colitem # call the data table of dates

    tree.const <- read.tree(out.fast.tree.file)

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


    write.tree(dater.tree, file = paste0(sub.dir.rename,"/calibrated.tree.nwk")) # for further manipulation in an easy format


    ## Count transmission events per each transmission network ##
    #############################################################


    trans.net <- simpact.trans.net

    transmissions.vec.i <- list()
    time.calendar.vec.i <- list()

    for (i in 1:length(IDs.transm)){

      p <- IDs.transm[i]

      trans.net.i <- simpact.trans.net[[p]]

      Infec.time.i <- trans.net.i$InfecTime + 1977

      min.val = 1987
      max.val = datalist$itable$population.simtime[1] + 1977 # round(max(trans.net.i$InfecTime)) + 1977

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
        if(is.null(numb.i)==TRUE){
          numb.tra <- c(numb.tra, 0)
        }else{
          numb.tra <- c(numb.tra, numb.i)
        }
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
    max.val = datalist$itable$population.simtime[1] + 1977 # round(max(trans.net.i$InfecTime)) + 1977

    step.int=1

    d <- (max.val-min.val)/step.int

    int.node.vec <- vector()
    for (i in 1:d) {
      inf <- 1986+i
      sup <- 1987+i
      int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt >= inf]
      if(is.null(int.node.age.i)==TRUE){
        int.node.vec <- c(int.node.vec,0)
      }else{
        int.node.vec <- c(int.node.vec,length(int.node.age.i))
      }

    }


    # Transmissison events and internal nodes within same time intervals

    numb.tra <- trans.sum
    int.node.vec <- int.node.vec


    # x <- i.vec
    # plot(x, int.node.vec, type="b", col="red", lwd=2,
    #      xlab = "Calendar time",
    #      ylab = "Count") # 1 > 1
    # lines(x, numb.tra, col='green3', type='b', lwd=2)
    # legend("topleft", legend = c("Internal nodes", "Transmission events"),
    #        col=c("red","green3"), pch=1)

    ## Transmission network plotting: Union of all transmission networks ##
    #######################################################################

    network.list <- list()

    for (i in 1:length(IDs.transm)){

      p <- IDs.transm[i]

      trans.net.i <- as.data.frame(simpact.trans.net[[p]])

      trans.net.i <- trans.net.i[-1,]

      trans.net.i$id <- paste("C.",p,".",trans.net.i$id, sep = "")
      trans.net.i$parent <- paste("C.",p,".",trans.net.i$parent, sep = "")

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


    # library(phytools)

    tree.cal <- read.tree(sub.dir.rename,"/calibrated.tree.nwk")

    H <- nodeHeights(dater.tree) # similar to node.depth.edgelength(tree)


    # It's clear from a casual inspection of the matrix that each parent node height (in the right column)
    # is represented twice and only twice. Thus, if we exclude the root node (zero height),
    # we can just take the mean of H[,1].

    feature3 <- mean(sort(H[,1])[3:nrow(H)]) # Statistics


    summary.df <- c(growthrate,
                    feature3) # <- mean(sort(H[,1])[3:nrow(H)]) # Statistics  )


    # Name the columns

    features.names <- c("Pop.growthrate","meanHeightTree")

    names(summary.df) <- c(features.names)

    return(summary.df)

  }

  # else{
  #
  # }




}

