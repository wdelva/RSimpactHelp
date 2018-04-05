
phylo.simpact.wrapper <- function(inputvector = inputvector){

  pacman::p_load(RSimpactCyan,RSimpactHelper,devtools,Rcpp,ape,expoTree,data.table,
                 readr,dplyr,adephylo,geiger,picante,igraph,ggplot2,
                 phyloTop,phytools,lme4,data.table, treedater, phangorn)

  # ## For sub-optimal sequence coverage

  # Kick-start the missing sequence data simulation study:

  # Let's start with verifying that we have parameter combinations that produce "useful" output.
  # i.e. output that is more or less in line with the demography, behaviour and epidemic that we want to consider.


  age.distr <- agedistr.creator(shape = 5, scale = 65)

  #source the input files var name <- cfg.list <- cfg
  source("R/Misc/David/Parallel_Simpact/min.input.creator.R")


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
  results <- simpact.run(configParams = cfg,
                         destDir = sub.dir.rename,
                         agedist = age.distr,
                         seed = seedid,
                         intervention = intervention)

  datalist <- readthedata(results)


  # Let's examine the output:
  #source the summary functions
  source("R/Misc/David/Parallel_Simpact/transmNetworkBuilder.diff3.R")
  source("R/Misc/David/Parallel_Simpact/00-Functions.R")

  growthrate <- pop.growth.calculator(datalist = datalist,
                                      timewindow = c(0, datalist$itable$population.simtime[1]))


  simpact.trans.net <- transmNetworkBuilder.diff3(datalist = datalist, endpoint = 40)


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
  file.copy(paste0("R/Misc/David/Parallel_Simpact/hiv.seq.B.pol.j.fasta"),paste0(sub.dir.rename,"/seed.seq.bis.nwk"))
  # add the number of tree in the file and
  write(n.tr,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)  # n.tr
  # the tree, to prepare the file to simulate the evolution of the virus across the tree

  write.tree(resolved.combined.tree,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)

  file.rename(from = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), to = paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk"))

  out.seq.gen.file <- paste0(sub.dir.rename,"/B.Epidemic_seed.seq.bis.sim.nwk.fasta")
  in.seq.gen.file <- paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk")

  system(paste0("R/Misc/David/Parallel_Simpact/seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230  -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< ",in.seq.gen.file," -z",seed," > ", out.seq.gen.file))

  # Function to tranform dates in named vector to be handled by treedater

  dates.Transform.NamedVector  <- function(dates=dates){

    dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # dates
    names(dates.val) <- as.character(dates$V1) # names are the names of the tips

    return(dates.val)
  }

  # 4.1. Construct phylogenetic tree

  # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree


  out.fast.tree.file <- paste0(sub.dir.rename,"/B.Epidemic_seed.seq.bis.sim.nwk.fasta.tree")

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


  # library(phytools)

  tree.cal <- read.tree(sub.dir.rename,"/calibrated.tree.nwk")

  H <- nodeHeights(dater.tree) # similar to node.depth.edgelength(tree)


  # It's clear from a casual inspection of the matrix that each parent node height (in the right column)
  # is represented twice and only twice. Thus, if we exclude the root node (zero height),
  # we can just take the mean of H[,1].

  feature3 <- mean(sort(H[,1])[3:nrow(H)]) # Statistics


  summary.df <- c(growthrate,
                   feature3) # <- mean(sort(H[,1])[3:nrow(H)]) # Statistics  )

  #delete all the file
  #unlink(paste0(sub.dir.rename,"/"), recursive = TRUE)

  # Name the columns

  features.names <- c("Pop.growthrate","meanHeightTree")

  names(summary.df) <- c(features.names)

  return(summary.df)

}

