
setwd("/home/david/RSimpactHelp/R/Projects_2017/Example3PaperSimpact/")

## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,
               phyclust, DECIPHER,treedater,geiger,picante)


#######################
# Step 1: Run Simpact #
#######################

# Set up the maodel

simulation.type <- "simpact-cyan"#"maxart"#

simpact.set.simulation(simulation.type)

agedist.data.frame <- agedistr.creator(shape = 5, scale = 70)


iv <- intervention.introduced(simulation.type = simulation.type)


# Set up parameters of the model

mastermodel.input <- input.params.creator(population.simtime = 40,
                                          population.numwomen = 4000,
                                          population.nummen = 3600,

                                          mortality.normal.weibull.shape = 5, # k for shape
                                          mortality.normal.weibull.scale = 70, # gamma for scale

                                          person.eagerness.man.dist.gamma.a = 0.85, #0.1
                                          person.eagerness.man.dist.gamma.b = 70, #100#3.5#5#10#20 #170
                                          person.eagerness.woman.dist.gamma.a = 0.1,
                                          person.eagerness.woman.dist.gamma.b = 70,#100#3.5#5#10#20#170

                                          #formation.hazard.agegapry.eagerness_diff = -0.110975,
                                          dissolution.alpha_0 = -0.6, #-0.1 # baseline
                                          dissolution.alpha_4 = -0.15,

                                          conception.alpha_base = -2, #c(-5, -1.5), #c(-4, -1.5)
                                          formation.hazard.agegapry.eagerness_diff = -0.1, # c(-0.5, 0),#-0.048 #-0.110975
                                          birth.boygirlratio = 0.501, # c(0.5,0.7), # 0.5024876 #101:100
                                          hivtransmission.param.a = -1.01, # c(-2,-1) # baseline of transmission: -1.0352239

                                          simulation.type = simulation.type)




# Run the model [Not to RUN]

#
# mastermodel.output <- simpact.run(configParams = mastermodel.input,
#                                   destDir = "temp",
#                                   agedist = agedist.data.frame,
#                                   intervention = iv,
#                                   seed = 123)


## Read the mode output
# mastermodel.datalist <- readthedata(mastermodel.output)

## Save the output
# save(mastermodel.datalist, file = "DummyMaster.datalist.RData")

# Read saved output data set
master.datalist <- get(load("DummyMaster.datalist.RData"))


#################################
# Step 2: Transmission networks #
#################################


# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder2.diff.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")

# Transmission network in a form of constructing a transmission network
simpact.trans.net <- transmNetworkBuilder.diff2(datalist = master.datalist, endpoint = 40)


###############################
# Step 3: Sequence simulation #
###############################


# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]] # transmission network for i^th seed

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
    num.i <- c(num.i,i)

    # Save the transmission tree
    write.tree(tree.i, file = paste("tree.model1.seed",i,".nwk", sep = ""))

    tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))

    # # count number of trees generated (normally one)
    # numb.tr <- function(tree=tree){
    #     if(length(tree) == 4){
    #         return(1)
    #     }else{
    #         return(length(tree))
    #     }
    # }
    #
    # # Simulate sequences, this require compiled toold seq-gen
    # # random sequence which will be simulated is choosen in the pool
    # seq.rand <- sample(1:30,1) # chose one sequence in first 30's in the seed sequence pool named "seed.seq.fasta"
    seq.rand <- 3
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("HIV.Pol.gene.fasta", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingdates_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds


#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1990+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}

# 4.1. Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]
  # External IQ-TREE
  # system(" ./iqtree-omp -s HIV.Env.gene.fasta -m GTR+R4 -nt AUTO -alrt 1000 -bb 1000")
  # Consensus tree written to HIV.Env.gene.fasta.contree
  # Reading input trees file HIV.Env.gene.fasta.contree
  # Log-likelihood of consensus tree: -10565.685
  #
  # Analysis results written to:
  #   IQ-TREE report:                HIV.Env.gene.fasta.iqtree
  # Maximum-likelihood tree:       HIV.Env.gene.fasta.treefile
  # Likelihood distances:          HIV.Env.gene.fasta.mldist
  #
  # Ultrafast bootstrap approximation results written to:
  #   Split support values:          HIV.Env.gene.fasta.splits.nex
  # Consensus tree:                HIV.Env.gene.fasta.contree
  # Screen log file:               HIV.Env.gene.fasta.log

  system(paste("./iqtree-omp -s", paste("sequence_all_seed_number_",id.trans,"_model1.fasta", sep = ""), "-m GTR+R4 -nt AUTO -alrt 1000 -bb 1000"))


  }



# 4.2. Internal node optimisation requires sampled dates
IDs.transm <- c(2,3,11)
for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]

  tree.const <- read.tree(paste("sequence_all_seed_number_",id.trans,"_model1.fasta.treefile", sep = ""))

  samp.dates <- read.csv(paste("samplingdates_seed_number_",id.trans,".csv", sep = ""))

  time.samp <- dates.Transform.NamedVector(dates=samp.dates)

  tree.tips <- as.numeric(tree.const$tip.label)

  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 1000) # s is the length of sequence

  # Save the tree
  write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))
  ###########################################################################################
  ###########################################################################################
#
#   # With changed label names, make sure the vector of sampling dates is renamed as wells
#
#   # Rename the tree by adding sampling date
  # g <- tree.const
  # a <- g$tip.label # label names
  # d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
  # g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # renamed.Ord.tree.dates <- Ord.tree.dates
  # names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # # Use of library(treedater) to calibrate internal nodes
  # dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 1000) # s is the length of sequence

#   # Save the tree
#   write.tree(dater.tree.g, file = paste("calibratedTree_",id.trans,"_pref.nwk", sep = ""))



}

