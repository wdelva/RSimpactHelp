## This code is a about complete work-flow to combine sexual behavioural data simulatin,
## phylodynamics with agent-based models, in a modular way.

## The code has mainly four components but in two steps: Epidemiology & Phylodynamics
##
## 1. Simulation of HIV epidemics from the bottom to the top using Simpact
##    dynamic sexual networks (relationships and some dependent variables)
##    are simulated, HIV transmission (with some dependent variables) and disease
##    progression monitoring (viral load & CD4 count).
##
## 2. From the output records of events which happened to individuals in the previous simulation,
##    with RSimpactHelp, data which are important to simulate virus evolution are processed:
##    the transmission networks, infection and sampling times, viral load, CD4 count,
##    time of treatment and disease progression
##
## 3. Simulation of molecular evolution of the virus across the transmission networks
##    Two models were used to simulate sequences of infected individual:
##
##    Model #1: inputs are a consensus sequence and a transmission tree and only
##              the virus variability between individuals is considered
##    Model #2: considering the within host dynamics, simulations of viruses
##              hisotry are done per each individual by Wrish-Fisher model by considering,
##              individual j is infected by one virus from a pool of individual i,
##              at sampling times, a consensus sequence for each individual is obtained from
##              a pool of viruses of that individual
##
## 4. Analysing the sequence with time-stamped tools and compute summary statistics from
##    phylogenetic trees

## Note: to run all the code, there are external compiled tools ms (for coalescent) and seq-gen (for sequence simulation)
##       iq-tree (for tree construction) which are not R packages, they need to be in working space and also the file
##       of seeds sequences (one file for consensus seeds and another one for isolated viruses)

## Clean the global environment
rm(list = ls())

## Set the working environment
setwd("/home/david/Documents/CombinedFrameworkHIVToy")


## Load required R packages
pacman::p_load(ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,
               phyclust, DECIPHER)

# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("DECIPHER")



## STEP 1: SIMPACT SIMULATIONS

## Note to run again, the output is already saved for further simulations and analysis

# # Set up the model
#
# simulation.type <- "simpact-cyan"#"maxart"#
#
# simpact.set.simulation(simulation.type)
#
# agedist.data.frame <- agedistr.creator(shape = 5, scale = 70)
#
#
# iv <- intervention.introduced(simulation.type = simulation.type)
#
#
# # Set up parameters of the model
#
# mastermodel.input <- input.params.creator(population.simtime = 40,
#                                           population.numwomen = 1000,
#                                           population.nummen = 1000,
#
#                                           mortality.normal.weibull.shape = 5, # k for shape
#                                           mortality.normal.weibull.scale = 70, # gamma for scale
#
#                                           person.eagerness.man.dist.gamma.a = 0.85, #0.1
#                                           person.eagerness.man.dist.gamma.b = 70, #100#3.5#5#10#20 #170
#                                           person.eagerness.woman.dist.gamma.a = 0.1,
#                                           person.eagerness.woman.dist.gamma.b = 70,#100#3.5#5#10#20#170
#
#                                           #formation.hazard.agegapry.eagerness_diff = -0.110975,
#                                           dissolution.alpha_0 = -0.6, #-0.1 # baseline
#                                           dissolution.alpha_4 = -0.15,
#
#                                           conception.alpha_base = -2, #c(-5, -1.5), #c(-4, -1.5)
#                                           formation.hazard.agegapry.eagerness_diff = -0.1, # c(-0.5, 0),#-0.048 #-0.110975
#                                           birth.boygirlratio = 0.501, # c(0.5,0.7), # 0.5024876 #101:100
#                                           hivtransmission.param.a = -1.01, # c(-2,-1) # baseline of transmission: -1.0352239
#
#                                           simulation.type = simulation.type)
#
#
#
#
# # Run the model
#
# mastermodel.output <- simpact.run(configParams = mastermodel.input,
#                                   destDir = "temp",
#                                   agedist = agedist.data.frame,
#                                   intervention = iv,
#                                   seed = 123)
#
#
# # Read the mode output
#
# mastermodel.datalist <- readthedata(mastermodel.output)
#
# save(mastermodel.datalist, file = "master.datalist.RData")
#


## STEP 2: PHYLODYNAMICS

# Loading the output of the previous step

master.datalist <- get(load("master.datalist.RData")) #, .GlobalEnv) #load(file="master.datalist.RData")

# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/Projects_2017/transNetworkNew.R")
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")

## 2.1 Model 1

## 2.1.1 Example with only one seed and its transmission network

# Transmission network in a form of constructing a transmission network
simpact.output.dt <- transmNetworkBuilder.diff(datalist = master.datalist, endpoint = 40)

# Choose only transmission network of seed 2
transNet2 <- simpact.output.dt[[2]]

# Transform the above transmission network in a transmission tree
transm.tree2 <- trans.network2tree(transNet2)

# Save the transmission tree
write.tree(transm.tree2, file = "tree.model1.seed2.nwk")


# Simulate sequences, this require compiled toold seq-gen
# random sequence which will be simulated is choosen in the pool
seq.rand <- sample(1:30,1) # chose one sequence in first 30's in the seed sequence pool named "seed.seq.fasta"
# These sequences must be consensus sequence here where we use transmission gtree

# Bit of exercise to automate this simulation:
# Put the tree below the sequences and between sequences and trees
# insert the number of trees below.
# In our settings, only one tree is required and the virus sequence is chosen among the sequences above the tree
# The trees can then be placed in this file at the end, after a line stating how many trees there are.
# Example:
# 4 50
# Taxon1 ATCTTTGTAGTCATCGCCGTATTAGCATTCTTAGATCTAA
# Taxon2 ATCCTAGTAGTCGCTTGCGCACTAGCCTTCCGAAATCTAG
# Taxon3 ACTTCTGTGTTTACTGAGCTACTAGCTTCCCTAAATCTAG
# Taxon4 ATTCCTATATTCGCTAATTTCTTAGCTTTCCTGAATCTGG
# 1
# (((Taxon1:0.2,Taxon2:0.2):0.1,Taxon3:0.3):0.1,Taxon4:0.4);
tr <- read.tree(file = "tree.model1.seed2.nwk")

# Function to get the number of trees in a file
numb.tr <- function(tree){
  if(length(tree) == 4){
    return(1)
  }else{
    return(length(tree))
  }
}

seed.id <- transNet2$id[1]
# count number of trees generated (normally one)
n.tr <- numb.tr(tr)
# call the seed sequences - pool of viruses and rename the file
file.copy(paste("seed.seq.fasta", sep = ""),paste("seed.seq.bis.nwk", sep = ""))
# add the number of tree in the file and
write(n.tr,file = "seed.seq.bis.nwk", append = TRUE)
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tr,file = "seed.seq.bis.nwk", append = TRUE)
file.rename(from = "seed.seq.bis.nwk", to = paste("seed.seq.bis",seed.id,".nwk", sep = ""))

system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",seed.id,".nwk > sequence_all_seed_",seed.id,"_model1.fasta",sep = ""))


## 2.1.2 Example with all seed individuals and their transmission networks

# Note: transmission network with less than 3 individuals will not be considered

trans.net <- simpact.output.dt # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]]

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
    num.i <- c(num.i,i)

    # Save the transmission tree
    write.tree(tree.i, file = paste("tree.model1.seed",i,".nwk", sep = ""))

    # count number of trees generated (normally one)
    numb.tr <- function(tree=tree){
      if(length(tree) == 4){
        return(1)
      }else{
        return(length(tree))
      }
    }

    # Simulate sequences, this require compiled toold seq-gen
    # random sequence which will be simulated is choosen in the pool
    seq.rand <- sample(1:30,1) # chose one sequence in first 30's in the seed sequence pool named "seed.seq.fasta"

    tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    n.tr <- numb.tr(tree = tr)

    seed.id <- tree.n$id[1]
    # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("seed.seq.fasta", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(1,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

num.i # vector of of seeds chosen in the list of seeds

## 2.2 Model 2

## 2.2.1 Example with only one seed and its transmission network

# Transmission network in a raw form with suppllement data that
# the previous transmission network, it includes viral load, CD4, etc.
simpact.output.raw <- transNetworkNew(datalist=master.datalist, endpoint=40)

seed2.dat <- as.data.frame(simpact.output.raw[[2]])

## Within host dynamics per individual: ms > produce coalescent trees
## When sampled, each individual, many copies are taken and these copies have different demographics events
## When an individual i transmit an infection to a=individual j, we assume one virus seed the infection in new host
for(i in 1:length(seed2.dat$id)){

  don <- seed2.dat$parent[i] # further we do not consider -1 (universal infector)
  rec <- seed2.dat$id[i]

  # Viral load at infection times and sampling times
  vl.inf <- 90 # viral load of the donors at time point of infecting someone else
  vl.samp <- 90 # viral load of the recipients at sampling time
  t.inf <- 2.0 # infection time for recipients, time when donor i transmits to recipient j
  t.samp <- 4.0 # sampling time for recipients, time when recipient j has been sampled (diadnosed or removed)

  # Compartments within the body: 3 islands
  subpo1.inf <- round(vl.inf/8)
  subpo2.inf <- round(vl.inf/6)
  subpo3.inf <- vl.inf - (subpo1.inf+subpo2.inf)

  subpo1.samp <- round(vl.samp/8)
  subpo2.samp <- round(vl.samp/6)
  subpo3.samp <- vl.samp - (subpo1.samp+subpo2.samp)

  topaste.inf <- paste("./msa", vl.inf," 1 -T -t",t.inf," -I 3",subpo1.inf, subpo2.inf, subpo3.inf, "-ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -n 1 0.5 -n 2 0.5 -n 3 0.5 -g 1 0.3 -g 2 0.2 -g 3 0.3 > ")
  topaste.samp <- paste("./msa", vl.samp," 1 -T -t",t.samp," -I 3",subpo1.samp, subpo2.samp, subpo3.samp, "-ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -n 1 0.5 -n 2 0.5 -n 3 0.5 -g 1 0.3 -g 2 0.2 -g 3 0.3 > ")
  system(paste(topaste.inf, "tree_at_from_",don,"_to_",rec,".nwk", sep = "")) # transmission
  system(paste(topaste.samp, "tree_at_for_",rec,"_samp.nwk", sep = "")) # sampling
}



## Generate sequences for each individual at sampling time ##
#############################################################

## The seed or (and) its fist recipients are responsible of the subsequent infections


# (i) Transmissions by the seed [1] in XXX$id

seed.id <- seed2.dat$id[1] # Seed individual [1] in XXX$id
# rec.seed <- seed2.dat$id[2] # First transmission by the seed

rec.ids <- seed2.dat$id
par.ids <- seed2.dat$parent

for(i in 2:length(seed2.dat$id)){

  rec.id <- rec.ids[i]
  par.id <- par.ids[i]
  rec.seed <- rec.id # recipient from seed donor

  if(par.id == seed.id){ # this means this individual is infected by the seed element

    # donor coalescent tree at transmission times by the seed
    seed.tree.trans <- read.tree(paste("tree_at_from_",seed.id,"_to_",rec.seed,".nwk", sep = ""))

    # compute number of trees in the file
    numb.tr <- function(tree=tree){
      if(length(tree) == 4){
        return(1)
      }else{
        return(length(tree))
      }
    }
    numb.trees.seed <- numb.tr(tree=seed.tree.trans) # count number of trees generated (normally one)

    file.copy(paste("seed.seq.fasta", sep = ""),paste("seed.seq.bis.nwk", sep = "")) # call the seed sequences - pool of viruses

    write(numb.trees.seed,file = "seed.seq.bis.nwk", append = TRUE) # add the number of tree in the file and
    write.tree(seed.tree.trans,file = "seed.seq.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
    file.rename(from = "seed.seq.bis.nwk", to = paste("seed.seq.bis",seed.id,".nwk", sep = ""))

    seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

    # Sequence of the seed at the first tramsmission event
    system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",seed.id,".nwk > sequence_from_",seed.id,"_to_",rec.seed,".fasta",sep = ""))
  }
}




# (ii) Sampling of the seed [1] in XXX$id

# donor coalescent tree at sampling time (by the seed)
seed.tree.samp <- read.tree(paste("tree_at_for_",seed.id,"_samp.nwk", sep = "")) # tree at sampling time of the seed
numb.tree.samp <- numb.tr(tree=seed.tree.samp) # count number of trees generated (normally one)

file.copy(paste("seed.seq.fasta", sep = ""),paste("seed.seq.bis.nwk", sep = "")) # call the seed sequences - pool of viruses

write(numb.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # add the number of tree in the file and
write.tree(seed.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
file.rename(from = "seed.seq.bis.nwk", to = paste("seed.seq.bis",seed.id,".nwk", sep = ""))

seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

# Sequence of the seed at the sampling (diagnosis/removal) event

system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",seed.id,".nwk > sequence_at_samp_",seed.id,".fasta",sep = ""))


# (iii) Sampling of recipients from the seed as donor in XXX$id[2]

rec.ids <- seed2.dat$id
par.ids <- seed2.dat$parent

for(i in 2:length(seed2.dat$id)){

  rec.id <- rec.ids[i]
  par.id <- par.ids[i]
  rec.seed <- rec.id

  if(par.id == seed.id){ # if the parent (donor) is the seed

    # donor coalescent tree at sampling time for the first recipient
    seed.tree.samp <- read.tree(paste("tree_at_for_",rec.seed,"_samp.nwk", sep = "")) # tree at sampling time of the seed
    numb.tree.samp <- numb.tr(tree=seed.tree.samp) # count number of trees generated (normally one)

    file.copy(paste("sequence_from_",seed.id,"_to_",rec.seed,".fasta", sep = ""),paste("seed.seq.bis.nwk", sep = "")) # call the seed sequences - pool of viruses
    write(numb.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # add the number of tree in the file and
    write.tree(seed.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
    file.rename(from = "seed.seq.bis.nwk", to = paste("seed.seq.bis",seed.id,".nwk", sep = ""))

    seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

    # Sequence of the first recipient at the sampling (diagnosis/removal) event
    system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",seed.id,".nwk > sequence_at_samp_",rec.seed,".fasta",sep = ""))
  }
}


# Others

# We know that individual j got infection from i which was transmitted by one virus from i
# having sequences of i at transmission time, choose only one from the pool which will
# be evolving in j



# (iv) Sequences at infetion time >> sequence of the donors

# Find parent of the current donor here "h": due to above mentionned argument
parent.find <- function(h){
  for(l in 1:length(seed2.dat$id)){
    if(seed2.dat$id[l] == h){
      par <- seed2.dat$parent[l]
    }
  }
  return(par)
}

for(i in 3:length(seed2.dat$id)){ # start at 3 bcz the seed and 1st transmission already simulated
  # random number for random sequence to infect new individual
  seq.rand <- sample(1:10,1)
  # recipients - k and donors - h
  k <- seed2.dat$id[i] # 219 recipients
  h <- seed2.dat$parent[i] # 1408 donors

  prev.par <- parent.find(h)
  if(prev.par!=-1){
    # tree of the recipient
    tr.ms <- read.tree(paste("tree_at_from_",h,"_to_",k,".nwk", sep = "")) # tree of donor at infecting time
    numb.trees <- numb.tr(tr.ms)
    file.copy(paste("sequence_from_",prev.par,"_to_",h,".fasta", sep = ""),paste("Sequence_",h,".bis.nwk", sep = "")) # h is going to give one of what (s)he got previously from prev.par
    write(numb.trees,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
    write.tree(tr.ms,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
    file.rename(from = paste("Sequence_",h,".bis.nwk", sep = ""), to = paste("seed.seq.bis",h,".nwk", sep = ""))
    # file.remove()
    system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",h,".nwk > sequence_from_",h,"_to_",k,".fasta",sep = ""))
  }
}

# (v) Sequences at sampling time >> sequence of recipients

# Find parent of the current donor here "h": due to above mentionned argument
parent.find <- function(h){
  for(l in 1:length(seed2.dat$id)){
    if(seed2.dat$id[l] == h){
      par <- seed2.dat$parent[l]
    }
  }
  return(par)
}


for(i in 3:length(seed2.dat$id)){ # start at 3 bcz the seed and 1st transmission already simulated
  # random number for random sequence to infect new individual
  seq.rand <- sample(1:10,1)
  # recipients - k and donors - h
  k <- seed2.dat$id[i] # [3]=219, recipients
  h <- seed2.dat$parent[i] # [3]=1408, donors

  prev.par <- parent.find(h)

  if(prev.par!=-1){
    # tree of the recipient
    tr.ms <- read.tree(paste("tree_at_for_",k,"_samp.nwk", sep = "")) # tree of the recipient at sampling time
    numb.trees <- numb.tr(tr.ms)
    file.copy(paste("sequence_from_",prev.par,"_to_",h,".fasta", sep = ""),paste("Sequence_",h,".bis.nwk", sep = "")) # pool of seq of the donor to the current donor
    write(numb.trees,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
    write.tree(tr.ms,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
    file.rename(from = paste("Sequence_",h,".bis.nwk", sep = ""), to = paste("seed.seq.bis",h,".nwk", sep = ""))
    # file.remove()
    system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",h,".nwk > sequence_at_samp_",k,".fasta",sep = ""))
  }
}


## Construct then consensus sequences ##
########################################
#
# library(DECIPHER)
# library(phyclust)

label.names <- vector()
for(i in 1:length(seed2.dat$id)){
  k <- seed2.dat$id[i]
  dna.samp <- read.phylip(file = paste("sequence_at_samp_",k,".fasta", sep = ""), sep = ",")
  matrix.dna <- dna.samp$org.code
  dna.vec <- vector()
  for(i in 1:nrow(matrix.dna)){
    dna.vec <- c(dna.vec,matrix.dna[i,])
  }
  consensus.seq <- ConsensusSequence(DNAStringSet(dna.vec), threshold=0.8)

  write.dna(consensus.seq,file = "consensus.seq.raw.fas", format = "fasta", nbcol=-1,
            append = TRUE) # the name of seq all are "1"

  label.names <- c(label.names,k) # Individuals IDs
}

# read the consensus sequence all with name "1" and rename them accoridng to individuals IDs
d <- read.FASTA(file = "consensus.seq.raw.fas")
names(d) <- label.names

write.dna(d,file = "consensus.seq.final.model2.fas", format = "fasta", nbcol=-1)


# 2.2.2 Example with all seeds and their transmission networks

# Transmission network in a raw form with suppllement data that
# the previous transmission network, it includes viral load, CD4, etc.
simpact.output.raw <- transNetworkNew(datalist=master.datalist, endpoint=40)

for(v in 1:length(simpact.output.raw)){

  seed2.dat <- as.data.frame(simpact.output.raw[[v]])

  if(nrow(seed2.dat) >= 3){


    ## Within host dynamics per individual: ms > produce coalescent trees
    ## When sampled, each individual, many copies are taken and these copies have different demographics events
    ## When an individual i transmit an infection to a=individual j, we assume one virus seed the infection in new host
    for(i in 1:length(seed2.dat$id)){

      don <- seed2.dat$parent[i] # further we do not consider -1 (universal infector)
      rec <- seed2.dat$id[i]

      # Viral load at infection times and sampling times
      vl.inf <- 90 # viral load of the donors at time point of infecting someone else
      vl.samp <- 90 # viral load of the recipients at sampling time
      t.inf <- 2.0 # infection time for recipients, time when donor i transmits to recipient j
      t.samp <- 4.0 # sampling time for recipients, time when recipient j has been sampled (diadnosed or removed)

      # Compartments within the body: 3 islands
      subpo1.inf <- round(vl.inf/8)
      subpo2.inf <- round(vl.inf/6)
      subpo3.inf <- vl.inf - (subpo1.inf+subpo2.inf)

      subpo1.samp <- round(vl.samp/8)
      subpo2.samp <- round(vl.samp/6)
      subpo3.samp <- vl.samp - (subpo1.samp+subpo2.samp)

      topaste.inf <- paste("./msa", vl.inf," 1 -T -t",t.inf," -I 3",subpo1.inf, subpo2.inf, subpo3.inf, "-ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -n 1 0.5 -n 2 0.5 -n 3 0.5 -g 1 0.3 -g 2 0.2 -g 3 0.3 > ")
      topaste.samp <- paste("./msa", vl.samp," 1 -T -t",t.samp," -I 3",subpo1.samp, subpo2.samp, subpo3.samp, "-ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -n 1 0.5 -n 2 0.5 -n 3 0.5 -g 1 0.3 -g 2 0.2 -g 3 0.3 > ")
      system(paste(topaste.inf, "tree_at_from_",don,"_to_",rec,"_net_",v,".nwk", sep = "")) # transmission
      system(paste(topaste.samp, "tree_at_for_",rec,"_samp_net_",v,".nwk", sep = "")) # sampling
    }



    ## Generate sequences for each individual at sampling time ##
    #############################################################

    ## The seed or (and) its fist recipients are responsible of the subsequent infections


    # (i) Transmissions by the seed [1] in XXX$id

    seed.id <- seed2.dat$id[1] # Seed individual [1] in XXX$id
    # rec.seed <- seed2.dat$id[2] # First transmission by the seed

    rec.ids <- seed2.dat$id
    par.ids <- seed2.dat$parent

    for(i in 2:length(seed2.dat$id)){

      rec.id <- rec.ids[i]
      par.id <- par.ids[i]
      rec.seed <- rec.id # recipient from seed donor

      if(par.id == seed.id){ # this means this individual is infected by the seed element

        # donor coalescent tree at transmission times by the seed
        seed.tree.trans <- read.tree(paste("tree_at_from_",seed.id,"_to_",rec.seed,"_net_",v,".nwk", sep = ""))

        # compute number of trees in the file
        numb.tr <- function(tree=tree){
          if(length(tree) == 4){
            return(1)
          }else{
            return(length(tree))
          }
        }
        numb.trees.seed <- numb.tr(tree=seed.tree.trans) # count number of trees generated (normally one)

        file.copy(paste("seed.seq.fasta", sep = ""),paste("seed.seq.bis.nwk", sep = "")) # call the seed sequences - pool of viruses

        write(numb.trees.seed,file = "seed.seq.bis.nwk", append = TRUE) # add the number of tree in the file and
        write.tree(seed.tree.trans,file = "seed.seq.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
        file.rename(from = "seed.seq.bis.nwk", to = paste("seed.seq.bis",seed.id,".nwk", sep = ""))

        seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

        # Sequence of the seed at the first tramsmission event
        system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",seed.id,".nwk > sequence_from_",seed.id,"_to_",rec.seed,"_net_",v,".fasta",sep = ""))
      }
    }




    # (ii) Sampling of the seed [1] in XXX$id

    # donor coalescent tree at sampling time (by the seed)
    seed.tree.samp <- read.tree(paste("tree_at_for_",seed.id,"_samp_net_",v,".nwk", sep = "")) # tree at sampling time of the seed
    numb.tree.samp <- numb.tr(tree=seed.tree.samp) # count number of trees generated (normally one)

    file.copy(paste("seed.seq.fasta", sep = ""),paste("seed.seq.bis.nwk", sep = "")) # call the seed sequences - pool of viruses

    write(numb.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # add the number of tree in the file and
    write.tree(seed.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
    file.rename(from = "seed.seq.bis.nwk", to = paste("seed.seq.bis",seed.id,".nwk", sep = ""))

    seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

    # Sequence of the seed at the sampling (diagnosis/removal) event

    system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",seed.id,".nwk > sequence_at_samp_",seed.id,"_net_",v,".fasta",sep = ""))


    # (iii) Sampling of recipients from the seed donor in XXX$id[2]

    rec.ids <- seed2.dat$id
    par.ids <- seed2.dat$parent

    for(i in 2:length(seed2.dat$id)){

      rec.id <- rec.ids[i]
      par.id <- par.ids[i]
      rec.seed <- rec.id

      if(par.id == seed.id){

        # donor coalescent tree at sampling time for the first recipient
        seed.tree.samp <- read.tree(paste("tree_at_for_",rec.seed,"_samp_net_",v,".nwk", sep = "")) # tree at sampling time of the seed
        numb.tree.samp <- numb.tr(tree=seed.tree.samp) # count number of trees generated (normally one)

        file.copy(paste("sequence_from_",seed.id,"_to_",rec.seed,"_net_",v,".fasta", sep = ""),paste("seed.seq.bis.nwk", sep = "")) # call the seed sequences - pool of viruses
        write(numb.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # add the number of tree in the file and
        write.tree(seed.tree.samp,file = "seed.seq.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
        file.rename(from = "seed.seq.bis.nwk", to = paste("seed.seq.bis",seed.id,".nwk", sep = ""))

        seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

        # Sequence of the first recipient at the sampling (diagnosis/removal) event
        system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",seed.id,".nwk > sequence_at_samp_",rec.seed,"_net_",v,".fasta",sep = ""))
      }
    }


    # Others

    # We know that individual j got infection from i which was transmitted by one virus from i
    # having sequences of i at transmission time, choose only one from the pool which will
    # be evolving in j



    # (iv) Sequences at infetion time >> sequence of the donors

    # Find parent of the current donor here "h": due to above mentionned argument
    parent.find <- function(h){
      for(l in 1:length(seed2.dat$id)){
        if(seed2.dat$id[l] == h){
          par <- seed2.dat$parent[l]
        }
      }
      return(par)
    }

    for(i in 3:length(seed2.dat$id)){ # start at 3 bcz the seed and 1st transmission already simulated
      # random number for random sequence to infect new individual
      seq.rand <- sample(1:10,1)
      # recipients - k and donors - h
      k <- seed2.dat$id[i] # 219 recipients
      h <- seed2.dat$parent[i] # 1408 donors

      prev.par <- parent.find(h)
      if(prev.par!=-1){
        # tree of the recipient
        tr.ms <- read.tree(paste("tree_at_from_",h,"_to_",k,"_net_",v,".nwk", sep = "")) # tree of donor at infecting time
        numb.trees <- numb.tr(tr.ms)
        file.copy(paste("sequence_from_",prev.par,"_to_",h,"_net_",v,".fasta", sep = ""),paste("Sequence_",h,".bis.nwk", sep = "")) # h is going to give one of what (s)he got previously from prev.par
        write(numb.trees,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        write.tree(tr.ms,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        file.rename(from = paste("Sequence_",h,".bis.nwk", sep = ""), to = paste("seed.seq.bis",h,".nwk", sep = ""))
        # file.remove()
        system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",h,".nwk > sequence_from_",h,"_to_",k,"_net_",v,".fasta",sep = ""))
      }
    }

    # (v) Sequences at sampling time >> sequence of recipients

    # Find parent of the current donor here "h": due to above mentionned argument
    parent.find <- function(h){
      for(l in 1:length(seed2.dat$id)){
        if(seed2.dat$id[l] == h){
          par <- seed2.dat$parent[l]
        }
      }
      return(par)
    }


    for(i in 3:length(seed2.dat$id)){ # start at 3 bcz the seed and 1st transmission already simulated
      # random number for random sequence to infect new individual
      seq.rand <- sample(1:10,1)
      # recipients - k and donors - h
      k <- seed2.dat$id[i] # [3]=219, recipients
      h <- seed2.dat$parent[i] # [3]=1408, donors

      prev.par <- parent.find(h)

      if(prev.par!=-1){
        # tree of the recipient
        tr.ms <- read.tree(paste("tree_at_for_",k,"_samp_net_",v,".nwk", sep = "")) # tree of the recipient at sampling time
        numb.trees <- numb.tr(tr.ms)
        file.copy(paste("sequence_from_",prev.par,"_to_",h,"_net_",v,".fasta", sep = ""),paste("Sequence_",h,".bis.nwk", sep = "")) # pool of seq of the donor to the current donor
        write(numb.trees,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        write.tree(tr.ms,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        file.rename(from = paste("Sequence_",h,".bis.nwk", sep = ""), to = paste("seed.seq.bis",h,".nwk", sep = ""))
        # file.remove()
        system(paste("./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k",seq.rand,"< seed.seq.bis",h,".nwk > sequence_at_samp_",k,"_net_",v,".fasta",sep = ""))
      }
    }


    ## Construct then consensus sequences ##
    ########################################
    #
    # library(DECIPHER)
    # library(phyclust)

    label.names <- vector()
    for(i in 1:length(seed2.dat$id)){
      k <- seed2.dat$id[i]
      dna.samp <- read.phylip(file = paste("sequence_at_samp_",k,"_net_",v,".fasta", sep = ""), sep = ",")
      matrix.dna <- dna.samp$org.code
      dna.vec <- vector()
      for(i in 1:nrow(matrix.dna)){
        dna.vec <- c(dna.vec,matrix.dna[i,])
      }
      consensus.seq <- ConsensusSequence(DNAStringSet(dna.vec), threshold=0.8)
      write.dna(consensus.seq,file = paste("consensus.seq.raw_net_",v,".fas", sep = ""), format = "fasta", nbcol=-1,
                append = TRUE) # the name of seq all are "1"

      label.names <- c(label.names,k) # Individuals IDs
    }

    # read the consensus sequence all with name "1" and rename them accoridng to individuals IDs
    d <- read.FASTA(file = paste("consensus.seq.raw_net_",v,".fas", sep = ""))
    names(d) <- label.names

    write.dna(d,file = paste("consensus.seq.final_net_",v,"_model2.fas", sep = ""), format = "fasta", nbcol=-1)


  } # end loop of chosen transmission networks (with at least 3 people)

} # end loop of transmission networks

