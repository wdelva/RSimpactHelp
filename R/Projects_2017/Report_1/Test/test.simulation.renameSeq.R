

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,dplyr,
               phyclust, DECIPHER,treedater,geiger,picante)


setwd("~/RSimpactHelp/R/Projects_2017/Report_1/")

# Read saved output data set
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalist.RData"))

# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff2.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")


simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)


setwd("~/RSimpactHelp/R/Projects_2017/Report_1/Test/")

# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # 1977+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}


# Transmission trees and Sequences

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

    tr.labels <- tr$

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

    # call the seed sequences and rename the file
    file.copy(paste("C.pool.gene.pol.fasta", sep = ""),paste("hiv.seq.C.pol.j.fasta",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("hiv.seq.C.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("hiv.seq.C.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("hiv.seq.C.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.C.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.393,0.172,0.223,0.212  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -n1 -k",seq.rand,"<hiv.seq.C.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >C.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

  }
}
