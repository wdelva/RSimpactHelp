#' Simulate evolutionary sequences accross transmission networks
#'
#' By considering transmission trees, the later guide the simulation of
#' HIV consensus sequences at sampling times
#' Within the working directory must be seq-gen compiled tool and a file of alignment of seeds sequences
#' Inputs is a transmission networks with viral load, rate of change of viral load, and time from infection to sampling
#'
seed=123

source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")


# Transmission network in a form of constructing a transmission network
simpact.output.dt <- transmNetworkBuilder.diff(datalist = master.datalist, endpoint = 40)

# Note: will consider transmission network with at least 3 individuals

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
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("seed.seq.fasta", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(1,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3,0.2,0.2,0.3 -a 0.9 -g 4 -r 1 1 1 1 1 1 -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed
  }
}
