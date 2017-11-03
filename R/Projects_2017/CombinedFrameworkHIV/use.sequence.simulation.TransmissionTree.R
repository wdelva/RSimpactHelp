#' Simulate evolutionary sequences accross transmission networks
#'
#' By considering transmission trees, the later guide the simulation of
#' HIV consensus sequences at sampling times
#' Within the working directory must be seq-gen compiled tool and a file of alignment of seeds sequences
#' Inputs is a transmission networks with viral load, rate of change of viral load, and time from infection to sampling
#'
#'
#'


pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings, dplyr, adephylo,
               phyclust, DECIPHER,treedater,geiger,picante)


setwd("/home/david/RSimpactHelp/R/Projects_2017/CombinedFrameworkHIV/")

# Read saved output data set
master.datalist <- get(load("master.datalist.RData"))

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



## Functions for external tips distances

extern.dist.tips <- function(phylo.tree){
  # Branch length of the tips
  tree <- phylo.tree
  # tree$edge.length<-round(tree$edge.length,3)
  tree$edge.length<-tree$edge.length # remove rounding
  n<-length(tree$tip.label)
  ee<-setNames(tree$edge.length[sapply(1:n,function(x,y)
    which(y==x),y=tree$edge[,2])],tree$tip.label)

  d <- as.data.frame(ee)
  branc.leng.tips <- as.data.frame(d$ee) # yes we can find them in phylo.tree$edge.length but here they are not rounded whereas
  # with ee are rounded with 3 digit
  tips.id <- as.data.frame(as.character(tree$tip.label))

  tips.id <- tree$tip.label
  tips.branch.leng.raw <- cbind(tips.id,branc.leng.tips)

  names(tips.branch.leng.raw) <- c("Ind", "branch.len")

  tips.branch.leng <- tips.branch.leng.raw
  return(tips.branch.leng)
}



num.i #  2  7 12 13

net.trans.2 <- trans.net[[2]] # 72
net.trans.7 <- trans.net[[7]] # 76
net.trans.12 <- trans.net[[12]] # 4
net.trans.13 <- trans.net[[13]] # 28

tr.2 <- read.tree(paste("tree.model1.seed",2,".nwk", sep = ""))
tr.7 <- read.tree(paste("tree.model1.seed",7,".nwk", sep = ""))
tr.12 <- read.tree(paste("tree.model1.seed",12,".nwk", sep = ""))
tr.13 <- read.tree(paste("tree.model1.seed",13,".nwk", sep = ""))



# Distance from tips to root
root.tr2 <- distRoot(tr.2)
root.tr7 <- distRoot(tr.7)
root.tr12 <- distRoot(tr.12)
root.tr13 <- distRoot(tr.13)

# External distance of tips
ext.tr2 <- extern.dist.tips(tr.2)
ext.tr7 <- extern.dist.tips(tr.7)
ext.tr12 <- extern.dist.tips(tr.12)
ext.tr13 <- extern.dist.tips(tr.13)

# Diff /b/ infec. and samp. times
dt.2 <- net.trans.2$itimes-net.trans.2$dtimes
dt.7 <- net.trans.7$itimes-net.trans.7$dtimes
dt.12 <- net.trans.12$itimes-net.trans.12$dtimes
dt.13 <- net.trans.13$itimes-net.trans.13$dtimes


# chose 12

# tips 1 and 3

rt1 <- root.tr12[[2]]
rt3 <- root.tr12[[3]]

d1 <- rt1 - rt3

ext1 <- ext.tr12$branch.len[2]
ext3 <- ext.tr12$branch.len[3]

d2 <- ext1 - ext3


x=as.data.table(net.trans.12)

# Graphic representation
x$itimes <- 30-x$itimes
x$dtimes <- 30-x$dtimes

a <- x$itimes
b <- x$dtimes

x$itimes <- b
x$dtimes <- a

trepi <- list()

trepi$itimes <- x$itimes
trepi$dtimes <- x$dtimes
trepi$id <- x$id
trepi$parent <- x$parent

t <- trans.network2tree(transnetwork = trepi)

source("/home/david/RSimpactHelp/R/epi2tree2.R")
source("/home/david/RSimpactHelp/R/epi2tree3.R")

t.r <- epi2tree3(net.trans.12)

library(adephylo)
distRoot(t.r)
