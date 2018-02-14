# Merging phylogenetic trees

# Statement: In Simpact, for each seed is associated a transmission network
# Sequence data are simulated through transmission trees which is correct since in real world
# we have different transmission chains.
# However, sequence analysis in real world does not capture this reality
# Thus, during  sequence simulation in Simpact we need to bind all transmission trees
# and have a set of sequence data for all, not subdivided data sets.


library(ape)
library(geiger)

tree1 <- rtree(4)
tree2 <- rtree(5)
tree3 <- rtree(20)

# bind.tree(x, y, where = "root")

combined.tree <- bind.tree(tree1, tree2, where=1) 
plot(combined.tree)


# For consistency checking: external branch length must be the same for combined tree and single trees

## OK: Validated


# Functions for external tips distances
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

t.1 <- extern.dist.tips(tree1)
t.2 <- extern.dist.tips(tree2)
t.comb <- extern.dist.tips(combined.tree)


setwd("~/Desktop/Test_Combined_Transm_Chains")

datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalist.RData"))

library(RSimpactCyan)
library(RSimpactHelper)
library(igraph)
library(dplyr)

# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff3.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
source("/home/david/RSimpactHelp/R/epi2tree2.R")

# simpact.trans.net <- transmNetworkBuilder.diff3(datalist = datalist, endpoint = 40)
# 
# save(simpact.trans.net, file = "simpact.trans.net.RData")

simpact.trans.net <- get(load("simpact.trans.net.RData"))



seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){
  
  tree.n <- trans.net[[i]] # transmission network for i^th seed
  
  if(nrow(as.data.frame(tree.n)) >= 3){
 
    tree.i <- trans.network2tree(transnetwork = tree.n)
    tree.j <- tree.i
    tree.j$tip.label <- paste(i,".", tree.j$tip.label, sep = "")
    
    
    # Save the transmission tree
    write.tree(tree.j, file = paste("tree.model1.seed",i,".nwk", sep = ""))
    
    tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    
    
    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(paste(i,".", tree.n$id, sep = ""), tree.n$SampTime)) # IDs and their samling times in the transmission network
    
    write.csv(id.samplingtime,file=paste("samplingtimes_seed_",i,".csv", sep = ""))
    
  }