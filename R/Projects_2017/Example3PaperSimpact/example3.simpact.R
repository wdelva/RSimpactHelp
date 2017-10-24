
setwd("/home/david/RSimpactHelp/R/Projects_2017/Example3PaperSimpact/")

## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,
               phyclust, DECIPHER,treedater,geiger,picante)


##############################
# Step 1: Simpact simulation #
##############################



# Read saved output data set
master.datalist <- get(load("DummyMaster.datalist.RData"))



#################################
# Step 2: Transmission networks #
#################################

# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")

simpact.trans.net <- transmNetworkBuilder.diff(datalist = master.datalist, endpoint = 40)

tran.net <- simpact.trans.net

################################
# Step 3: Sequence simulations #
################################

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds

IDs.transm <- c(2, 3, 11)

# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1990+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}


################################
# Step 4: Trees reconstruction #
################################


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
  
  # system(paste("./iqtree-omp -s", paste("sequence_all_seed_number_",id.trans,"_model1.fasta", sep = ""), "-m GTR+R4 -nt AUTO -alrt 1000 -bb 1000"))
  
  
}



# 4.2. Internal node optimisation requires sampled dates

IDs.transm <- c(2,3,11)

for (j in 1:length(IDs.transm)){
  
  id.trans <- IDs.transm[j]
  
  tree.const <- read.tree(paste("sequence_all_seed_number_",id.trans,"_model1.fasta.treefile", sep = ""))
  
  samp.dates <- read.csv(paste("samplingdates_seed_number_",id.trans,".csv", sep = ""))
  
  time.samp <- dates.Transform.NamedVector(dates=samp.dates)
  
  tree.tips <- as.numeric(tree.const$tip.label)
  
  # Sequence labels ordered in the same order of the tree tips labels
  
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

# Chose tree of seed 3 with  314 tips and 313 internal nodes.
tree.const <- read.tree(paste("sequence_all_seed_number_",3,"_model1.fasta.treefile", sep = ""))
samp.dates <- read.csv(paste("samplingdates_seed_number_",3,".csv", sep = ""))

time.samp <- dates.Transform.NamedVector(dates=samp.dates)

tree.tips <- as.numeric(tree.const$tip.label)

# Sequence labels ordered in the same order of the tree tips labels

Ord.tree.dates <- vector()
for(i in 1:length(tree.tips)){
  for(j in 1:length(time.samp)){
    if(tree.tips[i] == samp.dates$V1[j]){
      Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
    }
  }
}

g <- tree.const
a <- g$tip.label # label names
d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates

renamed.Ord.tree.dates <- Ord.tree.dates
names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates

# Use of library(treedater) to calibrate internal nodes
dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 1000) # s is the length of sequence


# Node age with picante package

N <- node.age(dater.tree.g)

# potential transmission times
int.node.age <- N$Ti # internal nodes ages
##########################################


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


## Transmission network of seed 3

tra.net.3 <- tran.net[[3]]

tra.net.3$dtimes <- abs(tra.net.3$dtimes-40)+1980 #(endpoint=40)
tra.net.3$itimes <- abs(tra.net.3$itimes-40)+1980 #(endpoint=40) 10>1990, -> +1980

min.val = 1990
max.val = round(max(tra.net.3$itimes))


step.int=1
d <- (max.val-min.val)/step.int

dat.f.trans <- as.data.frame(tra.net.3)
dt.node.age.dt <- int.node.age

numb.tra <- vector()
i.vec <- vector()
int.node.vec <- vector()
for (i in 1:d) {
  inf <- 1989+i
  sup <- 1990+i
  dat.f.trans.i <- dat.f.trans[which(dat.f.trans$itimes <= sup & dat.f.trans$itimes  >= inf),]
  numb.i <- nrow(dat.f.trans.i)
  numb.tra <- c(numb.tra, numb.i)
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt >= inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}



graph.build <- as.data.frame(trans.net[[3]])

graph.build[,4] <- as.character(graph.build$parent) # donors
graph.build[,3] <- as.character(graph.build$id) # recipients
gag = as.matrix(graph.build)
ga.graph = graph.edgelist(gag[,4:3])

V(ga.graph)$color <- "red"

transNet.yrs.Old <- delete.vertices(ga.graph, "-1")



### Plot figures 
#################

# 1. Transmission network from simpact                           # 1 #
plot.igraph(transNet.yrs.Old, edge.arrow.size=0.1, vertex.size=5,
            edge.color="black",
            asp = 1,
            xlim = c(-1, 2),
            ylim = c(-0.5,0.5),
            vertex.frame.color="black",
            vertex.label.color="black",
            #vertex.label = NULL,
            layout = layout_with_kk,
            edge.width = 1,
            vertex.label.cex=0.1,
            vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0
            #main = "True transmission network"
)

# 2. Phylogenetic tree

plot(dater.tree.g, show.tip.label=FALSE,
     edge.width=1,
     edge.color="blue") # Try a few different settings!
axisPhylo(backward = FALSE)


# 3. Transmission events and internal nodes
x <- i.vec
plot(x, int.node.vec, type="b", col="red", lwd=2,
     xlab = "Calendar time",
     ylab = "Count") # 1 > 1
lines(x, numb.tra, col='green3', type='b', lwd=2)
legend("topleft", legend = c("Internal nodes", "Transmission events"),
       col=c("red","green3"), pch=1)



