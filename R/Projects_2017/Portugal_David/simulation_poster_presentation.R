# Simulations for Portugal Poster presentation
rm(list = ls())
setwd("/home/david/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/")

pacman::p_load(ape,expoTree,data.table,phylosim,RSimpactHelper)


##### Section 1: Transmission tree/network and sequence simulation #####
########################################################################


# 1.1. Load data from master model output

# track records of dynamic sexual networks simulated using ABM
# with SImpact
master.datalist <- get(load("master.datalist.RData")) #, .GlobalEnv) #load(file="master.datalist.RData")

# head(master.datalist)

datalist <- master.datalist


# 1.2. Construct transmission epi object to be handled by epi2tree function
# to build a transmission tree

# With many seed IDs, each has its own transmission network
# let take, seed = 2, 7 & 13


# Transmission network with different sampling/removal times
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.baseline.R")
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")
transm.ls <- transmNetworkBuilder.diff(datalist = datalist,endpoint = 40)

# epi object
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
transnetwork <- transm.ls[[2]] # seed 2 considered
epi.tree <- trans.network2tree(transnetwork = transnetwork)


# 1.3. Simulate consensus sequence

source("/home/david/RSimpactHelp/R/sequence.simulation.R")

library(readr)

hiv.seq.raw <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
## Remove the break line in the string of DNA
clean.hiv.seq <-  gsub("\n", "", hiv.seq.raw)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position

hiv.seq.env <- substr(clean.hiv.seq, 6172,8742) # true c(6172,8742)

# Calculate the nucleotides frequencies
#
#       # library(Biostrings)
#       seq1 = DNAString(hiv.seq.env) # nulceotides
#
#       # Chech frequencies
#       freq <- letterFrequency(seq1, letters = c("A", "C", "G", "T"))/nchar(hiv.seq.env)
#       # > letterFrequency(gag, letters = c("A", "C", "T", "G"))/nchar(hiv.seq.gag)

freq <- c(0.3353293,0.2035928,0.2628077,0.1982701)
rate <- list("a"=0.2, "b"=0.6, "c"=0.12,"d"=0.001, "e"=0.25, "f"=0.24)

# Sequence simulation
sim <- sequence.simulation(transtree = epi.tree, seedSeq = hiv.seq.env, alpha = 0.90,
                           rate.list = rate, base.freq = freq)
saveAlignment.PhyloSim(sim,file = paste("HIVSeq_fullNetwork.fasta",sep=""), skip.internal = TRUE, paranoid = TRUE)
saveAlignment.PhyloSim(sim,file = paste("HIVSeq_fullNetworkNode.fasta",sep=""), paranoid = TRUE)


##### Section 2: True trees and transmission networks #####
###########################################################

# 2.1. True transmission network: transnetwork on 1.2

# Visualise and summarise
library(igraph)
library(network)


graph.net <- as.data.frame(cbind(transnetwork$parent, transnetwork$id))

graph.build <- graph.net

graph.build[,1] <- as.character(graph.net[,1]) # donors
graph.build[,2] <- as.character(graph.net[,2]) # recipients
gag = as.matrix(graph.build)
ga.graph = graph.edgelist(gag[,1:2])

V(ga.graph)$color <- "red"

# Transmission network from simpact
plot.igraph(ga.graph, edge.arrow.size=0, vertex.size=7,
            vertex.frame.color="black", vertex.label.color="black",
            vertex.label = V(ga.graph)$name,layout = layout_with_kk,
            vertex.label.cex=0.6, vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0,
            main = "Transmission network from Simpact")

# Histogram of the network degree
hist(degree(ga.graph))

# Diameter is essentially the longest path between two vertices
test.graph <- ga.graph

# How large is the network (I know I set this when we made the network,
# but what if I had not?)

test.graph      # Tells me that it is an IGRAPH object with 100 nodes and 197 links,
# made with the Barabasi algorithm
V(test.graph)   # gives the vertex sequence
E(test.graph)   # gives the edge sequence (edge list)

# Additionally I will use the NetIndices package,
# since its function "GenInd()" outputs several network properties
library(NetIndices)

# The "GenInd()" function requires an input of an adjacency matrix
test.graph.adj<-get.adjacency(test.graph,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

test.graph.properties<-GenInd(test.graph.adj)

# The function output consists of 10 network properties.
# I will consider five of them here:

test.graph.properties$N            #number of nodes

test.graph.properties$Ltot        #number of links

test.graph.properties$LD        #link density (average # of links per node)

test.graph.properties$C            #the connectance of the graph
# This function measures connectance as L/(N*(N-1)) where L is links, and N is nodes
# Connectance can also be calculated as L/(N^2)

# The degree of a node refers to the number of links associated with a node.
# Degree can be measured as the links going in ("in degree"), out ("out degree"), or both.
# The degree() function takes a graph input and gives the degree of specified nodes.
# With the argument "v=V(graph)" you tell the function to give the degree of all nodes in the graph,
# while the "mode" argument specifies in, out, or both.

in.deg.testgraph<-degree(test.graph,v=V(test.graph),mode="in")
out.deg.testgraph<-degree(test.graph,v=V(test.graph),mode="out")
all.deg.testgraph<-degree(test.graph,v=V(test.graph),mode="all")

# Degree distribution is the cumulative frequency of nodes with a given degree
# this, like degree() can be specified as "in", "out", or "all"
deg.distr<-degree.distribution(test.graph,cumulative=T,mode="all")

# Using the power.law.fit() function I can fit a power law to the degree distribution
power<-power.law.fit(all.deg.testgraph)

# The output of the power.law.fit() function tells me what the exponent of the power law is ($alpha)
# and the log-likelihood of the parameters used to fit the power law distribution ($logLik)
# Also, it performs a Kolmogov-Smirnov test to test whether the given degree distribution could have
# been drawn from the fitted power law distribution.
# The function thus gives me the test statistic ($KS.stat) and p-vaule ($KS.p) for that test

# Then I can plot the degree distribution
plot(deg.distr,log="xy",
     ylim=c(.01,10),
     bg="black",pch=21,
     xlab="Degree",
     ylab="Cumulative Frequency")

# And the expected power law distribution
lines(1:20,10*(1:20)^((-power$alpha)+1))


# Graphs typically have a Poisson distribution (if they are random),
# power law (preferential attachment), or truncated power law (many real networks) degree distribution



# Another way
##############

dd = degree.distribution(test.graph, mode = "all", cumulative = FALSE)
# Plot degree distribution

# write a function to plot the degree distribution
plot_degree_distribution = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "Degree Distribution")
}



plot_degree_distribution(test.graph)


# plot and fit the power law distribution
fit_power_law = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  beta = -cozf[[1]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("Beta =", round(beta, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}


fit_power_law(test.graph)

# Ro
# Diameter is essentially the longest path between two vertices
diameter(test.graph)
# Gives me the length of the diameter while

nodes.diameter<-get.diameter(test.graph)
# Gives me the labels for each node that participates in the diameter

# I can look at the diameter graphically also
# First I will define the node and edge attributes
V(test.graph)$color<-"skyblue"
# I want all the nodes to be skyblue
V(test.graph)$size<-7
# I want all the nodes to be size=7
V(test.graph)[nodes.diameter]$color<-"darkgreen"
V(test.graph)[nodes.diameter]$size<-10
V(test.graph)[nodes.diameter]$label.color<-"white"
# but the nodes in the diameter should be darkgreen and larger than the rest
# with a white label instead of black
# this will make the diameter pop out of the larger network
E(test.graph)$color<-"grey"
# all non-diameter edges will be grey
E(test.graph,path=nodes.diameter)$color<-"darkgreen"
E(test.graph,path=nodes.diameter)$width<-2
# Edges in the diameter will be darkgreen and a little extra wide

# If you do not set the attributes of all of the nodes and edges then it will
# default such that you only see what you have defined

# Now when I plot the diameter will be larger than everything else, and darkgreen instead
# of grey/blue
par(mar=c(.1,.1,.1,.1))
plot.igraph(test.graph,
            layout=layout.fruchterman.reingold,
            vertex.label.cex=.5,
            edge.arrow.size=.5)

# Clustering coefficient is the proportion of
# a nodes neighbors that can be reached by other neighbors
# in igraph this property is apparently called "transitivity"

transitivity(test.graph)
# gives the clustering coefficient of the whole network

transitivity(test.graph,type="local")
# gives the clustering coefficient of each node

# Betweenness is the number of shortest paths between two nodes that go through each node of interest

graph.betweenness<-betweenness(test.graph,v=V(test.graph))
graph.edge.betweenness<-edge.betweenness(test.graph,e=E(test.graph))

# Closeness refers to how connected a node is to its neighbors

graph.closeness<-closeness(test.graph,vids=V(test.graph))

# Clustering coefficient, betweenness, and closeness
# all describe the small world properties of the network.
# A network with small world properties is one in which
# it takes a relatively short path to get from one node to the next
# (e.g., six degrees of separation)


# 2.2. True tree: epi.tree on 1.2

plot(epi.tree)
axisPhylo() # add timescale

# phylogenetic tree characteristics

library(phangorn)
library(ape)
library(apTreeshape)

summary(epi.tree)

epi.tree.shape <- as.treeshape(epi.tree)

summary(epi.tree.shape)

##### Section 3: Phylogenetic tree reconstruction from sequence data #####
##########################################################################

seq.sim <- read.FASTA("HIVSeq_fullNetwork.fasta")
tree.dat <- phyDat(seq.sim, type = "DNA")
tree.ml <- dist.ml(tree.dat)
tree.sim <- upgma(tree.ml)

phylo.tree.full <- tree.sim


# Comparison 1: true tree & reconstructed full tree / transmission networks

# Comparison 1.1: Trees

# This function computes the topological distance between two phylogenetic trees
# or among trees in a list (if y = NULL using different methods.
# The trees are always considered as unrooted.
dist.topo(epi.tree, phylo.tree.full, method = "PH85")

# This function plots two trees face to face with the links if specified.
# It is possible to rotate the branches of each tree around the nodes by clicking.
#plot
cophyloplot(epi.tree, phylo.tree.full, assoc=association, length.line=4, space=28, gap=3)
#plot with rotations
cophyloplot(epi.tree, phylo.tree.full, assoc=association, length.line=4, space=28, gap=3, rotate=TRUE)

library(distory)
#
# phylo.diff plots two trees side by side,
# highlighting edges unique to each tree in red.
#
# distinct.edges finds the edges present in the first argument
# not in the second.
phylo.diff(epi.tree,phylo.tree.full)

library(phangorn)

# treedist computes different tree distance methods and RF.dist the Robinson-Foulds
# or symmetric distance. The Robinson-Foulds distance only depends on the toplogy of the trees.
# If edge weights should be considered wRF.dist calculates
# the weighted RF distance (Robinson & Foulds 1981). and KF.dist calculates
# the branch score distance (Kuhner & Felsenstein 1994).
# path.dist computes the path difference metric as described in Steel and Penny 1993).
# sprdist computes the approximate SPR distance
# (Oliveira Martins et al. 2008, de Oliveira Martins 2016).

tree1 <- epi.tree
tree2 <- phylo.tree.full

dist.fulltree <- treedist(tree1, tree2, check.labels = TRUE)

sprdist(tree1, tree2)

SPR.dist(tree1, tree2)

RF.dist(tree1, tree2, normalize = FALSE, check.labels = TRUE,
        rooted = FALSE)

wRF.dist(tree1, tree2, normalize = FALSE, check.labels = TRUE,
         rooted = FALSE)

KF.dist(tree1, tree2 , check.labels = TRUE, rooted = FALSE)

path.dist(tree1, tree2, check.labels = TRUE, use.weight = FALSE)


# Comparison 1.2: Transmission networks


full.true.graph <- ga.graph # full transmission network

full.true.graph.adj<-get.adjacency(full.true.graph,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

full.true.graph.properties<-GenInd(full.true.graph.adj)

fit_power_law(full.true.graph)

source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.full.graph <- ConnectNearBy(phylo.tree.full)


built.full.graph.adj<-get.adjacency(built.full.graph,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.full.graph.properties<-GenInd(full.true.graph.adj)

fit_power_law(built.full.graph)

# Comparion 2: pruned true tree & reconstructed tree (with sequences of individuals who form the pruned tree)
# and their respective transmission networks

# Subtree 1

# Prune the true epi.tree
keep1<-c("2","4","6","8","10")

pruned.epi.tree1<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep1, epi.tree$tip.label)])
# Or
# pruned.epi.tree1<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep1))
# write.tree(pruned.epi.tree1)

# Construct a new phylotree with selected individuals
seq.sim1 <- read.FASTA("HIVSeq_fullNetwork1.fasta")
tree.dat1 <- phyDat(seq.sim1, type = "DNA")
tree.ml1 <- dist.ml(tree.dat1)
tree.sim1 <- upgma(tree.ml1)
phylo.subtree1 <- tree.sim1

dist.subtree1 <- treedist(pruned.epi.tree1, phylo.subtree1, check.labels = TRUE)

# Subnetwork 1

# from the true transmission network
true.subgraph1 <- delete_vertices(ga.graph, keep1) # transmission subnetwork

true.subgraph.adj1<-get.adjacency(true.subgraph1,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties1<-GenInd(true.subgraph.adj1)

fit_power_law(true.subgraph1)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph1 <- ConnectNearBy(phylo.subtree1)


built.subraph.adj1<-get.adjacency(built.subgraph1,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties1<-GenInd(built.subraph.adj1)

fit_power_law(built.subgraph1)

# Subtree 2

# Prune the true epi.tree
keep2<-c("2","4","6","8","10")

pruned.epi.tree2<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep2, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim2 <- read.FASTA("HIVSeq_fullNetwork2.fasta")
tree.dat2 <- phyDat(seq.sim2, type = "DNA")
tree.ml2 <- dist.ml(tree.dat2)
tree.sim2 <- upgma(tree.ml2)
phylo.subtree2 <- tree.sim2

dist.subtree2 <- treedist(pruned.epi.tree2, phylo.subtree2, check.labels = TRUE)

# Subnetwork 2

# from the true transmission network
true.subgraph2 <- delete_vertices(ga.graph, keep2) # transmission subnetwork

true.subgraph.adj2<-get.adjacency(true.subgraph2,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties2<-GenInd(true.subgraph.adj2)

fit_power_law(true.subgraph2)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph2 <- ConnectNearBy(phylo.subtree2)


built.subraph.adj2<-get.adjacency(built.subgraph2,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties2<-GenInd(built.subraph.adj2)

fit_power_law(built.subgraph2)

# Subtree 3

# Prune the true epi.tree
keep3<-c("2","4","6","8","10")

pruned.epi.tree3<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep3, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim3 <- read.FASTA("HIVSeq_fullNetwork3.fasta")
tree.dat3 <- phyDat(seq.sim3, type = "DNA")
tree.ml3 <- dist.ml(tree.dat3)
tree.sim3 <- upgma(tree.ml3)
phylo.subtree3 <- tree.sim3

dist.subtree3 <- treedist(pruned.epi.tree3, phylo.subtree3, check.labels = TRUE)

# Subnetwork 3

# from the true transmission network
true.subgraph3 <- delete_vertices(ga.graph, keep3) # transmission subnetwork

true.subgraph.adj3<-get.adjacency(true.subgraph3,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties3<-GenInd(true.subgraph.adj3)

fit_power_law(true.subgraph3)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph3 <- ConnectNearBy(phylo.subtree3)


built.subraph.adj3<-get.adjacency(built.subgraph3,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties3<-GenInd(built.subraph.adj3)

fit_power_law(built.subgraph3)


# Subtree 4

# Prune the true epi.tree
keep4<-c("2","4","6","8","10")

pruned.epi.tree4<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep4, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim4 <- read.FASTA("HIVSeq_fullNetwork4.fasta")
tree.dat4 <- phyDat(seq.sim4, type = "DNA")
tree.ml4 <- dist.ml(tree.dat4)
tree.sim4 <- upgma(tree.ml4)
phylo.subtree4 <- tree.sim4

dist.subtree4 <- treedist(pruned.epi.tree4, phylo.subtree4, check.labels = TRUE)

# Subnetwork 4

# from the true transmission network
true.subgraph4 <- delete_vertices(ga.graph, keep4) # transmission subnetwork

true.subgraph.adj4<-get.adjacency(true.subgraph4,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties4<-GenInd(true.subgraph.adj4)

fit_power_law(true.subgraph4)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph4 <- ConnectNearBy(phylo.subtree4)


built.subraph.adj4<-get.adjacency(built.subgraph4,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties4<-GenInd(built.subraph.adj4)

fit_power_law(built.subgraph4)

# Subtree 5

# Prune the true epi.tree
keep5<-c("2","4","6","8","10")

pruned.epi.tree5<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep5, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim5 <- read.FASTA("HIVSeq_fullNetwork5.fasta")
tree.dat5 <- phyDat(seq.sim5, type = "DNA")
tree.ml5 <- dist.ml(tree.dat5)
tree.sim5 <- upgma(tree.ml5)
phylo.subtree5 <- tree.sim5

dist.subtree5 <- treedist(pruned.epi.tree5, phylo.subtree5, check.labels = TRUE)

# Subnetwork 5

# from the true transmission network
true.subgraph5 <- delete_vertices(ga.graph, keep5) # transmission subnetwork

true.subgraph.adj5<-get.adjacency(true.subgraph5,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties5<-GenInd(true.subgraph.adj5)

fit_power_law(true.subgraph5)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph5 <- ConnectNearBy(phylo.subtree5)


built.subraph.adj5<-get.adjacency(built.subgraph5,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties5<-GenInd(built.subraph.adj5)

fit_power_law(built.subgraph5)


# Subtree 6

# Prune the true epi.tree
keep6<-c("2","4","6","8","10")

pruned.epi.tree6<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep6, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim6 <- read.FASTA("HIVSeq_fullNetwork6.fasta")
tree.dat6 <- phyDat(seq.sim6, type = "DNA")
tree.ml6 <- dist.ml(tree.dat6)
tree.sim6 <- upgma(tree.ml6)
phylo.subtree6 <- tree.sim6

dist.subtree6 <- treedist(pruned.epi.tree6, phylo.subtree6, check.labels = TRUE)


# Subnetwork 6

# from the true transmission network
true.subgraph6 <- delete_vertices(ga.graph, keep6) # transmission subnetwork

true.subgraph.adj6<-get.adjacency(true.subgraph6,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties6<-GenInd(true.subgraph.adj6)

fit_power_law(true.subgraph6)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph6 <- ConnectNearBy(phylo.subtree6)


built.subraph.adj6<-get.adjacency(built.subgraph6,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties6<-GenInd(built.subraph.adj6)

fit_power_law(built.subgraph6)


# Subtree 7

# Prune the true epi.tree
keep7<-c("2","4","6","8","10")

pruned.epi.tree7<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep7, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim7 <- read.FASTA("HIVSeq_fullNetwork7.fasta")
tree.dat7 <- phyDat(seq.sim7, type = "DNA")
tree.ml7 <- dist.ml(tree.dat7)
tree.sim7 <- upgma(tree.ml7)
phylo.subtree7 <- tree.sim7

dist.subtree7 <- treedist(pruned.epi.tree7, phylo.subtree7, check.labels = TRUE)


# Subnetwork 7

# from the true transmission network
true.subgraph7 <- delete_vertices(ga.graph, keep7) # transmission subnetwork

true.subgraph.adj7<-get.adjacency(true.subgraph7,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties7<-GenInd(true.subgraph.adj7)

fit_power_law(true.subgraph7)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph7 <- ConnectNearBy(phylo.subtree7)


built.subraph.adj7<-get.adjacency(built.subgraph7,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties7<-GenInd(built.subraph.adj7)

fit_power_law(built.subgraph7)


# Subtree 8

# Prune the true epi.tree
keep8<-c("2","4","6","8","10")

pruned.epi.tree8<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep8, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim8 <- read.FASTA("HIVSeq_fullNetwork8.fasta")
tree.dat8 <- phyDat(seq.sim8, type = "DNA")
tree.ml8 <- dist.ml(tree.dat8)
tree.sim8 <- upgma(tree.ml8)
phylo.subtree8 <- tree.sim8

dist.subtree8 <- treedist(pruned.epi.tree8, phylo.subtree8, check.labels = TRUE)


# Subnetwork 8

# from the true transmission network
true.subgraph8 <- delete_vertices(ga.graph, keep8) # transmission subnetwork

true.subgraph.adj8<-get.adjacency(true.subgraph8,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties8<-GenInd(true.subgraph.adj8)

fit_power_law(true.subgraph8)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph8 <- ConnectNearBy(phylo.subtree8)


built.subraph.adj8<-get.adjacency(built.subgraph8,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties8<-GenInd(built.subraph.adj8)

fit_power_law(built.subgraph8)



# Subtree 9

# Prune the true epi.tree
keep9<-c("2","4","6","8","10")

pruned.epi.tree9<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep9, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim9 <- read.FASTA("HIVSeq_fullNetwork9.fasta")
tree.dat9 <- phyDat(seq.sim9, type = "DNA")
tree.ml9 <- dist.ml(tree.dat9)
tree.sim9 <- upgma(tree.ml9)
phylo.subtree9 <- tree.sim9

dist.subtree9 <- treedist(pruned.epi.tree9, phylo.subtree9, check.labels = TRUE)


# Subnetwork 9

# from the true transmission network
true.subgraph9 <- delete_vertices(ga.graph, keep9) # transmission subnetwork

true.subgraph.adj9<-get.adjacency(true.subgraph9,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties9<-GenInd(true.subgraph.adj9)

fit_power_law(true.subgraph9)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph9 <- ConnectNearBy(phylo.subtree9)


built.subraph.adj9<-get.adjacency(built.subgraph1,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties9<-GenInd(built.subraph.adj9)

fit_power_law(built.subgraph9)


# Subtree 10

# Prune the true epi.tree
keep10<-c("2","4","6","8","10")

pruned.epi.tree10<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep10, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim10 <- read.FASTA("HIVSeq_fullNetwork10.fasta")
tree.dat10 <- phyDat(seq.sim10, type = "DNA")
tree.ml10 <- dist.ml(tree.dat10)
tree.sim10 <- upgma(tree.ml10)
phylo.subtree10 <- tree.sim10

dist.subtree10 <- treedist(pruned.epi.tree10, phylo.subtree10, check.labels = TRUE)


# Subnetwork 10

# from the true transmission network
true.subgraph10 <- delete_vertices(ga.graph, keep10) # transmission subnetwork

true.subgraph.adj10<-get.adjacency(true.subgraph10,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties10<-GenInd(true.subgraph.adj10)

fit_power_law(true.subgraph10)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph10 <- ConnectNearBy(phylo.subtree10)


built.subraph.adj10<-get.adjacency(built.subgraph10,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties10<-GenInd(built.subraph.adj10)

fit_power_law(built.subgraph10)


# Subtree 11

# Prune the true epi.tree
keep11<-c("2","4","6","8","10")

pruned.epi.tree11<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep11, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim11 <- read.FASTA("HIVSeq_fullNetwork11.fasta")
tree.dat11 <- phyDat(seq.sim11, type = "DNA")
tree.ml11 <- dist.ml(tree.dat11)
tree.sim11 <- upgma(tree.ml11)
phylo.subtree11 <- tree.sim11

dist.subtree11 <- treedist(pruned.epi.tree11, phylo.subtree11, check.labels = TRUE)


# Subnetwork 11

# from the true transmission network
true.subgraph11 <- delete_vertices(ga.graph, keep11) # transmission subnetwork

true.subgraph.adj11<-get.adjacency(true.subgraph11,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties11<-GenInd(true.subgraph.adj11)

fit_power_law(true.subgraph11)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph11 <- ConnectNearBy(phylo.subtree11)


built.subraph.adj11<-get.adjacency(built.subgraph11,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties11<-GenInd(built.subraph.adj11)

fit_power_law(built.subgraph11)


# Subtree 12

# Prune the true epi.tree
keep12<-c("2","4","6","8","10")

pruned.epi.tree12<-drop.tip(epi.tree,epi.tree$tip.label[-match(keep12, epi.tree$tip.label)])

# Construct a new phylotree with selected individuals
seq.sim12 <- read.FASTA("HIVSeq_fullNetwork12.fasta")
tree.dat12 <- phyDat(seq.sim12, type = "DNA")
tree.ml12 <- dist.ml(tree.dat12)
tree.sim12 <- upgma(tree.ml12)
phylo.subtree12 <- tree.sim12

dist.subtree12 <- treedist(pruned.epi.tree12, phylo.subtree12, check.labels = TRUE)


# Subnetwork 12

# from the true transmission network
true.subgraph12 <- delete_vertices(ga.graph, keep12) # transmission subnetwork

true.subgraph.adj12<-get.adjacency(true.subgraph12,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

true.subgraph.properties1<-GenInd(true.subgraph.adj12)

fit_power_law(true.subgraph12)

# reconstructed from the phylogenetic tree
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

built.subgraph12 <- ConnectNearBy(phylo.subtree12)


built.subraph.adj12<-get.adjacency(built.subgraph12,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s

built.subgraph.properties12<-GenInd(built.subraph.adj12)

fit_power_law(built.subgraph12)
