# Simulations for Portugal Poster presentation

setwd("/home/david/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/")

pacman::p_load(ape,expoTree,data.table)

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
transm.ls <- transmNetworkBuilder.baseline(datalist = datalist,endpoint = 40)

# epi object
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
transnetwork <- transm.ls[[2]] # seed 2 considered
epi.tree <- trans.network2tree(transnetwork = transnetwork)


# 1.3. Simulate consensus sequence

source("/home/david/RSimpactHelp/R/sequence.simulation.R")
hiv.seq.raw <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
## Remove the break line in the string of DNA
clean.hiv.seq <-  gsub("\n", "", hiv.seq.raw)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position

hiv.seq.env <- substr(clean.hiv.seq, 6172,6200) # true c(6172,8742)

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

##### Section 2: Trees and transmission networks comparisons #####
##################################################################

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
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
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
# 2.1. True tree: epi.tree on 1.2

plot(epi.tree)
axisPhylo() # add timescale

