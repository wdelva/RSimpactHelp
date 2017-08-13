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

source("/home/david/RSimpactHelp/R/time.wind.network.R")
source("/home/david/RSimpactHelp/R/time.point.network.R")
tw.net <- time.wind.network(datalist = datalist, duration = c(1,15))
tp.net <- time.point.network(datalist = datalist, time = 30)


plot.igraph(tw.net, edge.arrow.size=0, vertex.size=7,
            vertex.label = V(tw.net)$name,layout = layout_with_kk,
            vertex.label.cex=0.6, vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0)

plot.igraph(tp.net, edge.arrow.size=0, vertex.size=7,
            vertex.label = V(tp.net)$name,layout = layout_with_kk,
            vertex.label.cex=0.6, vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0)

# 1.2. Construct transmission epi object to be handled by epi2tree function
# to build a transmission tree

# With many seed IDs, each has its own transmission network
# let take, seed = 2, other seeds 7 & 13 have big networks


# Transmission network with different sampling/removal times
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.baseline.R")
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")
transm.ls <- transmNetworkBuilder.diff(datalist = datalist,endpoint = 40)

# epi object
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
transnetwork <- transm.ls[[2]] # seed 2 considered
epi.tree <- trans.network2tree(transnetwork = transnetwork)

# Comapring timeMRCA
source("/home/david/RSimpactHelp/R/time.mrca.matrix.R")
time.mrca <- time.mrca.matrix(tree = epi.tree) # for the transmission tree

source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/infection_age.R")
infage <- infection_age(datalist = datalist, endpoint = 40)
# infage[[2]]

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

# # Sequence simulation: DONE, it takes long time
# sim <- sequence.simulation(transtree = epi.tree, seedSeq = hiv.seq.env, alpha = 0.90,
#                            rate.list = rate, base.freq = freq)
# saveAlignment.PhyloSim(sim,file = paste("HIVSeq_fullNetwork.fasta",sep=""), skip.internal = TRUE, paranoid = TRUE)
# saveAlignment.PhyloSim(sim,file = paste("HIVSeq_fullNetworkNode.fasta",sep=""), paranoid = TRUE)


##### Section 2: True trees and transmission networks #####
###########################################################

# 2.1. True transmission network: transnetwork on 1.2

# Visualise and summarise
library(igraph)
library(network)


source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/infection_age.R")
infage <- infection_age(datalist = datalist, endpoint = 40)

graph.net <- infage[[2]]

graph.build <- graph.net

graph.build[,1] <- as.character(graph.net[,1]) # donors
graph.build[,2] <- as.character(graph.net[,2]) # recipients
gag = as.matrix(graph.build)
ga.graph = graph.edgelist(gag[,1:2])
E(ga.graph)$weight <- infage[[2]]$infecage

V(ga.graph)$color <- "red"

# Transmission network from simpact
plot.igraph(ga.graph, edge.arrow.size=0, vertex.size=7,
            vertex.frame.color="black", vertex.label.color="black",
            vertex.label = V(ga.graph)$name,layout = layout_with_kk,
            vertex.label.cex=0.6, vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0,
            main = "Transmission network from Simpact")




# Network properties
test.graph <- ga.graph

source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/properties_network.R")
propeties.net <- properties_network(graph = test.graph)

# 2.2. True tree: epi.tree on 1.2

plot(epi.tree, main = "Transmission tree from Simpact")
axisPhylo() # add timescale

# Tree properties

source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/properties_tree.R")
properties.net <- properties_tree(tree = epi.tree)


##### Section 3: Phylogenetic tree reconstruction from sequence data #####
##########################################################################

seq.sim <- read.FASTA("HIVSeq_fullNetwork.fasta")
tree.dat <- phyDat(seq.sim, type = "DNA")
tree.ml <- dist.ml(tree.dat)
tree.sim <- upgma(tree.ml)

phylo.tree.full <- as.phylo(tree.sim)


# Comparison 1: true tree & reconstructed full tree / transmission networks

# Comparison 1.1: Trees


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

dist.fulltree <- treedist(tree1, tree2, check.labels = F) # remember label issue

# symmetric.difference = RF.dist()
RF.dist(tree1, tree2, normalize = FALSE, check.labels = F,
        rooted = FALSE)


## Individuals keepend per size group

# Transmission Old less or equal than 10 yrs
transNet.yrs.Old <- delete_edges(ga.graph, E(ga.graph)[weight>=10])

transNet.full <- ga.graph

seq.sim.size_full <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_72.fasta")
tree.dat.full <- phyDat(seq.sim.size_full, type = "DNA")
tree.ml.full <- dist.ml(tree.dat.full)
tree.sim.full <- upgma(tree.ml.full)

built.net.full <- ConnectNearBy(phylo.tree = epi.tree)# tree.sim.full)

prop.transNet.full <- properties_network(graph = transNet.full)
prop.built.net.full <- properties_network(graph = built.net.full)

time.mrca <- time.mrca.matrix(tree = tree.sim.full)

phylo.tree <- read.nexus("~/BEAST_COMPONENTS/beast.v.2.4.5/bin/envGenseqPosterDavidfull_consensus.nex")

matrix.time <- time.mrca.matrix(tree = phylo.tree)

# Construct transmission network from phylogeny
source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

# Size 5

seq.sim.size_5 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_5.fasta")
tree.dat5 <- phyDat(seq.sim.size_5, type = "DNA")
tree.ml5 <- dist.ml(tree.dat5)
tree.sim5 <- upgma(tree.ml5)
keep.rand.5 <- c("28", "19", "64", "46", "34")

pruned.epi.tree.size5<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.5))

diff.tree.size5 <- RF.dist(tree.sim5, pruned.epi.tree.size5, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size5 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.5))
built.net.size5 <- ConnectNearBy(phylo.tree = tree.sim5)

prop.pruned.net.size5 <-properties_network(graph = pruned.net.size5)
prop.built.net.size5 <-properties_network(graph = built.net.size5)

prop.tree.sim5 <- properties_tree(tree = tree.sim5)
prop.pruned.epi.tree.size5 <- properties_tree(tree = pruned.epi.tree.size5)

# tMRCA
# time.mrca.size5 <- time.mrca.matrix(tree = tree.sim5)

# Size 10

seq.sim.size_10 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_10.fasta")
tree.dat10 <- phyDat(seq.sim.size_10, type = "DNA")
tree.ml10 <- dist.ml(tree.dat10)
tree.sim10 <- upgma(tree.ml10)
keep.rand.10 <- c("28", "19", "64", "46",
                  "34", "13", "14", "17", "18", "30")

pruned.epi.tree.size10<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.10))

diff.tree.size10 <- RF.dist(tree.sim10, pruned.epi.tree.size10, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size10 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.10))
built.net.size10 <- ConnectNearBy(phylo.tree = tree.sim10)

prop.pruned.net.size10 <-properties_network(graph = pruned.net.size10)
prop.built.net.size10 <-properties_network(graph = built.net.size10)

prop.tree.sim10 <- properties_tree(tree = tree.sim10)
prop.pruned.epi.tree.size10 <- properties_tree(tree = pruned.epi.tree.size10)


# Size 15

seq.sim.size_15 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_15.fasta")
tree.dat15 <- phyDat(seq.sim.size_15, type = "DNA")
tree.ml15 <- dist.ml(tree.dat15)
tree.sim15 <- upgma(tree.ml15)
keep.rand.15 <- c("28", "19", "64", "46",
                  "34", "13", "14", "17", "18", "30",
                  "41", "51", "58", "39", "56")

pruned.epi.tree.size15<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.15))

diff.tree.size15 <- RF.dist(tree.sim15, pruned.epi.tree.size15, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size15 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.15))
built.net.size15 <- ConnectNearBy(phylo.tree = tree.sim15)

prop.pruned.net.size15 <-properties_network(graph = pruned.net.size15)
prop.built.net.size15 <-properties_network(graph = built.net.size15)

prop.tree.sim15 <- properties_tree(tree = tree.sim15)
prop.pruned.epi.tree.size15 <- properties_tree(tree = pruned.epi.tree.size15)

# Size 20

seq.sim.size_20 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_20.fasta")
tree.dat20 <- phyDat(seq.sim.size_20, type = "DNA")
tree.ml20 <- dist.ml(tree.dat20)
tree.sim20 <- upgma(tree.ml20)
keep.rand.20 <- c("28", "19", "64", "46", "34", "13", "14", "17",
                  "18", "30", "41", "51", "58", "39", "56", "60", "62",
                  "68", "71", "54")

pruned.epi.tree.size20<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.20))


diff.tree.size20 <- RF.dist(tree.sim20, pruned.epi.tree.size20, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size20 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.20))
built.net.size20 <- ConnectNearBy(phylo.tree = tree.sim20)

prop.pruned.net.size20 <-properties_network(graph = pruned.net.size20)
prop.built.net.size5 <-properties_network(graph = built.net.size20)

prop.tree.sim20 <- properties_tree(tree = tree.sim20)
prop.pruned.epi.tree.size20 <- properties_tree(tree = pruned.epi.tree.size20)

# Size 25

seq.sim.size_25 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_25.fasta")
tree.dat25 <- phyDat(seq.sim.size_25, type = "DNA")
tree.ml25 <- dist.ml(tree.dat25)
tree.sim25 <- upgma(tree.ml25)
keep.rand.25 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51", "58",
                  "39", "56", "60", "62",
                  "68", "71", "54", "31", "35", "43", "44", "45")

pruned.epi.tree.size25<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.25))

diff.tree.size25 <- RF.dist(tree.sim25, pruned.epi.tree.size25, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size25 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.25))
built.net.size25 <- ConnectNearBy(phylo.tree = tree.sim25)

prop.pruned.net.size25 <-properties_network(graph = pruned.net.size25)
prop.built.net.size25 <-properties_network(graph = built.net.size25)

prop.tree.sim25 <- properties_tree(tree = tree.sim25)
prop.pruned.epi.tree.size25 <- properties_tree(tree = pruned.epi.tree.size25)

# Size 30

seq.sim.size_30 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_30.fasta")
tree.dat30 <- phyDat(seq.sim.size_30, type = "DNA")
tree.ml30 <- dist.ml(tree.dat30)
tree.sim30 <- upgma(tree.ml30)
keep.rand.30 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51", "58",
                  "39", "56", "60", "62", "68", "71", "54",
                  "31", "35", "43", "44", "45", "49", "50",
                  "55", "9","16")

pruned.epi.tree.size30<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.30))

diff.tree.size30 <- RF.dist(tree.sim30, pruned.epi.tree.size30, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size30 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.30))
built.net.size30 <- ConnectNearBy(phylo.tree = tree.sim30)

prop.pruned.net.size30 <-properties_network(graph = pruned.net.size30)
prop.built.net.size30 <-properties_network(graph = built.net.size30)

prop.tree.sim30 <- properties_tree(tree = tree.sim30)
prop.pruned.epi.tree.size30 <- properties_tree(tree = pruned.epi.tree.size30)

# Size 35

seq.sim.size_35 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_35.fasta")
tree.dat35 <- phyDat(seq.sim.size_35, type = "DNA")
tree.ml35 <- dist.ml(tree.dat35)
tree.sim35 <- upgma(tree.ml35)
keep.rand.35 <- c("28", "19", "64", "46", "34",
                  "13", "14", "17", "18", "30",
                  "41", "51", "58", "39", "56",
                  "60", "62", "68", "71", "54",
                  "31", "35", "43", "44", "45",
                  "49", "50", "55", "9",  "16",
                  "26", "48", "3",  "5",  "7" )

pruned.epi.tree.size35<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.35))

diff.tree.size35 <- RF.dist(tree.sim35, pruned.epi.tree.size35, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size35 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.35))
built.net.size35 <- ConnectNearBy(phylo.tree = tree.sim35)

prop.pruned.net.size35 <-properties_network(graph = pruned.net.size35)
prop.built.net.size35 <-properties_network(graph = built.net.size35)

prop.tree.sim35 <- properties_tree(tree = tree.sim35)
prop.pruned.epi.tree.size35 <- properties_tree(tree = pruned.epi.tree.size35)

# Size 40

seq.sim.size_40 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_40.fasta")
tree.dat40 <- phyDat(seq.sim.size_40, type = "DNA")
tree.ml40 <- dist.ml(tree.dat40)
tree.sim40 <- upgma(tree.ml40)
keep.rand.40 <- c("28", "19", "64", "46", "34",
                  "13", "14", "17", "18", "30",
                  "41", "51", "58", "39", "56",
                  "60", "62", "68", "71", "54",
                  "31", "35", "43", "44", "45",
                  "49", "50", "55", "9",  "16",
                  "26", "48", "3",  "5",  "7",
                  "8",  "10", "11", "57", "61")

pruned.epi.tree.size40<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.40))

diff.tree.size40 <- RF.dist(tree.sim40, pruned.epi.tree.size40, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size40 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.40))
built.net.size40 <- ConnectNearBy(phylo.tree = tree.sim40)

prop.pruned.net.size40 <-properties_network(graph = pruned.net.size40)
prop.built.net.size40 <-properties_network(graph = built.net.size40)

prop.tree.sim40 <- properties_tree(tree = tree.sim40)
prop.pruned.epi.tree.size40 <- properties_tree(tree = pruned.epi.tree.size40)

# Size 45

seq.sim.size_45 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_45.fasta")
tree.dat45 <- phyDat(seq.sim.size_45, type = "DNA")
tree.ml45 <- dist.ml(tree.dat45)
tree.sim45 <- upgma(tree.ml45)
keep.rand.45 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51",
                  "58", "39", "56", "60", "62", "68",
                  "71", "54", "31", "35", "43", "44",
                  "45", "49", "50", "55", "9",  "16",
                  "26", "48", "3",  "5",  "7",  "8" ,
                  "10", "11", "57", "61", "59", "67",
                  "38", "70", "69")

pruned.epi.tree.size45<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.45))

diff.tree.size45 <- RF.dist(tree.sim45, pruned.epi.tree.size45, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size45 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.45))
built.net.size45 <- ConnectNearBy(phylo.tree = tree.sim45)

prop.pruned.net.size45 <-properties_network(graph = pruned.net.size45)
prop.built.net.size45 <-properties_network(graph = built.net.size45)

prop.tree.sim45 <- properties_tree(tree = tree.sim45)
prop.pruned.epi.tree.size45 <- properties_tree(tree = pruned.epi.tree.size45)

# Size 50

seq.sim.size_50 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_50.fasta")
tree.dat50 <- phyDat(seq.sim.size_50, type = "DNA")
tree.ml50 <- dist.ml(tree.dat50)
tree.sim50 <- upgma(tree.ml50)
keep.rand.50 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51",
                  "58", "39", "56", "60", "62", "68",
                  "71", "54", "31", "35", "43", "44",
                  "45", "49", "50", "55", "9",  "16",
                  "26", "48", "3",  "5",  "7",  "8",
                  "10", "11", "57", "61", "59", "67",
                  "38", "70", "69", "37", "42", "4",
                  "12", "15")

pruned.epi.tree.size50<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.50))

diff.tree.size50 <- RF.dist(tree.sim50, pruned.epi.tree.size50, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size50 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.50))
built.net.size50 <- ConnectNearBy(phylo.tree = tree.sim50)

prop.pruned.net.size50 <-properties_network(graph = pruned.net.size50)
prop.built.net.size50 <-properties_network(graph = built.net.size50)

prop.tree.sim50 <- properties_tree(tree = tree.sim50)
prop.pruned.epi.tree.size50 <- properties_tree(tree = pruned.epi.tree.size50)

# Size 55
seq.sim.size_55 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_55.fasta")
tree.dat55 <- phyDat(seq.sim.size_55, type = "DNA")
tree.ml55 <- dist.ml(tree.dat55)
tree.sim55 <- upgma(tree.ml55)
keep.rand.55 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51",
                  "58", "39", "56", "60", "62", "68",
                  "71", "54", "31", "35", "43", "44",
                  "45", "49", "50", "55", "9" , "16",
                  "26", "48", "3" , "5" , "7" , "8",
                  "10", "11", "57", "61","59", "67",
                  "38", "70", "69", "37", "42", "4",
                  "12", "15", "29", "40", "52", "24", "23")

pruned.epi.tree.size55<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.55))

diff.tree.size55 <- RF.dist(tree.sim55, pruned.epi.tree.size55, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size55 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.55))
built.net.size55 <- ConnectNearBy(phylo.tree = tree.sim55)

prop.pruned.net.size55 <-properties_network(graph = pruned.net.size55)
prop.built.net.size55 <-properties_network(graph = built.net.size55)

prop.tree.sim55 <- properties_tree(tree = tree.sim55)
prop.pruned.epi.tree.size55 <- properties_tree(tree = pruned.epi.tree.size55)

# Size 60

seq.sim.size_60 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_60.fasta")
tree.dat60 <- phyDat(seq.sim.size_60, type = "DNA")
tree.ml60 <- dist.ml(tree.dat60)
tree.sim60 <- upgma(tree.ml60)
keep.rand.60 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51",
                  "58", "39", "56", "60", "62", "68",
                  "71", "54", "31", "35", "43", "44",
                  "45", "49", "50", "55", "9" , "16",
                  "26", "48", "3" , "5" , "7" , "8",
                  "10", "11", "57", "61", "59", "67",
                  "38", "70", "69", "37", "42", "4",
                  "12", "15", "29", "40", "52", "24",
                  "23", "32", "33", "6" , "2" , "47")

pruned.epi.tree.size60<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.60))

diff.tree.size60 <- RF.dist(tree.sim60, pruned.epi.tree.size60, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size60 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.60))
built.net.size60 <- ConnectNearBy(phylo.tree = tree.sim60)

prop.pruned.net.size60 <-properties_network(graph = pruned.net.size60)
prop.built.net.size60 <-properties_network(graph = built.net.size60)

prop.tree.sim60 <- properties_tree(tree = tree.sim60)
prop.pruned.epi.tree.size60 <- properties_tree(tree = pruned.epi.tree.size60)


# Size 65

seq.sim.size_65 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_65.fasta")
tree.dat65 <- phyDat(seq.sim.size_65, type = "DNA")
tree.ml65 <- dist.ml(tree.dat65)
tree.sim65 <- upgma(tree.ml65)
keep.rand.65 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51",
                  "58", "39", "56", "60", "62", "68",
                  "71", "54", "31", "35", "43", "44",
                  "45", "49", "50", "55", "9",  "16",
                  "26", "48", "3" , "5", "7",  "8",
                  "10", "11", "57", "61", "59", "67",
                  "38", "70", "69", "37", "42", "4",
                  "12", "15", "29","40", "52", "24",
                  "23", "32", "33", "6",  "2",  "47",
                  "1" , "53", "63", "65", "66")

pruned.epi.tree.size65<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.65))


diff.tree.size65 <- RF.dist(tree.sim65, pruned.epi.tree.size65, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size65 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.65))
built.net.size65 <- ConnectNearBy(phylo.tree = tree.sim65)

prop.pruned.net.size65 <-properties_network(graph = pruned.net.size65)
prop.built.net.size65 <-properties_network(graph = built.net.size65)

prop.tree.sim65 <- properties_tree(tree = tree.sim65)
prop.pruned.epi.tree.size65 <- properties_tree(tree = pruned.epi.tree.size65)

# Size 70

seq.sim.size_70 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_70.fasta")
tree.dat70 <- phyDat(seq.sim.size_70, type = "DNA")
tree.ml70 <- dist.ml(tree.dat70)
tree.sim70 <- upgma(tree.ml70)
keep.rand.70 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51",
                  "58", "39", "56", "60", "62", "68",
                  "71", "54", "31", "35", "43", "44",
                  "45", "49", "50", "55", "9" , "16",
                  "26", "48", "3" , "5", "7",  "8",
                  "10", "11", "57", "61", "59", "67",
                  "38", "70", "69", "37", "42", "4",
                  "12", "15", "29", "40", "52", "24",
                  "23", "32", "33", "6" , "2",  "47",
                  "1" , "53", "63", "65", "66", "20",
                  "21", "22", "25", "27")

pruned.epi.tree.size70<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.70))

diff.tree.size70 <- RF.dist(tree.sim70, pruned.epi.tree.size70, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size70 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.70))
built.net.size70 <- ConnectNearBy(phylo.tree = tree.sim70)

prop.pruned.net.size70 <-properties_network(graph = pruned.net.size70)
prop.built.net.size70 <-properties_network(graph = built.net.size70)

prop.tree.sim70 <- properties_tree(tree = tree.sim70)
prop.pruned.epi.tree.size70 <- properties_tree(tree = pruned.epi.tree.size70)


# Size 72

seq.sim.size_72 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_72.fasta")
tree.dat72 <- phyDat(seq.sim.size_72, type = "DNA")
tree.ml72 <- dist.ml(tree.dat72)
tree.sim72 <- upgma(tree.ml72)
keep.rand.72 <- c("28", "19", "64", "46", "34", "13",
                  "14", "17", "18", "30", "41", "51",
                  "58", "39", "56", "60", "62", "68",
                  "71", "54", "31", "35", "43", "44",
                  "45", "49", "50", "55", "9" , "16",
                  "26", "48", "3" , "5" ,"7"  ,"8",
                  "10", "11", "57", "61", "59", "67",
                  "38", "70", "69", "37", "42", "4",
                  "12", "15", "29", "40", "52", "24",
                  "23", "32", "33", "6" , "2" , "47",
                  "1" , "53", "63", "65", "66", "20",
                  "21", "22", "25", "27", "36", "0" )


pruned.epi.tree.size72<-drop.tip(epi.tree, setdiff(epi.tree$tip.label, keep.rand.72))

diff.tree.size72 <- RF.dist(tree.sim72, pruned.epi.tree.size72, normalize = FALSE, check.labels = F,
                           rooted = T)

pruned.net.size72 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),keep.rand.72))
built.net.size72 <- ConnectNearBy(phylo.tree = tree.sim72)

prop.pruned.net.size72 <-properties_network(graph = pruned.net.size72)
prop.built.net.size72 <-properties_network(graph = built.net.size72)

prop.tree.sim72 <- properties_tree(tree = tree.sim72)
prop.pruned.epi.tree.size72 <- properties_tree(tree = pruned.epi.tree.size72)

#
#         XXXXXXXXXXXXXXXXXXXXXXXXXX 60 percent XXXXXXXXXXXXXXXXXXXX
#

# random, without clusters (elements with high degree) & bridges (elements which bridge two or more clusters)



seq.sim.random_60_1 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_1.fasta")
tree.datrandom_60_1 <- phyDat(seq.sim.random_60_1, type = "DNA")
tree.mlrandom_60_1 <- dist.ml(tree.datrandom_60_1)
tree.simrandom_60_1 <- upgma(tree.mlrandom_60_1)
keep.rand.random_60_1 <- c("9",  "26", "48", "8" , "10",
                           "14", "17", "18", "30", "41",
                           "51", "58", "60", "62", "68",
                           "71", "43", "44", "45", "49",
                           "50", "55", "59", "67", "70" ,
                           "69", "42", "53", "63", "65" ,
                           "46", "20", "21", "22", "25",
                           "27", "36", "12", "15", "29" ,
                           "40", "52", "32", "1" )


seq.sim.random_60_2 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_2.fasta")
tree.datrandom_60_2 <- phyDat(seq.sim.random_60_2, type = "DNA")
tree.mlrandom_60_2 <- dist.ml(tree.datrandom_60_2)
tree.simrandom_60_2 <- upgma(tree.mlrandom_60_2)
keep.rand.random_60_2 <-c("8",  "10", "11", "34", "13",
                          "17", "30", "41", "51", "58",
                          "60", "62", "68", "71", "43",
                          "44", "45", "49", "50", "55",
                          "57", "61", "59", "67", "70",
                          "69", "42", "53", "63", "65",
                          "64", "46", "21", "22", "25",
                          "27", "36", "12", "15", "32",
                          "33", "47", "1" , "0" )


seq.sim.random_60_3 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_3.fasta")
tree.datrandom_60_3 <- phyDat(seq.sim.random_60_3, type = "DNA")
tree.mlrandom_60_3 <- dist.ml(tree.datrandom_60_3)
tree.simrandom_60_3 <- upgma(tree.mlrandom_60_3)
keep.rand.random_60_3 <-c("9",  "16", "26", "48", "5",
                          "7",  "8" , "13", "14", "17",
                          "18", "30", "51", "58", "56",
                          "62", "71", "35", "43", "49",
                          "50", "55", "57", "61", "59",
                          "67", "69", "42", "53", "63",
                          "65", "66", "64", "46", "20",
                          "21", "27", "15", "29", "40",
                          "52", "47", "1" , "0" )


seq.sim.random_60_4 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_4.fasta")
tree.datrandom_60_4 <- phyDat(seq.sim.random_60_4, type = "DNA")
tree.mlrandom_60_4 <- dist.ml(tree.datrandom_60_4)
tree.simrandom_60_4 <- upgma(tree.mlrandom_60_4)
keep.rand.random_60_4 <-c("9",  "26", "48", "5",  "7",
                          "8",  "10", "11", "34", "13",
                          "14", "17", "18", "30", "41",
                          "51", "58", "56", "60", "62",
                          "68", "71", "35", "43", "44",
                          "45", "49", "61", "59", "67",
                          "70", "69", "42", "53", "63",
                          "65", "66", "64", "46", "20",
                          "21", "22", "12", "15")

seq.sim.random_60_5 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_5.fasta")
tree.datrandom_60_5 <- phyDat(seq.sim.random_60_5, type = "DNA")
tree.mlrandom_60_5 <- dist.ml(tree.datrandom_60_5)
tree.simrandom_60_5 <- upgma(tree.mlrandom_60_5)
keep.rand.random_60_5 <-c("14", "17", "18", "30", "41",
                          "51", "58", "56", "60", "62",
                          "68", "71", "35", "43", "44",
                          "45", "49", "50", "55", "57",
                          "61", "59", "67", "70", "69",
                          "42", "53", "63", "65", "66",
                          "64", "46", "20", "21", "22",
                          "25", "27", "36", "12", "15",
                          "29", "40", "52", "32")

seq.sim.random_60_6 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_6.fasta")
tree.datrandom_60_6 <- phyDat(seq.sim.random_60_6, type = "DNA")
tree.mlrandom_60_6 <- dist.ml(tree.datrandom_60_6)
tree.simrandom_60_6 <- upgma(tree.mlrandom_60_6)
keep.rand.random_60_6 <-c("41", "51", "58", "56", "60",
                          "62", "68", "71", "35", "43",
                          "44", "45", "49", "50", "55",
                          "57", "61", "59", "67", "70",
                          "69", "42", "53", "63", "65",
                          "66", "64", "46", "20", "21",
                          "22", "25", "27", "36", "12",
                          "15", "29", "40", "52", "32",
                          "33", "47", "1" , "0" )

seq.sim.random_60_7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_7.fasta")
tree.datrandom_60_7 <- phyDat(seq.sim.random_60_7, type = "DNA")
tree.mlrandom_60_7 <- dist.ml(tree.datrandom_60_7)
tree.simrandom_60_7 <- upgma(tree.mlrandom_60_7)
keep.rand.random_60_7 <-c("9",  "16", "26", "48", "5",
                          "7",  "8" , "10", "11", "34",
                          "13", "14", "17", "18", "30",
                          "41", "51", "58", "56", "60",
                          "62", "68", "71", "35", "43",
                          "44", "45", "49", "50", "55",
                          "57", "61", "59", "67", "70",
                          "69", "42", "53", "63", "65",
                          "66", "64", "46", "20")

seq.sim.random_60_8 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_8.fasta")
tree.datrandom_60_8 <- phyDat(seq.sim.random_60_8, type = "DNA")
tree.mlrandom_60_8 <- dist.ml(tree.datrandom_60_8)
tree.simrandom_60_8 <- upgma(tree.mlrandom_60_8)
keep.rand.random_60_8 <-c("9",  "16", "26", "48", "5",
                          "7",  "8" , "10", "11", "34",
                          "13", "58", "56", "60", "62",
                          "68", "71", "35", "43", "44",
                          "45", "49", "50", "55", "57",
                          "42", "53", "63", "65", "66",
                          "64", "46", "20", "21", "22",
                          "25", "27", "15", "29", "40",
                          "52", "32", "1" , "0" )


seq.sim.random_60_9 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_9.fasta")
tree.datrandom_60_9 <- phyDat(seq.sim.random_60_9, type = "DNA")
tree.mlrandom_60_9 <- dist.ml(tree.datrandom_60_9)
tree.simrandom_60_9 <- upgma(tree.mlrandom_60_9)
keep.rand.random_60_9 <-c("9",  "16", "10", "11", "34",
                          "13", "14", "17", "18", "30",
                          "41", "51", "58", "56", "60",
                          "62", "68", "71", "35", "43",
                          "44", "45", "49", "55", "57",
                          "61", "59", "67", "70", "42",
                          "53", "63", "65", "64", "46",
                          "20", "21", "22", "25", "15",
                          "29", "52", "32", "33")


seq.sim.random_60_10 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_10.fasta")
tree.datrandom_60_10 <- phyDat(seq.sim.random_60_10, type = "DNA")
tree.mlrandom_60_10 <- dist.ml(tree.datrandom_60_10)
tree.simrandom_60_10 <- upgma(tree.mlrandom_60_10)
keep.rand.random_60_10 <-c("22", "25", "27", "36", "12",
                           "9",  "7",  "8",  "14", "17",
                           "18", "30", "16", "56", "68",
                           "71", "35", "43", "45", "49",
                           "50", "55", "57", "61", "59",
                           "42", "53", "65", "15", "29",
                           "40", "66", "26", "48", "5",
                           "41", "33", "47", "1",  "0",
                           "51", "58", "64", "52")

seq.sim.random_60_11 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_11.fasta")
tree.datrandom_60_11 <- phyDat(seq.sim.random_60_11, type = "DNA")
tree.mlrandom_60_11 <- dist.ml(tree.datrandom_60_11)
tree.simrandom_60_11 <- upgma(tree.mlrandom_60_11)
keep.rand.random_60_11 <-c("42", "53", "63", "43", "44",
                           "45", "50", "55", "57", "61",
                           "59", "67", "65", "66", "64",
                           "46", "20", "25", "27", "9",
                           "16", "26", "48", "68", "71",
                           "36", "12", "15", "29", "40",
                           "52", "5",  "34", "33", "47",
                           "1",  "14", "17", "18", "30",
                           "41", "51", "58", "62")

seq.sim.random_60_12 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent/random_60_percent_12.fasta")
tree.datrandom_60_12 <- phyDat(seq.sim.random_60_12, type = "DNA")
tree.mlrandom_60_12 <- dist.ml(tree.datrandom_60_12)
tree.simrandom_60_12 <- upgma(tree.mlrandom_60_12)
keep.rand.random_60_12 <-c("12", "15", "48", "5" , "18",
                           "30", "41", "58", "56", "60",
                           "16", "62", "68", "35", "43",
                           "7",  "8",  "10", "14", "17",
                           "44", "45", "49", "50", "55",
                           "57", "61", "59", "70", "69",
                           "20", "21", "22", "25", "40",
                           "52", "32", "33", "47", "42",
                           "53", "63", "65", "66")



## XXXXXX


# random, without clusters (elements with high degree) & bridges (elements which bridge two or more clusters)



seq.sim.random_60_mix_1 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_1_mix.fasta")
tree.datrandom_60_mix_1 <- phyDat(seq.sim.random_60_mix_1, type = "DNA")
tree.mlrandom_60_mix_1 <- dist.ml(tree.datrandom_60_mix_1)
tree.simrandom_60_mix_1 <- upgma(tree.mlrandom_60_mix_1)
keep.rand.random_60_mix_1 <- c("23", "2",  "38", "28", "3",
                               "6",  "24", "54", "31", "37",
                               "4",  "39", "19", "62", "68",
                               "71", "43", "44", "45", "49",
                               "50", "55", "59", "67", "70",
                               "69", "42", "53", "63", "65",
                               "46", "20", "21", "22", "25",
                               "27", "36", "12", "15", "29",
                               "40", "52", "32", "1" )


seq.sim.random_60_mix_2 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_2_mix.fasta")
tree.datrandom_60_mix_2 <- phyDat(seq.sim.random_60_mix_2, type = "DNA")
tree.mlrandom_60_mix_2 <- dist.ml(tree.datrandom_60_mix_2)
tree.simrandom_60_mix_2 <- upgma(tree.mlrandom_60_mix_2)
keep.rand.random_60_mix_2 <-c("34", "13", "17", "30", "41",
                              "51", "58", "60", "62", "68",
                              "23", "2" , "38", "28", "3",
                              "6" , "24", "54", "31", "37",
                              "4" , "39", "19", "71", "43",
                              "44", "45", "49", "50", "55",
                              "57", "61", "59", "67", "70",
                              "69", "42", "46", "21", "22",
                              "25", "27", "36", "12")


seq.sim.random_60_mix_3 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_3_mix.fasta")
tree.datrandom_60_mix_3 <- phyDat(seq.sim.random_60_mix_3, type = "DNA")
tree.mlrandom_60_mix_3 <- dist.ml(tree.datrandom_60_mix_3)
tree.simrandom_60_mix_3 <- upgma(tree.mlrandom_60_mix_3)
keep.rand.random_60_mix_3 <-c("9",  "16", "26", "48", "5",
                              "7",  "8" , "13", "14", "17",
                              "18", "30", "51", "58", "56",
                              "62", "71", "35", "43", "49",
                              "50", "55", "57", "61", "59",
                              "67", "69", "42", "53", "23",
                              "2" , "38", "28", "3" , "6",
                              "24", "54", "31", "37", "4",
                              "39", "19", "40", "52")


seq.sim.random_60_mix_4 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_4_mix.fasta")
tree.datrandom_60_mix_4 <- phyDat(seq.sim.random_60_mix_4, type = "DNA")
tree.mlrandom_60_mix_4 <- dist.ml(tree.datrandom_60_mix_4)
tree.simrandom_60_mix_4 <- upgma(tree.mlrandom_60_mix_4)
keep.rand.random_60_mix_4 <-c("28", "3",  "6",  "24", "54",
                              "31", "37", "4",  "39", "19",
                              "13", "14", "17", "18", "30",
                              "41", "51", "58", "56", "60",
                              "62", "68", "71", "35", "43",
                              "44", "45", "49", "61", "59",
                              "67", "70", "69", "42", "53",
                              "63", "65", "66", "64", "46",
                              "20", "21", "22", "12")

seq.sim.random_60_mix_5 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_5_mix.fasta")
tree.datrandom_60_mix_5 <- phyDat(seq.sim.random_60_mix_5, type = "DNA")
tree.mlrandom_60_mix_5 <- dist.ml(tree.datrandom_60_mix_5)
tree.simrandom_60_mix_5 <- upgma(tree.mlrandom_60_mix_5)
keep.rand.random_60_mix_5 <-c("43", "44", "45", "49", "50", "55",
                              "57", "61", "59", "67", "70", "69",
                              "42", "53", "63", "65", "66", "64",
                              "46", "20", "21", "22", "23", "2" ,
                              "38", "28", "3",  "6" , "24", "27",
                              "36", "12", "15", "29", "40", "52",
                              "54", "31", "37", "4" , "39", "19",
                              "25", "32")

seq.sim.random_60_mix_6 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_6_mix.fasta")
tree.datrandom_60_mix_6 <- phyDat(seq.sim.random_60_mix_6, type = "DNA")
tree.mlrandom_60_mix_6 <- dist.ml(tree.datrandom_60_mix_6)
tree.simrandom_60_mix_6 <- upgma(tree.mlrandom_60_mix_6)
keep.rand.random_60_mix_6 <-c("23", "2" , "38", "28", "3",
                              "59", "67", "70", "69", "42",
                              "53", "63", "65", "66", "64",
                              "46", "20", "21", "22", "25",
                              "6" , "24", "54", "31", "37",
                              "4" , "39", "19", "50", "55",
                              "57", "61", "27", "36", "12",
                              "15", "29", "40", "52", "32",
                              "33", "47", "1" , "0")

seq.sim.random_60_mix_7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_7_mix.fasta")
tree.datrandom_60_mix_7 <- phyDat(seq.sim.random_60_mix_7, type = "DNA")
tree.mlrandom_60_mix_7 <- dist.ml(tree.datrandom_60_mix_7)
tree.simrandom_60_mix_7 <- upgma(tree.mlrandom_60_mix_7)
keep.rand.random_60_mix_7 <-c("23", "2",  "38", "28", "3",
                              "61", "59", "67", "70", "9",
                              "16", "26", "48", "5" , "7",
                              "8" , "10", "11", "34", "13",
                              "6" , "24", "54", "31", "37",
                              "4" , "39", "19", "50", "55",
                              "57", "14", "17", "18", "30",
                              "69", "42", "53", "63", "65",
                              "66", "64", "46", "20")

seq.sim.random_60_mix_8 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_8_mix.fasta")
tree.datrandom_60_mix_8 <- phyDat(seq.sim.random_60_mix_8, type = "DNA")
tree.mlrandom_60_mix_8 <- dist.ml(tree.datrandom_60_mix_8)
tree.simrandom_60_mix_8 <- upgma(tree.mlrandom_60_mix_8)
keep.rand.random_60_mix_8 <-c("23", "2",  "38", "28", "3",
                              "16", "26", "48", "5" , "7",
                              "8" , "10", "11", "34", "13",
                              "58", "56", "60", "62", "68",
                              "71", "35", "43", "44", "45",
                              "49", "6" , "24", "54", "31",
                              "37", "4" , "39", "19", "50",
                              "55", "57", "15", "29", "40",
                              "52", "32", "1" , "0" )


seq.sim.random_60_mix_9 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_9_mix.fasta")
tree.datrandom_60_mix_9 <- phyDat(seq.sim.random_60_mix_9, type = "DNA")
tree.mlrandom_60_mix_9 <- dist.ml(tree.datrandom_60_mix_9)
tree.simrandom_60_mix_9 <- upgma(tree.mlrandom_60_mix_9)
keep.rand.random_60_mix_9 <-c("23", "2" , "38", "28", "3",
                              "6" , "24", "54", "31", "37",
                              "4" , "39", "19", "16", "10",
                              "11", "34", "13", "14", "17",
                              "18", "30", "41", "51", "58",
                              "56", "60", "62", "68", "71",
                              "35", "43", "44", "45", "49",
                              "55", "57", "61", "59", "67",
                              "70", "42", "53", "63")


seq.sim.random_60_mix_10 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_10_mix.fasta")
tree.datrandom_60_mix_10 <- phyDat(seq.sim.random_60_mix_10, type = "DNA")
tree.mlrandom_60_mix_10 <- dist.ml(tree.datrandom_60_mix_10)
tree.simrandom_60_mix_10 <- upgma(tree.mlrandom_60_mix_10)
keep.rand.random_60_mix_10 <-c("23", "2" , "38", "28", "3",
                               "61", "59", "42", "7",  "8",
                               "14", "17", "5",  "41", "33",
                               "47", "1",  "53", "6",  "24",
                               "54", "31", "37", "4",  "39",
                               "19", "65", "15", "29", "40",
                               "66", "26", "48", "22", "25",
                               "27", "36", "12", "9",  "0",
                               "51", "58", "64", "52")

seq.sim.random_60_mix_11 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_11_mix.fasta")
tree.datrandom_60_mix_11 <- phyDat(seq.sim.random_60_mix_11, type = "DNA")
tree.mlrandom_60_mix_11 <- dist.ml(tree.datrandom_60_mix_11)
tree.simrandom_60_mix_11 <- upgma(tree.mlrandom_60_mix_11)
keep.rand.random_60_mix_11 <-c("68", "71", "36", "12", "15",
                               "29", "40", "52", "5",  "61",
                               "23", "2",  "38", "28", "3",
                               "6",  "24", "54", "31", "37",
                               "4",  "39", "19", "59", "67",
                               "65", "66", "64", "46", "20",
                               "25", "27", "34", "33", "47",
                               "1",  "14", "17", "18", "30",
                               "41", "51", "58", "62")

seq.sim.random_60_mix_12 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/random_60_percent_mix_bridges_clusters/random_60_percent_12_mix.fasta")
tree.datrandom_60_mix_12 <- phyDat(seq.sim.random_60_mix_12, type = "DNA")
tree.mlrandom_60_mix_12 <- dist.ml(tree.datrandom_60_mix_12)
tree.simrandom_60_mix_12 <- upgma(tree.mlrandom_60_mix_12)
keep.rand.random_60_mix_12 <-c("23", "2",  "31", "37", "4",
                               "39", "19", "7",  "8",  "10",
                               "14", "17", "44", "45", "49",
                               "50", "55", "57", "61", "59",
                               "70", "69", "38", "28", "3",
                               "6",  "24", "54", "18", "30",
                               "41", "58", "56", "60", "16",
                               "62", "32", "33", "47", "42",
                               "53", "63", "65", "66")

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
