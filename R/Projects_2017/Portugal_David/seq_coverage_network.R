      ##########  URGENCE

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

source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/infection_age.R")
infage <- infection_age(datalist = datalist, endpoint = 40)

graph.net <- infage[[7]]

graph.build <- graph.net

graph.build[,1] <- as.character(graph.net[,1]) # donors
graph.build[,2] <- as.character(graph.net[,2]) # recipients
gag = as.matrix(graph.build)
ga.graph = graph.edgelist(gag[,1:2])
E(ga.graph)$weight <- infage[[2]]$infecage

V(ga.graph)$color <- "red"

transNet.yrs.Old <- delete.vertices(ga.graph, "-1")


########################################
##### Baseline: same sampling time #####     SEED 2
########################################

source("/home/david/RSimpactHelp/R/ConnectNearBy.R")
source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/properties_network.R")
source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/properties_tree.R")
# Size 15

seq.sim.size_15 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/baseline_seq/size_15.fasta")
tree.dat15 <- phyDat(seq.sim.size_15, type = "DNA")
tree.ml15 <- dist.ml(tree.dat15)
tree.sim15 <- upgma(tree.ml15)
keep.rand.15 <- c("57 ", "61 ", "38 ", "44 ", "54 ", "14 ",
                  "32 ", "48 ", "3 ",  "5 " , "58 ", "39 ",
                  "56 ", "60 ", "59 ")



pruned.net.size15 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.15)))
built.net.size15 <- ConnectNearBy(phylo.tree = tree.sim15, epsilon = 0.0001)

prop.pruned.net.size15 <-properties_network(graph = pruned.net.size15)
prop.built.net.size15 <-properties_network(graph = built.net.size15)

# Size 30

seq.sim.size_30 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/baseline_seq/size_30.fasta")
tree.dat30 <- phyDat(seq.sim.size_30, type = "DNA")
tree.ml30 <- dist.ml(tree.dat30)
tree.sim30 <- upgma(tree.ml30)
keep.rand.30 <- c("57 ", "61 ", "38 ", "44 ", "54 ", "14 ",
                  "32 ", "48 ", "3 ",  "5 " , "58 ", "39 ",
                  "56 ", "60 ", "59 ", "67 ", "62 ", "33 ",
                  "26 ", "68 ", "71 ", "43 ", "49 ", "50 ",
                  "55 ", "12 ", "15 ", "31 ", "35 ", "69 ")


pruned.net.size30 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.30)))
built.net.size30 <- ConnectNearBy(phylo.tree = tree.sim30, epsilon = 0.0001)

prop.pruned.net.size30 <-properties_network(graph = pruned.net.size30)
prop.built.net.size30 <-properties_network(graph = built.net.size30)


# Size 45

seq.sim.size_45 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/baseline_seq/size_45.fasta")
tree.dat45 <- phyDat(seq.sim.size_45, type = "DNA")
tree.ml45 <- dist.ml(tree.dat45)
tree.sim45 <- upgma(tree.ml45)
keep.rand.45 <- c("57 ", "61 ", "38 ", "44 ", "54 ", "14 ", "32 ",
                  "48 ", "3 " , "5 ",  "58 ", "39 ", "56 ", "60 ",
                  "59 ", "67 ", "62 ", "33 ", "26 ", "68 ", "71 ",
                  "43 ", "49 ", "50 ", "55 ", "12 ", "15 ", "31 ",
                  "35 ", "69 ", "37 ", "45 ", "70 ", "42 ", "53 ",
                  "63 ", "36 ", "4 ",  "19 ", "64 ", "46 ", "20 ",
                  "65 ", "22 ", "29 ")

pruned.net.size45 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.45)))
built.net.size45 <- ConnectNearBy(phylo.tree = tree.sim45, epsilon = 0.0001)

prop.pruned.net.size45 <-properties_network(graph = pruned.net.size45)
prop.built.net.size45 <-properties_network(graph = built.net.size45)

# Size 60

seq.sim.size_60 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/baseline_seq/size_60.fasta")
tree.dat60 <- phyDat(seq.sim.size_60, type = "DNA")
tree.ml60 <- dist.ml(tree.dat60)
tree.sim60 <- upgma(tree.ml60)
keep.rand.60 <- c("57 ", "61 ", "38 ", "44 ", "54 ", "14 ",
                  "32 ", "48 ", "3 ",  "5 " , "58 ", "39 ",
                  "56 ", "60 ", "59 ", "67 ", "62 ", "33 ",
                  "26 ", "68 ", "71 ", "43 ", "49 ", "50 ",
                  "55 ", "12 ", "15 ", "31 ", "35 ", "69 ",
                  "37 ", "45 ", "70 ", "42 ", "53 ", "63 ",
                  "36 ", "4 " , "19 ", "64 ", "46 ", "20 ",
                  "65 ", "22 ", "29 ", "40 ", "52 ", "0 ",
                  "24 ", "23 ", "8 ",  "10 ", "11 ", "66 ",
                  "28 ", "25 ", "6 " , "2 " , "47 ", "1 " )

pruned.net.size60 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.60)))
built.net.size60 <- ConnectNearBy(phylo.tree = tree.sim60, epsilon = 0.0001)

prop.pruned.net.size60 <-properties_network(graph = pruned.net.size60)
prop.built.net.size60 <-properties_network(graph = built.net.size60)

# Size 72

seq.sim.size_72 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/baseline_seq/size_72.fasta")
tree.dat72 <- phyDat(seq.sim.size_72, type = "DNA")
tree.ml72 <- dist.ml(tree.dat72)
tree.sim72 <- upgma(tree.ml72)
keep.rand.72 <- c("57 ", "61 ", "38 ", "44 ", "54 ", "14 ", "32 ",
                  "48 ", "3 ",  "5 ",  "58 ", "39 ", "56 ", "60 ",
                  "59 ", "67 ", "62 ", "33 ", "26 ", "68 ", "71 ",
                  "43 ", "49 ", "50 ", "55 ", "12 ", "15 ", "31 ",
                  "35 ", "69 ", "37 ", "45 ", "70 ", "42 ", "53 ",
                  "63 ", "36 ", "4 ",  "19 ", "64 ", "46 ", "20 ",
                  "65 ", "22 ", "29 ", "40 ", "52 ", "0 " , "24 ",
                  "23 ", "8 ",  "10 ", "11 ", "66 ", "28 ", "25 ",
                  "6 ",  "2 ",  "47 ", "1 " , "7 ",  "34 ", "13 ",
                  "17 ", "9 ",  "27 ", "16 ", "18 ", "30 ", "41 ", "21 ", "51 ")


pruned.net.size72 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.72)))
built.net.size72 <- ConnectNearBy(phylo.tree = tree.sim72, epsilon = 0.0001)

prop.pruned.net.size72 <-properties_network(graph = pruned.net.size72)
prop.built.net.size72 <-properties_network(graph = built.net.size72)


#########################################
####### Different sampling time #########    SEED 2
#########################################

# Size 15

seq.sim.size_15 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_15.fasta")
tree.dat15 <- phyDat(seq.sim.size_15, type = "DNA")
tree.ml15 <- dist.ml(tree.dat15)
tree.sim15 <- upgma(tree.ml15)
keep.rand.15 <- c("28 ", "19 ", "64 ", "46 ",
                  "34 ", "13 ", "14 ", "17 ", "18 ", "30 ",
                  "41 ", "51 ", "58 ", "39 ", "56 ")


pruned.net.size15 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.15)))
built.net.size15 <- ConnectNearBy(phylo.tree = tree.sim15, epsilon = 0.0001)

prop.pruned.net.size15 <-properties_network(graph = pruned.net.size15)
prop.built.net.size15 <-properties_network(graph = built.net.size15)


# Size 30

seq.sim.size_30 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_30.fasta")
tree.dat30 <- phyDat(seq.sim.size_30, type = "DNA")
tree.ml30 <- dist.ml(tree.dat30)
tree.sim30 <- upgma(tree.ml30)
keep.rand.30 <- c("28 ", "19 ", "64 ", "46 ", "34 ", "13 ",
                  "14 ", "17 ", "18 ", "30 ", "41 ", "51 ", "58 ",
                  "39 ", "56 ", "60 ", "62 ", "68 ", "71 ", "54 ",
                  "31 ", "35 ", "43 ", "44 ", "45 ", "49 ", "50 ",
                  "55 ", "9 ","16 ")

pruned.net.size30 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.30)))
built.net.size30 <- ConnectNearBy(phylo.tree = tree.sim30, epsilon = 0.0001)

prop.pruned.net.size30 <-properties_network(graph = pruned.net.size30)
prop.built.net.size30 <-properties_network(graph = built.net.size30)

# Size 45

seq.sim.size_45 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_45.fasta")
tree.dat45 <- phyDat(seq.sim.size_45, type = "DNA")
tree.ml45 <- dist.ml(tree.dat45)
tree.sim45 <- upgma(tree.ml45)
keep.rand.45 <- c("28 ", "19 ", "64 ", "46 ", "34 ", "13 ",
                  "14 ", "17 ", "18 ", "30 ", "41 ", "51 ",
                  "58 ", "39 ", "56 ", "60 ", "62 ", "68 ",
                  "71 ", "54 ", "31 ", "35 ", "43 ", "44 ",
                  "45 ", "49 ", "50 ", "55 ", "9 ",  "16 ",
                  "26 ", "48 ", "3 ",  "5 ",  "7 ",  "8 " ,
                  "10 ", "11 ", "57 ", "61 ", "59 ", "67 ",
                  "38 ", "70 ", "69 ")

pruned.net.size45 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.45)))
built.net.size45 <- ConnectNearBy(phylo.tree = tree.sim45, epsilon = 0.0001)

prop.pruned.net.size45 <-properties_network(graph = pruned.net.size45)
prop.built.net.size45 <-properties_network(graph = built.net.size45)

# Size 60

seq.sim.size_60 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_60.fasta")
tree.dat60 <- phyDat(seq.sim.size_60, type = "DNA")
tree.ml60 <- dist.ml(tree.dat60)
tree.sim60 <- upgma(tree.ml60)
keep.rand.60 <- c("28 ", "19 ", "64 ", "46 ", "34 ", "13 ",
                  "14 ", "17 ", "18 ", "30 ", "41 ", "51 ",
                  "58 ", "39 ", "56 ", "60 ", "62 ", "68 ",
                  "71 ", "54 ", "31 ", "35 ", "43 ", "44 ",
                  "45 ", "49 ", "50 ", "55 ", "9 " , "16 ",
                  "26 ", "48 ", "3 " , "5 " , "7 " , "8 ",
                  "10 ", "11 ", "57 ", "61 ", "59 ", "67 ",
                  "38 ", "70 ", "69 ", "37 ", "42 ", "4 ",
                  "12 ", "15 ", "29 ", "40 ", "52 ", "24 ",
                  "23 ", "32 ", "33 ", "6 " , "2 " , "47 ")


pruned.net.size60 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.60)))
built.net.size60 <- ConnectNearBy(phylo.tree = tree.sim60, epsilon = 0.0001)

prop.pruned.net.size60 <-properties_network(graph = pruned.net.size60)
prop.built.net.size60 <-properties_network(graph = built.net.size60)


# Size 72

seq.sim.size_72 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Split_DNA_data/seq_sizes_without_considering_specific_seq/size_72.fasta")
tree.dat72 <- phyDat(seq.sim.size_72, type = "DNA")
tree.ml72 <- dist.ml(tree.dat72)
tree.sim72 <- upgma(tree.ml72)
keep.rand.72 <- c("28 ", "19 ", "64 ", "46 ", "34 ", "13 ",
                  "14 ", "17 ", "18 ", "30 ", "41 ", "51 ",
                  "58 ", "39 ", "56 ", "60 ", "62 ", "68 ",
                  "71 ", "54 ", "31 ", "35 ", "43 ", "44 ",
                  "45 ", "49 ", "50 ", "55 ", "9 " , "16 ",
                  "26 ", "48 ", "3 " , "5 " ,"7 "  ,"8 ",
                  "10 ", "11 ", "57 ", "61 ", "59 ", "67 ",
                  "38 ", "70 ", "69 ", "37 ", "42 ", "4 ",
                  "12 ", "15 ", "29 ", "40 ", "52 ", "24 ",
                  "23 ", "32 ", "33 ", "6 " , "2 " , "47 ",
                  "1 ", "53 ", "63 ", "65 ", "66 ", "20 ",
                  "21 ", "22 ", "25 ", "27 ", "36 ", "0 " )


pruned.net.size72 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.72)))
built.net.size72 <- ConnectNearBy(phylo.tree = tree.sim72, epsilon = 0.0001)

prop.pruned.net.size72 <-properties_network(graph = pruned.net.size72)
prop.built.net.size72 <-properties_network(graph = built.net.size72)

### FOR SEED 7

########################################
##### Baseline: same sampling time #####     SEED 7
########################################
library(phangorn)

# Size 15
seq.sim.size_15.seed7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Seq_coverage_base_SEED7/size_15_seed7.fasta")
tree.dat15.seed7 <- phyDat(seq.sim.size_15.seed7, type = "DNA")
tree.ml15.seed7 <- dist.ml(tree.dat15.seed7)
tree.sim15.seed7 <- upgma(tree.ml15.seed7)
keep.rand.15.seed7 <- c("10 ", "12 ", "19 ", "4 ",  "2 ",  "3 ",  "0 ",
                        "56 ", "72 ", "74 ", "71 ", "65 ", "45 ", "47 ", "14 ")
# why no "10 " and "0 "
pruned.net.size15.seed7 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.15.seed7)))
built.net.size15.seed7 <- ConnectNearBy(phylo.tree = tree.sim15.seed7, epsilon = 0.0001)

prop.pruned.net.size15.seed7 <-properties_network(graph = pruned.net.size15.seed7)
prop.built.net.size15.seed7 <-properties_network(graph = built.net.size15.seed7)

# Size 30

seq.sim.size_30.seed7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Seq_coverage_base_SEED7/size_30_seed7.fasta")
tree.dat30.seed7 <- phyDat(seq.sim.size_30.seed7, type = "DNA")
tree.ml30.seed7 <- dist.ml(tree.dat30.seed7)
tree.sim30.seed7 <- upgma(tree.ml30.seed7)
keep.rand.30.seed7 <- c("10 ", "12 ", "19 ", "4 ",  "2 " , "3 ",  "0 ",
                        "56 ", "72 ", "74 ", "71 ", "65 ", "45 ", "47 ",
                        "14 ", "11 ", "17 ", "20 ", "35 ", "15 ", "39 ",
                        "75 ", "73 ", "64 ", "41 ", "18 ", "13 ", "30 ",
                        "70 ", "34 ")

pruned.net.size30.seed7 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.30.seed7)))
built.net.size30.seed7 <- ConnectNearBy(phylo.tree = tree.sim30.seed7, epsilon = 0.0001)

prop.pruned.net.size30.seed7 <-properties_network(graph = pruned.net.size30.seed7)
prop.built.net.size30.seed7 <-properties_network(graph = built.net.size30.seed7)


# Size 45

seq.sim.size_45.seed7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Seq_coverage_base_SEED7/size_45_seed7.fasta")
tree.dat45.seed7 <- phyDat(seq.sim.size_45.seed7, type = "DNA")
tree.ml45.seed7 <- dist.ml(tree.dat45.seed7)
tree.sim45.seed7 <- upgma(tree.ml45.seed7)
keep.rand.45.seed7 <- c("10 ", "12 ", "19 ", "4 ",  "2 ",  "3 ",  "0 ",
                        "56 ", "72 ", "74 ", "71 ", "65 ", "45 ", "47 ",
                        "14 ", "11 ", "17 ", "20 ", "35 ", "15 ", "39 ",
                        "75 ", "73 ", "64 ", "41 ", "18 ", "13 ", "30 ",
                        "70 ", "34 ", "55 ", "58 ", "48 ", "36 ", "7 ",
                        "6 ", "66 ", "44 ", "9 ",  "61 " ,"22 " ,"24 ",
                        "67 ", "60 ", "32 ")

pruned.net.size45.seed7 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.45.seed7)))
built.net.size45.seed7 <- ConnectNearBy(phylo.tree = tree.sim45.seed7, epsilon = 0.0001)

prop.pruned.net.size45.seed7 <-properties_network(graph = pruned.net.size45.seed7)
prop.built.net.size45.seed7 <-properties_network(graph = built.net.size45.seed7)

# Size 60

seq.sim.size_60.seed7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Seq_coverage_base_SEED7/size_60_seed7.fasta")
tree.dat60.seed7 <- phyDat(seq.sim.size_60.seed7, type = "DNA")
tree.ml60.seed7 <- dist.ml(tree.dat60.seed7)
tree.sim60.seed7 <- upgma(tree.ml60.seed7)
keep.rand.60.seed7 <- c("10 ", "12 ", "19 ", "4 ",  "2 ",  "3 ",  "0 ",
                        "56 ", "72 ", "74 ", "71 ", "65 ", "45 ", "47 ",
                        "14 ", "11 ", "17 ", "20 ", "35 ", "15 ", "39 ",
                        "75 ", "73 ", "64 ", "41 ", "18 ", "13 ", "30 ",
                        "70 ", "34 ", "55 ", "58 ", "48 ", "36 ", "7 ",
                        "6 ", "66 ", "44 ", "9 ",  "61 ", "22 ", "24 ",
                        "67 ", "60 ", "32 ", "68 ", "28 ", "43 ", "51 ",
                        "62 ", "59 ", "63 ", "1 ",  "49 ", "52 ", "69 ",
                        "8 ",  "21 ", "42 ", "25 ")


pruned.net.size60.seed7 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.60.seed7)))
built.net.size60.seed7 <- ConnectNearBy(phylo.tree = tree.sim60.seed7, epsilon = 0.0001)

prop.pruned.net.size60.seed7 <-properties_network(graph = pruned.net.size60.seed7)
prop.built.net.size60.seed7 <-properties_network(graph = built.net.size60.seed7)

# Size 75

seq.sim.size_75.seed7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Seq_coverage_base_SEED7/size_75_seed7.fasta")
tree.dat75.seed7 <- phyDat(seq.sim.size_75.seed7, type = "DNA")
tree.ml75.seed7 <- dist.ml(tree.dat75.seed7)
tree.sim75.seed7 <- upgma(tree.ml75.seed7)
keep.rand.75.seed7 <- c( "10 ", "12 ", "19 ", "4 ",  "2 ",  "3 ",  "0 ",
                         "56 ", "72 ", "74 ", "71 ", "65 ", "45 ", "47 ",
                         "14 ", "11 ", "17 ", "20 ", "35 ", "15 ", "39 ",
                         "75 ", "73 ", "64 ", "41 ", "18 ", "13 ", "30 ",
                         "70 ", "34 ", "55 ", "58 ", "48 ", "36 ", "7 ",
                         "6 " , "66 ", "44 ", "9 ",  "61 ", "22 ", "24 ",
                         "67 ", "60 ", "32 ", "68 ", "28 ", "43 ", "51 ",
                         "62 ", "59 ", "63 ", "1 ",  "49 ", "52 ", "69 ",
                         "8 " , "21 ", "42 ", "25 ", "29 ", "38 ", "40 ",
                         "16 ", "23 ", "26 ", "53 ", "5 ",  "27 ", "54 ",
                         "46 ", "50 ", "57 ", "37 ", "31 ")


pruned.net.size75.seed7 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.75.seed7)))
built.net.size75.seed7 <- ConnectNearBy(phylo.tree = tree.sim75.seed7, epsilon = 0.0001)

prop.pruned.net.size75.seed7 <-properties_network(graph = pruned.net.size75.seed7)
prop.built.net.size75.seed7 <-properties_network(graph = built.net.size75.seed7)

# Size 76
seq.sim.size_76.seed7 <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/Seq_coverage_base_SEED7/size_76_seed7.fasta")
tree.dat76.seed7 <- phyDat(seq.sim.size_76.seed7, type = "DNA")
tree.ml76.seed7 <- dist.ml(tree.dat76.seed7)
tree.sim76.seed7 <- upgma(tree.ml76.seed7)
keep.rand.76.seed7 <- c( "10 ", "12 ", "19 ", "4 ",  "2 ",  "3 ",  "0 ",
                         "56 ", "72 ", "74 ", "71 ", "65 ", "45 ", "47 ",
                         "14 ", "11 ", "17 ", "20 ", "35 ", "15 ", "39 ",
                         "75 ", "73 ", "64 ", "41 ", "18 ", "13 ", "30 ",
                         "70 ", "34 ", "55 ", "58 ", "48 ", "36 ", "7 ",
                         "6 " , "66 ", "44 ", "9 ",  "61 ", "22 ", "24 ",
                         "67 ", "60 ", "32 ", "68 ", "28 ", "43 ", "51 ",
                         "62 ", "59 ", "63 ", "1 ",  "49 ", "52 ", "69 ",
                         "8 " , "21 ", "42 ", "25 ", "29 ", "38 ", "40 ",
                         "16 ", "23 ", "26 ", "53 ", "5 ",  "27 ", "54 ",
                         "46 ", "50 ", "57 ", "37 ", "31 ","33 ")


pruned.net.size76.seed7 <- delete_vertices(transNet.yrs.Old, setdiff(V(transNet.yrs.Old),as.numeric(keep.rand.76.seed7)))
built.net.size76.seed7 <- ConnectNearBy(phylo.tree = tree.sim76.seed7, epsilon = 0.0001)

prop.pruned.net.size76.seed7 <-properties_network(graph = pruned.net.size76.seed7)
prop.built.net.size76.seed7 <-properties_network(graph = built.net.size76.seed7)


transNet <- list()
transNet$itimes <- c(35.00000, 34.79618, 31.38415)
transNet$dtimes <- c(19.01449, 10.27123, 15.33656)
transNet$id <- c(0, 1, 2)
transNet$parent <- c(-1,  0,  0)
library(expoTree)
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
a = trans.network2tree(transnetwork = transNet)
plot(a) # error


