# Script for the inputs to treetime tool:

# sequences in fasta format, sampling dates in csv format, and
# a naive phylogenetic tree built with NJ method in nwk format


# Add dating to tips and work with GTR model in constructing phylogenetic trees

rm(list = ls())
setwd("/home/david/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/TreeTime")

pacman::p_load(ape, expoTree, data.table, phylosim,
               RSimpactHelper, readr, phangorn, Biostrings)

source("/home/david/RSimpactHelp/R/transmNetworkBuilder.baseline.R")

source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")

source("/home/david/RSimpactHelp/R/trans.network2tree.R")

source("/home/david/RSimpactHelp/R/sequence.simulation.R")


# 0. Load data from master model output

# track records of dynamic sexual networks simulated using ABM
# with SImpact
master.datalist <- get(load("master.datalist.RData")) #, .GlobalEnv) #load(file="master.datalist.RData")

# head(master.datalist)

datalist <- master.datalist

# 1. Transmission network
##########################

transm.ls <- transmNetworkBuilder.diff(datalist = datalist,endpoint = 40)

transnetwork <- transm.ls[[2]]

# 2. Transmission tree
#######################

epi.tree <- trans.network2tree(transnetwork = transnetwork)


# 3. Simulate sequences
########################

# Bacause it takes into account base frequencies

# get, clean and subset the gene of interest of HIV-
library(readr)
hiv.seq.raw <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
## Remove the break line in the string of DNA
clean.hiv.seq <-  gsub("\n", "", hiv.seq.raw)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position

hiv.seq.env <- substr(clean.hiv.seq, 1,100) # true c(6172,8742) for env gene,

# Calculate the nucleotides frequencies
#
# library(Biostrings)
seq1 = DNAString(hiv.seq.env) # nulceotides

# Chech frequencies
freq <- letterFrequency(seq1, letters = c("A", "C", "G", "T"))/nchar(hiv.seq.env)

# freq <- c(0.3353293,0.2035928,0.2628077,0.1982701)

# substitution rates
rate <- list("a"=0.2, "b"=0.6, "c"=0.12,"d"=0.001, "e"=0.25, "f"=0.24)

# Sequence simulation and saving
sim <- sequence.simulation(transtree = epi.tree, seedSeq = hiv.seq.env, alpha = 0.90,
                           rate.list = rate, base.freq = freq)
saveAlignment.PhyloSim(sim,file = paste("seq.inputTreeTime.fasta",sep=""), skip.internal = TRUE, paranoid = TRUE)

# 4. Construct phylogenetic tree
#################################

seq.sim <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/TreeTime/seq.inputTreeTime.fasta")
dna.dat <- phyDat(seq.sim, type = "DNA")
dna.dist <- dist.ml(dna.dat, model = "F81", bf = freq) # F81
tree.ini <- NJ(dna.dist)

write.tree(tree.ini, file = "NaiveTreeTime.nwk",tree.names = TRUE)

# read.tree.ini <- read.tree(file = "NaiveTree.nwk")

# 5. Getting sampling dates in same order as the sequences

id.samplingtime <- as.data.frame(cbind(transnetwork$id, transnetwork$dtimes)) # IDs and their samling times in the transmission network

id.sequence <-as.numeric(labels(seq.sim)) # IDs in the sequences file

# reorder the IDs as in the sequence file
# we need this when we call sampling dates to be in the same order as the sequences
samp.time.ord <- vector()
for(i in 1:length(id.sequence)){
  for(j in 1:length(id.sequence)){
    if(id.sequence[i] == id.samplingtime$V1[j]){
      samp.time <- id.samplingtime$V2[j]
    }
  }
  samp.time.ord <- c(samp.time.ord, samp.time)
}
id.samplingtime.ord <- cbind(id.sequence,samp.time.ord) # IDs and sampling time as in the same order as the sequences in fasta file

# Sampling dates: if we simulated from 1990
sampling.dates <- 1990+samp.time.ord # in the order tof the sequences in fasta file if we ruan the model from 1990

write.csv2(sampling.dates,file="samplingdates.csv", append = T)

dates <- read.csv2("samplingdates.csv")

# to get the dates

dates.dt <- dates[,2]


