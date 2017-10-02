rm(list = ls())
# setwd("/home/david/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/")


setwd("C:/Users/niyukuri/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/") # Windows

pacman::p_load(phylosim, dplyr, EasyABC, RSimpactCyan, RSimpactHelper,
               data.table, expoTree, readr, Biostrings)

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
# let take, seed = 2, other seeds 7 & 13 have big networks


# Transmission network with different sampling/removal times
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.baseline.R")
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")

source("C:/Users/niyukuri/Documents/New folder/RSimpactHelp/R/transmNetworkBuilder.baseline.R")
source("C:/Users/niyukuri/Documents/New folder/RSimpactHelp/R/transmNetworkBuilder.diff.R")
transm.ls.diff <- transmNetworkBuilder.diff(datalist = datalist,endpoint = 40)

# epi object
source("/home/david/RSimpactHelp/R/trans.network2tree.R")

source("C:/Users/niyukuri/Documents/New folder/RSimpactHelp/R/trans.network2tree.R")
transnetwork.diff <- transm.ls.diff[[2]] # seed 2 considered
epi.tree.diff <- trans.network2tree(transnetwork = transnetwork.diff)



transm.ls.base <- transmNetworkBuilder.baseline(datalist = datalist,endpoint = 40)

# epi object
transnetwork.base <- transm.ls.base[[2]] # seed 2 considered
epi.tree.base <- trans.network2tree(transnetwork = transnetwork.base)

# 1.3. Simulate consensus sequence

source("/home/david/RSimpactHelp/R/sequence.simulation.R")

source("C:/Users/niyukuri/Documents/New folder/RSimpactHelp/R/sequence.simulation.R")

hiv.seq <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")

hiv.seq <- read_file("C:/Users/niyukuri/Documents/New folder/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
## Remove the break line in the string of DNA
clean.hiv.seq <-  gsub("\n", "", hiv.seq)

clean.hiv <-  gsub("\r", "", clean.hiv.seq)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position
# nchar(clean.hiv.seq)
# https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html

# Genes: gag (790, 2086), pol (2085, 5096), env (6225, 8795), tat, rev, vif

# Mutation rates: average of whole genome: 3 × 10−5,

# gag:  (9.3 ± 2.3) × 10−5  (http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002251)
# In the above work, the average rate in Vivo was (4.1 ± 1.7) × 10−3

# https://retrovirology.biomedcentral.com/articles/10.1186/1742-4690-10-49
# Bwteen host mutation:  10-2 substitutions per site per year (subst ·site-1·year-1),
# Within host mutation: whereas the latter are closer to 10-3 subst ·site-1·year-1

hiv.seq.genome <- clean.hiv

# Calculate the nucleotides frequencies

seq = DNAString(hiv.seq.genome, start = 1) # nulceotides

# Check frequencies
freq <- letterFrequency(seq, letters = c("A", "C", "G", "T"))/nchar(hiv.seq.genome)


# Simulate the sequence by substitution process GTR only
# with nucleotide substitution rate = (0.8 - 1.7)*10^-3 per site/year >>> 1.25*10^-3 (average)
rate <- list("a"=0.00125, "b"=0.00125, "c"=0.00125,"d"=0.00125, "e"=0.00125, "f"=0.00125)

seq <- as.character(seq)

sequence.sim.diff <- sequence.simulation(transtree = epi.tree.diff, seedSeq = seq, alpha = 0.90,
                                    rate.list = rate, base.freq = freq)


# Save the resulting alignment:
saveAlignment(sequence.sim.diff,file="SeqNetworkdiff.fasta",skip.internal=TRUE)


sequence.sim.base <- sequence.simulation(transtree = epi.tree.base, seedSeq = seq, alpha = 0.90,
                                         rate.list = rate, base.freq = freq)


# Save the resulting alignment:
saveAlignment(sequence.sim.base,file="SeqNetworkbase.fasta",skip.internal=TRUE)
