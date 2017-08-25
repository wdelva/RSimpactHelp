pacman::p_load(ape, expoTree, data.table, phylosim,
               RSimpactHelper, readr, phangorn, Biostrings)

source("/home/david/RSimpactHelp/R/transmNetworkBuilder.baseline.R")

source("/home/david/RSimpactHelp/R/trans.network2tree.R")

source("/home/david/RSimpactHelp/R/sequence.simulation.R")

source("/home/david/RSimpactHelp/R/time.mrca.matrix.R")


# 1. Transmission network
transm.ls <- transmNetworkBuilder.baseline(datalist = datalist,endpoint = 40)

# 2. Transmission tree
epi.tree <- tree0 <- trans.network2tree(transnetwork = transnetwork)


# 3. Simulate sequences: Evolutionary model F81 (not JC69),

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

hiv.seq.env <- substr(clean.hiv.seq, 6172,8742) # true c(6172,8742)

# Calculate the nucleotides frequencies
#
# library(Biostrings)
seq1 = DNAString(hiv.seq.env) # nulceotides

# Chech frequencies
freq <- letterFrequency(seq1, letters = c("A", "C", "G", "T"))/nchar(hiv.seq.env)
# > letterFrequency(gag, letters = c("A", "C", "T", "G"))/nchar(hiv.seq.gag)

freq <- c(0.3353293,0.2035928,0.2628077,0.1982701)


# 4. Construct phylogenetic tree
seq.sim <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/baseline_seq/Baseline_HIVSeq_fullNetwork.fasta")
tree.dat <- phyDat(seq.sim, type = "DNA")
tree.ml <- dist.ml(tree.dat, model = "F81", bf = freq)
tree.sim <- upgma(tree.ml)
