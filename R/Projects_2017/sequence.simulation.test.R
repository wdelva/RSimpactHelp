
pacman::p_load(phylosim, dplyr, EasyABC, RSimpactCyan, RSimpactHelper, expoTree, readr, Biostrings)

# 1. Get a transmission network
transm.ls <- transmNetworkBuilder.baseline(datalist = datalist,endpoint = 40) # different trandmission networks
            # or use transmNetworkBuilder.diff(datalist = datalist,endpoint = 40)

transm.ls1 <- transm.ls[[1]] # choose one, here the first



# 2. Get a transmission tree
tree0 <- trans.network2tree(transnetwork = transm.ls1)

# random tree e.g tree0 <- rtree(n=12, rooted = TRUE)



# 3. Simulate consensus sequences for all individuals in the transmission network


# Get input sequence

library(readr)
hiv.seq <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
## Remove the break line in the string of DNA
clean.hiv.seq <-  gsub("\n", "", hiv.seq)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position

hiv.seq.env <- substr(clean.hiv.seq, 6172,8742)

# Calculate the nucleotides frequencies
library(Biostrings)
seq1 = DNAString(hiv.seq.env) # nulceotides

# Check frequencies
freq <- letterFrequency(seq1, letters = c("A", "C", "G", "T"))/nchar(hiv.seq.env)

hivSeq <- hiv.seq.env

# Simulate the sequence by substitution process GTR only
# with nucleotide substitution rate = (0.8 - 1.7)*10^-3 per site/year >>> 1.25*10^-3 (average)
rate <- list("a"=0.00125, "b"=0.00125, "c"=0.00125,"d"=0.00125, "e"=0.00125, "f"=0.00125)

sequence.sim <- sequence.simulation(transtree = tree0, seedSeq = hivSeq, alpha = 0.90,
                           rate.list = rate, base.freq = freq)



saveAlignment.PhyloSim(sequence.sim, file = paste("file_name",sep="")) # alignment saved with internal nodes
saveAlignment.PhyloSim(sequence.sim, file = paste("file_name",sep=""), skip.internal = TRUE, paranoid = TRUE) # alignment saved without internal nodes

