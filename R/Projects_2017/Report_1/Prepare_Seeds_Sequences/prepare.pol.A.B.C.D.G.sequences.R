
setwd("~/RSimpactHelp/R/Projects_2017/Report_1/Prepare_Seeds_Sequences/")

## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,
               phyclust, DECIPHER,treedater,geiger,picante)

# From unique sequence (pol gene) in the files below per each strain,
# simulate other 30 sequences under a coalescent tree

# hiv.seq.A.pol.j.fasta
# hiv.seq.B.pol.j.fasta
# hiv.seq.C.pol.j.fasta
# hiv.seq.D.pol.j.fasta
# hiv.seq.G.pol.j.fasta

gnrTime <- 4.0*365/1.2 # infection time for recipients, time interval from when the current donor j
# got the infection from a donor i (or seeding event) until when current donor j transmits to recipient k


system(paste("./msa 30 1 -T -t", gnrTime," > tree.sequence.seed.nwk"))

tree0 <- read.tree("tree.sequence.seed.nwk")


## When using seq-gen to simulate sequence, we need to change a bit the structure of
## input sequence file by adding the tree as follow in the example below:

# Bit of exercise to automate this simulation:
# Put the tree below the sequences and between sequences and trees
# insert the number of trees below.
# In our settings, only one tree is required and the virus sequence is chosen among the sequences above the tree
# The trees can then be placed in this file at the end, after a line stating how many trees there are.
# Example:
# 4 50
# Taxon1 ATCTTTGTAGTCATCGCCGTATTAGCATTCTTAGATCTAA
# Taxon2 ATCCTAGTAGTCGCTTGCGCACTAGCCTTCCGAAATCTAG
# Taxon3 ACTTCTGTGTTTACTGAGCTACTAGCTTCCCTAAATCTAG
# Taxon4 ATTCCTATATTCGCTAATTTCTTAGCTTTCCTGAATCTGG
# 1
# (((Taxon1:0.2,Taxon2:0.2):0.1,Taxon3:0.3):0.1,Taxon4:0.4);

# Pol gene A
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.A.pol.j.fasta", sep = ""),paste("hiv.seq.A.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.A.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.A.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.A.pol.j.fasta.nwk", to = paste("hiv.seq.A.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.3935, 0.1708, 0.2060, 0.2297 -a 0.9 -g 4 -i 0.5230 -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00507 -n 1 -k 1 <hiv.seq.A.pol.j.nwk >A.pool.gene.pol.fasta")
# -s evolutionary rate of pol gene within-host is 5.07 * 10^-1 susb per year

# Pol gene B
###########


# call the seed sequences and rename the file
file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("hiv.seq.B.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.B.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.B.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.B.pol.j.fasta.nwk", to = paste("hiv.seq.B.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -g 4 -i 0.5 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -s 0.00507 -n 1 -k 1 <hiv.seq.B.pol.j.nwk >B.pool.gene.pol.fasta")

#
#
# # Pol gene C
# ###########
#
# # call the seed sequences and rename the file
# file.copy(paste("hiv.seq.C.pol.j.fasta", sep = ""),paste("hiv.seq.C.pol.j.fasta.nwk", sep = ""))
# # add the number of tree in the file and
# write(1,file = "hiv.seq.C.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# # the tree, to prepare the file to simulate the evolution of the virus across the tree
# write.tree(tree0,file = "hiv.seq.C.pol.j.fasta.nwk", append = TRUE)
# file.rename(from = "hiv.seq.C.pol.j.fasta.nwk", to = paste("hiv.seq.C.pol.j.nwk", sep = ""))
#
# system("./seq-gen -mGTR -f 0.393,0.172,0.223,0.212  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -s 0.0045 -n 1 -k 1 <hiv.seq.C.pol.j.nwk >C.pool.gene.pol.fasta")
#
#
#
# # Pol gene D
# ###########
#
# # call the seed sequences and rename the file
# file.copy(paste("hiv.seq.D.pol.j.fasta", sep = ""),paste("hiv.seq.D.pol.j.fasta.nwk", sep = ""))
# # add the number of tree in the file and
# write(1,file = "hiv.seq.D.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# # the tree, to prepare the file to simulate the evolution of the virus across the tree
# write.tree(tree0,file = "hiv.seq.D.pol.j.fasta.nwk", append = TRUE)
# file.rename(from = "hiv.seq.D.pol.j.fasta.nwk", to = paste("hiv.seq.D.pol.j.nwk", sep = ""))
#
# system("./seq-gen -mGTR -f 0.394,0.162,0.225,0.218  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -s 0.0045 -n 1 -k 1 <hiv.seq.D.pol.j.nwk >D.pool.gene.pol.fasta")
#


# Pol gene G
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.G.pol.j.fasta", sep = ""),paste("hiv.seq.G.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.G.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.G.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.G.pol.j.fasta.nwk", to = paste("hiv.seq.G.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.3987, 0.1563, 0.2202, 0.2249 -a 0.9460 -g 4 -i 0.5120 -r 1.4520, 9.9166, 1.3332, 1.2652, 14.9356, 1.0000 -s 0.00507 -n 1 -k 1 <hiv.seq.G.pol.j.nwk >G.pool.gene.pol.fasta")

