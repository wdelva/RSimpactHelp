
setwd("~/RSimpactHelp/R/Projects_2017/Report_1/Prepare_Seeds_Sequences/")

## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,
               phyclust, DECIPHER,treedater,geiger,picante)

# From unique sequence (pol gene) in the files below per each strain,
# simulate other 30 sequences under a coalescent tree

# hiv.seq.A.pol.i.fasta
# hiv.seq.B.pol.i.fasta
# hiv.seq.C.pol.i.fasta
# hiv.seq.D.pol.i.fasta
# hiv.seq.G.pol.i.fasta


system('./msa 30 1 -T -t 10.0 > tree.sequence.seed.nwk')

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
file.copy(paste("hiv.seq.A.pol.i.fasta", sep = ""),paste("hiv.seq.A.pol.i.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.A.pol.i.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.A.pol.i.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.A.pol.i.fasta.nwk", to = paste("hiv.seq.A.pol.i.nwk", sep = ""))

#
# sequence.sim.pol <- sequence.simulation(transtree = tree0, seedSeq = seq.pol, alpha = 0.90,
#                                     rate.list = rate.pol, base.freq = freq.pol)
# saveAlignment(sequence.sim.pol, file=paste("HIV.Pol.gene.fasta", sep = ""),
#               skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3887782 0.1646746 0.2277556 0.2187915  -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <hiv.seq.A.pol.i.nwk >A.pool.gene.pol.fasta")


# Pol gene B
###########


# call the seed sequences and rename the file
file.copy(paste("hiv.seq.B.pol.i.fasta", sep = ""),paste("hiv.seq.B.pol.i.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.B.pol.i.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.B.pol.i.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.B.pol.i.fasta.nwk", to = paste("hiv.seq.B.pol.i.nwk", sep = ""))

#
# sequence.sim.pol <- sequence.simulation(transtree = tree0, seedSeq = seq.pol, alpha = 0.90,
#                                     rate.list = rate.pol, base.freq = freq.pol)
# saveAlignment(sequence.sim.pol, file=paste("HIV.Pol.gene.fasta", sep = ""),
#               skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3887782 0.1646746 0.2277556 0.2187915  -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <hiv.seq.B.pol.i.nwk >B.pool.gene.pol.fasta")



# Pol gene C
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.C.pol.i.fasta", sep = ""),paste("hiv.seq.C.pol.i.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.C.pol.i.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.C.pol.i.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.C.pol.i.fasta.nwk", to = paste("hiv.seq.C.pol.i.nwk", sep = ""))

#
# sequence.sim.pol <- sequence.simulation(transtree = tree0, seedSeq = seq.pol, alpha = 0.90,
#                                     rate.list = rate.pol, base.freq = freq.pol)
# saveAlignment(sequence.sim.pol, file=paste("HIV.Pol.gene.fasta", sep = ""),
#               skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3887782 0.1646746 0.2277556 0.2187915  -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <hiv.seq.C.pol.i.nwk >C.pool.gene.pol.fasta")



# Pol gene D
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.D.pol.i.fasta", sep = ""),paste("hiv.seq.D.pol.i.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.D.pol.i.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.D.pol.i.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.D.pol.i.fasta.nwk", to = paste("hiv.seq.D.pol.i.nwk", sep = ""))

#
# sequence.sim.pol <- sequence.simulation(transtree = tree0, seedSeq = seq.pol, alpha = 0.90,
#                                     rate.list = rate.pol, base.freq = freq.pol)
# saveAlignment(sequence.sim.pol, file=paste("HIV.Pol.gene.fasta", sep = ""),
#               skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3887782 0.1646746 0.2277556 0.2187915  -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <hiv.seq.D.pol.i.nwk >D.pool.gene.pol.fasta")



# Pol gene G
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.G.pol.i.fasta", sep = ""),paste("hiv.seq.G.pol.i.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.G.pol.i.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.G.pol.i.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.G.pol.i.fasta.nwk", to = paste("hiv.seq.G.pol.i.nwk", sep = ""))

#
# sequence.sim.pol <- sequence.simulation(transtree = tree0, seedSeq = seq.pol, alpha = 0.90,
#                                     rate.list = rate.pol, base.freq = freq.pol)
# saveAlignment(sequence.sim.pol, file=paste("HIV.Pol.gene.fasta", sep = ""),
#               skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3887782 0.1646746 0.2277556 0.2187915  -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <hiv.seq.G.pol.i.nwk >G.pool.gene.pol.fasta")


