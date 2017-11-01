
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
file.copy(paste("hiv.seq.A.pol.j.fasta", sep = ""),paste("hiv.seq.A.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.A.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.A.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.A.pol.j.fasta.nwk", to = paste("hiv.seq.A.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.386,0.169,0.231,0.214  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -n1 <hiv.seq.A.pol.j.nwk >A.pool.gene.pol.fasta")


# Pol gene B
###########


# call the seed sequences and rename the file
file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("hiv.seq.B.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.B.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.B.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.B.pol.j.fasta.nwk", to = paste("hiv.seq.B.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.389,0.165,0.228,0.219  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -n1 <hiv.seq.B.pol.j.nwk >B.pool.gene.pol.fasta")



# Pol gene C
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.C.pol.j.fasta", sep = ""),paste("hiv.seq.C.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.C.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.C.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.C.pol.j.fasta.nwk", to = paste("hiv.seq.C.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.393,0.172,0.223,0.212  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -n1 <hiv.seq.C.pol.j.nwk >C.pool.gene.pol.fasta")



# Pol gene D
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.D.pol.j.fasta", sep = ""),paste("hiv.seq.D.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.D.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.D.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.D.pol.j.fasta.nwk", to = paste("hiv.seq.D.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.394,0.162,0.225,0.218  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -n1 <hiv.seq.D.pol.j.nwk >D.pool.gene.pol.fasta")



# Pol gene G
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.G.pol.j.fasta", sep = ""),paste("hiv.seq.G.pol.j.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.G.pol.j.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.G.pol.j.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.G.pol.j.fasta.nwk", to = paste("hiv.seq.G.pol.j.nwk", sep = ""))

system("./seq-gen -mGTR -f 0.390,0.170,0.228,0.213  -a 0.9 -g 4 -r 3.37,14.50,1.44,1.21,14.50,1.00 -n1 <hiv.seq.G.pol.j.nwk >G.pool.gene.pol.fasta")

