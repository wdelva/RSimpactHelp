
setwd("~/RSimpactHelp/R/Projects_2017/Report_1/Prepare_Seeds_Sequences/")

# Scripts to prepare seeds sequences
# A pool of sequences evolved from reference sequence from LOS ALAMOS
# We assume the HIV-1 was introduced by someone for 10 years and we simulate
# the virus evolution and take 30 virions to be the seeds of HIV epidemic in a given population.

## Load required R packages
pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,
               phyclust, DECIPHER)
source("/home/david/RSimpactHelp/R/sequence.simulation.R")
# Take the hole genome of HIV-1

hiv.seq.B <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt") # Subtypes B
## Remove the break line in the string of DNA
clean.hiv.seq.B <-  gsub("\n", "", hiv.seq.B)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position
# nchar(clean.hiv.seq.B)
# https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html

# Genes: gag (790, 2086), pol (2085, 5096), env (6225, 8795), tat, rev, vif

# Mutation rates: average of whole genome: 3 × 10−5,

# gag:  (9.3 ± 2.3) × 10−5  (http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002251)
# In the above work, the average rate in Vivo was (4.1 ± 1.7) × 10−3 >> take mean 2.9 × 10−3

# https://retrovirology.biomedcentral.com/articles/10.1186/1742-4690-10-49
# Bwteen host mutation:  10-2 substitutions per site per year (subst ·site-1·year-1),
# Within host mutation: whereas the latter are closer to 10-3 subst ·site-1·year-1

hiv.seq.B.env <- substr(clean.hiv.seq.B, 6172,8742) # 2571
hiv.seq.B.gag <- substr(clean.hiv.seq.B, 790,2086) # 1297
hiv.seq.B.pol <- substr(clean.hiv.seq.B, 2085,5096) # 3012

# Detect recombinations with REGA, it doesn't take above 1000 ncls
hiv.seq.B.env.part1.800 <- substr(clean.hiv.seq.B, 6172,6971)
hiv.seq.B.env.part2.800 <- substr(clean.hiv.seq.B, 6972,7771)
hiv.seq.B.env.part3.971 <- substr(clean.hiv.seq.B, 7772, 8742)
write(hiv.seq.B.env.part1.800, file = "hiv.seq.B.env.part1.800.fasta")
write(hiv.seq.B.env.part2.800, file = "hiv.seq.B.env.part2.800.fasta")
write(hiv.seq.B.env.part3.971, file = "hiv.seq.B.env.part3.971.fasta")

hiv.seq.B.gag.part1.800 <- substr(clean.hiv.seq.B, 790,1589)
hiv.seq.B.gag.part2.497 <- substr(clean.hiv.seq.B, 1590,2086)
write(hiv.seq.B.gag.part1.800, file = "hiv.seq.B.gag.part1.800.fasta")
write(hiv.seq.B.gag.part2.497, file = "hiv.seq.B.gag.part2.497.fasta")

hiv.seq.B.pol.part1.800 <- substr(clean.hiv.seq.B, 2085,2884)
hiv.seq.B.pol.part2.800 <- substr(clean.hiv.seq.B, 2885,3684)
hiv.seq.B.pol.part3.800 <- substr(clean.hiv.seq.B, 3685,4484)
hiv.seq.B.pol.part4.612 <- substr(clean.hiv.seq.B, 4485,5096)
write(hiv.seq.B.pol.part1.800, file = "hiv.seq.B.pol.part1.800.fasta")
write(hiv.seq.B.pol.part2.800, file = "hiv.seq.B.pol.part2.800.fasta")
write(hiv.seq.B.pol.part3.800, file = "hiv.seq.B.pol.part3.800.fasta")
write(hiv.seq.B.pol.part4.612, file = "hiv.seq.B.pol.part4.612.fasta")


# Genes

write(hiv.seq.B.env, file = "hiv.seq.B.env.i.fasta")
write(hiv.seq.B.gag, file = "hiv.seq.B.gag.i.fasta")
write(hiv.seq.B.pol, file = "hiv.seq.B.pol.i.fasta")

# Calculate the nucleotides frequencies

seq.env = DNAString(hiv.seq.B.env) # nulceotides
seq.gag = DNAString(hiv.seq.B.gag) # nulceotides
seq.pol = DNAString(hiv.seq.B.pol) # nulceotides
# Check frequencies

freq.env <- letterFrequency(seq.env, letters = c("A", "C", "G", "T"))/nchar(seq.env)
# A         C         G         T
# 0.3492804 0.1695838 0.2399844 0.2411513
freq.gag <- letterFrequency(seq.gag, letters = c("A", "C", "G", "T"))/nchar(seq.gag)
# A         C         G         T
# 0.3808790 0.1819584 0.2467232 0.1904395

freq.pol <- letterFrequency(seq.pol, letters = c("A", "C", "G", "T"))/nchar(seq.pol)
# A         C         G         T
# 0.3887782 0.1646746 0.2277556 0.2187915

system('./msa 30 1 -T -t 10.0 > tree.sequence.seed.nwk')

tree0 <- read.tree("tree.sequence.seed.nwk")

# Simulate the sequence by substitution process GTR only
# with nucleotide substitution rate = (0.8 - 1.7)*10^-3 per site/year >>> 1.25*10^-3 (average)

# But it is gene and dynamics (between-host or within-host) specific

rate.env <- list("a"=0.00125, "b"=0.00125, "c"=0.00125,"d"=0.00125, "e"=0.00125, "f"=0.00125)
rate.gag <- list("a"=0.00125, "b"=0.00125, "c"=0.00125,"d"=0.00125, "e"=0.00125, "f"=0.00125)
rate.pol <- list("a"=0.00125, "b"=0.00125, "c"=0.00125,"d"=0.00125, "e"=0.00125, "f"=0.00125)

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


# Env gene B
##########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.B.env.i.fasta", sep = ""),paste("hiv.seq.B.env.i.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.B.env.i.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.B.env.i.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.B.env.i.fasta.nwk", to = paste("hiv.seq.B.env.i.nwk", sep = ""))


# sequence.sim.env <- sequence.simulation(transtree = tree0, seedSeq = seq.env, alpha = 0.90,
#                                     rate.list = rate.env, base.freq = freq.env)
# saveAlignment(sequence.sim.env, file=paste("HIV.Env.gene.fasta", sep = ""),
#               skip.internal = TRUE, paranoid = TRUE)

system("./seq-gen -mGTR -f 0.3492804 0.1695838 0.2399844 0.2411513 -a 0.9 -g 4 -r 0.10125,0.00125,0.00125,0.60125,0.00125,0.00125  -n1 <hiv.seq.B.env.i.nwk >HIV.Env.gene.fasta")



# Gag gene B
###########

# call the seed sequences and rename the file
file.copy(paste("hiv.seq.B.gag.i.fasta", sep = ""),paste("hiv.seq.B.gag.i.fasta.nwk", sep = ""))
# add the number of tree in the file and
write(1,file = "hiv.seq.B.gag.i.fasta.nwk", append = TRUE) # 1 is the number of tree across which we simulate the sequences
# the tree, to prepare the file to simulate the evolution of the virus across the tree
write.tree(tree0,file = "hiv.seq.B.gag.i.fasta.nwk", append = TRUE)
file.rename(from = "hiv.seq.B.gag.i.fasta.nwk", to = paste("hiv.seq.B.gag.i.nwk", sep = ""))

#
# sequence.sim.gag <- sequence.simulation(transtree = tree0, seedSeq = seq.gag, alpha = 0.90,
#                                     rate.list = rate.gag, base.freq = freq.gag)
# saveAlignment(sequence.sim.gag, file=paste("HIV.Gag.gene.fasta", sep = ""),
#               skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3808790 0.1819584 0.2467232 0.1904395 -a 0.9 -g 4 -r 0.00125,0.00125,0.00125,0.00125,0.00125,0.00125 -n1 <hiv.seq.B.gag.i.nwk >HIV.Gag.gene.fasta")

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


### Treat indels

setwd("~/RSimpactHelp/R/Projects_2017/Report_1/Sequences/Untitled Folder/")

dna.A <- read.FASTA("hiv.seq.A.pol.i.fasta")
dna.B <- read.FASTA("hiv.seq.B.pol.i.fasta")
dna.C <- read.FASTA("hiv.seq.C.pol.i.fasta")
dna.D <- read.FASTA("hiv.seq.D.pol.i.fasta")
dna.G <- read.FASTA("hiv.seq.G.pol.i.fasta")

A <- del.gaps(dna.A)
B <- del.gaps(dna.B)
C <- del.gaps(dna.C)
D <- del.gaps(dna.D)
G <- del.gaps(dna.G)

write.dna(A, file = "hiv.seq.A.pol.j.fasta",nbcol = 1, colsep = " ", colw = 4000)
write.dna(B, file = "hiv.seq.B.pol.j.fasta",nbcol = 1, colsep = " ", colw = 4000)
write.dna(C, file = "hiv.seq.C.pol.j.fasta",nbcol = 1, colsep = " ", colw = 4000)
write.dna(D, file = "hiv.seq.D.pol.j.fasta",nbcol = 1, colsep = " ", colw = 4000)
write.dna(G, file = "hiv.seq.G.pol.j.fasta",nbcol = 1, colsep = " ", colw = 4000)

#
# > dna.A <- read.FASTA("hiv.seq.A.pol.i.fasta")
# > dna.B <- read.FASTA("hiv.seq.B.pol.i.fasta")
# > dna.C <- read.FASTA("hiv.seq.C.pol.i.fasta")
# > dna.D <- read.FASTA("hiv.seq.D.pol.i.fasta")
# > dna.G <- read.FASTA("hiv.seq.G.pol.i.fasta")
# > dna.A
# 1 DNA sequence in binary format stored in a list.
#
# Sequence length: 3012
#
# Label:
#
#
#   Base composition:
#   a     c     g     t
# 0.386 0.169 0.231 0.214
# > dna.B
# 1 DNA sequence in binary format stored in a list.
#
# Sequence length: 3012
#
# Label:
#
#
#   Base composition:
#   a     c     g     t
# 0.389 0.165 0.228 0.219
# > dna.C
# 1 DNA sequence in binary format stored in a list.
#
# Sequence length: 3012
#
# Label:
#
#
#   Base composition:
#   a     c     g     t
# 0.393 0.172 0.223 0.212
# > dna.D
# 1 DNA sequence in binary format stored in a list.
#
# Sequence length: 3012
#
# Label:
#
#
#   Base composition:
#   a     c     g     t
# 0.394 0.162 0.225 0.218
# > dna.G
# 1 DNA sequence in binary format stored in a list.
#
# Sequence length: 3012
#
# Label:
#
#
#   Base composition:
#   a     c     g     t
# 0.390 0.170 0.228 0.213
# >

# Monday 6 Nov. 2017 > Subtyping

# Subtype A
# read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
clean.hiv.seq.A <- read_file("hiv.seq.A.pol.j.fasta")
clean.hiv.seq.B <- read_file("hiv.seq.B.pol.j.fasta")
clean.hiv.seq.C <- read_file("hiv.seq.C.pol.j.fasta")
clean.hiv.seq.D <- read_file("hiv.seq.D.pol.j.fasta")
clean.hiv.seq.G <- read_file("hiv.seq.G.pol.j.fasta")


hiv.seq.A.pol.part1.800 <- substr(clean.hiv.seq.A, 1,800)
hiv.seq.A.pol.part2.800 <- substr(clean.hiv.seq.A, 801,1600)
hiv.seq.A.pol.part3.800 <- substr(clean.hiv.seq.A, 1601,2400)
hiv.seq.A.pol.part4.612 <- substr(clean.hiv.seq.A, 2401,3012)
write(hiv.seq.A.pol.part1.800, file = "hiv.seq.A.pol.part1.800.fasta")
write(hiv.seq.A.pol.part2.800, file = "hiv.seq.A.pol.part2.800.fasta")
write(hiv.seq.A.pol.part3.800, file = "hiv.seq.A.pol.part3.800.fasta")
write(hiv.seq.A.pol.part4.612, file = "hiv.seq.A.pol.part4.612.fasta")

# Subtype B

hiv.seq.B.pol.part1.800 <- substr(clean.hiv.seq.B, 1,800)
hiv.seq.B.pol.part2.800 <- substr(clean.hiv.seq.B, 801,1600)
hiv.seq.B.pol.part3.800 <- substr(clean.hiv.seq.B, 1601,2400)
hiv.seq.B.pol.part4.612 <- substr(clean.hiv.seq.B, 2401,3012)
write(hiv.seq.B.pol.part1.800, file = "hiv.seq.B.pol.part1.800.fasta")
write(hiv.seq.B.pol.part2.800, file = "hiv.seq.B.pol.part2.800.fasta")
write(hiv.seq.B.pol.part3.800, file = "hiv.seq.B.pol.part3.800.fasta")
write(hiv.seq.B.pol.part4.612, file = "hiv.seq.B.pol.part4.612.fasta")

# Subtype C

hiv.seq.C.pol.part1.800 <- substr(clean.hiv.seq.C, 1,800)
hiv.seq.C.pol.part2.800 <- substr(clean.hiv.seq.C, 801,1600)
hiv.seq.C.pol.part3.800 <- substr(clean.hiv.seq.C, 1601,2400)
hiv.seq.C.pol.part4.612 <- substr(clean.hiv.seq.C, 2401,3012)
write(hiv.seq.C.pol.part1.800, file = "hiv.seq.C.pol.part1.800.fasta")
write(hiv.seq.C.pol.part2.800, file = "hiv.seq.C.pol.part2.800.fasta")
write(hiv.seq.C.pol.part3.800, file = "hiv.seq.C.pol.part3.800.fasta")
write(hiv.seq.C.pol.part4.612, file = "hiv.seq.C.pol.part4.612.fasta")

# Subtype D

hiv.seq.D.pol.part1.800 <- substr(clean.hiv.seq.D, 1,800)
hiv.seq.D.pol.part2.800 <- substr(clean.hiv.seq.D, 801,1600)
hiv.seq.D.pol.part3.800 <- substr(clean.hiv.seq.D, 1601,2400)
hiv.seq.D.pol.part4.612 <- substr(clean.hiv.seq.D, 2401,3012)
write(hiv.seq.D.pol.part1.800, file = "hiv.seq.D.pol.part1.800.fasta")
write(hiv.seq.D.pol.part2.800, file = "hiv.seq.D.pol.part2.800.fasta")
write(hiv.seq.D.pol.part3.800, file = "hiv.seq.D.pol.part3.800.fasta")
write(hiv.seq.D.pol.part4.612, file = "hiv.seq.D.pol.part4.612.fasta")

# Subtype G

hiv.seq.G.pol.part1.800 <- substr(clean.hiv.seq.G, 1,800)
hiv.seq.G.pol.part2.800 <- substr(clean.hiv.seq.G, 801,1600)
hiv.seq.G.pol.part3.800 <- substr(clean.hiv.seq.G, 1601,2400)
hiv.seq.G.pol.part4.612 <- substr(clean.hiv.seq.G, 2401,3012)
write(hiv.seq.G.pol.part1.800, file = "hiv.seq.G.pol.part1.800.fasta")
write(hiv.seq.G.pol.part2.800, file = "hiv.seq.G.pol.part2.800.fasta")
write(hiv.seq.G.pol.part3.800, file = "hiv.seq.G.pol.part3.800.fasta")
write(hiv.seq.G.pol.part4.612, file = "hiv.seq.G.pol.part4.612.fasta")

# All splitted seq will be subtyped
