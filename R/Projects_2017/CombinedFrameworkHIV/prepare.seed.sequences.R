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

hiv.seq <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
## Remove the break line in the string of DNA
clean.hiv.seq <-  gsub("\n", "", hiv.seq)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position
# nchar(clean.hiv.seq)
# https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html

# Genes: gag (790, 2086), pol (2085, 5096), env (6225, 8795), tat, rev, vif

# Mutation rates: average of whole genome: 3 × 10−5,

# gag:  (9.3 ± 2.3) × 10−5  (http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002251)
# In the above work, the average rate in Vivo was (4.1 ± 1.7) × 10−3 >> take mean 2.9 × 10−3

# https://retrovirology.biomedcentral.com/articles/10.1186/1742-4690-10-49
# Bwteen host mutation:  10-2 substitutions per site per year (subst ·site-1·year-1),
# Within host mutation: whereas the latter are closer to 10-3 subst ·site-1·year-1

hiv.seq.env <- substr(clean.hiv.seq, 6172,8742)
hiv.seq.gag <- substr(clean.hiv.seq, 790,2086)
hiv.seq.pol <- substr(clean.hiv.seq, 2085,5096)

# Calculate the nucleotides frequencies

seq.env = DNAString(hiv.seq.env) # nulceotides
seq.gag = DNAString(hiv.seq.gag) # nulceotides
seq.pol = DNAString(hiv.seq.pol) # nulceotides
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


sequence.sim.env <- sequence.simulation(transtree = tree0, seedSeq = seq.env, alpha = 0.90,
                                    rate.list = rate.env, base.freq = freq.env)
saveAlignment(sequence.sim.env, file=paste("HIV.Env.gene.fasta", sep = ""),
              skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3492804 0.1695838 0.2399844 0.2411513 -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <tree.sequence.seed.nwk >HIV.Env.gene.fasta")


sequence.sim.gag <- sequence.simulation(transtree = tree0, seedSeq = seq.gag, alpha = 0.90,
                                    rate.list = rate.gag, base.freq = freq.gag)
saveAlignment(sequence.sim.gag, file=paste("HIV.Gag.gene.fasta", sep = ""),
              skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3808790 0.1819584 0.2467232 0.1904395 -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <tree.sequence.seed.nwk >HIV.Gag.gene.fasta")


sequence.sim.pol <- sequence.simulation(transtree = tree0, seedSeq = seq.pol, alpha = 0.90,
                                    rate.list = rate.pol, base.freq = freq.pol)
saveAlignment(sequence.sim.pol, file=paste("HIV.Pol.gene.fasta", sep = ""),
              skip.internal = TRUE, paranoid = TRUE)
system("./seq-gen -mGTR -f 0.3887782 0.1646746 0.2277556 0.2187915  -a 0.9 -g 4 -r 0.00125 0.00125 0.00125 0.00125 0.00125 0.00125 -n1 <tree.sequence.seed.nwk >HIV.Pol.gene.fasta")


