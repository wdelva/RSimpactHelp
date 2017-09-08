# Add dating to tips and work with GTR model in constructing phylogenetic trees

rm(list = ls())
setwd("/home/david/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/")

pacman::p_load(ape, expoTree, data.table, phylosim,
               RSimpactHelper, readr, phangorn, Biostrings)

source("/home/david/RSimpactHelp/R/transmNetworkBuilder.baseline.R")

source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")

source("/home/david/RSimpactHelp/R/trans.network2tree.R")

source("/home/david/RSimpactHelp/R/sequence.simulation.R")

source("/home/david/RSimpactHelp/R/time.mrca.matrix.R")


# 0. Load data from master model output

# track records of dynamic sexual networks simulated using ABM
# with SImpact
master.datalist <- get(load("master.datalist.RData")) #, .GlobalEnv) #load(file="master.datalist.RData")

# head(master.datalist)

datalist <- master.datalist

# 1. Transmission network
transm.ls <- transmNetworkBuilder.diff(datalist = datalist,endpoint = 40)

transnetwork <- transm.ls[[2]]

# 2. Transmission tree
epi.tree <- trans.network2tree(transnetwork = transnetwork)


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


# 4. Construct phylogenetic tree with any model (specificaaly to me GTR)

seq.sim <- read.FASTA("~/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/baseline_seq/Baseline_HIVSeq_fullNetwork.fasta")
tree.dat <- phyDat(seq.sim, type = "DNA")
tree.ml <- dist.ml(tree.dat, model = "F81", bf = freq) # F81
tree.sim <- upgma(tree.ml)


fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"), format = "phylip")

dm <- dist.ml(primates)
treeUPGMA <- upgma(dm)
treeNJ <- NJ(dm)

# Parsimony
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)

# Maximum likelihood
fit = pml(treeNJ, data=primates)

methods(class="pml")

fitJC <- optim.pml(fit, TRUE)
logLik(fitJC)

fitGTR <- update(fit, k=4, inv=0.2)

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 0))

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

# Model selection
anova(fitJC, fitGTR)

SH.test(fitGTR, fitJC)

AIC(fitJC)


mt = modelTest(primates)


# Estimate the dates of a rooted phylogenetic tree from the tip dates.

estimate.mu(t, node.dates, p.tol = 0.05)
estimate.dates(t, node.dates, mu = estimate.mu(t, node.dates),
               min.date = -.Machine$double.xmax, show.steps = 0,
               opt.tol = 1e-8, nsteps = 1000, lik.tol = 0,
               is.binary = is.binary.phylo(t))

# node.dates
# a numeric vector of dates for the tips, in the same order as 't$tip.label' or a vector of dates for all of the nodes.

# nodes and their infection & sampling times
id.info <- as.data.frame(cbind(transnetwork$id, (40-transnetwork$itimes),
                               (40-transnetwork$dtimes)))

# tips as appear on the tree
v.num <- as.numeric(tree.sim$tip.label)

# sampling dates in same order of tips as on the tree
v.dat <- as.vector(id.info$V3[v.num])
d = as.phylo(tree.sim)
f = estimate.mu(d, v.dat, p.tol = 0.05)


estimate.dates(d, v.dat, mu = estimate.mu(d, v.dat, p.tol = 0.05, df = 1),
               min.date = -.Machine$double.xmax, show.steps = 0,
               opt.tol = 1e-8, nsteps = 1000, lik.tol = 0,
               is.binary = is.binary.phylo(d))


estimate.mu <- function (t, node.dates, p.tol = 0.05, df = 1)
{
  g <- glm(node.depth.edgelength(t)[1:length(node.dates)] ~
             node.dates, na.action = na.omit)
  null.g <- glm(node.depth.edgelength(t)[1:length(node.dates)] ~
                  1, na.action = na.omit)
  if ((1 - pchisq(AIC(null.g) - AIC(g) + 2, df = 1)) > p.tol) {
    warning(paste("Cannot reject null hypothesis (p=", (1 -
                                                          pchisq(AIC(null.g) - AIC(g))), ")"))
  }
  coef(g)[[2]]
}


#### Renaming automatically sequences adding their sampling times

id.samplingtime <- as.data.frame(cbind(transnetwork$id, transnetwork$dtimes))

id.sequence <-as.numeric(labels(seq.sim))

# reorder
samp.time <- vector()
for(i in 1:length(id.sequence)){
  for(j in 1:length(id.sequence)){
    if(id.sequence[i] == id.samplingtime$V1[j]){
      sampltime <- id.samplingtime$V2[j]
    }
  }
  samp.time <- c(samp.time, sampltime)
}
id.samplingtime.ord <- cbind(id.sequence,samp.time) # IDs and sampling time as in the same order as the sequences in fasta file

samplinng.time <- 1990+samp.time # in the order tof the sequences in fasta file if we ruan the model from 1990

write.csv2(samplinng.time,file="samplingdates.csv", append = F)

dates <- read.csv2("samplingdates.csv")

# to get the dates

dates.dt <- dates[,2]


seq.sim
v <- matrix(nrow = length(id.sequence), ncol = 1)
for(i in 1:length(id.sequence)){
  v[i,1] <- paste(id.samplingtime.ord[,1][i],"_",id.samplingtime.ord[,2][i], sep = "")
}
w <- as.matrix(v)
write.csv2(w,file="samplingtimes.csv", sep = "", append = F)

n <- read.csv2("samplingtimes.csv")

# rename lables

labels(seq.sim) <- n[2]
