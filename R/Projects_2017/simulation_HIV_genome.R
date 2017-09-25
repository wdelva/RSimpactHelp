rm(list = ls())
setwd("/home/david/Dropbox/Niyukuri/Abstract_SACEMA_Research_Days_2017/check/")

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

source("/home/david/RSimpactHelp/R/time.wind.network.R")
source("/home/david/RSimpactHelp/R/time.point.network.R")
tw.net <- time.wind.network(datalist = datalist, duration = c(1,15))
tp.net <- time.point.network(datalist = datalist, time = 30)


plot.igraph(tw.net, edge.arrow.size=0, vertex.size=7,
            vertex.label = V(tw.net)$name,layout = layout_with_kk,
            vertex.label.cex=0.6, vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0)

plot.igraph(tp.net, edge.arrow.size=0, vertex.size=7,
            vertex.label = V(tp.net)$name,layout = layout_with_kk,
            vertex.label.cex=0.6, vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0)

# 1.2. Construct transmission epi object to be handled by epi2tree function
# to build a transmission tree

# With many seed IDs, each has its own transmission network
# let take, seed = 2, other seeds 7 & 13 have big networks


# Transmission network with different sampling/removal times
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.baseline.R")
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")
transm.ls <- transmNetworkBuilder.diff(datalist = datalist,endpoint = 40)

# epi object
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
transnetwork <- transm.ls[[2]] # seed 2 considered
epi.tree <- trans.network2tree(transnetwork = transnetwork)

# Comapring timeMRCA
source("/home/david/RSimpactHelp/R/time.mrca.matrix.R")
time.mrca <- time.mrca.matrix(tree = epi.tree) # for the transmission tree

source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/infection_age.R")
infage <- infection_age(datalist = datalist, endpoint = 40)
# infage[[2]]

# 1.3. Simulate consensus sequence

source("/home/david/RSimpactHelp/R/sequence.simulation.R")


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
# In the above work, the average rate in Vivo was (4.1 ± 1.7) × 10−3

# https://retrovirology.biomedcentral.com/articles/10.1186/1742-4690-10-49
# Bwteen host mutation:  10-2 substitutions per site per year (subst ·site-1·year-1),
# Within host mutation: whereas the latter are closer to 10-3 subst ·site-1·year-1

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




#________________________________
# Simulating partitions

# The following example demonstrates how to use the processes and site-
# and process-specific parameters to simulate "partitions" with different properties.
# We will simulate four partitions:

# . Partition 1: sites in range 1:25 evolving by JC + d Gamma with a shape parameter alpha = 1
# . Partition 2: sites in range 26:50 evolving by JC + d Gamma with a shape parameter alpha = 0.5
# . Partition 3: sites in range 51:75 evolving by HKY + d Gamma with a shape parameter alpha = 1
# . Partition 4: sites in range 76:100 evolving by HKY + d Gamma with a shape parameter alpha = 0.5

# First construct two substitution process objects:
jc69 <- JC69()
hky <- HKY(rate.params = list("Alpha"=5, "Beta"=2),
           base.freqs = c(4,3,2,1)/10)

# Construct a root sequence object of length 100:
root.seq <- NucleotideSequence(length = 100)

# Attach process jc69 to range 1:50:
attachProcess(root.seq, jc69,1:50)

# Attach process hky to range 51:100:
attachProcess(root.seq, hky, 51:100)

# Sample rate multipliers in the four partitions:
plusGamma(root.seq, jc69, 1, 1:25)
plusGamma(root.seq, jc69, 0.5, 26:50)
plusGamma(root.seq, hky, 1, 51:75)
plusGamma(root.seq, hky, 0.5, 76:100)


# Construct the PhyloSim object, sample states, set root sequence,
# set the phylo object (random coalescent tree for three taxa) and run the simulation:
sim <- Simulate(PhyloSim(root.seq = sampleStates(root.seq),
                         phylo = rcoal(3)))

# Plot the alignment alongside the tree, skip sequences at ancestral nodes:
plot(sim, num.pages = 1, plot.ancestors = FALSE)


#______________________________________________________________________
# Simulating many replicates

# Constructing Sequence objects with a large number of sites is expensive,
# so it is a good idea to do that outside the cycle when simulating many
# replicates with the same root sequence length.
# Reusing the root sequence object is easy, but do not forget to do the
# modifications needed to get independent simulations (e.g. clearing the states of the root sequence,
# resampling the rate multipliers).
# The following code illustrates how to simulate many replicates under the JC69 + d Gamma model.


# Construct the root sequence object and attach the substitution process:
p <- JC69()
root.seq <- NucleotideSequence(length = 10)
attachProcess(root.seq,p)

# Read the required phylogeny from file (this will remain fixed in the simulated replicates):
tree <- read.tree("3taxa.nwk")


# Simulate three replicates. Note that the states are cleared and resampled; the rate multipliers are
# resampled as well. The resulting alignments are stored in files aln_1.fas, aln_2.fas, aln_3.fas.
for(i in 1:3){cat(paste("\n\n Simulating replication ", i, "\n\n", sep=""))
  clearStates(root.seq)
  plusGamma(root.seq,p,0.25)
  sampleStates(root.seq)

  sim <- Simulate(PhyloSim(root.seq = root.seq,
                           phylo = tree
  ))

  saveAlignment(sim,file=paste("aln_",i,".fas",sep=""))
}

#--------------
# Simulating many replicates in parallel

# The the speed of the above method for simulating replicates can be improved on a multicore
# machine by running many replicates in parallel by using the mclapply method from the parallel
# package (currently not available on Windows operating systems).

# Under default settings, the mclapply method launches one replication per core and
# this approach needs enough memory to run all of them in parallel.

# The following code illustrates how to simulate many replicates in parallel under the JC69 + d Gamma model.

# Construct the root sequence object and attach the substitution process:
p <- JC69()
root.seq <- NucleotideSequence(length = 50)
attachProcess(root.seq,p)

# Read the required phylogeny from file (this will remain fixed in the simulated replicates):
tree <- read.tree("3taxa.nwk")

# Function to simulate a single replication:

sim.replicate <- function(i){
  name <- paste("replication_",i, sep="")
  clearStates(root.seq)
  plusGamma(root.seq, p, 0.25)
  sampleStates(root.seq)

  sim<-Simulate(PhyloSim(
    name=name,
    root.seq=root.seq,
    phylo=tree
  ),
  quiet = TRUE)
  saveAlignment(sim, file=paste("aln_",i,".fas",sep=""))
  return(sim)
  # return(TRUE)
}


# Note that the states are cleared and resampled; the rate multipliers are resampled as well.
# The resulting alignments are stored in files aln_1.fas, aln_2.fas, alni_3.fas.
# Memory can be saved by throwing away the objects generated by the replication
# by returning TRUE (or any other fixed value) from the sim.replicate function.

# Load the parallel package if available:

have.mcore <- is.element("parallel", installed.packages()[,1])
if(have.mcore){
  library(parallel)
}

# Run replicates in parallel, print the resulting PhyloSimobjects:
if(have.mcore){
  nr.replicates <- 3
  res.objects <- mclapply(1:nr.replicates, sim.replicate)
  print(res.objects)
}

# Print one of the resulting alignments:
if(have.mcore){
  print(res.objects[[1]]$alignment)
}


## Exercises:  Simulating partitions
## See also the package vignette (vignette("PhyloSim",package="phylosim")).
##

# load PhyloSim
library("phylosim")

# The following example demonstrates how to use the processes and site- and process-specific
# parameters to simulate ``partitions'' with different properties.

# We will simulate four partitions:


#	* Partition 1: sites in range \code{1:25} evolving by JC+d$\Gamma$ with a shape parameter
#	  alpha=1
#
#	* Partition 2: sites in range \code{26:50} evolving by JC+d$\Gamma$ with a shape parameter
#	  alpha=0.5
#
#	* Partition 3: sites in range \code{51:75} evolving by HKY+d$\Gamma$ with a shape parameter
#	  alpha=1
#
#	* Partition 4: sites in range \code{76:100} evolving by HKY+d$\Gamma$ with a shape parameter
#	  alpha=0.5

# First construct two substitution process objects:
jc69<-JC69()
hky<-HKY(rate.params=list( "Alpha"=5,"Beta"=2),
         base.freqs=c(4,3,2,1)/10
)

# Construct a root sequence object of length 100:
root.seq<-NucleotideSequence(length=100)

# Attach process jc69 to range 1:50:
attachProcess(root.seq,jc69,1:50)

# Attach process hky to range 51:100:
attachProcess(root.seq,hky,51:100)

# Sample rate multipliers in the four partitions:
plusGamma(root.seq,jc69,1,1:25)
plusGamma(root.seq,jc69,0.5,26:50)

plusGamma(root.seq,hky,1,51:75)
plusGamma(root.seq,hky,0.5,76:100)

# Construct the PhyloSim object, sample states, set root sequence, set the phylo object (random
# coalescent tree for three taxa) and run the simulation:
sim<-Simulate(PhyloSim(
  root.seq=sampleStates(root.seq),
  phylo=rcoal(3)
))

# Plot the alignment alongside the tree, skip sequences at ancestral nodes:
plot(sim,plot.ancestors=FALSE)



##
## Simulating rate variation,
## insertions and deletions.
##

# load the package
library(phylosim)

#	# Construct a phylo object for the
#	# simulation, scale total tree length to 2:
#
tmp<-PhyloSim(phylo=rcoal(3))
scaleTree(tmp,2/tmp$treeLength)
tmp$treeLength
t<-tmp$phylo

# construct a GTR process object
gtr<-GTR(
  name="MyGTR",
  rate.params=list(
    "a"=1, "b"=2, "c"=3,
    "d"=1, "e"=2, "f"=3
  ),
  base.freqs=c(2,2,1,1)/6
)

# get object summary
summary(gtr)

# get a bubble plot
plot(gtr)

# construct root sequence object
s<-NucleotideSequence(length=30)

# attach process via virtual field
s$processes<-list(list(gtr))

# sample states from the equilibrium
# distribution of the attached processes

sampleStates(s)

# create among-sites rate variation by sampling
# the "rate.multiplier" site-process specific parameter
# from a discrete gamma distribution (GTR+G).
plusGamma(s,gtr,shape=0.5)

# make the range 11:20 invariable
setRateMultipliers(s,gtr,0,11:20)

# get the rate multipliers for s and gtr
getRateMultipliers(s,gtr)

# construct deletion process object
# proposing length in range 1:3
d<-DiscreteDeletor(
  rate=1,
  name="MyDel",
  sizes=c(1:3),
  probs=c(3/6,2/6,1/6)
)
# get object
summary(d)

# plot deletion length distribution
plot(d)

# attach d to s
attachProcess(s,d)

# create a region rejecting all deletions
setDeletionTolerance(s,d,0,11:20)

# construct insertion process object
# proposing length in range 1:3
i<-DiscreteInsertor(
  rate=1,
  name="MyDel",
  sizes=c(1:2),
  probs=c(1/2,1/2),
  template.seq=NucleotideSequence(length=1,processes=list(list(JC69())))
)

# states will be sampled from the JC69 equilibrium distribution
# get object
summary(i)

# plot insertion length distribution
plot(i)

# attach i to s
attachProcess(s,i)

# create a region rejecting all insertions
setInsertionTolerance(s,i,0,11:20)

# plot total site rates
plot(s)

# construct simulation object
sim<-PhyloSim(root.seq=s, phylo=t)

# get object summary
summary(sim)

# plot tree
plot(sim)

# run simulation
Simulate(sim)

# get the list of recorded per-branch event counts
getBranchEvents(sim)

# export the number of subtitions as a phylo object
subst<-exportStatTree(sim,"substitution")
subst

# plot the exported phylo object
plot(subst)

# plot tree and alignment
plot(sim)
# save alingment
file<-paste("PhyloSim_fasta_",Sys.getpid(),".fas",sep="");
saveAlignment(sim,file=file,paranoid=TRUE);


##
## Evolving a genomic region containing a small "gene"
#
# The following code demonstrates how to simulate a genomic region
# containing a small "gene" with two exons and the flanking noncoding regions.

# Load the package:
library(phylosim)
# Enable fast mode:
PSIM_FAST<-TRUE

# Construct the alphabet objects:
nuc.alph    <- NucleotideAlphabet();
codon.alph  <- CodonAlphabet()

# Construct the root sequence object:
root.sequence<-Sequence(length=150)

# Define coding and noncoding regions:
coding.pos      <-c(31:50, 91:110)
noncoding.pos   <- (1:150)[ -coding.pos ]

# Set alphabets:
setAlphabets(root.sequence, list(nuc.alph), noncoding.pos)
setAlphabets(root.sequence, list(codon.alph), coding.pos)

# Construct the substitution processes:
k80     <-K80(rate.params=list("Alpha"=2,"Beta"=1), base.freqs=c(2, 1, 2, 1)/4)
gy94    <-GY94(kappa=2, omega.default=0.1, scale.nuc=TRUE)

# Set up indel length distribution:
id.dist<-exp(6:1)/sum(exp(6:1))

# Construct the deletion processes:
del.nc   <- DiscreteDeletor(rate=0.1, sizes=1:6, probs=id.dist )
del.c    <- DiscreteDeletor(rate=0.03, sizes=1:6, probs=id.dist  )

# Construct insertion processes:
ins.nc  <- DiscreteInsertor(rate=0.1, sizes=1:6, probs=id.dist )
ins.c   <- DiscreteInsertor(rate=0.03, sizes=1:6, probs=id.dist )

# Set the template sequences:
ins.nc$templateSeq  <- NucleotideSequence(length=1,processes=list(list( k80, del.nc, ins.nc ) ))
ins.c$templateSeq   <- CodonSequence(length=1,processes=list(list( gy94, del.c, ins.c ) ))

# Attach processes to root sequence:
setProcesses(root.sequence, list(list(k80, del.nc, ins.nc)), noncoding.pos)
setProcesses(root.sequence, list(list(gy94, del.c, ins.c)), coding.pos)

# Fix the stop codon:
start.pos   <- coding.pos[1]
setStates(root.sequence, "ATG", start.pos);                       # Set the state.
setRateMultipliers(root.sequence, gy94, 0, c(start.pos) )         # Make the site invariable.
setDeletionTolerance(root.sequence, del.c, 0, c(start.pos));      # Make the site reject deletions.
setInsertionTolerance(root.sequence,ins.c ,0, c(start.pos) );     # Make the site reject neighboring insertions.

# Construct a substitution process acting on stop codons only:
stop.alphabet   <- Alphabet(symbols=c("TAG", "TAA", "TGA"))
stop.subst      <-GeneralSubstitution(
  alphabet=stop.alphabet,
  rate.list=list("TAG->TAA"=1,
                 "TAG->TGA"=2,
                 "TAA->TAG"=3,
                 "TAA->TGA"=1,
                 "TGA->TAG"=2,
                 "TGA->TAA"=3
  )
)
stop.pos    <- tail(coding.pos, 1)
root.sequence$sites[[stop.pos]]$alphabet    <- stop.alphabet    # Set alphabet for stop codon site.
setProcesses(root.sequence,list(list(stop.subst)), stop.pos)    # Set substitution process for stop codon site.

# Fix splicing sites:
splicing.sites<-c(51, 52, 89, 90);
setStates(root.sequence, c("G", "T", "A", "G"), splicing.sites) # Set site states.
setRateMultipliers(root.sequence, k80, 0, splicing.sites)       # Make sites invariable.
setDeletionTolerance(root.sequence, del.nc , 0, splicing.sites) # Make sites reject deletions.
setInsertionTolerance(root.sequence, ins.nc , 0, splicing.sites)# Make sites reject neighboring insertions.

# Sample site states:
sampleStates(root.sequence)

# Construct simulation object:
sim<-PhyloSim(
  phylo=read.tree("data/mammals.nwk"),
  root.seq=root.sequence
)
# Run the simulation:
Simulate(sim)
# Save the resulting alignment:
saveAlignment(sim,file="example_A1.fas",skip.internal=TRUE)


#
## Simulating many replicates
## See also the package vignette (vignette("PhyloSim",package="phylosim")).
##

# load PhyloSim
library("phylosim")

# Constructing Sequence objects with a large number of sites is expensive, so it is a good idea
# to do that outside the cycle when simulating many replicates with the same root sequence length.

# Reusing the root sequence object is easy, but do not forget to do the modifications needed to get
# independent simulations (e.g. clearing the states of the root sequence, resampling the rate multipliers).

# The following code illustrates how to simulate many replicates under the JC69+dG model.

# Construct the root sequence object and attach the substitution process:
p<-JC69();
root.seq<-NucleotideSequence(length=50)
attachProcess(root.seq,p)

# Read the required phylogeny from file (this will remain fixed in the simulated replicates):
tree<-read.tree("data/3taxa.nwk");

# Simulate three replicates. Note that the states are cleared and resampled; the rate multipliers
# are resampled as well. The resulting aligments are stored in files aln_1.fas, aln_2.fas, aln_3.fas.
for(i in 1:3){
  cat(paste("\n\nSimulating replication ",i,"\n\n",sep=""))

  clearStates(root.seq)
  plusGamma(root.seq,p,0.25)
  sampleStates(root.seq)

  sim<-Simulate(PhyloSim(
    root.seq=root.seq,
    phylo=tree
  ))

  saveAlignment(sim,file=paste("aln_",i,".fas",sep=""))
}


##
## Simulating many replicates in parallel
## See also the package vignette (vignette("PhyloSim",package="phylosim")).
##
# load PhyloSim
library("phylosim")

# The the speed of the above method for simulating replicates can be improved on a multicore machine by running
# many replicates in parallel by using the mclapply method from the multicore package.

# Under default settings, the mclapply method launches one replication per core and this approach needs enough memory to run all
# of them in parallel.

# The following code illustrates how to simulate many replicates in parallel under the JC69+dG model.

# Construct the root sequence object and attach the substitution process:
p<-JC69();
root.seq<-NucleotideSequence(length=50)
attachProcess(root.seq,p)

# Read the required phylogeny from file (this will remain fixed in the simulated replicates):
tree<-read.tree("data/3taxa.nwk");

# Function to simulate a single replication:
sim.replicate<-function(i){
  name<-paste("replication_",i,sep="")
  clearStates(root.seq)
  plusGamma(root.seq,p,0.25)
  sampleStates(root.seq)

  sim<-Simulate(PhyloSim(
    name=name,
    root.seq=root.seq,
    phylo=tree,
  ),
  quiet=TRUE
  )

  saveAlignment(sim,file=paste("aln_",i,".fas",sep=""))
  return(sim)
  # return(TRUE)
}

# Note that the states are cleared and resampled; the rate multipliers
# are resampled as well. The resulting alignments are stored in files aln_1.fas, aln_2.fas, aln_3.fas.

# Memory can be saved by throwing away the objects generated by the replication by returning TRUE from the
# sim.replicate function.

# Load the multicore package:
library(multicore)

# Run replicates in parallel:
nr.replicates <-3
res.objects<-mclapply(1:nr.replicates, sim.replicate)

# Print the resulting PhyloSim objects:
print(res.objects)

# Plot one of the resulting alignments:
plot(res.objects[[1]])


#
## Simulating amino acid sequences with "domains" and heterogeneous evolution.
## See also the package vignette (vignette("PhyloSim",package="phylosim")).
#
# Setting up the substitution processes:

# load PhyloSim
library("phylosim");

# Use the ll() method to list the methods and virtual fields implemented in the Sequence class:
ll(Sequence())

# Enable the "fast & careless mode":
PSIM_FAST <- TRUE

## The following: WAG, JTT and LG are proteins replacement matrices
wag<-WAG();	# Create a WAG substitution process.
jtt<-JTT();	# Create a JTT substitution process.
lg<-LG();	# Create a LG substitution process.
pam<-PAM();	# Create a PAM substitution process.

# Get an object summary for wag:
summary(wag);

# Get a bubble plot of wag:
plot(wag);

# Create a continous deletor process:
cont.del<-ContinuousDeletor(
  rate=0.5,	# global rate for this deletion process
  max.length=10,	# the maximum allowed deletion length
  dist=expression(rnorm(1,mean=5,sd=3))	# length sampling expression
);

# Creating the template sequence for the insertion process:
templ.seq.wag<-AminoAcidSequence(length=10); # this is just a sequence with length 10.

# Note that the template sequence state is undefined, so the states
# will be sampled from the equlibrium distribution of the substitution process(es).

# Clone the template sequence:
templ.seq.lg<-clone(templ.seq.wag);

# Create a continous insertor process object:
cont.ins.wag<-ContinuousInsertor(
  rate=0.5,	# global rate for this insertion process
  max.length=10,	# the maximum allowed insertion length
  dist=expression(rnorm(1,mean=5,sd=3))	# length sampling expression
);

# Create a continous insertor process object:
cont.ins.lg<-ContinuousInsertor(
  rate=0.005,	# global rate for this insertion process
  max.length=10,	# the maximum allowed insertion length
  dist=expression(rnorm(1,mean=5,sd=3)) # length sampling expression
);

# Setting up the template sequences for the insertion processes:
templ.seq.wag$processes<-list(list(wag,cont.ins.wag,cont.del));
templ.seq.lg$processes<-list(list(lg,cont.ins.lg,cont.del));

# Disabling write protection for the insertion processes:
cont.ins.wag$writeProtected<-FALSE;
cont.ins.lg$writeProtected<-FALSE;

# Setting the template sequence for the insertion processes:
cont.ins.wag$templateSeq<-templ.seq.wag;
cont.ins.lg$templateSeq<-templ.seq.lg;

# Setting up the insert hook for the insertion processes:
cont.ins.wag$insertHook<-function(seq,target.seq,event.pos,insert.pos){
  # Create rate variation among the sites of seq by the +I+G model.
  plusInvGamma(
    seq,		# the sequence object
    process=wag, 	# the substitution process
    pinv=0.4,    	# the proportion of invariant sites.
    shape=0.6	# gamma shape parameter
  );
  return(seq);

}

cont.ins.lg$insertHook<-function(seq,target.seq,event.pos,insert.pos){
  # Create rate variation among the sites of seq by the +I+G model.
  plusInvGamma(
    seq,		# the sequence
    process=lg, 	# the substitution process
    pinv=0.4,	# the proportion of invariant sites.
    shape=0.6	# gamma shape parameter
  );
  return(seq);

}

#
# Setting up the root sequence:
#

seq<-AminoAcidSequence(length=60); # Create a sequence of length 200.

# Create the process pattern:
process.pattern<-c(
  rep(list(list(wag,cont.del, cont.ins.wag)), times=20),	# Left linker model: WAG
  rep(list(list(jtt)), times=20),				# "Core" model: JTT
  rep(list(list(lg,cont.del, cont.ins.lg)), times=20)	# Right linker model: LG
);

# Apply the process pattern to the root sequence:
seq$processes<-process.pattern;

# Set up site specific rates:

# Iterate over sites:
for (i in 1:seq$length){
  # Set a low rate for the core sites:
  if(isAttached(seq$sites[[i]],jtt)){
    # Sample rate from a truncated normal distribution.
    while( (site.rate<-rnorm(1,mean=0.001,sd=0.01)) < 0 ){}
    # Set the rate multiplier.
    setRateMultipliers(seq,jtt,site.rate,index=i);
  }
  else if(isAttached(seq$sites[[i]],wag)){

    plusInvGamma(
      seq,          	# the sequence
      process=wag,   	# the substitution process
      pinv=0.4,     	# the proportion of invariant sites.
      shape=0.6,     	# gamma shape parameter.
      index=i		# index vector
    );

  }
  else if(isAttached(seq$sites[[i]],lg)){

    plusInvGamma(
      seq,		# the sequence
      process=lg,	# the substitution process
      pinv=0.4,	# the proportion of invariant sites.
      shape=0.6,	# gamma shape parameter.
      index=i		# index vector
    );

  }
}

# Sample the states from the attached substitution process(es):
sampleStates(seq);

print(seq); # Print the actual sequence.

plot(seq);  # Plot the "rate landscape".

# Read in the tree using APE:
tree<-read.tree(
  file="data/smalldemotree.nwk"	# the path to the tree file
);

# Create the simulation object:
sim<-PhyloSim(
  phylo=tree,	# the tree as an APE phylo object
  root.seq=seq	# the root sequence.
);

# Create a node hook function:
node.hook<-function(seq){
  for (site in seq$sites){
    if(isAttached(site,jtt)){
      attachProcess(site,pam); # Attach the PAM process to the core sites.
    }
  }
  return(seq);
}

# Attach the hook to node 8:
attachHookToNode(
  sim,		# PhyloSim object.
  node=8,		# the node
  fun=node.hook	# the node hook function
);

# Run the simulation:
Simulate(sim)

# Plot the resulting alingment alongside the tree:
plot(sim)

# Save the resulting alignment, skip internal nodes:
saveAlignment(
  sim,				# the phylo object
  file="example_V3.1_aln.fas",
  skip.internal=TRUE		# filename for alignment
);

# Disable fast mode:
rm(PSIM_FAST)


## Evolving codon sequences
## See also the package vignette (vignette("PhyloSim",package="phylosim")).
##

library(phylosim)

# Enable "fast & careless" mode:
PSIM_FAST<-TRUE;

# Construct a GY94 codon substitution model:
p<-GY94();

# Set the transition/transverion rate ratio:
p$kappa=2

# Sample codon frequencies from a normal distribution:
p$equDist<-abs(rnorm(61,mean=10,sd=3))

# Get object summary for p:
summary(p)

# Get a bubble plot of p:
plot(p,scale=0.8)

# Construct a discrete deletor process:
d<-DiscreteDeletor(
  rate=1,
  sizes=1:4,
  probs=c(4,3,2,1)/10
);

# Construct a discrete insertor process inserting neutrally evolving sites:
i<-DiscreteInsertor(
  rate=1.5,
  sizes=1:4,
  probs=c(4,3,2,1)/10,
  template.seq=CodonSequence(length=4,processes=list(list(p)))
);

# Construct root sequence and attach process p:

s<-CodonSequence(length=30,processes=list(list(p)))

# Sample omegas from a discrete model:
omegaVarM3(s,p,omegas=c(0,1,2),probs=c(2/4,1/4,1/4))

# Plot the omega values across sites:
plotParametersAtSites(s,p,"omega");

# Sample states:

sampleStates(s)

# Construct the simulation object:
sim<-PhyloSim(
  root.seq=s,
  phylo=read.tree("data/smalldemotree.nwk")
);

# Create a node hook function and attach to node 9:
node.hook<-function(seq){

  # Set all omegas to 1 (neutral):
  setOmegas(seq,p,1);
  # attach the deletion process:
  attachProcess(seq,d)
  # attach the insertion process:
  attachProcess(seq,i)

  return(seq);
}

attachHookToNode(
  sim,                    # PhyloSim object.
  node=9,                 # the node
  fun=node.hook           # the node hook function
);

# Disable fast mode just before simulation in order to preserve branch statistics:
rm(PSIM_FAST)

# Run the simulation:
Simulate(sim)

# Plot the resulting alingment alongside the tree:
plot(sim);

# Export the nonsynonymous substitution counts as a phylo object:
nsyn.subst<-exportStatTree(sim,"nr.nsyn.subst")

# Plot the exported phylo object:
plot(nsyn.subst)
nodelabels()

# Save the resulting alignment:
saveAlignment(sim,                            # the phylo object
              file="example_V3.2_aln.fas",      # filename for alignment
);


##
## Implementing an inverted duplication process
##
##

# load the package
library(phylosim)

# enable fast & careless mode
PSIM_FAST<-TRUE;

# construct a DiscreteInsertor process
ivd<-DiscreteInsertor(rate=0.04,sizes=c(4,6),probs=c(2/3,1/3));

# set template sequence just to make the process object happy:
ivd$templateSeq<-NucleotideSequence(length=1);

# Replace the function object stored in the
# generateBy virtual field. See the documentation of the
# GeneralInsertor class.
ivd$generateBy<-function(process=NA,length=NA,target.seq=NA,event.pos=NA,insert.pos=NA){
  # get the target sequence length
  target.length<-target.seq$length;
  # construct a vector with the positions to copy:
  positions<-(insert.pos+1):(insert.pos + length)
  # discard illegal positions:
  positions<-positions[ positions > 0 & positions <= target.length];
  # copy subsequence
  insert<-copySubSequence(target.seq,positions,process);
  # reverse complement sequence,
  # take care, the class of this objects is "Sequence":
  revComp.NucleotideSequence(insert);
  # do not allow nested insertions:
  setRateMultipliers(insert,ivd,0);
  # return insert
  return(insert);
}

# Now we have a process which performs inverted duplications.

# construct a JC69 process object
p<-JC69();

# construct root sequence object
s<-NucleotideSequence(length=50)

# attach processes via virtual field
s$processes<-list(list(p,ivd))

# sample states from the equilibrium
# distribution of the attached processes

sampleStates(s)
# detach the substitution process:
detachProcess(s,p)

# create among-sites rate variation for the inverted duplication
# process by sampling rate multipliers from an I+G model:
plusGamma(s,ivd,pinv=0.5,shape=0.5)

# construct simulation object
sim<-PhyloSim(root.seq=s, phylo=read.tree("data/smalldemotree.nwk"));

# run simulation
Simulate(sim)

# plot tree and alignment
plot(sim)
# save alingment
saveAlignment(sim,file="example_V3.3.fas");

# disable fast & careless mode
rm(PSIM_FAST)




