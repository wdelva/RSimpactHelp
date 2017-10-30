rm(list=ls())
setwd("/home/david/RSimpactHelp/R/Projects_2017/Example3PaperSimpact2/")

## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings, dplyr, adephylo,
               phyclust, DECIPHER,treedater,geiger,picante)

#######################
# Step 1: Run Simpact #
#######################

# Set up the maodel

age.distr <- agedistr.creator(shape = 5, scale = 65)
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                 population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
                                 population.nummen = 2000, #2500,
                                 population.numwomen = 2000, #2500,
                                 hivseed.time = 10,
                                 hivseed.type = "amount",
                                 hivseed.amount = 10, #30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 hivtransmission.param.a = -1,
                                 hivtransmission.param.b = -90,
                                 hivtransmission.param.c = 0.5,
                                 hivtransmission.param.f1 = log(2), #log(inputvector[2]) , #log(2),
                                 hivtransmission.param.f2 = log(log(1.4) / log(2)) / 5, #log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5, #log(log(1.4) / log(2)) / 5,
                                 formation.hazard.agegapry.gap_factor_man_age = -0.01, #-0.01472653928518528523251061,
                                 formation.hazard.agegapry.gap_factor_woman_age = -0.01, #-0.0726539285185285232510561,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 formation.hazard.agegapry.gap_factor_man_exp = -1, #-6,#-1.5,
                                 formation.hazard.agegapry.gap_factor_woman_exp = -1, #-6,#-1.5,
                                 formation.hazard.agegapry.gap_agescale_man = 0.25, #inputvector[3], # 0.25,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.25, #inputvector[3], # 0.25,#-0.30000007,#-0.03,
                                 debut.debutage = 15,
                                 conception.alpha_base = -2.5#inputvector[14]#-2.5#,
                                 #person.art.accept.threshold.dist.fixed.value = 0
)


cfg.list["formation.hazard.agegapry.baseline"] <- 2
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.fraction.log_viralload"] <- 0.3
cfg.list["dropout.interval.dist.uniform.min"] <- 100
cfg.list["dropout.interval.dist.uniform.max"] <- 200

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.cd4.threshold"] <- 0






# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
art.intro <- list()
art.intro["time"] <- 25
art.intro["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[4] ######### 0.5
art.intro["diagnosis.baseline"] <- 0#100
art.intro["monitoring.cd4.threshold"] <- 100 # 1200

### add something about diagnosis
art.intro["diagnosis.agefactor"] <- 0
art.intro["diagnosis.genderfactor"] <- 0
art.intro["diagnosis.diagpartnersfactor"] <- 0
art.intro["diagnosis.isdiagnosedfactor"] <- 0
### end of add-on about diagnosis



#art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

art.intro2 <- list()
art.intro2["time"] <- 25 + 5 # inputvector[5] ######### 30
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 3 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 38
art.intro5["monitoring.cd4.threshold"] <- 5000 # This is equivalent to immediate access
art.intro5["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[8] ########### 0.75

# tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status

interventionlist <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist # scenario(interventionlist, tasp.indicator)




inputvector <- c(123,1.1, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.5, 2.8, -0.2, -0.2, -2.5, -0.52, -0.05)


cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[3]
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[3]
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]




cfg <- cfg.list

cfg["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
cfg["monitoring.fraction.log_viralload"] <- 0.3
cfg["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

seedid <- inputvector[1]
#cfg["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
#cfg["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10] ######### -0.5
cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10] ######### -0.5
cfg["formation.hazard.agegapry.baseline"] <- inputvector[11]

cfg["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
cfg["conception.alpha_base"] <- inputvector[14] #is conception.alpha.base (higher up)
cfg["dissolution.alpha_0"] <- inputvector[15]
cfg["dissolution.alpha_4"] <- inputvector[16]


# Run Simpact
results <- simpact.run(configParams = cfg,
                       destDir = "temp",
                       agedist = age.distr,
                       seed = seedid,
                       intervention = intervention)

datalist <- readthedata(results)


## Save the output
# save(datalist, file = "TESTMaster.datalist.RData")

# Read saved output data set
datalist <- get(load("TESTMaster.datalist.RData"))

# Calculate pervalence
prevalence.df <- prevalence.calculator(datalist = datalist,
                                       agegroup = c(15, 30),
                                       timepoint = 10)
# prevalence.df

# Calculate incidence
incidence.df <- incidence.calculator(datalist = datalist,
                                     agegroup = c(15, 30), timewindow = c(5, 15))
# incidence.df

# Plot prevalence
prevalence.df.plot <-prevalence.plotter(datalist = datalist, agegroup = c(15, 50))



#################################
# Step 2: Transmission networks #
#################################


# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff2.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")

# Transmission networks in a form of constructing transmission trees
simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = datalist$itable$population.simtime[1])


###############################
# Step 3: Sequence simulation #
###############################


# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]] # transmission network for i^th seed

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)

    sharedbranch.lengths <- vcv(tree.i)
    root.to.leave.distances <- diag(sharedbranch.lengths)

    num.trees <- c(num.trees,tree.n$id[1])
    num.i <- c(num.i,i)

    # Save the transmission tree
    write.tree(tree.i, file = paste("tree.model1.seed",i,".nwk", sep = ""))

    tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))

    # # count number of trees generated (normally one)
    # numb.tr <- function(tree=tree){
    #     if(length(tree) == 4){
    #         return(1)
    #     }else{
    #         return(length(tree))
    #     }
    # }
    #
    # # Simulate sequences, this require compiled toold seq-gen
    # # random sequence which will be simulated is choosen in the pool
    # seq.rand <- sample(1:30,1) # chose one sequence in first 30's in the seed sequence pool named "seed.seq.fasta"
    seq.rand <- 3
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("HIV.Pol.gene.fasta", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 3.1426, 1.6107, 6.0991, 9.1822, 6.8009, 1.0000  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # # z: seed
    #
    # Model = GTR+I+G
    # partition = 012345
    # -lnL = 15886.7332
    # K = 16
    # freqA = 0.3905
    # freqC = 0.1681
    # freqG = 0.2271
    # freqT = 0.2144
    # R(a) [AC] =  3.1426
    # R(b) [AG] =  1.6107
    # R(c) [AT] =  6.0991
    # R(d) [CG] =  9.1822
    # R(e) [CT] =  6.8009
    # R(f) [GT] =  1.0000
    # p-inv = 0.0000
    # gamma shape = 1.4240
    #

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingdates_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds
# 2 - 5 - 6
# 307 - 507 - 2519

#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # 1977+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}

# 4.1. Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]
  # External IQ-TREE
  # system(" ./iqtree-omp -s HIV.Env.gene.fasta -m GTR+R4 -nt AUTO -alrt 1000 -bb 1000")
  # Consensus tree written to HIV.Env.gene.fasta.contree
  # Reading input trees file HIV.Env.gene.fasta.contree
  # Log-likelihood of consensus tree: -10565.685
  #
  # Analysis results written to:
  #   IQ-TREE report:                HIV.Env.gene.fasta.iqtree
  # Maximum-likelihood tree:       HIV.Env.gene.fasta.treefile
  # Likelihood distances:          HIV.Env.gene.fasta.mldist
  #
  # Ultrafast bootstrap approximation results written to:
  #   Split support values:          HIV.Env.gene.fasta.splits.nex
  # Consensus tree:                HIV.Env.gene.fasta.contree
  # Screen log file:               HIV.Env.gene.fasta.log

  system(paste("./iqtree-omp -s", paste("sequence_all_seed_number_",id.trans,"_model1.fasta", sep = ""), "-m GTR+R4 -nt AUTO -alrt 1000 -bb 1000"))


}



# 4.2. Internal node optimisation requires sampled dates
IDs.transm <- c(2,5,6)
for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]

  tree.const <- read.tree(paste("sequence_all_seed_number_",id.trans,"_model1.fasta.treefile", sep = ""))

  samp.dates <- read.csv(paste("samplingdates_seed_number_",id.trans,".csv", sep = ""))

  time.samp <- dates.Transform.NamedVector(dates=samp.dates)

  tree.tips <- as.numeric(tree.const$tip.label)

  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 3012, omega0 = 3.2) # s is the length of sequence

  # dater(tre, sts, s = 1000, omega0 = NA, minblen = NA, maxit = 100, abstol = .001, quiet = TRUE, searchRoot=5, temporalConstraints=TRUE, strictClock=FALSE, estimateSampleTimes=NULL, estimateSampleTimes_densities = list())

  otip <- outlier.tips(dater.tree)

  # Warning messages:
  #   1: In dater(t, sts, s = s, omega0 = omega0, minblen = minblen, maxit = maxit,  :
  #                 Root to tip regression predicts a substition rate less than zero. Tree may be poorly rooted or there may be small temporal signal.
  #               2: In dater(t, sts, s = s, omega0 = omega0, minblen = minblen, maxit = maxit,  :
  #                             Root to tip regression predicts a substition rate less than zero. Tree may be poorly rooted or there may be small temporal signal.
  #                           3: In dater(t, sts, s = s, omega0 = omega0, minblen = minblen, maxit = maxit,  :
  #                                         Root to tip regression predicts a substition rate less than zero. Tree may be poorly rooted or there may be small temporal signal.
  #                                       4: In dater(t, sts, s = s, omega0 = omega0, minblen = minblen, maxit = maxit,  :
  #                                                     Root to tip regression predicts a substition rate less than zero. Tree may be poorly rooted or there may be small temporal signal.
  #                                                   5: In dater(t, sts, s = s, omega0 = omega0, minblen = minblen, maxit = maxit,  :
  #                                                                 Root to tip regression predicts a substition rate less than zero. Tree may be poorly rooted or there may be small temporal signal.
  #

  # Save the tree
  write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))
  ###########################################################################################
  ###########################################################################################
  #
  #   # With changed label names, make sure the vector of sampling dates is renamed as wells
  #
  #   # Rename the tree by adding sampling date
  # g <- tree.const
  # a <- g$tip.label # label names
  # d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
  # g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # renamed.Ord.tree.dates <- Ord.tree.dates
  # names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # # Use of library(treedater) to calibrate internal nodes
  # dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 1000) # s is the length of sequence

  #   # Save the tree
  #   write.tree(dater.tree.g, file = paste("calibratedTree_",id.trans,"_pref.nwk", sep = ""))



}

########################################### ------------------------------------ ########################

# Chose tree of seed 3 with  314 tips and 313 internal nodes.

tree.const <- read.tree(paste("sequence_all_seed_number_",5,"_model1.fasta.treefile", sep = ""))
samp.dates <- read.csv(paste("samplingdates_seed_number_",5,".csv", sep = ""))

time.samp <- dates.Transform.NamedVector(dates=samp.dates)

tree.tips <- as.numeric(tree.const$tip.label)

# Sequence labels ordered in the same order of the tree tips labels

Ord.tree.dates <- vector()
for(i in 1:length(tree.tips)){
  for(j in 1:length(time.samp)){
    if(tree.tips[i] == samp.dates$V1[j]){
      Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
    }
  }
}

# paste("tree.model1.seed",i,".nwk", sep = "")
g <- read.tree("tree.model1.seed5.nwk")  #tree.const

a <- g$tip.label # label names
d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates

renamed.Ord.tree.dates <- Ord.tree.dates
names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates

# Use of library(treedater) to calibrate internal nodes
dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 3012)#,  omega0 = 3.2, minblen = 3.2, maxit = 1000) # s is the length of sequence

otip.g <- outlier.tips(dater.tree.g)

# Node age with picante package

N <- node.age(dater.tree.g)

# potential transmission times
int.node.age <- N$Ti # internal nodes ages
##########################################


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


## Transmission network of seed 3

tra.net.3 <- trans.net[[5]]

tra.net.3$dtimes <- datalist$itable$population.simtime[1] - tra.net.3$dtimes + 1977 #(endpoint=40)
tra.net.3$itimes <- datalist$itable$population.simtime[1] - tra.net.3$itimes + 1977 #(endpoint=40) 10>1990, -> +1980

min.val = 1977
max.val = 1977 + datalist$itable$population.simtime[1]


step.int=1
d <- (max.val-min.val)/step.int

dat.f.trans <- as.data.frame(tra.net.3)
dt.node.age.dt <- int.node.age

numb.tra <- vector()
i.vec <- vector()
int.node.vec <- vector()
for (i in 1:d) {
  inf <- 1976+i
  sup <- 1977+i
  dat.f.trans.i <- dat.f.trans[which(dat.f.trans$itimes < sup & dat.f.trans$itimes  >= inf),]
  numb.i <- nrow(dat.f.trans.i)
  numb.tra <- c(numb.tra, numb.i)
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age < sup & dt.node.age.dt >= inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}



graph.build <- as.data.frame(trans.net[[5]])

graph.build[,4] <- as.character(graph.build$parent) # donors
graph.build[,3] <- as.character(graph.build$id) # recipients
gag = as.matrix(graph.build)
ga.graph = graph.edgelist(gag[,4:3])

V(ga.graph)$color <- "red"

transNet.yrs.Old <- delete.vertices(ga.graph, "-1")


## Use phangorn for inputTree Seed5
seq.sim.size_full <- read.dna("~/RSimpactHelp/R/Projects_2017/Example3PaperSimpact/sequence_all_seed_number_5_model1.fasta")
tree.dat.full <- phyDat(seq.sim.size_full, type = "DNA")
tree.ml.full <- dist.ml(tree.dat.full)
tree.sim.full <- nj(tree.ml.full)

### Plot figures
#################

# 1. Transmission network from simpact                           # 1 #
plot.igraph(transNet.yrs.Old, edge.arrow.size=0.1, vertex.size=5,
            edge.color="black",
            asp = 1,
            xlim = c(-1, 2),
            ylim = c(-0.5,0.5),
            vertex.frame.color="black",
            vertex.label.color="black",
            #vertex.label = NULL,
            layout = layout_with_kk,
            edge.width = 1,
            vertex.label.cex=0.1,
            vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0
            #main = "True transmission network"
)

# 2. Phylogenetic tree

plot(dater.tree.g, show.tip.label=FALSE,
     #edge.width=1,
     edge.color="blue") # Try a few different settings!
axisPhylo(backward = FALSE)


# 3. Transmission events and internal nodes
x <- i.vec
plot(x, int.node.vec, type="b", col="red", lwd=2,
     xlab = "Calendar time",
     ylab = "Count") # 1 > 1
lines(x, numb.tra, col='green3', type='b', lwd=2)
legend("topleft", legend = c("Internal nodes", "Transmission events"),
       col=c("red","green3"), pch=1)

### Removing non diagnosed individuals

seq.rm.5 <- read.dna("sequence_all_seed_number_5_model1.fasta")

# Access labels
labels(seq.rm.5)

# Non diagnosed individuals
dtframe <- as.data.frame(tra.net.3)
diag.ind <- dtframe[which(dtframe$dtypes==0),]
non.diag.died.ind <- dtframe[which(dtframe$dtypes==-1),]
non.diag.live.ind <- dtframe[which(dtframe$dtypes==1),]

ids.non.diag <- c(non.diag.died.ind$id, non.diag.live.ind$id)

w <- seq.rm.5[ids.non.diag, drop=TRUE]

## tips to root

library(adephylo)

rt <- distRoot(g)


## distance between tips
distTip <- distTips(g)

## branch lengths of tips
# Branch length of the tips
tree <- phylo.tree
# tree$edge.length<-round(tree$edge.length,3)
tree$edge.length<-tree$edge.length # remove rounding
n<-length(tree$tip.label)
ee<-setNames(tree$edge.length[sapply(1:n,function(x,y)
  which(y==x),y=tree$edge[,2])],tree$tip.label)

d <- as.data.frame(ee)
branc.leng.tips <- as.data.frame(d$ee) # yes we can find them in phylo.tree$edge.length but here they are not rounded whereas
# with ee are rounded with 3 digit
tips.id <- as.data.frame(as.character(tree$tip.label))

tips.id <- tree$tip.label
tips.branch.leng.raw <- cbind(tips.id,branc.leng.tips)

names(tips.branch.leng.raw) <- c("Ind", "branch.len")

tips.branch.leng <- tips.branch.leng.raw

# branch lengths of tips == dtimes-itimes
sort(tips.branch.leng[,2])
sort(tra.net.3$dtimes-tra.net.3$itimes)

