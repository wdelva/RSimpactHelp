# ## For sub-optimal sequence coverage
#

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/")


## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings, dplyr, adephylo,
               phyclust, DECIPHER,treedater,geiger,picante)



#######################
# Step 1: Run Simpact # # for A: 4000 individuals(1800M&2200W) and B: 8000 individuals(3800M&4200W)
#######################

## Run Simpact for specific parameter combination

age.distr <- agedistr.creator(shape = 5, scale = 65)
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                 population.msm = "no",
                                 population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
                                 population.nummen = 3800, # 3800, #2500,
                                 population.numwomen = 4200, #4200, #2500,
                                 hivseed.time = 10, # 10,
                                 hivseed.type = "amount",
                                 hivseed.amount = 40, #30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 hivtransmission.param.a = -1, # -1,
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
cfg.list["monitoring.fraction.log_viralload"] <- 0 # 0.3
cfg.list["dropout.interval.dist.uniform.min"] <- 1000
cfg.list["dropout.interval.dist.uniform.max"] <- 2000

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.cd4.threshold"] <- 10000
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 1
cfg.list["diagnosis.baseline"] <- -2






# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
art.intro <- list()
art.intro["time"] <- 0.0001 # 25
art.intro["person.art.accept.threshold.dist.fixed.value"] <- 1 # 0.5 # inputvector[4] ######### 0.5
art.intro["diagnosis.baseline"] <- -2 # 0#100
art.intro["monitoring.cd4.threshold"] <- 10000 # 1200

### add something about diagnosis
art.intro["diagnosis.agefactor"] <- 0
art.intro["diagnosis.genderfactor"] <- 0
art.intro["diagnosis.diagpartnersfactor"] <- 0
art.intro["diagnosis.isdiagnosedfactor"] <- 0
### end of add-on about diagnosis



#art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

# art.intro2 <- list()
# art.intro2["time"] <- 25 + 5 # inputvector[5] ######### 30
# art.intro2["monitoring.cd4.threshold"] <- 200
#
# art.intro3 <- list()
# art.intro3["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
# art.intro3["monitoring.cd4.threshold"] <- 350
#
# art.intro4 <- list()
# art.intro4["time"] <- 3 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
# art.intro4["monitoring.cd4.threshold"] <- 500
#
# art.intro5 <- list()
# art.intro5["time"] <- 38
# art.intro5["monitoring.cd4.threshold"] <- 5000 # This is equivalent to immediate access
# art.intro5["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[8] ########### 0.75

# tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status

interventionlist <- list(art.intro) #, art.intro2, art.intro3, art.intro4, art.intro5)

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
# cfg["monitoring.fraction.log_viralload"] <- 0.3
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

# #
# # # # Run Simpact
# results <- simpact.run(configParams = cfg,
#                        destDir = "temp",
#                        agedist = age.distr,
#                        seed = seedid,
#                        intervention = intervention)

datalist <- readthedata(results)
#
# table(datalist$etable$eventname)
#
# ## Save the output
# save(datalist, file = "MasterModelSubOptimalSeqCovearge.datalistA.RData")


# Read saved output data set
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistB.RData"))


table(datalist$etable$eventname) # check events

# Prevlence per age group and time at a given time-point
prevalence.df <- prevalence.calculator(datalist = datalist,
                                       agegroup = c(15, 30),
                                       timepoint = 10)

prevalence.df.plot <-prevalence.plotter(datalist = datalist, agegroup = c(15, 50))


# Incidence per age group in a given time-window
incidence.df <- incidence.calculator(datalist = datalist,
                                     agegroup = c(15, 30), timewindow = c(10, 20))



###########################################
# Step 2: Construct transmission networks #
###########################################


# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff2.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
source("/home/david/RSimpactHelp/R/epi2tree2.R")


simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)

smallest.branches <- rep(NA, times = length(simpact.trans.net))
for (list.element in 1:length(simpact.trans.net)){
  net.list <- simpact.trans.net[[list.element]]
  if(length(net.list$id) > 2){
    tree.tr <- epi2tree2(net.list)
    smallest.branch <- min(tree.tr$edge.length)
    smallest.branches[list.element] <- smallest.branch
  }
}
min(smallest.branches, na.rm = TRUE) #
## seeds and transmission network sizes:
# A: 10 >>6, 18 >>20, 19 >>3, 21 >>3, 22 >>1403, 24 >>716, 27 >>370, 32 >>20, 37>>5
# B: 4 >>3, 12 >>1079, 19 >>4, 21 >>3, 22 >>3032, 24 >>809, 28 >>3, 34 >>3 39 >>5
which(smallest.branches!="NA")

###################################################################################################

# Scenario 1 - A
#
# one subtype of the virus (HIV-1-A) for all seeds
# complete sampling for a transmission network of one seed
# same sampling time interval (e.g.: five or three years) for a transmission network of all seeds

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario1_B/")



## (i) Sequence simulation

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
    seq.rand <- 1
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]

    # call the seed sequences and rename the file
    file.copy(paste("hiv.seq.A.pol.j.fasta", sep = ""),paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -s 0.00475 -n1 -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds


#######################
# Scenario 1 Analysis #
#######################

setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/Complete_Transmission_Seed_22/")
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistB.RData"))



# 1.1. Full phylogenetic tree for seed 22 - 100% of sequences

# Sampling dates in calender time
dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1977+40-as.numeric(dates$V2) # dates datalist$itable$population.simtime[1] - dates$V2 + 1977
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.100.3032.fasta", sep = ""), paste(">A.seed.22.seq.100.3032.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22 <- read.tree(paste("A.seed.22.seq.100.3032.fasta.tree", sep = ""))

samp.dates.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22 <- dates.Transform.NamedVector(dates=samp.dates.22)

tree.tips.22 <- as.numeric(tree.fasttree.22$tip.label)

Ord.tree.dates.22 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22)){
  for(j in 1:length(time.samp.22)){
    if(tree.tips.22[i] == samp.dates.22$V1[j]){
      Ord.tree.dates.22 <- c(Ord.tree.dates.22, time.samp.22[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22 <- dater(tree.fasttree.22,
                       Ord.tree.dates.22,
                       s = 3012,
                       omega0 = 0.00475) # s is the length of sequence
d=node.age(dater.tree.22)

save(dater.tree.22, file = paste("dated.tree.A.seed.22.seq.100.3032.Rdata", sep = ""))

# 1.2. Full network for seed 2 - 100% of individuals

simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)
simpact.trans.net.22 <- simpact.trans.net[[22]]


graph.net.22 <- as.data.frame(simpact.trans.net.22)

graph.build.22 <- graph.net.22[,3:4]

graph.build.22[,2] <- as.character(graph.build.22[,2]) # donors
graph.build.22[,1] <- as.character(graph.build.22[,1]) # recipients
graph.22 = as.matrix(graph.build.22)
graph.f.22 = graph.edgelist(graph.22[,1:2])
# E(graph.f.22)$weight <- infage[[2]]$infecage
# V(graph.f.22)$color <- "red"
source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/properties_network.R")
properties.trans.net.22 <- properties_network(graph = graph.f.22)
save(properties.trans.net.22, file = paste("trans.net.A.seed.22.seq.100.3032.Rdata", sep = ""))


# 2.1. Phylogenetic tree for seed 22 - 90% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.90.2729.fasta", sep = ""), paste(">A.seed.22.seq.90.2729.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.90 <- read.tree(paste("A.seed.22.seq.90.2729.fasta.tree", sep = ""))

samp.dates.22.90 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.90 <- dates.Transform.NamedVector(dates=samp.dates.22.90)

tree.tips.22.90 <- as.numeric(tree.fasttree.22.90$tip.label)

Ord.tree.dates.22.90 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.90)){
  for(j in 1:length(time.samp.22.90)){
    if(tree.tips.22.90[i] == samp.dates.22.90$V1[j]){
      Ord.tree.dates.22.90 <- c(Ord.tree.dates.22.90, time.samp.22.90[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.90 <- dater(tree.fasttree.22.90,
                          Ord.tree.dates.22.90,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.90, file = paste("dated.tree.A.seed.22.seq.90.2729.Rdata", sep = ""))



# 3.1. Phylogenetic tree for seed 22 - 80% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.80.2426.fasta", sep = ""), paste(">A.seed.22.seq.80.2426.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.80 <- read.tree(paste("A.seed.22.seq.80.2426.fasta.tree", sep = ""))

samp.dates.22.80 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.80 <- dates.Transform.NamedVector(dates=samp.dates.22.80)

tree.tips.22.80 <- as.numeric(tree.fasttree.22.80$tip.label)

Ord.tree.dates.22.80 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.80)){
  for(j in 1:length(time.samp.22.80)){
    if(tree.tips.22.80[i] == samp.dates.22.80$V1[j]){
      Ord.tree.dates.22.80 <- c(Ord.tree.dates.22.80, time.samp.22.80[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.80 <- dater(tree.fasttree.22.80,
                          Ord.tree.dates.22.80,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.80, file = paste("dated.tree.A.seed.22.seq.80.2426.Rdata", sep = ""))



# 4.1. Phylogenetic tree for seed 22 - 70% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.70.2122.fasta", sep = ""), paste(">A.seed.22.seq.70.2122.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.70 <- read.tree(paste("A.seed.22.seq.70.2122.fasta.tree", sep = ""))

samp.dates.22.70 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.70 <- dates.Transform.NamedVector(dates=samp.dates.22.70)

tree.tips.22.70 <- as.numeric(tree.fasttree.22.70$tip.label)

Ord.tree.dates.22.70 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.70)){
  for(j in 1:length(time.samp.22.70)){
    if(tree.tips.22.70[i] == samp.dates.22.70$V1[j]){
      Ord.tree.dates.22.70 <- c(Ord.tree.dates.22.70, time.samp.22.70[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.70 <- dater(tree.fasttree.22.70,
                          Ord.tree.dates.22.70,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.70, file = paste("dated.tree.A.seed.22.seq.70.2122.Rdata", sep = ""))


# 5.1. Phylogenetic tree for seed 22 - 60% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.60.1819.fasta", sep = ""), paste(">A.seed.22.seq.60.1819.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.60 <- read.tree(paste("A.seed.22.seq.60.1819.fasta.tree", sep = ""))

samp.dates.22.60 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.60 <- dates.Transform.NamedVector(dates=samp.dates.22.60)

tree.tips.22.60 <- as.numeric(tree.fasttree.22.60$tip.label)

Ord.tree.dates.22.60 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.60)){
  for(j in 1:length(time.samp.22.60)){
    if(tree.tips.22.60[i] == samp.dates.22.60$V1[j]){
      Ord.tree.dates.22.60 <- c(Ord.tree.dates.22.60, time.samp.22.60[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.60 <- dater(tree.fasttree.22.60,
                          Ord.tree.dates.22.60,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.60, file = paste("dated.tree.A.seed.22.seq.60.1819.Rdata", sep = ""))


# 6.1. Phylogenetic tree for seed 22 - 50% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.50.1516.fasta", sep = ""), paste(">A.seed.22.seq.50.1516.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.50 <- read.tree(paste("A.seed.22.seq.50.1516.fasta.tree", sep = ""))

samp.dates.22.50 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.50 <- dates.Transform.NamedVector(dates=samp.dates.22.50)

tree.tips.22.50 <- as.numeric(tree.fasttree.22.50$tip.label)

Ord.tree.dates.22.50 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.50)){
  for(j in 1:length(time.samp.22.50)){
    if(tree.tips.22.50[i] == samp.dates.22.50$V1[j]){
      Ord.tree.dates.22.50 <- c(Ord.tree.dates.22.50, time.samp.22.50[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.50 <- dater(tree.fasttree.22.50,
                          Ord.tree.dates.22.50,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.50, file = paste("dated.tree.A.seed.22.seq.50.1516.Rdata", sep = ""))


# 7.1. Phylogenetic tree for seed 22 - 40% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.40.1213.fasta", sep = ""), paste(">A.seed.22.seq.40.1213.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.40 <- read.tree(paste("A.seed.22.seq.40.1213.fasta.tree", sep = ""))

samp.dates.22.40 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.40 <- dates.Transform.NamedVector(dates=samp.dates.22.40)

tree.tips.22.40 <- as.numeric(tree.fasttree.22.40$tip.label)

Ord.tree.dates.22.40 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.40)){
  for(j in 1:length(time.samp.22.40)){
    if(tree.tips.22.40[i] == samp.dates.22.40$V1[j]){
      Ord.tree.dates.22.40 <- c(Ord.tree.dates.22.40, time.samp.22.40[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.40 <- dater(tree.fasttree.22.40,
                          Ord.tree.dates.22.40,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.40, file = paste("dated.tree.A.seed.22.seq.40.1213.Rdata", sep = ""))


# 8.1. Phylogenetic tree for seed 22 - 30% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.30.910.fasta", sep = ""), paste(">A.seed.22.seq.30.910.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.30 <- read.tree(paste("A.seed.22.seq.30.910.fasta.tree", sep = ""))

samp.dates.22.30 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.30 <- dates.Transform.NamedVector(dates=samp.dates.22.30)

tree.tips.22.30 <- as.numeric(tree.fasttree.22.30$tip.label)

Ord.tree.dates.22.30 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.30)){
  for(j in 1:length(time.samp.22.30)){
    if(tree.tips.22.30[i] == samp.dates.22.30$V1[j]){
      Ord.tree.dates.22.30 <- c(Ord.tree.dates.22.30, time.samp.22.30[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.30 <- dater(tree.fasttree.22.30,
                          Ord.tree.dates.22.30,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.30, file = paste("dated.tree.A.seed.22.seq.30.910.Rdata", sep = ""))


# 9.1. Phylogenetic tree for seed 22 - 20% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.20.606.fasta", sep = ""), paste(">A.seed.22.seq.20.606.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.20 <- read.tree(paste("A.seed.22.seq.20.606.fasta.tree", sep = ""))

samp.dates.22.20 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/samplingtimes_seed_number_22.csv")

time.samp.22.20 <- dates.Transform.NamedVector(dates=samp.dates.22.20)

tree.tips.22.20 <- as.numeric(tree.fasttree.22.20$tip.label)

Ord.tree.dates.22.20 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.20)){
  for(j in 1:length(time.samp.22.20)){
    if(tree.tips.22.20[i] == samp.dates.22.20$V1[j]){
      Ord.tree.dates.22.20 <- c(Ord.tree.dates.22.20, time.samp.22.20[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.20 <- dater(tree.fasttree.22.20,
                          Ord.tree.dates.22.20,
                          s = 2012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.20, file = paste("dated.tree.A.seed.22.seq.20.606.Rdata", sep = ""))

# Scenario 2 - A
#
# one subtype of the virus (HIV-1-A) for all seeds
# complete sampling for a transmission network of all seeds
# same sampling time interval (e.g.: five or three years) for a transmission network of all seeds

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario2_B/")


## (i) Sequence simulation

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
    seq.rand <- 1
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]

    # call the seed sequences and rename the file
    file.copy(paste("hiv.seq.A.pol.j.fasta", sep = ""),paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -s 0.00475 -n1 -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds


#######################
# Scenario 2 Analysis #
#######################

setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/")

datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistB.RData"))
simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)
save(simpact.trans.net, file="simpact.trans.net.B.RData")
seeds.great3 <- c(4, 12, 19, 21, 22, 24, 28, 34, 39) # seeds with at least 3 transmission events

# Remain with Seeds which have sampled individuals in past 7 years
seed.7yrs.samp <- vector()
for (i in 1:length(seeds.great3)){
  k <- seeds.great3[i]
  net.i <- as.data.frame(simpact.trans.net[[k]])
  net.filtered.i <- filter(net.i, dtimes<=7)
  dtimes <- net.i$dtimes<=7
  dtimes.elements <- which(dtimes)
  if(length(dtimes.elements)!=0){
    seed.7yrs.samp <- c(seed.7yrs.samp, k)
    save(net.filtered.i, file=paste("filtered.trans.net.seed",k,".Rdata", sep = ""))
  }
}

net.seed12 = get(load("filtered.trans.net.seed12.Rdata")) # length(net.seed12$id) = 314
net.seed22 = get(load("filtered.trans.net.seed22.Rdata")) # length(net.seed22$id) = 1579
net.seed24 = get(load("filtered.trans.net.seed24.Rdata")) # length(net.seed24$id) = 11

# IDs to exclude in the sequences
id.12 <- setdiff(simpact.trans.net[[12]]$id, net.seed12$id) # Done! - 765 Ok
id.22 <- setdiff(simpact.trans.net[[22]]$id, net.seed22$id) # Done! - 1453 Ok
id.24 <- setdiff(simpact.trans.net[[24]]$id, net.seed24$id) # Done! - 798 Ok


#
# seq.sim <- read.dna("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/A.Epidemic12.Sequences.gene.pol.314.fasta")
# tree.dat <- phyDat(seq.sim, type = "DNA")
# tree.ml <- dist.ml(tree.dat,model = "JC69")
# tree.sim <- upgma(tree.ml)
#
# d <- sort(tree.sim$tip.label)


# # rename sampling times automatically and make a single file for sampling dates of seeds 12, 22, 24
# First load and rename samplingtimes_seed_number_12.csv, .. > Tsamplingtimes_seed_number_12.csv, ..
# After load and rename Tsamplingtimes_seed_number_12.csv, .. in dd

#
# transfSampT.12 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.csv")
# transfSampT.12$V1 <- paste("A.12.",as.character(transfSampT.12$V1), sep = "")
# write.csv(transfSampT.12, file = "Tsamplingtimes_seed_number_12.csv", quote = TRUE)
#
#
# transfSampT.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.csv")
# transfSampT.22$V1 <- paste("A.22.",as.character(transfSampT.12$V1), sep = "")
# write.csv(transfSampT.22, file = "Tsamplingtimes_seed_number_22.csv", quote = TRUE)
#
# transfSampT.24 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_24.csv")
# transfSampT.24$V1 <- paste("A.24.",as.character(transfSampT.12$V1), sep = "")
# write.csv(transfSampT.24, file = "Tsamplingtimes_seed_number_24.csv", quote = TRUE)
#
#
# transfSampTimes.12 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_12.csv")
# transfSampTimes.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_22.csv")
# transfSampTimes.24 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_24.csv")
#
# d <- rbind(transfSampTimes.12[,3:4], transfSampTimes.22[,3:4], transfSampTimes.24[,3:4])
# dd <- as.data.frame(d)
# dd$V1 <- as.character(dd$V1)
# dd$V2 <- as.numeric(as.character(dd$V2)) # solve factors issue
#
# write.csv(dd, file = "samplingtimes_seed_number_22.24.27.csv")

# 1.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 100% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.1904.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.1904.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.1904.fasta.tree", sep = ""))

samp.dates.12.22.24 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24$V1 <- as.character(samp.dates.12.22.24$V1)
samp.dates.12.22.24$V2 <- as.numeric(as.character(samp.dates.12.22.24$V2))

time.samp.12.22.24 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24)

tree.tips.12.22.24 <- (tree.fasttree.12.22.24$tip.label) # length==1904

c.samp.dates.12.22.24<- as.character(samp.dates.12.22.24$V1)

Ord.tree.dates.12.22.24 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24)){
  for(j in 1:length(time.samp.12.22.24)){
    if(tree.tips.12.22.24[i] == c.samp.dates.12.22.24[j]){
      Ord.tree.dates.12.22.24 <- c(Ord.tree.dates.12.22.24, time.samp.12.22.24[j])
      t.tip <- c(t.tip,tree.tips.12.22.24[i])
      t.date <- c(t.date, c.samp.dates.12.22.24[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24 <- dater(tree.fasttree.12.22.24,
                             Ord.tree.dates.12.22.24,
                             s = 3012,
                             omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24, file = paste("dated.tree.A.seed.12.22.24Seq.p100.1904.Rdata", sep = ""))
#e <- get(load("dated.tree.A.seed.12.22.24seq.1904.Rdata"))

tr.scenario2B <- get(load("dated.tree.A.seed.12.22.24Seq.p100.1904.Rdata"))
write.tree(tr.scenario2B, file = "tr.scenario2B.12.22.24.1904.HIVA.tree")

# 2.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 90% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.90.1714.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.90.1714.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p90 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.90.1714.fasta.tree", sep = ""))

samp.dates.12.22.24.p90 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.p90$V1 <- as.character(samp.dates.12.22.24.p90$V1)
samp.dates.12.22.24.p90$V2 <- as.numeric(as.character(samp.dates.12.22.24.p90$V2))

time.samp.12.22.24.p90 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p90)

tree.tips.12.22.24.p90 <- (tree.fasttree.12.22.24.p90$tip.label) # length==1904

c.samp.dates.12.22.24.p90<- as.character(samp.dates.12.22.24.p90$V1)

Ord.tree.dates.12.22.24.p90 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p90)){
  for(j in 1:length(time.samp.12.22.24.p90)){
    if(tree.tips.12.22.24.p90[i] == c.samp.dates.12.22.24.p90[j]){
      Ord.tree.dates.12.22.24.p90 <- c(Ord.tree.dates.12.22.24.p90, time.samp.12.22.24.p90[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p90[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p90[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p90 <- dater(tree.fasttree.12.22.24.p90,
                                 Ord.tree.dates.12.22.24.p90,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p90, file = paste("dated.tree.A.seed.12.22.24Seq.p90.1714.Rdata", sep = ""))



# 3.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 80% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.80.1523.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.80.1523.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p80 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.80.1523.fasta.tree", sep = ""))

samp.dates.12.22.24.p80 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")


samp.dates.12.22.24.p80$V1 <- as.character(samp.dates.12.22.24.p80$V1)
samp.dates.12.22.24.p80$V2 <- as.numeric(as.character(samp.dates.12.22.24.p80$V2))

time.samp.12.22.24.p80 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p80)

tree.tips.12.22.24.p80 <- (tree.fasttree.12.22.24.p80$tip.label) # length==1904

c.samp.dates.12.22.24.p80<- as.character(samp.dates.12.22.24.p80$V1)

Ord.tree.dates.12.22.24.p80 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p80)){
  for(j in 1:length(time.samp.12.22.24.p80)){
    if(tree.tips.12.22.24.p80[i] == c.samp.dates.12.22.24.p80[j]){
      Ord.tree.dates.12.22.24.p80 <- c(Ord.tree.dates.12.22.24.p80, time.samp.12.22.24.p80[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p80[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p80[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p80 <- dater(tree.fasttree.12.22.24.p80,
                                 Ord.tree.dates.12.22.24.p80,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p80, file = paste("dated.tree.A.seed.12.22.24seq.p80.1523.Rdata", sep = ""))


# 4.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 70% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.70.1333.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.70.1333.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p70 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.70.1333.fasta.tree", sep = ""))

samp.dates.12.22.24.p70 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.p70$V1 <- as.character(samp.dates.12.22.24.p70$V1)
samp.dates.12.22.24.p70$V2 <- as.numeric(as.character(samp.dates.12.22.24.p70$V2))

time.samp.12.22.24.p70 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p70)

tree.tips.12.22.24.p70 <- (tree.fasttree.12.22.24.p70$tip.label) # length==1904

c.samp.dates.12.22.24.p70<- as.character(samp.dates.12.22.24.p70$V1)

Ord.tree.dates.12.22.24.p70 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p70)){
  for(j in 1:length(time.samp.12.22.24.p70)){
    if(tree.tips.12.22.24.p70[i] == c.samp.dates.12.22.24.p70[j]){
      Ord.tree.dates.12.22.24.p70 <- c(Ord.tree.dates.12.22.24.p70, time.samp.12.22.24.p70[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p70[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p70[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p70 <- dater(tree.fasttree.12.22.24.p70,
                                 Ord.tree.dates.12.22.24.p70,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p70, file = paste("dated.tree.A.seed.12.22.24seq.p70.1333.Rdata", sep = ""))


# 5.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 60% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.60.1142.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.60.1142.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p60 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.60.1142.fasta.tree", sep = ""))

samp.dates.12.22.24.p60 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.p60$V1 <- as.character(samp.dates.12.22.24.p60$V1)
samp.dates.12.22.24.p60$V2 <- as.numeric(as.character(samp.dates.12.22.24.p60$V2))

time.samp.12.22.24.p60 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p60)

tree.tips.12.22.24.p60 <- (tree.fasttree.12.22.24.p60$tip.label) # length==1904

c.samp.dates.12.22.24.p60<- as.character(samp.dates.12.22.24.p60$V1)

Ord.tree.dates.12.22.24.p60 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p60)){
  for(j in 1:length(time.samp.12.22.24.p60)){
    if(tree.tips.12.22.24.p60[i] == c.samp.dates.12.22.24.p60[j]){
      Ord.tree.dates.12.22.24.p60 <- c(Ord.tree.dates.12.22.24.p60, time.samp.12.22.24.p60[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p60[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p60[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p60 <- dater(tree.fasttree.12.22.24.p60,
                                 Ord.tree.dates.12.22.24.p60,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p60, file = paste("dated.tree.A.seed.12.22.24seq.p60.1142.Rdata", sep = ""))



# 6.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 50% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.50.952.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.50.952.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p50 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.50.952.fasta.tree", sep = ""))

samp.dates.12.22.24.p50 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.p50$V1 <- as.character(samp.dates.12.22.24.p50$V1)
samp.dates.12.22.24.p50$V2 <- as.numeric(as.character(samp.dates.12.22.24.p50$V2))

time.samp.12.22.24.p50 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p50)

tree.tips.12.22.24.p50 <- (tree.fasttree.12.22.24.p50$tip.label) # length==1904

c.samp.dates.12.22.24.p50<- as.character(samp.dates.12.22.24.p50$V1)

Ord.tree.dates.12.22.24.p50 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p50)){
  for(j in 1:length(time.samp.12.22.24.p50)){
    if(tree.tips.12.22.24.p50[i] == c.samp.dates.12.22.24.p50[j]){
      Ord.tree.dates.12.22.24.p50 <- c(Ord.tree.dates.12.22.24.p50, time.samp.12.22.24.p50[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p50[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p50[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p50 <- dater(tree.fasttree.12.22.24.p50,
                                 Ord.tree.dates.12.22.24.p50,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p50, file = paste("dated.tree.A.seed.12.22.24seq.p50.952.Rdata", sep = ""))



# 7.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 40% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.40.762.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.40.762.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p40 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.40.762.fasta.tree", sep = ""))

samp.dates.12.22.24.p40 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.p40$V1 <- as.character(samp.dates.12.22.24.p40$V1)
samp.dates.12.22.24.p40$V2 <- as.numeric(as.character(samp.dates.12.22.24.p40$V2))

time.samp.12.22.24.p40 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p40)

tree.tips.12.22.24.p40 <- (tree.fasttree.12.22.24.p40$tip.label) # length==1904

c.samp.dates.12.22.24.p40<- as.character(samp.dates.12.22.24.p40$V1)

Ord.tree.dates.12.22.24.p40 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p40)){
  for(j in 1:length(time.samp.12.22.24.p40)){
    if(tree.tips.12.22.24.p40[i] == c.samp.dates.12.22.24.p40[j]){
      Ord.tree.dates.12.22.24.p40 <- c(Ord.tree.dates.12.22.24.p40, time.samp.12.22.24.p40[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p40[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p40[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p40 <- dater(tree.fasttree.12.22.24.p40,
                                 Ord.tree.dates.12.22.24.p40,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p40, file = paste("dated.tree.A.seed.12.22.24seq.p40.762.Rdata", sep = ""))




# 8.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 30% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.30.571.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.30.571.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p30 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.30.571.fasta.tree", sep = ""))

samp.dates.12.22.24.p30 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.p30$V1 <- as.character(samp.dates.12.22.24.p30$V1)
samp.dates.12.22.24.p30$V2 <- as.numeric(as.character(samp.dates.12.22.24.p30$V2))

time.samp.12.22.24.p30 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p30)

tree.tips.12.22.24.p30 <- (tree.fasttree.12.22.24.p30$tip.label) # length==1904

c.samp.dates.12.22.24.p30<- as.character(samp.dates.12.22.24.p30$V1)

Ord.tree.dates.12.22.24.p30 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p30)){
  for(j in 1:length(time.samp.12.22.24.p30)){
    if(tree.tips.12.22.24.p30[i] == c.samp.dates.12.22.24.p30[j]){
      Ord.tree.dates.12.22.24.p30 <- c(Ord.tree.dates.12.22.24.p30, time.samp.12.22.24.p30[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p30[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p30[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p30 <- dater(tree.fasttree.12.22.24.p30,
                                 Ord.tree.dates.12.22.24.p30,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p30, file = paste("dated.tree.A.seed.12.22.24seq.p30.571.Rdata", sep = ""))



# 9.1 Phylogenetic tree for all seeds 12-22-24 in past 7 years - 20% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.20.381.fasta", sep = ""), paste(">A.Epidemic12.22.24Sequences.gene.pol.7yrs.20.381.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.p20 <- read.tree(paste("A.Epidemic12.22.24Sequences.gene.pol.7yrs.20.381.fasta.tree", sep = ""))

samp.dates.12.22.24.p20 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.p20$V1 <- as.character(samp.dates.12.22.24.p20$V1)
samp.dates.12.22.24.p20$V2 <- as.numeric(as.character(samp.dates.12.22.24.p20$V2))

time.samp.12.22.24.p20 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.p20)

tree.tips.12.22.24.p20 <- (tree.fasttree.12.22.24.p20$tip.label) # length==1904

c.samp.dates.12.22.24.p20<- as.character(samp.dates.12.22.24.p20$V1)

Ord.tree.dates.12.22.24.p20 <- vector() # order the dates according to tips order in the tree
t.tip <- vector()
t.date <- vector()
for(i in 1:length(tree.tips.12.22.24.p20)){
  for(j in 1:length(time.samp.12.22.24.p20)){
    if(tree.tips.12.22.24.p20[i] == c.samp.dates.12.22.24.p20[j]){
      Ord.tree.dates.12.22.24.p20 <- c(Ord.tree.dates.12.22.24.p20, time.samp.12.22.24.p20[j])
      t.tip <- c(t.tip,tree.tips.12.22.24.p20[i])
      t.date <- c(t.date, c.samp.dates.12.22.24.p20[j])
    }
  }
}

# calibrate internal nodes
dater.tree.12.22.24.p20 <- dater(tree.fasttree.12.22.24.p20,
                                 Ord.tree.dates.12.22.24.p20,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.p20, file = paste("dated.tree.A.seed.12.22.24seq.p20.381.Rdata", sep = ""))




# Scenario 3 - A- B- G
#
# different subtypes of the virus (HIV-1-A-B-G) for all seeds
# complete sampling for a transmission network of one seed
# same sampling time interval of seven years

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario3_B/")

## (i) Sequence simulation

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
    seq.rand <- 1
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]

    if(i==12 || i==24){ # subtype G
      # call the seed sequences and rename the file
      file.copy(paste("hiv.seq.G.pol.j.fasta", sep = ""),paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""))
      # add the number of tree in the file and
      write(n.tr,file = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
      # the tree, to prepare the file to simulate the evolution of the virus across the tree
      write.tree(tr,file = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
      file.rename(from = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

      system(paste("./seq-gen -mGTR -f 0.3987, 0.1563, 0.2202, 0.2249 -a 0.9 -i 0.80 -g 4 -r 1.4520, 9.9166, 1.3332, 1.2652, 14.9356, 1.0000 -s 0.00475 -n1 -k 1 <hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >G.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


      # a: shape parameter of Gamma > Gamma Rate Heterogeneity
      # g: category of Gamma > Discrete Gamma Rate Heterogeneity
      # r: rate matrix
      # z: seed

      # Keep sampling dates
      id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

      write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

    }

    else if(i==22){ # subtype A
      # call the seed sequences and rename the file
      file.copy(paste("hiv.seq.A.pol.j.fasta", sep = ""),paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""))
      # add the number of tree in the file and
      write(n.tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
      # the tree, to prepare the file to simulate the evolution of the virus across the tree
      write.tree(tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
      file.rename(from = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

      system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -s 0.00475 -n1  -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


      # a: shape parameter of Gamma > Gamma Rate Heterogeneity
      # g: category of Gamma > Discrete Gamma Rate Heterogeneity
      # r: rate matrix
      # z: seed

      # Keep sampling dates
      id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

      write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

    }
    else{ # subtype B
      # call the seed sequences and rename the file
      file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""))
      # add the number of tree in the file and
      write(n.tr,file = paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
      # the tree, to prepare the file to simulate the evolution of the virus across the tree
      write.tree(tr,file = paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
      file.rename(from = paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.B.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

      system(paste("./seq-gen -mGTR -f 0.3935, 0.1708, 0.2060, 0.2297 -a 0.9 -i 0.80 -g 4 -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475 -n1 -k 1 <hiv.seq.B.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >B.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


      # a: shape parameter of Gamma > Gamma Rate Heterogeneity
      # g: category of Gamma > Discrete Gamma Rate Heterogeneity
      # r: rate matrix
      # z: seed

      # Keep sampling dates
      id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network
      write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))
    }
  }

}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds


#######################
# Scenario 3 Analysis #
#######################


setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/")
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistB.RData"))
simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)
save(simpact.trans.net, file="simpact.trans.net.B.RData")
seeds.great3 <- c(4, 12, 19, 21, 22, 24, 28, 34, 39) # seeds with at least 3 transmission events

# Remain with Seeds which have sampled individuals in past 7 years
seed.7yrs.samp <- vector()
for (i in 1:length(seeds.great3)){
  k <- seeds.great3[i]
  net.i <- as.data.frame(simpact.trans.net[[k]])
  net.filtered.i <- filter(net.i, dtimes<=7)
  dtimes <- net.i$dtimes<=7
  dtimes.elements <- which(dtimes)
  if(length(dtimes.elements)!=0){
    seed.7yrs.samp <- c(seed.7yrs.samp, k)
    save(net.filtered.i, file=paste("filtered.trans.net.seed",k,".Rdata", sep = ""))
  }
}

net.seed12 = get(load("filtered.trans.net.seed12.Rdata")) # length(net.seed12$id) = 314, subtype G
net.seed22 = get(load("filtered.trans.net.seed22.Rdata")) # length(net.seed22$id) = 1579, subtype A
net.seed24 = get(load("filtered.trans.net.seed24.Rdata")) # length(net.seed24$id) = 11, subtype G

# IDs to exclude in the sequences
id.12 <- setdiff(simpact.trans.net[[12]]$id, net.seed12$id) # Done! - 765 Ok
id.22 <- setdiff(simpact.trans.net[[22]]$id, net.seed22$id) # Done! - 1453 Ok
id.24 <- setdiff(simpact.trans.net[[24]]$id, net.seed24$id) # Done! - 798 Ok


# # rename sampling times automatically and make a single file for sampling dates of seeds 12, 22, 24
#
#
# transfSampT.12 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.csv")
# transfSampT.12$V1 <- paste("A.12.",as.character(transfSampT.12$V1), sep = "")
# write.csv(transfSampT.12, file = "Tsamplingtimes_seed_number_12.csv", quote = TRUE)
#
#
# transfSampT.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.csv")
# transfSampT.22$V1 <- paste("A.22.",as.character(transfSampT.12$V1), sep = "")
# write.csv(transfSampT.22, file = "Tsamplingtimes_seed_number_22.csv", quote = TRUE)
#
# transfSampT.24 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_24.csv")
# transfSampT.24$V1 <- paste("A.24.",as.character(transfSampT.12$V1), sep = "")
# write.csv(transfSampT.24, file = "Tsamplingtimes_seed_number_24.csv", quote = TRUE)
#
#
# transfSampTimes.12 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_12.csv")
# transfSampTimes.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_22.csv")
# transfSampTimes.24 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_24.csv")
#
# d <- rbind(transfSampTimes.12[,3:4], transfSampTimes.22[,3:4], transfSampTimes.24[,3:4])
# dd <- as.data.frame(d)
# dd$V1 <- as.character(dd$V1)
# dd$V2 <- as.numeric(as.character(dd$V2)) # solve factors issue
#
# write.csv(dd, file = "samplingtimes_seed_number_22.24.27.csv")


# 1.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 100% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.314.1579.11.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.314.1579.11.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.314.1579.11.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G$V1 <- as.character(samp.dates.12.22.24.G.A.G$V1)
samp.dates.12.22.24.G.A.G$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G$V2))

time.samp.12.22.24.G.A.G <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G)

tree.tips.12.22.24.G.A.G <- (tree.fasttree.12.22.24.G.A.G$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G<- as.character(samp.dates.12.22.24.G.A.G$V1)

Ord.tree.dates.12.22.24.G.A.G <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G)){
  for(j in 1:length(time.samp.12.22.24.G.A.G)){
    if(tree.tips.12.22.24.G.A.G[i] == c.samp.dates.12.22.24.G.A.G[j]){
      Ord.tree.dates.12.22.24.G.A.G <- c(Ord.tree.dates.12.22.24.G.A.G, time.samp.12.22.24.G.A.G[j])
      t <- c(t, tree.tips.12.22.24.G.A.G[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G <- dater(tree.fasttree.12.22.24.G.A.G,
                                   Ord.tree.dates.12.22.24.G.A.G,
                                   s = 3012,
                                   omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p100.1904.Rdata", sep = ""))

tr.scenario3B <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p100.1904.Rdata"))
write.tree(tr.scenario3B, file = "tr.scenario3B.12.22.24.1904.HIVGAG.tree")

# 2.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 90% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.90.1714.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.90.1714.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p90 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.90.1714.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p90  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p90$V1 <- as.character(samp.dates.12.22.24.G.A.G.p90$V1)
samp.dates.12.22.24.G.A.G.p90$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p90$V2))

time.samp.12.22.24.G.A.G.p90 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p90)

tree.tips.12.22.24.G.A.G.p90 <- (tree.fasttree.12.22.24.G.A.G.p90$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p90<- as.character(samp.dates.12.22.24.G.A.G.p90$V1)

Ord.tree.dates.12.22.24.G.A.G.p90 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p90)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p90)){
    if(tree.tips.12.22.24.G.A.G.p90[i] == c.samp.dates.12.22.24.G.A.G.p90[j]){
      Ord.tree.dates.12.22.24.G.A.G.p90 <- c(Ord.tree.dates.12.22.24.G.A.G.p90, time.samp.12.22.24.G.A.G.p90[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p90[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p90 <- dater(tree.fasttree.12.22.24.G.A.G.p90,
                                       Ord.tree.dates.12.22.24.G.A.G.p90,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p90, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p90.1714.Rdata", sep = ""))




# 3.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 80% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.80.1523.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.80.1523.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p80 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.80.1523.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p80  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p80$V1 <- as.character(samp.dates.12.22.24.G.A.G.p80$V1)
samp.dates.12.22.24.G.A.G.p80$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p80$V2))

time.samp.12.22.24.G.A.G.p80 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p80)

tree.tips.12.22.24.G.A.G.p80 <- (tree.fasttree.12.22.24.G.A.G.p80$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p80<- as.character(samp.dates.12.22.24.G.A.G.p80$V1)

Ord.tree.dates.12.22.24.G.A.G.p80 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p80)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p80)){
    if(tree.tips.12.22.24.G.A.G.p80[i] == c.samp.dates.12.22.24.G.A.G.p80[j]){
      Ord.tree.dates.12.22.24.G.A.G.p80 <- c(Ord.tree.dates.12.22.24.G.A.G.p80, time.samp.12.22.24.G.A.G.p80[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p80[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p80 <- dater(tree.fasttree.12.22.24.G.A.G.p80,
                                       Ord.tree.dates.12.22.24.G.A.G.p80,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p80, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p80.1523.Rdata", sep = ""))


# 4.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 70% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.70.1333.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.70.1333.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p70 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.70.1333.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p70  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p70$V1 <- as.character(samp.dates.12.22.24.G.A.G.p70$V1)
samp.dates.12.22.24.G.A.G.p70$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p70$V2))

time.samp.12.22.24.G.A.G.p70 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p70)

tree.tips.12.22.24.G.A.G.p70 <- (tree.fasttree.12.22.24.G.A.G.p70$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p70<- as.character(samp.dates.12.22.24.G.A.G.p70$V1)

Ord.tree.dates.12.22.24.G.A.G.p70 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p70)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p70)){
    if(tree.tips.12.22.24.G.A.G.p70[i] == c.samp.dates.12.22.24.G.A.G.p70[j]){
      Ord.tree.dates.12.22.24.G.A.G.p70 <- c(Ord.tree.dates.12.22.24.G.A.G.p70, time.samp.12.22.24.G.A.G.p70[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p70[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p70 <- dater(tree.fasttree.12.22.24.G.A.G.p70,
                                       Ord.tree.dates.12.22.24.G.A.G.p70,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p70, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p70.1333.Rdata", sep = ""))


# 5.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 60% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.60.1142.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.60.1142.fasta.tree", sep = "")))
# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p60 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.70.1333.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p60  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p60$V1 <- as.character(samp.dates.12.22.24.G.A.G.p60$V1)
samp.dates.12.22.24.G.A.G.p60$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p60$V2))

time.samp.12.22.24.G.A.G.p60 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p60)

tree.tips.12.22.24.G.A.G.p60 <- (tree.fasttree.12.22.24.G.A.G.p60$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p60<- as.character(samp.dates.12.22.24.G.A.G.p60$V1)

Ord.tree.dates.12.22.24.G.A.G.p60 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p60)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p60)){
    if(tree.tips.12.22.24.G.A.G.p60[i] == c.samp.dates.12.22.24.G.A.G.p60[j]){
      Ord.tree.dates.12.22.24.G.A.G.p60 <- c(Ord.tree.dates.12.22.24.G.A.G.p60, time.samp.12.22.24.G.A.G.p60[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p60[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p60 <- dater(tree.fasttree.12.22.24.G.A.G.p60,
                                       Ord.tree.dates.12.22.24.G.A.G.p60,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p60, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p60.1142.Rdata", sep = ""))



# 6.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 50% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.50.952.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.50.952.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p50 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.50.952.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p50  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p50$V1 <- as.character(samp.dates.12.22.24.G.A.G.p50$V1)
samp.dates.12.22.24.G.A.G.p50$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p50$V2))

time.samp.12.22.24.G.A.G.p50 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p50)

tree.tips.12.22.24.G.A.G.p50 <- (tree.fasttree.12.22.24.G.A.G.p50$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p50<- as.character(samp.dates.12.22.24.G.A.G.p50$V1)

Ord.tree.dates.12.22.24.G.A.G.p50 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p50)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p50)){
    if(tree.tips.12.22.24.G.A.G.p50[i] == c.samp.dates.12.22.24.G.A.G.p50[j]){
      Ord.tree.dates.12.22.24.G.A.G.p50 <- c(Ord.tree.dates.12.22.24.G.A.G.p50, time.samp.12.22.24.G.A.G.p50[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p50[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p50 <- dater(tree.fasttree.12.22.24.G.A.G.p50,
                                       Ord.tree.dates.12.22.24.G.A.G.p50,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p50, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p50.952.Rdata", sep = ""))


# 7.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 40% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.40.762.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.40.762.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p40 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.40.762.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p40  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p40$V1 <- as.character(samp.dates.12.22.24.G.A.G.p40$V1)
samp.dates.12.22.24.G.A.G.p40$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p40$V2))

time.samp.12.22.24.G.A.G.p40 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p40)

tree.tips.12.22.24.G.A.G.p40 <- (tree.fasttree.12.22.24.G.A.G.p40$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p40<- as.character(samp.dates.12.22.24.G.A.G.p40$V1)

Ord.tree.dates.12.22.24.G.A.G.p40 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p40)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p40)){
    if(tree.tips.12.22.24.G.A.G.p40[i] == c.samp.dates.12.22.24.G.A.G.p40[j]){
      Ord.tree.dates.12.22.24.G.A.G.p40 <- c(Ord.tree.dates.12.22.24.G.A.G.p40, time.samp.12.22.24.G.A.G.p40[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p40[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p40 <- dater(tree.fasttree.12.22.24.G.A.G.p40,
                                       Ord.tree.dates.12.22.24.G.A.G.p40,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p40, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p40.762.Rdata", sep = ""))



# 8.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 30% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.30.571.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.30.571.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p30 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.30.571.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p30  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p30$V1 <- as.character(samp.dates.12.22.24.G.A.G.p30$V1)
samp.dates.12.22.24.G.A.G.p30$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p30$V2))

time.samp.12.22.24.G.A.G.p30 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p30)

tree.tips.12.22.24.G.A.G.p30 <- (tree.fasttree.12.22.24.G.A.G.p30$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p30<- as.character(samp.dates.12.22.24.G.A.G.p30$V1)

Ord.tree.dates.12.22.24.G.A.G.p30 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p30)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p30)){
    if(tree.tips.12.22.24.G.A.G.p30[i] == c.samp.dates.12.22.24.G.A.G.p30[j]){
      Ord.tree.dates.12.22.24.G.A.G.p30 <- c(Ord.tree.dates.12.22.24.G.A.G.p30, time.samp.12.22.24.G.A.G.p30[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p30[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p30 <- dater(tree.fasttree.12.22.24.G.A.G.p30,
                                       Ord.tree.dates.12.22.24.G.A.G.p30,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p30, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p30.571.Rdata", sep = ""))




# 8.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 20% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.20.381.fasta", sep = ""), paste(">G.A.GEpidemic12.22.24.Sequences.gene.pol.p.20.381.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.12.22.24.G.A.G.p20 <- read.tree(paste("G.A.GEpidemic12.22.24.Sequences.gene.pol.p.20.381.fasta.tree", sep = ""))

samp.dates.12.22.24.G.A.G.p20  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_12.22.24.csv", dec = ",")

samp.dates.12.22.24.G.A.G.p20$V1 <- as.character(samp.dates.12.22.24.G.A.G.p20$V1)
samp.dates.12.22.24.G.A.G.p20$V2 <- as.numeric(as.character(samp.dates.12.22.24.G.A.G.p20$V2))

time.samp.12.22.24.G.A.G.p20 <- dates.Transform.NamedVector(dates=samp.dates.12.22.24.G.A.G.p20)

tree.tips.12.22.24.G.A.G.p20 <- (tree.fasttree.12.22.24.G.A.G.p20$tip.label) # length==1098

c.samp.dates.12.22.24.G.A.G.p20<- as.character(samp.dates.12.22.24.G.A.G.p20$V1)

Ord.tree.dates.12.22.24.G.A.G.p20 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.12.22.24.G.A.G.p20)){
  for(j in 1:length(time.samp.12.22.24.G.A.G.p20)){
    if(tree.tips.12.22.24.G.A.G.p20[i] == c.samp.dates.12.22.24.G.A.G.p20[j]){
      Ord.tree.dates.12.22.24.G.A.G.p20 <- c(Ord.tree.dates.12.22.24.G.A.G.p20, time.samp.12.22.24.G.A.G.p20[j])
      t <- c(t, tree.tips.12.22.24.G.A.G.p20[i])
    }
  }
}


# calibrate internal nodes
dater.tree.12.22.24.G.A.G.p20 <- dater(tree.fasttree.12.22.24.G.A.G.p20,
                                       Ord.tree.dates.12.22.24.G.A.G.p20,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.12.22.24.G.A.G.p20, file = paste("dated.tree.G.A.G.seed.12.22.24seq.p20.381.Rdata", sep = ""))




# Scenario 4 - A- B- G
#
# different subtypes of the virus (HIV-1-A-B-G) for all seeds
# complete sampling in pasy seven years, and study the effect of gender and age on sampling

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario4_B/")


## (i) Sequence simulation

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
    seq.rand <- 1
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]

    if(i==12 || i==24){ # subtype G
      # call the seed sequences and rename the file
      file.copy(paste("hiv.seq.G.pol.j.fasta", sep = ""),paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""))
      # add the number of tree in the file and
      write(n.tr,file = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
      # the tree, to prepare the file to simulate the evolution of the virus across the tree
      write.tree(tr,file = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
      file.rename(from = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

      system(paste("./seq-gen -mGTR -f 0.3987, 0.1563, 0.2202, 0.2249 -a 0.9 -i 0.80 -g 4 -r 1.4520, 9.9166, 1.3332, 1.2652, 14.9356, 1.0000 -s 0.00475 -n1 -k 1 <hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >G.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


      # a: shape parameter of Gamma > Gamma Rate Heterogeneity
      # g: category of Gamma > Discrete Gamma Rate Heterogeneity
      # r: rate matrix
      # z: seed

      # Keep sampling dates
      id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

      write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

    }

    else if(i==22){ # subtype A
      # call the seed sequences and rename the file
      file.copy(paste("hiv.seq.A.pol.j.fasta", sep = ""),paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""))
      # add the number of tree in the file and
      write(n.tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
      # the tree, to prepare the file to simulate the evolution of the virus across the tree
      write.tree(tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
      file.rename(from = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

      system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -s 0.00475 -n1  -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


      # a: shape parameter of Gamma > Gamma Rate Heterogeneity
      # g: category of Gamma > Discrete Gamma Rate Heterogeneity
      # r: rate matrix
      # z: seed

      # Keep sampling dates
      id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

      write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

    }
    else{ # subtype B
      # call the seed sequences and rename the file
      file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""))
      # add the number of tree in the file and
      write(n.tr,file = paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
      # the tree, to prepare the file to simulate the evolution of the virus across the tree
      write.tree(tr,file = paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
      file.rename(from = paste("hiv.seq.B.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.B.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

      system(paste("./seq-gen -mGTR -f 0.3935, 0.1708, 0.2060, 0.2297 -a 0.9 -i 0.80 -g 4 -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475 -n1 -k 1 <hiv.seq.B.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >B.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


      # a: shape parameter of Gamma > Gamma Rate Heterogeneity
      # g: category of Gamma > Discrete Gamma Rate Heterogeneity
      # r: rate matrix
      # z: seed

      # Keep sampling dates
      id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network
      write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))
    }
  }

}
# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds


#######################
# Scenario 4 Analysis #
#######################


setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis4_B")
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistB.RData"))
simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)
save(simpact.trans.net, file="simpact.trans.net.B.RData")

seeds.great3 <- c(4, 12, 19, 21, 22, 24, 28, 34, 39) # seeds with at least 3 transmission events

# Remain with Seeds which have sampled individuals in past 7 years
seed.7yrs.samp <- vector()
for (i in 1:length(seeds.great3)){
  k <- seeds.great3[i]
  net.i <- as.data.frame(simpact.trans.net[[k]])
  net.filtered.i <- filter(net.i, dtimes<=7)
  dtimes <- net.i$dtimes<=7
  dtimes.elements <- which(dtimes)
  if(length(dtimes.elements)!=0){
    seed.7yrs.samp <- c(seed.7yrs.samp, k)
    save(net.filtered.i, file=paste("filtered.trans.net.seed",k,".Rdata", sep = ""))
  }
}

net.seed12 = as.data.frame(get(load("filtered.trans.net.seed12.Rdata"))) # length(net.seed12$id) = 314, subtype G
net.seed22 = as.data.frame(get(load("filtered.trans.net.seed22.Rdata"))) # length(net.seed22$id) = 1579, subtype A
net.seed24 = as.data.frame(get(load("filtered.trans.net.seed24.Rdata"))) # length(net.seed24$id) = 11, subtype G

# IDs to exclude in the sequences
id.12 <- setdiff(simpact.trans.net[[12]]$id, net.seed12$id) # Done! - 765 Ok
id.22 <- setdiff(simpact.trans.net[[22]]$id, net.seed22$id) # Done! - 1453 Ok
id.24 <- setdiff(simpact.trans.net[[24]]$id, net.seed24$id) # Done! - 798 Ok

## Let say we want to have sample 70% among these sampled in past year, and we want
## (i) Gender: 60% of them must be female (1) individuals, 40% males (0)
## (ii) Age: 90% of them must be below 35 years  old
## (iii) Gender + Age: 60% women where 90% being under 35 years

net.seed12$TOB <- (40 - net.seed12$TOB) # G
net.seed22$TOB <- (40 - net.seed22$TOB) # A
net.seed24$TOB <- (40 - net.seed24$TOB) # G

net.seed12$id <- paste("G.12.", net.seed12$id, sep = "")
net.seed22$id <- paste("A.22.", net.seed22$id, sep = "")
net.seed24$id <- paste("G.24.", net.seed24$id, sep = "")

# (i) 70% of 1904 = 1333, and 60% being women (Gender=1) > 800, 40% being men (Gender=0) > 533

net.seeds.12.22.24 <- as.data.table(rbind(net.seed12, net.seed22, net.seed24))

net.seeds.12.22.24.W.60 <- net.seeds.12.22.24[net.seeds.12.22.24$Gender==1]
net.seeds.12.22.24.M.40 <- net.seeds.12.22.24[net.seeds.12.22.24$Gender==0]

samp.women.60 <- sample(net.seeds.12.22.24.W.60$id, 800)
samp.men.40 <- sample(net.seeds.12.22.24.M.40$id, 533)

samp.women.60.men.40 <- c(samp.women.60, samp.men.40) # to be used for tree

remove.seq.i <- setdiff(net.seeds.12.22.24$id, samp.women.60.men.40) # 571

# (ii)  70% of 1904 = 1333, and 90% (1200) being below 35 years old

net.seeds.12.22.24.below.35yrs <- net.seeds.12.22.24[net.seeds.12.22.24$TOB <= 35] # 1675
id.net.seeds.12.22.24.remain <- setdiff(net.seeds.12.22.24$id, net.seeds.12.22.24.below.35yrs$id) # 229

samp.below.35yrs <- sample(net.seeds.12.22.24.below.35yrs$id, 1200) # r
samp.above.or.less.35yrs <- sample(id.net.seeds.12.22.24.remain, 133) #

samp.80.below.35yrs <- c(samp.below.35yrs, samp.above.or.less.35yrs) # to be used for tree

remove.seq.ii <- setdiff(net.seeds.12.22.24$id, samp.80.below.35yrs) # 571


# (iii) 70% of 1904 = 1333, and 50% being women (Gender=1) > 666, 50% being men (Gender=0) > 667
# and where 90% of women being being under 35 years (600), 67 women being of random age
net.seeds.12.22.24.W <- net.seeds.12.22.24[net.seeds.12.22.24$Gender==1] # 836
net.seeds.12.22.24.M <- net.seeds.12.22.24[net.seeds.12.22.24$Gender==0]

samp.men.50 <- sample(net.seeds.12.22.24.M$id, 667)

net.seeds.12.22.24.W.50.90.below.35yrs <- net.seeds.12.22.24.W[net.seeds.12.22.24.W$TOB <= 35] # 770

samp.50.90.below.35yrs <- sample(net.seeds.12.22.24.W.50.90.below.35yrs$id, 600)

dif.samp.women <- setdiff(net.seeds.12.22.24.W$id, net.seeds.12.22.24.W.50.90.below.35yrs$id)

samp.50.10.above.or.less.35yrs <-sample(dif.samp.women, 66)

samp.50.50.W.M.90W.below.35yrs <- c(samp.men.50, samp.60.90.below.35yrs, samp.50.10.above.or.less.35yrs) # to be used for tree

remove.seq.iii <- setdiff(net.seeds.12.22.24$id, samp.50.50.W.M.90W.below.35yrs) # 571


# > remove.seq.i
# [8] "G.12.658"  "G.12.700"  "G.12.727"  "G.12.747"  "G.12.755"  "G.12.759"  "G.12.766"
# [15] "G.12.769"  "G.12.770"  "G.12.772"  "G.12.775"  "G.12.789"  "G.12.792"  "G.12.806"
# [22] "G.12.808"  "G.12.813"  "G.12.817"  "G.12.825"  "G.12.828"  "G.12.830"  "G.12.834"
# [29] "G.12.841"  "G.12.846"  "G.12.848"  "G.12.851"  "G.12.856"  "G.12.862"  "G.12.864"
# [36] "G.12.868"  "G.12.873"  "G.12.880"  "G.12.882"  "G.12.891"  "G.12.893"  "G.12.894"
# [43] "G.12.899"  "G.12.901"  "G.12.912"  "G.12.918"  "G.12.919"  "G.12.922"  "G.12.924"
# [50] "G.12.931"  "G.12.932"  "G.12.933"  "G.12.934"  "G.12.941"  "G.12.947"  "G.12.952"
# [57] "G.12.953"  "G.12.955"  "G.12.959"  "G.12.960"  "G.12.966"  "G.12.968"  "G.12.970"
# [64] "G.12.974"  "G.12.975"  "G.12.976"  "G.12.980"  "G.12.981"  "G.12.990"  "G.12.991"
# [71] "G.12.993"  "G.12.995"  "G.12.997"  "G.12.999"  "G.12.1003" "G.12.1008" "G.12.1009"
# [78] "G.12.1013" "G.12.1025" "G.12.1029" "G.12.1031" "G.12.1033" "G.12.1035" "G.12.1051"
# [85] "G.12.1063" "G.12.1064" "G.12.1065" "G.12.1068" "G.12.1071" "G.12.1072" "A.22.682"
# [92] "A.22.756"  "A.22.814"  "A.22.888"  "A.22.912"  "A.22.956"  "A.22.966"  "A.22.975"
# [99] "A.22.981"  "A.22.1026" "A.22.1046" "A.22.1064" "A.22.1139" "A.22.1148" "A.22.1167"
# [106] "A.22.1182" "A.22.1193" "A.22.1205" "A.22.1219" "A.22.1225" "A.22.1228" "A.22.1244"
# [113] "A.22.1249" "A.22.1266" "A.22.1268" "A.22.1276" "A.22.1298" "A.22.1302" "A.22.1312"
# [120] "A.22.1324" "A.22.1347" "A.22.1367" "A.22.1368" "A.22.1374" "A.22.1388" "A.22.1391"
# [127] "A.22.1392" "A.22.1396" "A.22.1397" "A.22.1398" "A.22.1401" "A.22.1407" "A.22.1410"
# [134] "A.22.1416" "A.22.1423" "A.22.1442" "A.22.1467" "A.22.1468" "A.22.1477" "A.22.1481"
# [141] "A.22.1485" "A.22.1487" "A.22.1489" "A.22.1493" "A.22.1501" "A.22.1503" "A.22.1505"
# [148] "A.22.1512" "A.22.1538" "A.22.1555" "A.22.1566" "A.22.1567" "A.22.1585" "A.22.1591"
# [155] "A.22.1597" "A.22.1606" "A.22.1607" "A.22.1620" "A.22.1622" "A.22.1626" "A.22.1628"
# [162] "A.22.1633" "A.22.1647" "A.22.1649" "A.22.1661" "A.22.1666" "A.22.1673" "A.22.1684"
# [169] "A.22.1685" "A.22.1686" "A.22.1689" "A.22.1694" "A.22.1703" "A.22.1708" "A.22.1714"
# [176] "A.22.1723" "A.22.1733" "A.22.1740" "A.22.1748" "A.22.1751" "A.22.1755" "A.22.1757"
# [183] "A.22.1769" "A.22.1772" "A.22.1773" "A.22.1776" "A.22.1779" "A.22.1782" "A.22.1783"
# [190] "A.22.1787" "A.22.1802" "A.22.1803" "A.22.1813" "A.22.1821" "A.22.1822" "A.22.1824"
# [197] "A.22.1829" "A.22.1831" "A.22.1852" "A.22.1854" "A.22.1857" "A.22.1860" "A.22.1861"
# [204] "A.22.1864" "A.22.1867" "A.22.1869" "A.22.1873" "A.22.1876" "A.22.1877" "A.22.1882"
# [211] "A.22.1893" "A.22.1897" "A.22.1900" "A.22.1903" "A.22.1907" "A.22.1910" "A.22.1911"
# [218] "A.22.1912" "A.22.1913" "A.22.1914" "A.22.1917" "A.22.1918" "A.22.1923" "A.22.1926"
# [225] "A.22.1932" "A.22.1933" "A.22.1935" "A.22.1939" "A.22.1943" "A.22.1945" "A.22.1946"
# [232] "A.22.1948" "A.22.1950" "A.22.1952" "A.22.1954" "A.22.1956" "A.22.1962" "A.22.1968"
# [239] "A.22.1970" "A.22.1972" "A.22.1974" "A.22.1979" "A.22.1987" "A.22.1989" "A.22.1991"
# [246] "A.22.1992" "A.22.1995" "A.22.1996" "A.22.1998" "A.22.1999" "A.22.2001" "A.22.2011"
# [253] "A.22.2012" "A.22.2014" "A.22.2017" "A.22.2022" "A.22.2024" "A.22.2027" "A.22.2030"
# [260] "A.22.2034" "A.22.2036" "A.22.2041" "A.22.2054" "A.22.2059" "A.22.2063" "A.22.2066"
# [267] "A.22.2067" "A.22.2068" "A.22.2074" "A.22.2077" "A.22.2078" "A.22.2079" "A.22.2080"
# [274] "A.22.2081" "A.22.2093" "A.22.2101" "A.22.2104" "A.22.2108" "A.22.2109" "A.22.2110"
# [281] "A.22.2111" "A.22.2115" "A.22.2122" "A.22.2129" "A.22.2131" "A.22.2133" "A.22.2134"
# [288] "A.22.2135" "A.22.2136" "A.22.2140" "A.22.2142" "A.22.2143" "A.22.2144" "A.22.2148"
# [295] "A.22.2150" "A.22.2153" "A.22.2154" "A.22.2156" "A.22.2161" "A.22.2168" "A.22.2171"
# [302] "A.22.2173" "A.22.2175" "A.22.2179" "A.22.2180" "A.22.2186" "A.22.2187" "A.22.2189"
# [309] "A.22.2191" "A.22.2193" "A.22.2201" "A.22.2204" "A.22.2208" "A.22.2210" "A.22.2211"
# [316] "A.22.2215" "A.22.2223" "A.22.2224" "A.22.2225" "A.22.2229" "A.22.2230" "A.22.2233"
# [323] "A.22.2242" "A.22.2243" "A.22.2246" "A.22.2248" "A.22.2252" "A.22.2256" "A.22.2257"
# [330] "A.22.2260" "A.22.2261" "A.22.2269" "A.22.2270" "A.22.2273" "A.22.2277" "A.22.2278"
# [337] "A.22.2280" "A.22.2282" "A.22.2285" "A.22.2289" "A.22.2291" "A.22.2292" "A.22.2293"
# [344] "A.22.2295" "A.22.2305" "A.22.2306" "A.22.2307" "A.22.2308" "A.22.2323" "A.22.2324"
# [351] "A.22.2327" "A.22.2328" "A.22.2333" "A.22.2334" "A.22.2340" "A.22.2342" "A.22.2346"
# [358] "A.22.2347" "A.22.2351" "A.22.2352" "A.22.2357" "A.22.2359" "A.22.2360" "A.22.2362"
# [365] "A.22.2364" "A.22.2369" "A.22.2380" "A.22.2382" "A.22.2386" "A.22.2387" "A.22.2389"
# [372] "A.22.2394" "A.22.2395" "A.22.2398" "A.22.2406" "A.22.2407" "A.22.2408" "A.22.2416"
# [379] "A.22.2419" "A.22.2424" "A.22.2426" "A.22.2428" "A.22.2429" "A.22.2432" "A.22.2442"
# [386] "A.22.2447" "A.22.2448" "A.22.2449" "A.22.2451" "A.22.2454" "A.22.2455" "A.22.2457"
# [393] "A.22.2462" "A.22.2464" "A.22.2465" "A.22.2466" "A.22.2468" "A.22.2469" "A.22.2471"
# [400] "A.22.2473" "A.22.2476" "A.22.2478" "A.22.2480" "A.22.2481" "A.22.2483" "A.22.2485"
# [407] "A.22.2486" "A.22.2493" "A.22.2494" "A.22.2497" "A.22.2498" "A.22.2506" "A.22.2507"
# [414] "A.22.2508" "A.22.2515" "A.22.2522" "A.22.2527" "A.22.2528" "A.22.2529" "A.22.2533"
# [421] "A.22.2534" "A.22.2540" "A.22.2541" "A.22.2542" "A.22.2545" "A.22.2546" "A.22.2550"
# [428] "A.22.2555" "A.22.2557" "A.22.2562" "A.22.2563" "A.22.2568" "A.22.2572" "A.22.2574"
# [435] "A.22.2575" "A.22.2582" "A.22.2585" "A.22.2586" "A.22.2588" "A.22.2591" "A.22.2593"
# [442] "A.22.2594" "A.22.2595" "A.22.2597" "A.22.2598" "A.22.2610" "A.22.2617" "A.22.2622"
# [449] "A.22.2623" "A.22.2625" "A.22.2628" "A.22.2629" "A.22.2641" "A.22.2642" "A.22.2644"
# [456] "A.22.2649" "A.22.2652" "A.22.2654" "A.22.2655" "A.22.2666" "A.22.2667" "A.22.2668"
# [463] "A.22.2669" "A.22.2672" "A.22.2673" "A.22.2681" "A.22.2687" "A.22.2691" "A.22.2692"
# [470] "A.22.2693" "A.22.2695" "A.22.2697" "A.22.2698" "A.22.2705" "A.22.2706" "A.22.2709"
# [477] "A.22.2712" "A.22.2714" "A.22.2726" "A.22.2729" "A.22.2732" "A.22.2736" "A.22.2741"
# [484] "A.22.2742" "A.22.2744" "A.22.2746" "A.22.2747" "A.22.2749" "A.22.2751" "A.22.2759"
# [491] "A.22.2760" "A.22.2764" "A.22.2769" "A.22.2780" "A.22.2783" "A.22.2785" "A.22.2787"
# [498] "A.22.2788" "A.22.2789" "A.22.2796" "A.22.2799" "A.22.2802" "A.22.2810" "A.22.2816"
# [505] "A.22.2821" "A.22.2822" "A.22.2823" "A.22.2828" "A.22.2829" "A.22.2833" "A.22.2834"
# [512] "A.22.2836" "A.22.2839" "A.22.2854" "A.22.2859" "A.22.2860" "A.22.2862" "A.22.2864"
# [519] "A.22.2865" "A.22.2867" "A.22.2869" "A.22.2873" "A.22.2876" "A.22.2879" "A.22.2885"
# [526] "A.22.2888" "A.22.2889" "A.22.2892" "A.22.2895" "A.22.2898" "A.22.2899" "A.22.2900"
# [533] "A.22.2904" "A.22.2905" "A.22.2907" "A.22.2910" "A.22.2911" "A.22.2914" "A.22.2916"
# [540] "A.22.2918" "A.22.2919" "A.22.2921" "A.22.2922" "A.22.2923" "A.22.2924" "A.22.2931"
# [547] "A.22.2933" "A.22.2937" "A.22.2943" "A.22.2945" "A.22.2948" "A.22.2950" "A.22.2956"
# [554] "A.22.2965" "A.22.2969" "A.22.2971" "A.22.2984" "A.22.2985" "A.22.2986" "A.22.2990"
# [561] "A.22.2999" "A.22.3001" "A.22.3002" "A.22.3012" "A.22.3018" "A.22.3022" "A.22.3024"
# [568] "A.22.3027" "G.24.706"  "G.24.787"  "G.24.795"
# >

