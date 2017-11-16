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
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistA.RData"))


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

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario1_A/")



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

setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/Seed_22/")
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistA.RData"))


# 1.1. Full phylogenetic tree for seed 22 - 100% of sequences

# Sampling dates in calender time
dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1977+40-as.numeric(dates$V2) # dates datalist$itable$population.simtime[1] - dates$V2 + 1977
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.100.1403.fasta", sep = ""), paste(">A.seed.22.seq.100.1403.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22 <- read.tree(paste("A.seed.22.seq.100.1403.fasta.tree", sep = ""))

samp.dates.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22, file = paste("dated.tree.A.seed.22.seq.100.1403.Rdata", sep = ""))

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
save(properties.trans.net.22, file = paste("trans.net.A.seed.22.seq.100.1403.Rdata", sep = ""))


# 2.1. Phylogenetic tree for seed 22 - 90% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.90.1263.fasta", sep = ""), paste(">A.seed.22.seq.90.1263.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.90 <- read.tree(paste("A.seed.22.seq.90.1263.fasta.tree", sep = ""))

samp.dates.22.90 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22.90, file = paste("dated.tree.A.seed.22.seq.90.1263.Rdata", sep = ""))




# 3.1. Phylogenetic tree for seed 22 - 80% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.80.1122.fasta", sep = ""), paste(">A.seed.22.seq.80.1122.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.80 <- read.tree(paste("A.seed.22.seq.80.1122.fasta.tree", sep = ""))

samp.dates.22.80 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22.80, file = paste("dated.tree.A.seed.22.seq.80.1122.Rdata", sep = ""))



# 4.1. Phylogenetic tree for seed 22 - 70% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.70.982.fasta", sep = ""), paste(">A.seed.22.seq.70.982.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.70 <- read.tree(paste("A.seed.22.seq.70.982.fasta.tree", sep = ""))

samp.dates.22.70 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22.70, file = paste("dated.tree.A.seed.22.seq.70.982.Rdata", sep = ""))



# 5.1. Phylogenetic tree for seed 22 - 60% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.60.842.fasta", sep = ""), paste(">A.seed.22.seq.60.842.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.60 <- read.tree(paste("A.seed.22.seq.60.842.fasta.tree", sep = ""))

samp.dates.22.60 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22.60, file = paste("dated.tree.A.seed.22.seq.60.842.Rdata", sep = ""))


# 6.1. Phylogenetic tree for seed 22 - 50% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.50.701.fasta", sep = ""), paste(">A.seed.22.seq.50.701.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.50 <- read.tree(paste("A.seed.22.seq.50.701.fasta.tree", sep = ""))

samp.dates.22.50 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22.50, file = paste("dated.tree.A.seed.22.seq.50.701.Rdata", sep = ""))



# 7.1. Phylogenetic tree for seed 22 - 40% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.40.561.fasta", sep = ""), paste(">A.seed.22.seq.40.561.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.40 <- read.tree(paste("A.seed.22.seq.40.561.fasta.tree", sep = ""))

samp.dates.22.40 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22.40, file = paste("dated.tree.A.seed.22.seq.40.561.Rdata", sep = ""))




# 8.1. Phylogenetic tree for seed 22 - 30% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.30.421.fasta", sep = ""), paste(">A.seed.22.seq.30.421.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.30 <- read.tree(paste("A.seed.22.seq.30.421.fasta.tree", sep = ""))

samp.dates.22.30 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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

save(dater.tree.22.30, file = paste("dated.tree.A.seed.22.seq.30.421.Rdata", sep = ""))




# 9.1. Phylogenetic tree for seed 22 - 20% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.seed.22.seq.20.281.fasta", sep = ""), paste(">A.seed.22.seq.20.281.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.20 <- read.tree(paste("A.seed.22.seq.20.281.fasta.tree", sep = ""))

samp.dates.22.20 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_A/samplingtimes_seed_number_22.csv")

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
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.20, file = paste("dated.tree.A.seed.22.seq.20.281.Rdata", sep = ""))



# Scenario 2 - A
#
# one subtype of the virus (HIV-1-A) for all seeds
# complete sampling for a transmission network of all seeds
# same sampling time interval of 7 years for a transmission network of all seeds

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario2_A/")


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


setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/")
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistA.RData"))
simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)

seeds.great3 <- c(10, 18, 19, 21, 22, 24, 27, 32, 37) # seeds with at least 3 transmission events

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

net.seed22 = get(load("filtered.trans.net.seed22.Rdata")) # length(net.seed22$id) = 895
net.seed24 = get(load("filtered.trans.net.seed24.Rdata")) # length(net.seed24$id) = 24
net.seed27 = get(load("filtered.trans.net.seed27.Rdata")) # length(net.seed27$id) = 179

# IDs to exclude in the sequences
id.22 <- setdiff(simpact.trans.net[[22]]$id, net.seed22$id) # Done! - 508
id.24 <- setdiff(simpact.trans.net[[24]]$id, net.seed24$id) # Done! - 692
id.27 <- setdiff(simpact.trans.net[[27]]$id, net.seed27$id) # Done! - 191


# Tip to see forgotten IDs when deleting
# length(sort(net.seed22$id))
#
# seq.sim <- read.dna("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/A.Epidemic27.Sequences.gene.pol.7yrs.fasta")
# tree.dat <- phyDat(seq.sim, type = "DNA")
# tree.ml <- dist.ml(tree.dat,model = "JC69")
# tree.sim <- upgma(tree.ml)

# sort(as.numeric(tree.sim$tip.label))==sort(net.seed27$id)
# table(sort(as.numeric(tree.sim$tip.label))==sort(net.seed27$id)) # number of TRUE's must equal to length(net.seed27$id)

#
# length(sort(as.numeric(tree.sim$tip.label)))
#
# setdiff(sort(as.numeric(tree.sim$tip.label)), sort(net.seed22$id))

# # rename sampling times automatically and make a single file for sampling dates of seeds 22, 24, and 27

# transfSampTimes.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/Tsamplingtimes_seed_number_22.csv", dec = ",")
# transfSampTimes.22$V1 <- as.character(transfSampTimes.22$V1)
# z <- transfSampTimes.22$V1
# transfSampTimes.22$V1 <- paste("22.",z,sep = "")
# write.csv(transfSampTimes.22, file = "Tsamplingtimes_seed_number_22.csv")
#
# transfSampTimes.24 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/Tsamplingtimes_seed_number_24.csv", dec = ",")
# transfSampTimes.24$V1 <- as.character(transfSampTimes.24$V1)
# x <- transfSampTimes.24$V1
# transfSampTimes.24$V1 <- paste("24.",x,sep = "")
# write.csv(transfSampTimes.24, file = "Tsamplingtimes_seed_number_24.csv")
#
# transfSampTimes.27 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/Tsamplingtimes_seed_number_27.csv", dec = ",")
# transfSampTimes.27$V1 <- as.character(transfSampTimes.27$V1)
# y <- transfSampTimes.27$V1
# transfSampTimes.27$V1 <- paste("27.",y,sep = "")
# write.csv(transfSampTimes.27, file = "Tsamplingtimes_seed_number_27.csv")
#
# d <- rbind(transfSampTimes.22[,3:4], transfSampTimes.24[,3:4], transfSampTimes.27[,3:4])
# dd <- as.data.frame(d)
# dd$V1 <- as.character(dd$V1)
# dd$V2 <- as.numeric(as.character(dd$V2)) # solve factors issue
#
# write.csv(dd, file = "samplingtimes_seed_number_22.24.27.csv")

setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/")

# 1.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 100% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.1098.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.1098.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.24.27 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.1098.fasta.tree", sep = ""))

samp.dates.22.24.27 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27$V1 <- as.character(samp.dates.22.24.27$V1)
samp.dates.22.24.27$V2 <- as.numeric(as.character(samp.dates.22.24.27$V2))

time.samp.22.24.27 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27)

tree.tips.22.24.27 <- (tree.fasttree.22.24.27$tip.label) # length==1098

c.samp.dates.22.24.27<- as.character(samp.dates.22.24.27$V1)

Ord.tree.dates.22.24.27 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27)){
  for(j in 1:length(time.samp.22.24.27)){
    if(tree.tips.22.24.27[i] == c.samp.dates.22.24.27[j]){
      Ord.tree.dates.22.24.27 <- c(Ord.tree.dates.22.24.27, time.samp.22.24.27[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27 <- dater(tree.fasttree.22.24.27,
                          Ord.tree.dates.22.24.27,
                          s = 3012,
                          omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27, file = paste("dated.tree.A.seed.22.24.27seq.1098.Rdata", sep = ""))
#e <- get(load("dated.tree.A.seed.22.24.27seq.1098.Rdata"))



# 2.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 90% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.90.988.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.90.988.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p90 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.90.988.fasta.tree", sep = ""))

samp.dates.22.24.27.p90 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p90$V1 <- as.character(samp.dates.22.24.27.p90$V1)
samp.dates.22.24.27.p90$V2 <- as.numeric(as.character(samp.dates.22.24.27.p90$V2))

time.samp.22.24.27.p90 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p90)

tree.tips.22.24.27.p90 <- (tree.fasttree.22.24.27.p90$tip.label) # length==1098

c.samp.dates.22.24.27.p90<- as.character(samp.dates.22.24.27.p90$V1)

Ord.tree.dates.22.24.27.p90 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p90)){
  for(j in 1:length(time.samp.22.24.27.p90)){
    if(tree.tips.22.24.27.p90[i] == c.samp.dates.22.24.27.p90[j]){
      Ord.tree.dates.22.24.27.p90 <- c(Ord.tree.dates.22.24.27.p90, time.samp.22.24.27.p90[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p90 <- dater(tree.fasttree.22.24.27.p90,
                             Ord.tree.dates.22.24.27.p90,
                             s = 3012,
                             omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p90, file = paste("dated.tree.A.seed.22.24.27seq.p90.988.Rdata", sep = ""))

# W <- node.age(e)


# 3.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 80% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.80.878.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.80.878.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p80 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.80.878.fasta.tree", sep = ""))

samp.dates.22.24.27.p80 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p80$V1 <- as.character(samp.dates.22.24.27.p80$V1)
samp.dates.22.24.27.p80$V2 <- as.numeric(as.character(samp.dates.22.24.27.p80$V2))

time.samp.22.24.27.p80 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p80)

tree.tips.22.24.27.p80 <- (tree.fasttree.22.24.27.p80$tip.label) # length==1098

c.samp.dates.22.24.27.p80<- as.character(samp.dates.22.24.27.p80$V1)

Ord.tree.dates.22.24.27.p80 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p80)){
  for(j in 1:length(time.samp.22.24.27.p80)){
    if(tree.tips.22.24.27.p80[i] == c.samp.dates.22.24.27.p80[j]){
      Ord.tree.dates.22.24.27.p80 <- c(Ord.tree.dates.22.24.27.p80, time.samp.22.24.27.p80[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p80 <- dater(tree.fasttree.22.24.27.p80,
                                 Ord.tree.dates.22.24.27.p80,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p80, file = paste("dated.tree.A.seed.22.24.27seq.p80.878.Rdata", sep = ""))



# 4.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 70% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.70.769.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.70.769.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p70 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.70.769.fasta.tree", sep = ""))

samp.dates.22.24.27.p70 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p70$V1 <- as.character(samp.dates.22.24.27.p70$V1)
samp.dates.22.24.27.p70$V2 <- as.numeric(as.character(samp.dates.22.24.27.p70$V2))

time.samp.22.24.27.p70 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p70)

tree.tips.22.24.27.p70 <- (tree.fasttree.22.24.27.p70$tip.label) # length==1098

c.samp.dates.22.24.27.p70<- as.character(samp.dates.22.24.27.p70$V1)

Ord.tree.dates.22.24.27.p70 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p70)){
  for(j in 1:length(time.samp.22.24.27.p70)){
    if(tree.tips.22.24.27.p70[i] == c.samp.dates.22.24.27.p70[j]){
      Ord.tree.dates.22.24.27.p70 <- c(Ord.tree.dates.22.24.27.p70, time.samp.22.24.27.p70[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p70 <- dater(tree.fasttree.22.24.27.p70,
                                 Ord.tree.dates.22.24.27.p70,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p70, file = paste("dated.tree.A.seed.22.24.27seq.p70.769.Rdata", sep = ""))


# 5.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 60% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.60.659.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.60.659.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p60 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.60.659.fasta.tree", sep = ""))

samp.dates.22.24.27.p60 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p60$V1 <- as.character(samp.dates.22.24.27.p60$V1)
samp.dates.22.24.27.p60$V2 <- as.numeric(as.character(samp.dates.22.24.27.p60$V2))

time.samp.22.24.27.p60 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p60)

tree.tips.22.24.27.p60 <- (tree.fasttree.22.24.27.p60$tip.label) # length==1098

c.samp.dates.22.24.27.p60<- as.character(samp.dates.22.24.27.p60$V1)

Ord.tree.dates.22.24.27.p60 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p60)){
  for(j in 1:length(time.samp.22.24.27.p60)){
    if(tree.tips.22.24.27.p60[i] == c.samp.dates.22.24.27.p60[j]){
      Ord.tree.dates.22.24.27.p60 <- c(Ord.tree.dates.22.24.27.p60, time.samp.22.24.27.p60[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p60 <- dater(tree.fasttree.22.24.27.p60,
                                 Ord.tree.dates.22.24.27.p60,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p60, file = paste("dated.tree.A.seed.22.24.27seq.p60.659.Rdata", sep = ""))



# 6.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 50% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.50.549.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.50.549.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p50 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.50.549.fasta.tree", sep = ""))

samp.dates.22.24.27.p50 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p50$V1 <- as.character(samp.dates.22.24.27.p50$V1)
samp.dates.22.24.27.p50$V2 <- as.numeric(as.character(samp.dates.22.24.27.p50$V2))

time.samp.22.24.27.p50 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p50)

tree.tips.22.24.27.p50 <- (tree.fasttree.22.24.27.p50$tip.label) # length==1098

c.samp.dates.22.24.27.p50<- as.character(samp.dates.22.24.27.p50$V1)

Ord.tree.dates.22.24.27.p50 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p50)){
  for(j in 1:length(time.samp.22.24.27.p50)){
    if(tree.tips.22.24.27.p50[i] == c.samp.dates.22.24.27.p50[j]){
      Ord.tree.dates.22.24.27.p50 <- c(Ord.tree.dates.22.24.27.p50, time.samp.22.24.27.p50[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p50 <- dater(tree.fasttree.22.24.27.p50,
                                 Ord.tree.dates.22.24.27.p50,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p50, file = paste("dated.tree.A.seed.22.24.27seq.p50.549.Rdata", sep = ""))


# 7.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 40% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.40.440.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.40.440.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p40 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.40.440.fasta.tree", sep = ""))

samp.dates.22.24.27.p40 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p40$V1 <- as.character(samp.dates.22.24.27.p40$V1)
samp.dates.22.24.27.p40$V2 <- as.numeric(as.character(samp.dates.22.24.27.p40$V2))

time.samp.22.24.27.p40 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p40)

tree.tips.22.24.27.p40 <- (tree.fasttree.22.24.27.p40$tip.label) # length==1098

c.samp.dates.22.24.27.p40<- as.character(samp.dates.22.24.27.p40$V1)

Ord.tree.dates.22.24.27.p40 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p40)){
  for(j in 1:length(time.samp.22.24.27.p40)){
    if(tree.tips.22.24.27.p40[i] == c.samp.dates.22.24.27.p40[j]){
      Ord.tree.dates.22.24.27.p40 <- c(Ord.tree.dates.22.24.27.p40, time.samp.22.24.27.p40[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p40 <- dater(tree.fasttree.22.24.27.p40,
                                 Ord.tree.dates.22.24.27.p40,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p40, file = paste("dated.tree.A.seed.22.24.27seq.p40.440.Rdata", sep = ""))


# 8.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 30% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.30.330.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.30.330.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p30 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.30.330.fasta.tree", sep = ""))

samp.dates.22.24.27.p30 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p30$V1 <- as.character(samp.dates.22.24.27.p30$V1)
samp.dates.22.24.27.p30$V2 <- as.numeric(as.character(samp.dates.22.24.27.p30$V2))

time.samp.22.24.27.p30 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p30)

tree.tips.22.24.27.p30 <- (tree.fasttree.22.24.27.p30$tip.label) # length==1098

c.samp.dates.22.24.27.p30<- as.character(samp.dates.22.24.27.p30$V1)

Ord.tree.dates.22.24.27.p30 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p30)){
  for(j in 1:length(time.samp.22.24.27.p30)){
    if(tree.tips.22.24.27.p30[i] == c.samp.dates.22.24.27.p30[j]){
      Ord.tree.dates.22.24.27.p30 <- c(Ord.tree.dates.22.24.27.p30, time.samp.22.24.27.p30[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p30 <- dater(tree.fasttree.22.24.27.p30,
                                 Ord.tree.dates.22.24.27.p30,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p30, file = paste("dated.tree.A.seed.22.24.27seq.p30.330.Rdata", sep = ""))



# 9.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years - 20% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.20.220.fasta", sep = ""), paste(">A.Epidemic22.24.27Sequences.gene.pol.7yrs.20.220.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.p20 <- read.tree(paste("A.Epidemic22.24.27Sequences.gene.pol.7yrs.20.220.fasta.tree", sep = ""))

samp.dates.22.24.27.p20 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_A/7_Years_back_all_seeds/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.p20$V1 <- as.character(samp.dates.22.24.27.p20$V1)
samp.dates.22.24.27.p20$V2 <- as.numeric(as.character(samp.dates.22.24.27.p20$V2))

time.samp.22.24.27.p20 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.p20)

tree.tips.22.24.27.p20 <- (tree.fasttree.22.24.27.p20$tip.label) # length==1098

c.samp.dates.22.24.27.p20<- as.character(samp.dates.22.24.27.p20$V1)

Ord.tree.dates.22.24.27.p20 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.22.24.27.p20)){
  for(j in 1:length(time.samp.22.24.27.p20)){
    if(tree.tips.22.24.27.p20[i] == c.samp.dates.22.24.27.p20[j]){
      Ord.tree.dates.22.24.27.p20 <- c(Ord.tree.dates.22.24.27.p20, time.samp.22.24.27.p20[j])
    }
  }
}

# calibrate internal nodes
dater.tree.22.24.27.p20 <- dater(tree.fasttree.22.24.27.p20,
                                 Ord.tree.dates.22.24.27.p20,
                                 s = 3012,
                                 omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.p20, file = paste("dated.tree.A.seed.22.24.27seq.p20.220.Rdata", sep = ""))


# Scenario 3 - A- B- G
#
# different subtypes of the virus (HIV-1-A-B-G) for all seeds
# complete sampling for a transmission network of one seed
# same sampling time interval of seven years

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario3_A/")

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

    if(i==22 || i==27){ # subtype A
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

    else if(i==24){ # subtype G
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


setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/")
datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistA.RData"))
simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)

seeds.great3 <- c(10, 18, 19, 21, 22, 24, 27, 32, 37) # seeds with at least 3 transmission events

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

net.seed22 = get(load("filtered.trans.net.seed22.Rdata")) # length(net.seed22$id) = 895 # Subtype A
net.seed24 = get(load("filtered.trans.net.seed24.Rdata")) # length(net.seed24$id) = 24 # Subtype G
net.seed27 = get(load("filtered.trans.net.seed27.Rdata")) # length(net.seed27$id) = 179 # Subtype A

# IDs to exclude in the sequences
id.22 <- setdiff(simpact.trans.net[[22]]$id, net.seed22$id) # Done! - 508
id.24 <- setdiff(simpact.trans.net[[24]]$id, net.seed24$id) # Done! - 692
id.27 <- setdiff(simpact.trans.net[[27]]$id, net.seed27$id) # Done! - 191

# # Transform sampling times > calender times
#
# transfSampTimes.22 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_22.csv")
# transfSampTimes.22$V1 <- paste("A.22.",as.character(transfSampTimes.22$V1), sep = "")
# write.csv(transfSampTimes.22, file = "Tsamplingtimes_seed_number_22.csv", quote = TRUE)
#
# transfSampTimes.24 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_24.csv")
# transfSampTimes.24$V1 <- paste("G.24.",as.character(transfSampTimes.24$V1), sep = "")
# write.csv(transfSampTimes.24, file = "Tsamplingtimes_seed_number_24.csv", quote = TRUE)
#
# transfSampTimes.27 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/Tsamplingtimes_seed_number_27.csv")
# transfSampTimes.27$V1 <- paste("A.27.",as.character(transfSampTimes.27$V1), sep = "")
# write.csv(transfSampTimes.27, file = "Tsamplingtimes_seed_number_27.csv", quote = TRUE)
#
#
# d <- rbind(transfSampTimes.22[,3:4], transfSampTimes.24[,3:4], transfSampTimes.27[,3:4])
# dd <- as.data.frame(d)
# dd$V1 <- as.character(dd$V1)
# dd$V2 <- as.numeric(as.character(dd$V2)) # solve factors issue
#
# write.csv(dd, file = "samplingtimes_seed_number_22.24.27.csv")
#
#


setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/")

# 1.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 100% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G$V1 <- as.character(samp.dates.22.24.27.A.A.G$V1)
samp.dates.22.24.27.A.A.G$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G$V2))

time.samp.22.24.27.A.A.G <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G)

tree.tips.22.24.27.A.A.G <- (tree.fasttree.22.24.27.A.A.G$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G<- as.character(samp.dates.22.24.27.A.A.G$V1)

Ord.tree.dates.22.24.27.A.A.G <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G)){
  for(j in 1:length(time.samp.22.24.27.A.A.G)){
    if(tree.tips.22.24.27.A.A.G[i] == c.samp.dates.22.24.27.A.A.G[j]){
      Ord.tree.dates.22.24.27.A.A.G <- c(Ord.tree.dates.22.24.27.A.A.G, time.samp.22.24.27.A.A.G[j])
      t <- c(t, tree.tips.22.24.27.A.A.G[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G <- dater(tree.fasttree.22.24.27.A.A.G,
                             Ord.tree.dates.22.24.27.A.A.G,
                             s = 3012,
                             omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.Rdata", sep = ""))
# V <- node.age(dater.tree.22.24.27.A.A.G) # with the two clades I got one internal node of age 1365.782



# 2.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 90% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.90.988.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.90.988.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p90 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.90.988.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p90  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p90$V1 <- as.character(samp.dates.22.24.27.A.A.G.p90$V1)
samp.dates.22.24.27.A.A.G.p90$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p90$V2))

time.samp.22.24.27.A.A.G.p90 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p90)

tree.tips.22.24.27.A.A.G.p90 <- (tree.fasttree.22.24.27.A.A.G.p90$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p90<- as.character(samp.dates.22.24.27.A.A.G.p90$V1)

Ord.tree.dates.22.24.27.A.A.G.p90 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p90)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p90)){
    if(tree.tips.22.24.27.A.A.G.p90[i] == c.samp.dates.22.24.27.A.A.G.p90[j]){
      Ord.tree.dates.22.24.27.A.A.G.p90 <- c(Ord.tree.dates.22.24.27.A.A.G.p90, time.samp.22.24.27.A.A.G.p90[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p90[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p90 <- dater(tree.fasttree.22.24.27.A.A.G.p90,
                                       Ord.tree.dates.22.24.27.A.A.G.p90,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p90, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p90.988.Rdata", sep = ""))


# 3.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 80% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.80.878.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.80.878.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p80 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.80.878.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p80  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p80$V1 <- as.character(samp.dates.22.24.27.A.A.G.p80$V1)
samp.dates.22.24.27.A.A.G.p80$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p80$V2))

time.samp.22.24.27.A.A.G.p80 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p80)

tree.tips.22.24.27.A.A.G.p80 <- (tree.fasttree.22.24.27.A.A.G.p80$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p80<- as.character(samp.dates.22.24.27.A.A.G.p80$V1)

Ord.tree.dates.22.24.27.A.A.G.p80 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p80)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p80)){
    if(tree.tips.22.24.27.A.A.G.p80[i] == c.samp.dates.22.24.27.A.A.G.p80[j]){
      Ord.tree.dates.22.24.27.A.A.G.p80 <- c(Ord.tree.dates.22.24.27.A.A.G.p80, time.samp.22.24.27.A.A.G.p80[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p80[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p80 <- dater(tree.fasttree.22.24.27.A.A.G.p80,
                                       Ord.tree.dates.22.24.27.A.A.G.p80,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p80, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p80.878.Rdata", sep = ""))



# 4.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 70% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.70.769.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.70.769.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p70 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.70.769.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p70  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p70$V1 <- as.character(samp.dates.22.24.27.A.A.G.p70$V1)
samp.dates.22.24.27.A.A.G.p70$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p70$V2))

time.samp.22.24.27.A.A.G.p70 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p70)

tree.tips.22.24.27.A.A.G.p70 <- (tree.fasttree.22.24.27.A.A.G.p70$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p70<- as.character(samp.dates.22.24.27.A.A.G.p70$V1)

Ord.tree.dates.22.24.27.A.A.G.p70 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p70)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p70)){
    if(tree.tips.22.24.27.A.A.G.p70[i] == c.samp.dates.22.24.27.A.A.G.p70[j]){
      Ord.tree.dates.22.24.27.A.A.G.p70 <- c(Ord.tree.dates.22.24.27.A.A.G.p70, time.samp.22.24.27.A.A.G.p70[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p70[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p70 <- dater(tree.fasttree.22.24.27.A.A.G.p70,
                                       Ord.tree.dates.22.24.27.A.A.G.p70,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p70, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p70.769.Rdata", sep = ""))


# 5.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 60% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.60.659.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.60.659.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p60 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.60.659.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p60  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p60$V1 <- as.character(samp.dates.22.24.27.A.A.G.p60$V1)
samp.dates.22.24.27.A.A.G.p60$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p60$V2))

time.samp.22.24.27.A.A.G.p60 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p60)

tree.tips.22.24.27.A.A.G.p60 <- (tree.fasttree.22.24.27.A.A.G.p60$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p60<- as.character(samp.dates.22.24.27.A.A.G.p60$V1)

Ord.tree.dates.22.24.27.A.A.G.p60 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p60)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p60)){
    if(tree.tips.22.24.27.A.A.G.p60[i] == c.samp.dates.22.24.27.A.A.G.p60[j]){
      Ord.tree.dates.22.24.27.A.A.G.p60 <- c(Ord.tree.dates.22.24.27.A.A.G.p60, time.samp.22.24.27.A.A.G.p60[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p60[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p60 <- dater(tree.fasttree.22.24.27.A.A.G.p60,
                                       Ord.tree.dates.22.24.27.A.A.G.p60,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p60, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p60.659.Rdata", sep = ""))


# 6.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 50% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.50.549.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.50.549.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p50 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.50.549.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p50  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p50$V1 <- as.character(samp.dates.22.24.27.A.A.G.p50$V1)
samp.dates.22.24.27.A.A.G.p50$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p50$V2))

time.samp.22.24.27.A.A.G.p50 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p50)

tree.tips.22.24.27.A.A.G.p50 <- (tree.fasttree.22.24.27.A.A.G.p50$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p50<- as.character(samp.dates.22.24.27.A.A.G.p50$V1)

Ord.tree.dates.22.24.27.A.A.G.p50 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p50)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p50)){
    if(tree.tips.22.24.27.A.A.G.p50[i] == c.samp.dates.22.24.27.A.A.G.p50[j]){
      Ord.tree.dates.22.24.27.A.A.G.p50 <- c(Ord.tree.dates.22.24.27.A.A.G.p50, time.samp.22.24.27.A.A.G.p50[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p50[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p50 <- dater(tree.fasttree.22.24.27.A.A.G.p50,
                                       Ord.tree.dates.22.24.27.A.A.G.p50,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p50, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p50.549.Rdata", sep = ""))


# 7.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 40% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.40.440.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.40.440.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p40 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.40.440.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p40  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p40$V1 <- as.character(samp.dates.22.24.27.A.A.G.p40$V1)
samp.dates.22.24.27.A.A.G.p40$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p40$V2))

time.samp.22.24.27.A.A.G.p40 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p40)

tree.tips.22.24.27.A.A.G.p40 <- (tree.fasttree.22.24.27.A.A.G.p40$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p40<- as.character(samp.dates.22.24.27.A.A.G.p40$V1)

Ord.tree.dates.22.24.27.A.A.G.p40 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p40)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p40)){
    if(tree.tips.22.24.27.A.A.G.p40[i] == c.samp.dates.22.24.27.A.A.G.p40[j]){
      Ord.tree.dates.22.24.27.A.A.G.p40 <- c(Ord.tree.dates.22.24.27.A.A.G.p40, time.samp.22.24.27.A.A.G.p40[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p40[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p40 <- dater(tree.fasttree.22.24.27.A.A.G.p40,
                                       Ord.tree.dates.22.24.27.A.A.G.p40,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p40, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p40.440.Rdata", sep = ""))



# 8.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 30% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.30.330.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.30.330.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p30 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.30.330.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p30  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p30$V1 <- as.character(samp.dates.22.24.27.A.A.G.p30$V1)
samp.dates.22.24.27.A.A.G.p30$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p30$V2))

time.samp.22.24.27.A.A.G.p30 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p30)

tree.tips.22.24.27.A.A.G.p30 <- (tree.fasttree.22.24.27.A.A.G.p30$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p30<- as.character(samp.dates.22.24.27.A.A.G.p30$V1)

Ord.tree.dates.22.24.27.A.A.G.p30 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p30)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p30)){
    if(tree.tips.22.24.27.A.A.G.p30[i] == c.samp.dates.22.24.27.A.A.G.p30[j]){
      Ord.tree.dates.22.24.27.A.A.G.p30 <- c(Ord.tree.dates.22.24.27.A.A.G.p30, time.samp.22.24.27.A.A.G.p30[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p30[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p30 <- dater(tree.fasttree.22.24.27.A.A.G.p30,
                                       Ord.tree.dates.22.24.27.A.A.G.p30,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p30, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p30.330.Rdata", sep = ""))


# 9.1 Phylogenetic tree for all seeds 22-24-27 in past 7 years for subtypes A-A and G - 20% of individuals

# (i) tree construction
system(paste("./FastTree  <", paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.30.330.fasta", sep = ""), paste(">A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.30.330.fasta.tree", sep = "")))

# (ii) internal node calibration
tree.fasttree.22.24.27.A.A.G.p20 <- read.tree(paste("A.A.GEpidemic22.27.24Sequences.gene.pol.7yrs.895.179.24.p.30.330.fasta.tree", sep = ""))

samp.dates.22.24.27.A.A.G.p20  <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_A/7_yrs_all_seeds_diff_subtype/samplingtimes_seed_number_22.24.27.csv", dec = ",")

samp.dates.22.24.27.A.A.G.p20$V1 <- as.character(samp.dates.22.24.27.A.A.G.p20$V1)
samp.dates.22.24.27.A.A.G.p20$V2 <- as.numeric(as.character(samp.dates.22.24.27.A.A.G.p20$V2))

time.samp.22.24.27.A.A.G.p20 <- dates.Transform.NamedVector(dates=samp.dates.22.24.27.A.A.G.p20)

tree.tips.22.24.27.A.A.G.p20 <- (tree.fasttree.22.24.27.A.A.G.p20$tip.label) # length==1098

c.samp.dates.22.24.27.A.A.G.p20<- as.character(samp.dates.22.24.27.A.A.G.p20$V1)

Ord.tree.dates.22.24.27.A.A.G.p20 <- vector() # order the dates according to tips order in the tree
t <- vector() # tips
for(i in 1:length(tree.tips.22.24.27.A.A.G.p20)){
  for(j in 1:length(time.samp.22.24.27.A.A.G.p20)){
    if(tree.tips.22.24.27.A.A.G.p20[i] == c.samp.dates.22.24.27.A.A.G.p20[j]){
      Ord.tree.dates.22.24.27.A.A.G.p20 <- c(Ord.tree.dates.22.24.27.A.A.G.p20, time.samp.22.24.27.A.A.G.p20[j])
      t <- c(t, tree.tips.22.24.27.A.A.G.p20[i])
    }
  }
}


# calibrate internal nodes
dater.tree.22.24.27.A.A.G.p20 <- dater(tree.fasttree.22.24.27.A.A.G.p20,
                                       Ord.tree.dates.22.24.27.A.A.G.p20,
                                       s = 3012,
                                       omega0 = 0.00475) # s is the length of sequence

save(dater.tree.22.24.27.A.A.G.p20, file = paste("dated.tree.A.A.G.seed.22.24.27seq.895.179.24.p20.330.Rdata", sep = ""))



# Scenario 4 - A- B- G
#
# different subtypes of the virus (HIV-1-A-B-G) for all seeds
# complete sampling for a transmission network of all seeds
# same sampling time interval (e.g.: five or three years) for a transmission network of all seeds

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/Scenario4_A/")


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

    if(i==22 || i==27){ # subtype A
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

    else if(i==24){ # subtype G
      # call the seed sequences and rename the file
      file.copy(paste("hiv.seq.G.pol.j.fasta", sep = ""),paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""))
      # add the number of tree in the file and
      write(n.tr,file = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
      # the tree, to prepare the file to simulate the evolution of the virus across the tree
      write.tree(tr,file = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
      file.rename(from = paste("hiv.seq.G.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

      system(paste("./seq-gen -mGTR -f 0.3987, 0.1563, 0.2202, 0.2249 -a 0.9 -i 0.80 -g 4 -r 1.4520, 9.9166, 1.3332, 1.2652, 14.9356, 1.0000 -n1 -k 1 <hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >G.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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
