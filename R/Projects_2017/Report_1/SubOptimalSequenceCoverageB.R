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
# same sampling time interval (e.g.: five or three years) for a transmission network of one seed

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

    system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -n1 -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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

    system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -n1 -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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



# Scenario 3 - A- B- G
#
# different subtypes of the virus (HIV-1-A-B-G) for all seeds
# complete sampling for a transmission network of one seed
# same sampling time interval (e.g.: five or three years) for a transmission network of one seed

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

      system(paste("./seq-gen -mGTR -f 0.3987, 0.1563, 0.2202, 0.2249 -a 0.9 -i 0.80 -g 4 -r 1.4520, 9.9166, 1.3332, 1.2652, 14.9356, 1.0000 -n1 -k 1 <hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >G.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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

      system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -n1  -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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

      system(paste("./seq-gen -mGTR -f 0.3935, 0.1708, 0.2060, 0.2297 -a 0.9 -i 0.80 -g 4 -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -n1 -k 1 <hiv.seq.B.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >B.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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



# Scenario 4 - A- B- G
#
# different subtypes of the virus (HIV-1-A-B-G) for all seeds
# complete sampling for a transmission network of all seeds
# same sampling time interval (e.g.: five or three years) for a transmission network of all seeds

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

      system(paste("./seq-gen -mGTR -f 0.3987, 0.1563, 0.2202, 0.2249 -a 0.9 -i 0.80 -g 4 -r 1.4520, 9.9166, 1.3332, 1.2652, 14.9356, 1.0000 -n1 -k 1 <hiv.seq.G.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >G.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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

      system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -n1  -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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

      system(paste("./seq-gen -mGTR -f 0.3935, 0.1708, 0.2060, 0.2297 -a 0.9 -i 0.80 -g 4 -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -n1 -k 1 <hiv.seq.B.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >B.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


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

