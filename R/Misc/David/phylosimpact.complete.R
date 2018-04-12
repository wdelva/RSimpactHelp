# Within TEST_19_1_2018 folder we have two folders:

# seqgen : containing seq-gen and seed sequence
# Fasttree : containing FastTree


# The exercise is about:

# (i) running simpact
# (ii) simulate sequences
# (iii) construct and calibrate the phylogenetic tree from these sequences


setwd("/home/david/Desktop/TEST_19_1_2018/")

# ## For sub-optimal sequence coverage


## Load required packages

pacman::p_load(RSimpactCyan, RSimpactHelper, devtools, Rcpp, ape, expoTree,
               data.table, tidyr, phylosim, readr, dplyr, adephylo, geiger,
               picante, igraph, ggplot2, magrittr, lmtest, gsubfn, utils, pcaPP,
               phyloTop, phytools, lme4, data.table, treedater, phangorn, nlme,
               fitdistrplus, apTreeshape, DECIPHER)

#######################
# Step 1: Run Simpact #
#######################

## Run Simpact for specific parameter combination


destDir <- "/home/david/Desktop/TEST_19_1_2018/temp" # laptop



age.distr <- agedistr.creator(shape = 5, scale = 65)

#source the input files var name <- cfg.list <- cfg
cfg.list <- input.params.creator( # population.eyecap.fraction = 0.2, #0.21,#1,
                                 # population.msm = "no",
                                 population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
                                 population.nummen = 600, #3000, #600, # 3800, #2500,
                                 population.numwomen = 600, # 3000, #600, #4200, #2500,
                                 hivseed.time = 10, # 20,
                                 hivseed.type = "amount",
                                 hivseed.amount = 20, #30,
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
cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
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
cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.4
cfg.list["diagnosis.baseline"] <- -2


# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
art.intro <- list()
art.intro["time"] <- 20
art.intro["diagnosis.baseline"] <- -2 # 0#100
art.intro["monitoring.cd4.threshold"] <- 100 # 1200

### add something about diagnosis
art.intro["diagnosis.agefactor"] <- 0
art.intro["diagnosis.genderfactor"] <- 0
art.intro["diagnosis.diagpartnersfactor"] <- 0
art.intro["diagnosis.isdiagnosedfactor"] <- 0
### end of add-on about diagnosis
#art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
art.intro1 <- list()
art.intro1["time"] <- 22
art.intro1["diagnosis.baseline"] <- -2 # 0#100
art.intro1["monitoring.cd4.threshold"] <- 150 # 1200

art.intro2 <- list()
art.intro2["time"] <- 25 # inputvector[5] ######### 30
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 36
art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access

inputvector <- c(151,1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)


# tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist # scenario(interventionlist, tasp.indicator)

cfg.list["hivtransmission.param.f1"] <- log(inputvector[2])
cfg.list["hivtransmission.param.f2"] <- log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[3]
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[3]
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]

#cfg <- cfg.list


cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
# cfg["monitoring.fraction.log_viralload"] <- 0.3
cfg.list["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

seedid <- inputvector[1]
#cfg.list["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
cfg.list["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10] ######### -0.5
cfg.list["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10] ######### -0.5
cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[11]

cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
cfg.list["conception.alpha_base"] <- inputvector[14] #is conception.alpha.base (higher up)
cfg.list["dissolution.alpha_0"] <- inputvector[15]
cfg.list["dissolution.alpha_4"] <- inputvector[16]


results <- simpact.run(configParams = cfg.list,
                       destDir = destDir,
                       agedist = age.distr,
                       seed = seedid, #, Introducing ART has helped to keep the prevalence high
                       intervention = intervention)


datalist <- readthedata(results)



###########################################
# Step 2: Construct transmission networks #
###########################################


# Resource required RSimpactHelp function in my branch
# source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff3.R")
# source("/home/david/RSimpactHelp/R/epi2tree2.R")
# source("/home/david/RSimpactHelp/R/trans.network2tree.R")



simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)

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

which(smallest.branches!="NA")




###############################
# Step 3: Sequence simulation #
###############################


# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123

trans.net <- simpact.trans.net # all transmission networks

# source("/home/david/RSimpactHelp/R/sequence.simulation.segen.R")

dirseqgen <- "/home/david/Desktop/TEST_19_1_2018/seqgen"

sequence.simulation.seqgen(dir.seq = dirseqgen,
                           seq.gen.tool = "seq-gen",
                           datalist = datalist,
                           seeds.num = 123,
                           endpoint = 40,
                           limitTransmEvents = 3,
                           seed.file = "hiv.seq.C.pol.j.fasta")


#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################


# source("/home/david/RSimpactHelp/R/phylogenetic.tree.fasttree.R")


dirfasttree <- "/home/david/Desktop/TEST_19_1_2018/Fasttree"


tree.calib <- phylogenetic.tree.fasttree(dir.tree = dirfasttree,
                                         dir.seq = dirseqgen,
                                         fasttree.tool = "FastTree",
                                         dir.tree = dirfasttree,
                                         calendar.dates = "dates.csv",
                                         simseqfile = "C.Epidemic.sequence.fasta",
                                         count.start = 1977,
                                         endsim = 40)



