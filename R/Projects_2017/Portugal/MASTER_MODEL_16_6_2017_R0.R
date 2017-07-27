
setwd("/home/david/Dropbox/Fitting_Simpact/MasterModel/")

## Load required packages

pacman::p_load(RSimpactCyan, RSimpactHelper)


# Set up the maodel

simulation.type <- "simpact-cyan"#"maxart"#

simpact.set.simulation(simulation.type)

agedist.data.frame <- agedistr.creator(shape = 5, scale = 70)


iv <- intervention.introduced(simulation.type = simulation.type)


# Set up parameters of the model

mastermodel.input <- input.params.creator(population.simtime = 40,
                                  population.numwomen = 1000,
                                  population.nummen = 1000,

                                  mortality.normal.weibull.shape = 5, # k for shape
                                  mortality.normal.weibull.scale = 70, # gamma for scale

                                  person.eagerness.man.dist.gamma.a = 0.85, #0.1
                                  person.eagerness.man.dist.gamma.b = 70, #100#3.5#5#10#20 #170
                                  person.eagerness.woman.dist.gamma.a = 0.1,
                                  person.eagerness.woman.dist.gamma.b = 70,#100#3.5#5#10#20#170

                                  #formation.hazard.agegapry.eagerness_diff = -0.110975,
                                  dissolution.alpha_0 = -0.6, #-0.1 # baseline
                                  dissolution.alpha_4 = -0.15,

                                  conception.alpha_base = -2, #c(-5, -1.5), #c(-4, -1.5)
                                  formation.hazard.agegapry.eagerness_diff = -0.1, # c(-0.5, 0),#-0.048 #-0.110975
                                  birth.boygirlratio = 0.501, # c(0.5,0.7), # 0.5024876 #101:100
                                  hivtransmission.param.a = -1.01, # c(-2,-1) # baseline of transmission: -1.0352239

                                  simulation.type = simulation.type)




# Run the model

mastermodel.output <- simpact.run(configParams = mastermodel.input,
                          destDir = "temp",
                          agedist = agedist.data.frame,
                          intervention = iv,
                          seed = 123)


# Read the mode outpu

mastermodel.datalist <- readthedata(mastermodel.output)

save(mastermodel.datalist, file = "master.datalist.RData")




# Get the summary statistics as you wish

pop.growth <- pop.growth.calculator(datalist = mastermodel.datalist,
                      timewindow = c(0,
                                     timewindow.max=unique(mastermodel.datalist$itable$population.simtime)))

# incidence.calculator(datalist = mastermodel.datalist, agegroup = c(20, 25),
#                      timewindow = c(32, 34), only.active = "No")
#
# prevalence.calculator(datalist = mastermodel.datalist, agegroup = c(18, 20),
#                       timepoint = 34)


rels.rate <- relationship.rate.calculator(datalist = mastermodel.datalist,
                                          timewindow = c(0, 40),
                                          int = FALSE, by = 1)

transm.rate <- transmission.rate.calculator(datalist = mastermodel.datalist,
                                            timewindow = c(0, 40),
                                            int = FALSE, by = 1)


R0 <- basicnumber.calculator(datalist = mastermodel.datalist,
                                                   beta = 0.1508, trans.rate.int = trans.rate.int)

# Phylogenetic component

trans.net <- transmNetworkBuilder.diff(datalist = mastermodel.datalist,
                                       endpoint = 40, population.simtime=40)

hiv.seq.raw <- read_file("~/RSimpactHelp/R/Projects_2017/HIV_Seq_K03455.txt")
## Remove the break line in the string of DNA
clean.hiv.seq <-  gsub("\n", "", hiv.seq.raw)

## For any part of the DNA you want to study its evolution,
# retrieve the range of interest.
## We downloaded the K03455, HIV-1 sequence with 9719 nucleotides

## Choose the gene of interest: env from 6172 to 8742 nucl position

hiv.seq.env <- substr(clean.hiv.seq, 6172,6372) # true c(6172,8742)

# Calculate the nucleotides frequencies
#
#       # library(Biostrings)
#       seq1 = DNAString(hiv.seq.env) # nulceotides
#
#       # Chech frequencies
#       freq <- letterFrequency(seq1, letters = c("A", "C", "G", "T"))/nchar(hiv.seq.env)
#       # > letterFrequency(gag, letters = c("A", "C", "T", "G"))/nchar(hiv.seq.gag)

freq <- c(0.3353293,0.2035928,0.2628077,0.1982701)
rate <- list("a"=0.2, "b"=0.6, "c"=0.12,"d"=0.001, "e"=0.25, "f"=0.24)


phylo.alpha.1 <- vector()
phylo.beta.1 <- vector()

for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]]

  if(nrow(as.data.frame(tree.n)) >= 5){
    tree.i <- trans.network2tree(transnetwork = tree.n)


    # Sequence simulation
    sim <- sequence.simulation(transtree = tree0, seedSeq = hiv.seq.env, alpha = 0.90,
                               rate.list = rate, base.freq = freq)
    saveAlignment.PhyloSim(sim,file = paste("HIVSeq_name",i,".fasta",sep=""), skip.internal = TRUE, paranoid = TRUE)

    # read the sequences
    seq.sim <- read.FASTA(paste("HIVSeq_name",i,".fasta",sep=""))
    tree.dat <- phyDat(seq.sim, type = "DNA")
    tree.ml <- dist.ml(tree.dat)
    tree.sim <- upgma(tree.ml)

    # Tree statistics
    xy <- phylogenetictree.trend(tree = tree.sim)
    x <- xy$num.tree
    y <- xy$size.tree
    reg <- lm(log(y) ~ log(x))
    cozf <- coef(reg)

    phylo.alpha.1 <- c(phylo.alpha.1, cozf[[1]]) # intercept
    phylo.beta.1 <- c(phylo.beta.1, cozf[[2]]) # slope

  }

}

phylo.alpha <- sum(phylo.alpha.1)/length(phylo.alpha.1)

phylo.beta <- sum(phylo.beta.1)/length(phylo.beta.1)

summary.statistics.df <- as.data.frame(cbind(pop.growth, rels.rate, transm.rate, phylo.alpha, phylo.beta, R0))

write.csv(summary.statistics.df, file = "summary.statistics.df.Master.csv", sep = "")
