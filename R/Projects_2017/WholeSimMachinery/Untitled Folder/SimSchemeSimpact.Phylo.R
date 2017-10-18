
setwd("/home/david/RSimpactHelp/R/Projects_2017/WholeSimMachinery/")

## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings,
               phyclust, DECIPHER)


#######################
# Step 1: Run Simpact #
#######################

# Set up the maodel

simulation.type <- "simpact-cyan"#"maxart"#

simpact.set.simulation(simulation.type)

agedist.data.frame <- agedistr.creator(shape = 5, scale = 70)


iv <- intervention.introduced(simulation.type = simulation.type)


# Set up parameters of the model

mastermodel.input <- input.params.creator(population.simtime = 40,
                                          population.numwomen = 4000,
                                          population.nummen = 3600,

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


# Read the mode output
mastermodel.datalist <- readthedata(mastermodel.output)

# Save the output
save(mastermodel.datalist, file = "DummyMaster.datalist.RData")

# Read saved output data set
master.datalist <- get(load("DummyMaster.datalist.RData"))


#################################
# Step 2: Transmission networks #
#################################


# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")

# Transmission network in a form of constructing a transmission network
simpact.trans.net <- transmNetworkBuilder.diff(datalist = master.datalist, endpoint = 40)


###############################
# Step 3: Sequence simulation #
###############################


# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered


trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]]

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
    num.i <- c(num.i,i)

    # Save the transmission tree
    write.tree(tree.i, file = paste("tree.model1.seed",i,".nwk", sep = ""))

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
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("HIV.Pol.gene.fasta", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(1,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed
  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds


#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################


# 4.1. Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree



# 4.2. Internal node optimisation requires sampled dates
