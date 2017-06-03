# Fitting ABM with ABC

# 1. Run a master model with well estbalished parameters

pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper, phylosim, ape, lhs)

# run input.params.creator() and simpact.config.inputs()

all.sim.start <- proc.time()

set.new.seed <- 1
init.design.points <- 10 #set the initial design points
design.points.total <- 20 #argument init design point to this value
rep.sample <- ceiling(design.points.total/init.design.points) - 1


###### Generate the input parameters for the simulation ###############################################
inPUT.df.complete <- simpact.config.inputs(design.points = init.design.points, resample.count = 2,


                                           formation.hazard.agegapry.numrel_diff = c(-1,0), #-0.1
                                           formation.hazard.agegapry.eagerness_diff = c(-0.5, 0),#-0.048 #-0.110975
                                           conception.alpha_base = c(-5,0),# -2.35
                                           birth.boygirlratio = c(0.5,0.7) # 0.5024876 #101:100
                                           )
#


for (i in 1:rep.sample){
  inPUT.df.complete.new <- simpact.config.inputs.add.sample(datalist = inPUT.df.complete,

                                                            resample.points = init.design.points, set.seed.new = set.new.seed,
                                                            conception.alpha_base = c(-5, -0), #c(-4, -1.5)
                                                            formation.hazard.agegapry.numrel_diff = c(-1,0), #-0.1
                                                            formation.hazard.agegapry.eagerness_diff = c(-0.5, 0),#-0.048 #-0.110975
                                                            birth.boygirlratio = c(0.5,0.7) # 0.5024876 #101:100
                                                            )

  inPUT.df.complete <- rbind(inPUT.df.complete, inPUT.df.complete.new)
}


############### Argument sample above ##############################

#Select a chunk to send to process
min.chunk <- 2227
max.chunk <- 2227

if(max.chunk > nrow(inPUT.df.complete)){max.chunk <- nrow(inPUT.df.complete)}
if(min.chunk > nrow(inPUT.df.complete) || min.chunk < 1){min.chunk <- max.chunk}

inANDout.df.chunk <- inPUT.df.complete[min.chunk:max.chunk,]

#make sure there are no empty rows
inANDout.df.chunk <- inANDout.df.chunk[!is.na(inANDout.df.chunk$sim.id),]

#set how many time the single row will be repeated
sim_repeat <- 24

# number of cores per node
ncluster.use <- 2


#indicate the target statitics that you want to hit

# target statistics: average relationships for men and women,
# standard deviation  relationships for men and women, and average population growth rate

target.variables <- c("aver.rels.men", "aver.rels.women", "sd.rels.men", "sd.rels.women", "growth.rate")

##Each of these should be calculated after each run, else we give an NA

#set the prior names - varied parameters
preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

#rbind all the results for this chunk to be merged after
#Create a dataframe with NA for the summary statistics Will collect all the chunks with the sim.id to link back
chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 0, ncol = length(target.variables)+2))
names(chunk.summary.stats.df) <- c(target.variables, "sim.id")


############   MAIN Simulation is here without sequence simulation #######################

simpact4ABC.chunk.wrapper <- function(simpact.chunk.prior){

}
