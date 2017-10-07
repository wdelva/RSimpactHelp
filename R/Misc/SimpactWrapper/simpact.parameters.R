#I presume this will be called through source()
#Idea is to give the simpact parameters and their ranges
# to be calibrated during the simpact run

init.design.points <- 100 #set the initial design points

#choose if select ids or just a range
sel.list <- "list" #min.max

if(sel.list == "list"){
  #set the list of ids that you want to simulate
  sel.id.list <- c(75, 83)
  }else{
    #Select a chunk to process
    min.chunk <- 11
    max.chunk <- 12
}

#set how many time the single row will be repeated
sim_repeat <- 2
#number of cores per node
ncluster.use <- 4

###### Generate the input parameters for the simulation #################################
simpact.config.inputs(design.points = init.design.points, resample.count = 1,
                       conception.alpha_base = c(-5, -0.1),
                       person.art.accept.threshold.dist.fixed.value = c(0.4, 1),
                       person.agegap.man.dist.normal.sigma = c(0.5, 3.5),
                       person.agegap.woman.dist.normal.sigma = c(0.5, 3.5),
                       hivtransmission.param.f1 = c(log(2), log(3.5))
)
############### Argument sample above ##################################################

