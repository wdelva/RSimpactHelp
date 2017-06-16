
setwd("/home/david/Dropbox/Fitting_Simpact/")

## Load required packages

pacman::p_load(RSimpactCyan, RSimpactHelper)


# Set up the maodel

simulation.type <- "simpact-cyan"#"maxart"#

simpact.set.simulation(simulation.type)

agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)


iv <- intervention.introduced(simulation.type = simulation.type)


# Set up parameters of the model

mastermodel.input <- input.params.creator(population.simtime = 40,
                                  population.numwomen = 1000,
                                  population.nummen = 1000,
                                  simulation.type = simulation.type)


# Run the model

mastermodel.output <- simpact.run(configParams = mastermodel.input,
                          destDir = "temp",
                          agedist = agedist.data.frame,
                          intervention = iv,
                          seed = 2)


# Read the mode outpu

datalist.test <- readthedata(mastermodel.output)




# Get the summary statistics as you wish

pop.growth.calculator(datalist = datalist.test,
                      timewindow = c(0,
                                     timewindow.max=unique(datalist.test$itable$population.simtime)))

incidence.calculator(datalist = datalist.test, agegroup = c(20, 25),
                     timewindow = c(32, 34), only.active = "No")

prevalence.calculator(datalist = datalist.test, agegroup = c(18, 20),
                      timepoint = 34)
