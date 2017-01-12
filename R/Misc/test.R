## Testing the results from a simply simpact.run
pacman::p_load(RSimpactCyan, RSimpactHelper)

simulation.type <- "maxart"#"simpact-cyan"#

simpact.set.simulation(simulation.type)

agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)


iv <- intervention.introduced(simulation.type = simulation.type)


testinput <- input.params.creator(simulation.type = simulation.type)

testoutput <- simpact.run(configParams = testinput,
                          destDir = "temp",
                          agedist = agedist.data.frame,
                          intervention = iv,
                          seed = 20)

datalist.test <- readthedata(testoutput)


