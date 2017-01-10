### NOT sure who wrote this (but fixing the parameter initialisation for simpact)
pacman::p_load(dplyr, RSimpactCyan, RSimpactHelper)

#test the number range to use with the formation.hazard.agegapry.baseline (-5,5)
#calculate growth rate and prevelance15-49yo

testinput <- input.params.creator()

res.agegapry <- data.frame(list(NA,NA,NA))
names(res.agegapry) <- c("agegapry.baseline","pop.growth","prev15.49")
seed.id = 0
for (j in seq(0,3,0.1)){
  seed.id = seed.id + 1
  testinput$formation.hazard.agegapry.baseline <- j
  testoutput <- simpact.run(configParams = testinput,
                            destDir = "temp",
                            agedist = agedist.data.frame,
                            intervention = iv,
                            seed = 1)
  datalist.test <- readthedata(testoutput)
  c.pop.growth <- pop.growth.calculator(datalist = datalist.test)
  c.prev15.49 <- prevalence.calculator(datalist = datalist.test, agegroup = c(15, 50), timepoint = 34)
  c.prev15.49 <- c.prev15.49$pointprevalence[3]

  res.agegapry <- rbind(res.agegapry,c(j,c.pop.growth,c.prev15.49))

}

plot(res.agegapry$agegapry.baseline, res.agegapry$pop.growth,
     xlab = "agegapry.value", ylab = "Population Growth Summmary")

plot(res.agegapry$agegapry.baseline, res.agegapry$prev15.49,
     xlab = "agegapry.value", ylab = "Prevalence 15-49yo Summmary")



