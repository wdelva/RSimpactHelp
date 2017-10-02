## Testing the results from a simply simpact.run
pacman::p_load(RSimpactCyan, RSimpactHelper, lhs)

#Set up simulation parameters and initiate simulation
sim.start.full <- as.Date("1970-03-31")
maxart.starttime <- as.Date("2014-09-01")
maxart.endtime <- as.Date("2017-08-31")
sim.end.full <- as.Date("2019-03-31")
seed.hiv.time <- round(as.numeric(difftime(as.Date("1986-03-31"), sim.start.full, units = "days")/365.242),0)
sim.duration <- round(as.numeric(difftime(sim.end.full,sim.start.full, units = "days")/365.242),0)
simulation.type <- "maxart"#"simpact-cyan"#
simpact.set.simulation(simulation.type)


agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)
#agedist.data.frame <- agedistr.creator(shape = 2, scale = 25)
#This matches the 1970 UN population

iv <- intervention.introduced(simulation.type = simulation.type)

#initial population

init.population.total <- 3000
women.frac <- 0.5253
num.women.prop <- round(women.frac * init.population.total, 0)
num.men.prop <- init.population.total - num.women.prop

testinput <- input.params.creator(population.simtime = sim.duration,
                                  population.numwomen = num.women.prop, # 1051, #1000,
                                  population.nummen = num.men.prop, # 949, #1000,
                                  person.art.accept.threshold.dist.fixed.value = 0.9,
                                  debut.debutage = 10,  #maxart
                                  hivseed.time = seed.hiv.time, #1986 seed HIV
                                  simulation.type = simulation.type)

if(simulation.type == "maxart"){
  testinput$facilities.randomization <- "${SIMPACT_DATA_DIR}maxart-randomization.csv"
  testinput$maxart.starttime <- round(as.numeric(difftime(maxart.starttime ,sim.start.full, units = "days")/365),0)
  testinput$person.geo.dist2d.discrete.maskfile <-  ""
}

testoutput <- simpact.run(configParams = testinput,
                          destDir = "temp",
                          agedist = agedist.data.frame,
                          intervention = iv,
                          seed = 8)

datalist.test <- readthedata(testoutput)

datalist.test$ptable <- client.facility(datalist = datalist.test, site = "MaxART")
datalist.test$ptable$pfacility[datalist.test$ptable$pfacility.value > 15] <- "Not Hhohho"

table(datalist.test$ptable$pfacility)

pop.growth.calculator(datalist = datalist.test,
                      timewindow = c(0,
                          timewindow.max=datalist.test$itable$population.simtime[1]))

#check incidence within the study window

inc.study.time.u <- as.numeric(difftime(maxart.starttime, sim.start.full, units = "days"))
inc.study.time.d <- as.numeric(difftime(maxart.endtime, sim.start.full, units = "days"))

inc.study.time.u <- round(inc.study.time.u/365.242,0)
inc.study.time.d <- round(inc.study.time.d/365.242,0)


incidence.calculator(datalist = datalist.test, agegroup = c(10, 49),
                     timewindow = c(inc.study.time.u, inc.study.time.d),
                     only.active = "No")

prevalence.calculator(datalist = datalist.test, agegroup = c(10, 50),
                    timepoint = inc.study.time.d)


###################################################################################################


#Example of Latin hypercube sample (each row and column is filled
#with one point)
#With an addittion to re-sample after removing the rows that do not meet
#the success criterion

#One sample in each row and each column
#Sampling is random in each grid
#Generate higher dimmension sampling similarly
#Avoid the propabaility that all sampling points come from the same local region
set.seed(1)
n.points <- 20
rlhs.test <- randomLHS(n.points, 2)
break.point = 1/n.points
axis.range <- seq(0,1,break.point)
plot(rlhs.test, type = "p", pch = 19, xaxt = "n", yaxt = "n",
     col = "red", xlim = c(0,1), ylim = c(0,1),
     xlab = "Probability of x1", ylab = "Probability of x2",
     main = paste0(" Two dimension Augmented Latin Hypercube with 10 + 10 Points"),
     cex.lab = 2, cex.axis = 2, cex.main = 2, mgp=c(2.5,2.5,0))
abline(h = axis.range, v = axis.range, col = "blue", lty = 3, lwd = 1)
#set the axis
axis(1, at=axis.range)
axis(2, at=axis.range)

#pause to view the sample
Sys.sleep(3)

#do a resample and plot
rlhs.test.reg <- augmentLHS(rlhs.test,n.points)
rlhs.test.2 <- tail(rlhs.test.reg, n.points)
points(rlhs.test.2, type = "p", col = "black", pch = 19)

Sys.sleep(3)
#show how the resample happens 10 times more
for(i in 1:5){
  #set.seed(2)
  rlhs.test.reg <- augmentLHS(rlhs.test.reg,n.points)
  rlhs.test.3 <- tail(rlhs.test.reg, n.points)

  points(rlhs.test.3, type = "p", col = i+10 , pch = 19)
  Sys.sleep(5)
}





