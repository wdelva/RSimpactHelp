## Testing the results from a simply simpact.run
pacman::p_load(RSimpactCyan, RSimpactHelper, lhs)

simulation.type <- "simpact-cyan"#"maxart"#

simpact.set.simulation(simulation.type)

agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)


iv <- intervention.introduced(simulation.type = simulation.type)


testinput <- input.params.creator(population.simtime = 40,
                                  population.numwomen = 1000,
                                  population.nummen = 1000,
                                  simulation.type = simulation.type)


testoutput <- simpact.run(configParams = testinput,
                          destDir = "temp",
                          agedist = agedist.data.frame,
                          intervention = iv,
                          seed = 2)

datalist <- readthedata(testoutput)

pop.growth.calculator(datalist = datalist.test,
                      timewindow = c(0,
                          timewindow.max=unique(datalist.test$itable$population.simtime)))

incidence.calculator(datalist = datalist.test, agegroup = c(20, 25),
                     timewindow = c(32, 35), only.active = "No")

prevalence.calculator(datalist = datalist.test, agegroup = c(15, 50),
                    timepoint = 30)

#Example of Latin hypercube sample (each row and column is filled
#with one point)
#With an addittion to re-sample after removing the rows that do not meet
#the success criterion

#One sample in each row and each column
#Sampling is random in each grid
#Generate higher dimmension sampling similarly
#Avoid the propabaility that all sampling points come from the same local region
set.seed(50)
n.points <- 10
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





