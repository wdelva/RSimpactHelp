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

datalist.test <- readthedata(testoutput)

pop.growth.calculator(datalist = datalist.test,
                      timewindow = c(0,
                          timewindow.max=unique(datalist.test$itable$population.simtime)))

incidence.calculator(datalist = datalist.test, agegroup = c(20, 25),
                     timewindow = c(32, 34), only.active = "No")

prevalence.calculator(datalist = datalist.test, agegroup = c(18, 20),
                    timepoint = 34)

#Example of Latin hypercube sample (each row and column is filled
#with one point)
#With an addittion to re-sample after removing the rows that do not meet
#the success criterion

#One sample in each row and each column
#Sampling is random in each grid
#Generate higher dimmension sampling similarly
#Avoid the propabaility that all sampling points come from the same local region
set.seed(1)
n.points <- 20000
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


# Tree imbalance

id1 <- c(seq(from = 0, to = 14, by = 1))

par1 <- c(-1, 0,0,1,1,2,2,3,3,4,4,5,5,6,6)

itime1 <- c(seq(from = 1, to = 15, by = 1))

id2 <- c(seq(from = 0, to = 14, by = 1))

par2 <- c(seq(from = -1, to = 13, by =1))

itime2 <- c(seq(from = 1, to = 15, by = 1))

epi1 <- list()
epi1$itimes <- rev(itime1*0.1)
epi1$dtimes <- rep(0,length(itime1))
epi1$id <- id1
epi1$parent <- par1

tree1 <- epi2tree(epi1)

epi2 <- list()
epi2$itimes <- rev(itime2*0.1)
epi2$dtimes <- rep(0,length(itime2))
epi2$id <- id2
epi2$parent <- par2

tree2 <- epi2tree(epi2)

plot(tree1) # synetric tree >> slope is around -2

plot(tree2) # ladder like tree >> slope is not defined NA

trend1 <- phylogenetictree.trend(tree = tree1)
x1 = trend1$num.tree
y1 = trend1$size.tree
reg1 = lm(y1 ~ x1)

trend2 <- phylogenetictree.trend(tree = tree2)
x2 = trend2$num.tree
y2 = trend2$size.tree
reg2 = lm(y2 ~ x2)

