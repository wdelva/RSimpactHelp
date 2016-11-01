##neatly load all the packages listed in p_load()
pacman::p_load(lhs)

#Set working directory
dirname <- getwd()
#Set the number of design points
design.points <- 10


#Set the number of input simpact params to vary (these will the parameters that we seek to calibrate).
input.varied.params <- c("person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b", "conception.alpha_base",
                       "formation.hazard.agegapry.numrel_man", "formation.hazard.agegapry.eagerness_diff",
                       "formation.hazard.agegapry.gap_factor_man_exp", "person.agegap.man.dist.normal.mu",
                       "person.agegap.woman.dist.normal.mu","person.agegap.man.dist.normal.sigma",
                       "person.agegap.woman.dist.normal.sigma")


#Upper and lower bounds of the varies parameters
input.varied.params.boundaries <- list(person.eagerness.man.dist.gamma.a.min =0.1,
                                       person.eagerness.man.dist.gamma.a.max =2,
                                       person.eagerness.man.dist.gamma.b.min = 5,
                                       person.eagerness.man.dist.gamma.b.max = 60,
                                       conception.alpha_base.min = -3.6,
                                       conception.alpha_base.max = -1.2,
                                       formation.hazard.agegapry.numrel_man.min = -1.5,
                                       formation.hazard.agegapry.numrel_man.max = -0.1,
                                       formation.hazard.agegapry.eagerness_diff.min = -0.1,
                                       formation.hazard.agegapry.eagerness_diff.max = 0,
                                       formation.hazard.agegapry.gap_factor_man_exp.min = -1.5,
                                       formation.hazard.agegapry.gap_factor_man_exp.max =-0.4,
                                       person.agegap.man.dist.normal.mu.min = 0,
                                       person.agegap.man.dist.normal.mu.max = 4,
                                       person.agegap.woman.dist.normal.mu.min =0,
                                       person.agegap.woman.dist.normal.mu.max = 4,
                                       person.agegap.man.dist.normal.sigma.min = 0.5,
                                       person.agegap.man.dist.normal.sigma.max =2,
                                       person.agegap.woman.dist.normal.sigma.min =0.5,
                                       person.agegap.woman.dist.normal.sigma.max =2)


# Creating the LHS over the 0-1 uniform parameter space for the parameters to be estimated
variables <- length(input.varied.params)
set.seed(1)
rlhs <- randomLHS(design.points, variables)

#Select the config parameters that will be varied from the input config
lhs.df <- data.frame(matrix(NA, nrow = 1, ncol = length(input.varied.params)))
names(lhs.df) <- input.varied.params

#Repeat have to meet the design.points
lhs.df <- as.data.frame(lapply(lhs.df, rep, design.points))

#Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)
x.index <- 0
for (j in input.varied.params){
  x.index <- x.index + 1
  min.var <- input.varied.params.boundaries[paste(j,".min",sep = "")][[1]]
  max.var <- input.varied.params.boundaries[paste(j,".max",sep = "")][[1]]
  col.index <- which(colnames(lhs.df)==j)
  lhs.df[col.index] <- qunif(rlhs[ , x.index], min = as.numeric(min.var), max = as.numeric(max.var))
}

#Set those that need to use the same values  (NEED to make this auto as well)
lhs.df$person.eagerness.woman.dist.gamma.a <- lhs.df$person.eagerness.man.dist.gamma.a
lhs.df$person.eagerness.woman.dist.gamma.b <- lhs.df$person.eagerness.man.dist.gamma.b
lhs.df$formation.hazard.agegapry.numrel_woman <- lhs.df$formation.hazard.agegapry.numrel_man
lhs.df$formation.hazard.agegapry.gap_factor_woman_exp <- lhs.df$formation.hazard.agegapry.gap_factor_man_exp
#lhs.df$success.rows <- NA

##This will create the input file for the simpact
inPUT.df <- cbind.data.frame(sim.id = 1:design.points, rlhs, lhs.df)

##save(inPUT.df, file = paste0(dirname,"inputANDoutput-",design.points,"Points",variables,"Par",Sys.Date(),"-RData"))

write.csv(inPUT.df, file =paste0("INPUT.df","-",design.points,"Points",variables,"Par",Sys.Date(), ".csv"), row.names = FALSE)


