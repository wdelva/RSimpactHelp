##neatly load all the packages listed in p_load()
pacman::p_load(lhs)

#Set working directory
dirname <- getwd()
#Set the number of design points
design.points <- 100


#Set the number of input simpact params to vary (these will be the parameters that we seek to calibrate).
input.varied.params <- c("conception.alpha_base", "formation.hazard.agegapry.baseline")


#Upper and lower bounds of the varies parameters
input.varied.params.boundaries <- list(conception.alpha_base.min = -3.6,
                                       conception.alpha_base.max = -1.2,
                                       formation.hazard.agegapry.baseline.min = 1.5,
                                       formation.hazard.agegapry.baseline.max = 3)


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


#Name the xdesign dataframe
rlhs <- data.frame(rlhs)
names(rlhs) <- paste0("xdesign",1:length(rlhs))

##This will create the input file for the simpact
inPUT.df <- cbind.data.frame(sim.id = 1:design.points, rlhs, lhs.df)

##save(inPUT.df, file = paste0(dirname,"inputANDoutput-",design.points,"Points",variables,"Par",Sys.Date(),"-RData"))

write.csv(inPUT.df, file =paste0("SWAZIINPUT.df","-",design.points,"Points",variables,"Par",Sys.Date(), ".csv"), row.names = FALSE)


