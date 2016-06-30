library(fBasics) 
library(VarianceGamma) 
library(stats4) 
#Likelihood Estimation library
#library(MASS) 
library("fitdistrplus")
#estimate beta parameters using Max Likelihood Est ,le
library("EnvStats")
degree.distr.fitter <- function(dataframe = degree.df){
  
  # Fit various parametric distribution to the degree data
  # and evaluate / compare their goodness-of-fit.
  # Normal, Poisson, Negative Bionomial, Gamma, Beta, Log-normal, Weibull, ...
  # And return the parameters (e.g. mean and variance or scale and shape) of the
  # best fitting distribution.
  mean.data<-mean(dataframe$degree)
  #Normal distribution
  nrows<-length(dataframe$degree)
  #estimate parameters
  normal.est<-fitdistr(dataframe$degree,"normal") 
  x.normal= rnorm(n=nrows,m=normal.est$estimate[1], sd=normal.est$estimate[2]) 
  
  goodness.normal<-mean(abs(dataframe$degree-trunc(x.normal)))

  #Poisson distribution
  poisson.est<-fitdistr(dataframe$degree,"poisson") 
  x.poisson= rpois(n=nrows,lambda=poisson.est$estimate) 
  goodness.poisson<-mean(abs(dataframe$degree-trunc(x.poisson)))
  #Negative Bionomial distribution
  nbinom.est<-fitdist(dataframe$degree, "nbinom")
  x.nbinom=dnbinom(dataframe$degree, size=nbinom.est$estimate[1], mu=nbinom.est$estimate[2], log = FALSE)
  goodness.nbinom<-mean(abs(dataframe$degree-trunc(x.nbinom)))
  #' 
  #' #Gamma distribution
  #' #Error in stats::optim(x = c(97, 63, 24, 84, 0, 2, 8, 100, 1, 0, 53, 1,  : 
  #' #initial value in 'vmmin' is not finite
  #' fitdistr(dataframe$degree,"gamma") 
  #' x.gamma = rgamma(n=nrows,scale=0.83,shape=10.59) 
  #' hist(x.gamma) 
  #' qqplot(dataframe$degree,x.gamma) 
  #' #Beta distribution
  #' #Error in ebeta(dataframe$degree, method = "mle") : 
  #' #'x' must contain at least 2 non-missing distinct values, 
  #' #'and all non-missing values of 'x' must be between 0 and 1.
  #' ebeta(dataframe$degree, method = "mle")
  #' fitdistr(dataframe$degree,"beta",start=list(shape1=0.5,shape2=0.5)) 
  #' x.beta<-dbeta(dataframe$degree, shape1, shape2, ncp = 0, log = FALSE)
  #' 
  #' hist(x.beta) 
  #' qqplot(dataframe$degree,x.beta) 
  #' #Log-normal distribution
  #' #error: unsupported distribution
  #' fitdistr(dataframe$degree,"lnorm") 
  #' x.lognormal<-dlnorm(dataframe$degree, meanlog = 0, sdlog = 1, log = FALSE)
  #' #Weibul distribution
  #' #error, weibull values must be > 0
  #' fitdistr(dataframe$degree,"weibull") 
  #' x.weibull = rweibull(n=nrows,scale=3.5,shape=14.1) 
  #' #hist(x.weibull) 
  #' #qqplot(dataframe$degree,x.weibull) 
  tmin<-match(min(distr.goodness), distr.goodness)
  {if(tmin==1){return(x.normal)}
    else if(tmin==2){return(x.poisson)}
    else
    {return(x.nbinom)}
  }
}
degree.good <- degree.distr.fitter(dataframe = degree.df$degree.df)