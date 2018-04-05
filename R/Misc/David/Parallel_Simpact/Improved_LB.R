# Load balancing in function 'dynamicClusterApply' is achieved by the line 'if (j <= n) submit(d$node, j)'
# Insert a tracer to check whether j is ever smaller or equal than n
trace(parallel:::dynamicClusterApply, tracer = quote(print(paste(j, n, d$node))), at = list(c(4, 3, 5, 4, 4)))


# My suggested fix to actually use load balancing
#	- if chunksize = length(X) / length(cl), then equivalent to parallel::parLapplyLB
parLapplyLBc <- function (cl = NULL, X, fun, ..., chunksize = 1L) {
	# modified from parallel::parLapplyLB  (v3.2.4)
    cl <- parallel:::defaultCluster(cl)

    # x <- splitList(X, length(cl)) # this is the original code, but if length(x) == length(cl), then 'dynamicClusterApply' does not load balance the tasks
    # instead the modified version creates a list of indices longer than the number of workers so that load balancing can occur
    x <- if (chunksize == 1L) X else parallel:::splitList(X, max(length(cl), round(length(X) / chunksize)))

    do.call(c, parallel::clusterApplyLB(cl, x = x, fun = lapply, fun, ...), quote = TRUE)
}



cl <- parallel::makeCluster(2, type = "SOCK", outfile = "testLB.txt")

system.time(parallel::parLapply(cl, 1:10, Sys.sleep))
	#   user  system elapsed
	#  0.490   0.109  40.049

system.time(parallel::parLapplyLB(cl, 1:10, Sys.sleep))
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "3 2 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "4 2 2"
	#   user  system elapsed
	#  0.444   0.097  40.043
# ==> load balancing functionality is never used; speed of 'parLapplyLB' is equal to 'parLapply'

system.time(parLapplyLBc(cl, 1:10, Sys.sleep, chunksize = 1))
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "3 10 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "4 10 2"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "5 10 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "6 10 2"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "7 10 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "8 10 2"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "9 10 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "10 10 2"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "11 10 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "12 10 2"
	#   user  system elapsed
	#  0.362   0.080  30.046
# ==> load balancing is used; speed of 'parLapplyLBc' is in this example 25% faster than 'parLapplyLB'

system.time(parLapplyLBc(cl, 1:10, Sys.sleep, chunksize = 2))
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "3 5 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "4 5 2"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "5 5 1"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "6 5 2"
	#Tracing dynamicClusterApply(cl, fun, length(x), argfun) step 4,3,5,4,4
	#[1] "7 5 1"
	#   user  system elapsed
	#  0.389   0.083  33.057


# Clean-up
parallel::stopCluster(cl)
untrace(parallel:::dynamicClusterApply)

# System information
# ==> I am aware that I am not running R version 3.2.4, but I checked that this issue hasn't been addressed with R-devel_2016-03-30.tar.gz
sessionInfo()
	#R version 3.2.3 (2015-12-10)
	#Platform: x86_64-apple-darwin14.5.0 (64-bit)
	#Running under: OS X 10.10.5 (Yosemite)
	#
	#locale:
	#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
	#
	#attached base packages:
	#[1] stats     graphics  grDevices
	#[4] utils     datasets  methods
	#[7] base
	#
	#loaded via a namespace (and not attached):
	#[1] snow_0.4-1     tools_3.2.3
	#[3] parallel_3.2.3
