#####################
#
# Wrapper to run the model for a given parameter set
# Dominique Gravel
# June 2nd, 2016
# 
####################
rm(list = ls())

run_metaco <- function(pars, E, XY) {

	source("functions.R")
	source("landscape.R")
	source("main.R")

	load(pars)
	E <- read.table(E)
	XY <- read.table(XY)

	# Initial conditions
	N <- nrow(E)
	R <- ncol(pars$A)
	Y0 = matrix(0, nr = N, nc = R)
	rand = matrix(runif(N*R), nr = N, nc = R)
	Y0[rand < 0.5] = 1

	# Run the model
	nsteps = 1000
	set.seed(1)
	res = main(XY,E,Y0,pars,A,nsteps)
	return(res[[nsteps]])
}

source("parameters.R")
Y1 <- run_metaco("metapop.Rds", "E.txt", "XY.txt")
Y2 <- run_metaco("sorting.Rds", "E.txt", "XY.txt")
Y3 <- run_metaco("sorting_interactions.Rds", "E.txt", "XY.txt")
Y4 <- run_metaco("neutral.Rds", "E.txt", "XY.txt")



