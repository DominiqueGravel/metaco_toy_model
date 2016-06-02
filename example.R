#####################
#
# Example of a single run of the model
# Dominique Gravel
# June 2nd, 2016
# 
####################
rm(list = ls())

source("functions.R")
source("parameters.R")
source("landscape.R")
source("main.R")

# Landscape characteristics
XY = get_XY(N)
E = get_E(D,N)

# Initial conditions
Y0 = matrix(0, nr = N, nc = R)
rand = matrix(runif(N*R), nr = N, nc = R)
Y0[rand < 0.5] = 1

# Run the model
test = main(XY,E,Y0,pars,A,nsteps)

# Compute the number of species per time step
gamma = numeric(nsteps)

for(i in 1:nsteps) {
	p = apply(test[[i]],2,sum) / N
	ispresent = numeric(R)
	ispresent[p!=0] = 1
	gamma[i] = sum(ispresent)
}




