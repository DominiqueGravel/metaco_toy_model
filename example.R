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
nsteps = 1000
set.seed(1)
test = main(XY,E,Y0,pars,A,nsteps)

# Compute the number of species per time step
alpha_div = numeric(nsteps)
gamma_div = numeric(nsteps)

for(i in 1:nsteps) {
	alpha_div[i] = mean(apply(test[[i]],1,sum))
	p = apply(test[[i]],2,sum) / N
	persist = numeric(R)
	persist[p != 0] = 1
	gamma_div[i] = sum(persist)
}

dev.new(width = 8, height = 6)
par(mar = c(5,6,2,1))
plot(c(1:nsteps), alpha_div, type = "l", xlab = "Time", ylab = "Diversity", 
cex.lab = 1.25, cex.axis = 1.25)
title(main = "Local diversity")

dev.new(width = 8, height = 6)
par(mar = c(5,6,2,1))
plot(c(1:nsteps), gamma_div, type = "l", xlab = "Time", ylab = "Diversity", 
cex.lab = 1.25, cex.axis = 1.25)
title(main = "Regional diversity")

