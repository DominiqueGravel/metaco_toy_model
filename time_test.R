#####################
#
# Tests to check computation performance
# Dominique Gravel
# June 2nd, 2016
# 
####################
rm(list = ls())

source("functions.R")
source("parameters.R")
source("landscape.R")
source("main.R")

Nseq = 10^seq(1,4,length = 10)
res = numeric(10)

for(z in 1:10) {

	N = Nseq[z]

	# Landscape characteristics
	XY = get_XY(N)
	E = get_E(D,N)

	# Initial conditions
	Y0 = matrix(0, nr = N, nc = R)
	rand = matrix(runif(N*R), nr = N, nc = R)
	Y0[rand < 0.5] = 1

	nsteps = 1000
	res[10] = system.time(main(XY,E,Y0,pars,A,nsteps))[3]

	cat(N,'\n')
}
dev.new(width = 8, height = 6)
par(mar = c(5,6,2,1))
plot(Nseq, res, xlab = "Number of spatial units", ylab = "Time required", 
cex.lab = 1.25, cex.axis = 1.25)

dev.copy2pdf(file = "Performance_report.pdf")