diagnos = function(Y, E, XY) {	
	library(vegan)

	Z = pcnm(dist(XY))$vectors
	E2 = E^2
	X = cbind(E,E2)

	# E fraction
	rda_x = rda(Y~X, Z = Z, transfo = "hel" )
	Efrac = RsquareAdj(rda_x)$adj.r.squared

	# S fraction
	rda_s = rda(Y~Z, transfo = "hel" )
	Sfrac = RsquareAdj(rda_s)$adj.r.squared

	# Average absolute correlation
#	mean_cor = mean(abs(cor(Y)),na.rm=T)
	S = ncol(Y)
	cooc = t(Y)%*%Y / N
	P = apply(Y, 2, sum) / 100
	exp_cooc = matrix(P,nr=S,nc=S)*t(matrix(P,nr=S,nc=S))
	rel_cooc = log(cooc/exp_cooc)
	rel_cooc[rel_cooc==-Inf] = NA
	mean_cooc = mean(abs(rel_cooc),na.rm=T)

	c(Efrac,Sfrac,mean_cooc)

}

XY = read.table("input/XY.txt")[1:N,]
E = as.matrix(read.table("input/E.txt")[1:N,])

# Metapop
load("simulations/Metapop_run.Rds")
Y = run[[1000]]
diagnos(Y, E, XY)

load("simulations/Sorting_run.Rds")
Y = run[[1000]]
diagnos(Y, E, XY)

# NEUTRAL
load("simulations/Neutral_run.Rds")
Y = run[[1000]]
diagnos(Y, E, XY)

# COMPET
load("simulations/Compet_run.Rds")
Y = run[[1000]]
diagnos(Y, E, XY)

