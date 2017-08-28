#########################
# Example 2
# Species sorting with dispersal limitations
#########################

setwd("~/Bureau/Documents/Workshops/metaco2.0/metaco_toy_model")
rm(list = ls())

# Load libraries
library(HMSC)
library(scatterplot3d)
library(corrplot)
library(viridis)
library(vcd)

scenario <- "Compet"

# Number of patches
N = 100

# Number of environmental variables
D = 1

# Regional species richness
R = 10

# Effect of the environment on colonization
u_c = matrix(nr = D, nc = R)
u_c[1,] = seq(0.1,0.9,length=R)
s_c = matrix(Inf, nr = D, nc = R)

# Effect of the environment on extinction
u_e = matrix(nr = D, nc = R)
u_e[1,] = rep(0.5, R)
s_e = matrix(Inf, nr = D, nc = R)

# Mean dispersal
alpha = Inf

# Immigration
m = 0.001

# Colonization function
c_0 = rep(0.1, R) # Colonization at 0 interactions
c_max = rep(1, R) # Colonization at max interactions

# Extinction function
e_0 = rep(0.05, R) # Extinction at 0 interactions
e_min = rep(0, R) # Exinction at max interactions

# Sensitivity to interactions
d_c = 1
d_e = 0

# Interaction matrix
A = matrix(0,nr=R,nc=R)
d = as.matrix(dist(c(1:R),upper=TRUE,diag=T))
A[d==1] = -1 
diag(A) = 0

# Collect all parameters into a single
pars = list(u_c = u_c, u_e = u_e, s_c = s_c, s_e = s_e, alpha = alpha, m = m, 
c_0 = c_0, e_0 = e_0, c_max = c_max, e_min = e_min, d_c = d_c, d_e = d_e, A = A)

save(pars, file = paste("pars/",scenario,"_pars.Rds", sep=""))

#=======================
# Run the model 
#=======================

source("core_model/functions.R")
source("core_model/landscape.R")
source("core_model/main.R")

# Landscape characteristics
#XY = get_XY(N)
XY = read.table("input/XY.txt")[1:N,]

#E = get_E(D,N)
E = as.matrix(read.table("input/E.txt")[1:N,])

# Initial conditions
Y0 = matrix(0, nr = N, nc = R)
rand = matrix(runif(N*R), nr = N, nc = R)
Y0[rand < 0.5] = 1

# Run the model
nsteps = 1000
set.seed(1)
run = main(XY,E,Y0,pars,A,nsteps)
save(run, file = paste("simulations/",scenario,"_run.Rds", sep=""))

# Compute the number of species per time step
alpha_div = numeric(nsteps)
gamma_div = numeric(nsteps)

for(i in 1:nsteps) {
	alpha_div[i] = mean(apply(run[[i]],1,sum))
	p = apply(run[[i]],2,sum) / N
	persist = numeric(R)
	persist[p != 0] = 1
	gamma_div[i] = sum(persist)
}

#dev.new(width = 8, height = 6)
#par(mar = c(5,6,2,1))
#plot(c(1:nsteps), alpha_div, type = "l", xlab = "Time", ylab = "Diversity", 
#cex.lab = 1.25, cex.axis = 1.25)
#title(main = "Local diversity")

#dev.new(width = 8, height = 6)
#par(mar = c(5,6,2,1))
#plot(c(1:nsteps), gamma_div, type = "l", xlab = "Time", ylab = "Diversity", 
#cex.lab = 1.25, cex.axis = 1.25)
#title(main = "Regional diversity")

Y = run[[nsteps]]
S = ncol(Y)
cooc = t(Y)%*%Y / N
P = apply(Y, 2, sum) / 100
P
exp_cooc = matrix(P,nr=S,nc=S)*t(matrix(P,nr=S,nc=S))
rel_cooc = log(cooc/exp_cooc)
rel_cooc[rel_cooc==-Inf] = NA
diag(rel_cooc) = 1
mean(rel_cooc[A==-1], na.rm = T)
mean(rel_cooc[A==0], na.rm = T)
mean(rel_cooc,na.rm=T)
sd(rel_cooc,na.rm=T)

#=======================
# Run the HMSC ANALYSIS
#=======================

# Prepare the data
siteID = as.factor(c(1:nrow(E)))
coord <- data.frame(siteID = siteID, coordX = XY[,1], coorY  = XY[,2])

# Fit the model with space as an autocorrelation and site as a random effect
#data_full <- as.HMSCdata(Y = Y, X = cbind(E,E^2), Random = siteID, Auto = coord)
data_full <- as.HMSCdata(Y = Y, X = cbind(E,E^2), Random = siteID)
model_full <- hmsc(data_full, family = "probit", niter = 10000, nburn = 1000, thin = 10)
 
# Test for convergence
mixingParamX <- as.mcmc(model_full, parameters = "paramX")
tmp<-geweke.diag(mixingParamX)
pvalue2sided<-2*pnorm(-abs(tmp$z))
#mixingParamLatent <- as.mcmc(model_full, parameters = "paramLatent")
#mixingMeansParamX <- as.mcmc(model_full, parameters = "meansParamX")
#plot(mixingParamLatent, col = "blue")

# Variance partitioning
VP <- variPart(model_full, groupX = c("Intercept","E", "E2"))
Efract = apply(VP[,2:3],1,sum)
Jfract = VP[,4]
Sfract = VP[,5]

# Compute the McFadden R2
pred = predict(model_full, type = "response")
ll = matrix(nr = N, nc = R)
ll[Y==1] = log(pred[Y==1])
ll[Y==0] = log(1-pred[Y==0])

occup = apply(Y,2,sum)
null = matrix(occup, nr = N, nc = S, byrow = TRUE)/N
null_ll = matrix(nr = N, nc = R)
null_ll[Y==1] = log(null[Y==1])
null_ll[Y==0] = log(1-null[Y==0])
R2 = 1 - sum(ll[,occup>0])/sum(null_ll[,occup>0])

save(model_full, file = paste("models/",scenario,"_hmsc_full.Rds", sep=""))

#=======================
# Plot the results of variation partitioning
#=======================

dev.new(width = 10, height = 10)
ternaryplot(
	cbind(Efract,Jfract,Sfract), 
	pch = 19, 
	bg = "lightgray",
	grid_color = "white",
	col = "darkblue",
	main = paste(scenario, " - R2 = ",R2, dec = ""),
	dimnames = c("Environment", "Spatial autocorrelation", "Co-distribution")
)
dev.copy2pdf(file = paste("figures/",scenario,".pdf", dec = ""))

# Plot the heat maps
pal <- viridis(200)
dev.new(width = 10, height = 10)
codistr_full <- apply(corRandomEff(model_full)$Random, 1:2, mean)
corrplot(codistr_full, method = "color", col = pal, type = "lower",
diag = FALSE, tl.srt = 45)
title(scenario, cex = 2)
dev.copy2pdf(file = paste("figures/",scenario,".pdf", dec = ""))

#dev.new(width = 10, height = 10)
#codistr_E0 <- apply(corRandomEff(model_E0)$Random, 1:2, mean)
#codistr_E0 <- cor(Y)
#corrplot(codistr_E0, method = "color", col = pal, type = "lower",
#diag = FALSE, tl.srt = 45)
#title(paste(scenario, "Not controlling for E"), cex = 2)
#dev.copy2pdf(file = paste("figures/",scenario,"_withoutE.pdf", dec = ""))