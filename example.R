#-----------------------------------
### Example ---------------------
#-----------------------------------
# Parameter initialisation
N = 100
S = 5

# Landscape
A = numeric(N) + 1
land <- make.land(N=N,aggr=1,r=1000,model="Gradient", plot=T)
XY = land[[1]] # spatial coordinates
Env = land[[2]]
distMat = land[[3]]
adjMat = land[[4]]

# Species presence matrix
presMat = matrix(1,nc = S, nr = N)
presMat <- t(apply(presMat,1,pick))

# Species parameters (here following the SS or ME paradigms)
pars <- matrix(NA,S,10)
colnames(pars) <- c("uF","sF","maxF","uS","sS","maxS","uD","sD","maxD","gamma")
pars[,1] <- sample(1:100,S)
pars[,3] <- rep(0.9,S)
pars[,2] <- rep(30,S)
pars[,4] <- pars[,1]
pars[,6] <- pars[,3]
pars[,5] <- pars[,2]
pars[,7] <- rep(30,S)
pars[,9] <- rep(0.9,S)
pars[,8] <- pars[,2]
pars[,10] <- rep(0.1,S)
(pars=as.data.frame(pars))

resMat <- metaco_model(model="Full",  pars=pars, Env=Env, adjMat=adjMat, presMat=presMat, nsteps=1000)