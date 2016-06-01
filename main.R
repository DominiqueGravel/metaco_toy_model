# model: metacommunity paradigm. Options: Full, N, CC, SS, ME 
# pars: matrix of parameters
# adjMat: adjacency matrix ()
# Env: matrix of environmental variables
# presMat: initial presence-absence matrix
# nsteps: number of steps to run the model

# metaco model
metaco_model = function(model,pars,adjMat,Env,presMat,nsteps) {
  
  dat <- makeMat(as.data.frame(pars),Env,model=model)
  FMat=as.matrix(dat[[1]])
  SMat=as.matrix(dat[[2]])
  DMat=as.matrix(dat[[3]])
  
  ##################
  # Matrix in which we record the regional abundance over time  
  Series = matrix(nr=nsteps,nc=S+1)
  
  ############################################  
  # Loop over all time steps
  for(time in 1:nsteps) {   

#    paste(time)
    ######################
    # Test if there is colonization
    randCol = matrix(runif(N*S,0,1),nr=N,nc=S)
    ColMat = matrix(0,nr=N,nc=S)
    
    dat.compcol <- colcompMat(presMat,adjMat,SMat,FMat,DMat,model=model)
    IMat=as.matrix(dat.compcol[[1]])
    nMat=as.matrix(dat.compcol[[2]])
       
    # Perform the test
    ColMat[rowSums(presMat) == 0 & randCol < IMat] = 1    
    
    # Select the best colonizer among the colonists
    if(model=="ME"){
      best = apply(ColMat,1,max)
      ColMat[ColMat == 1 & ColMat != best] = 0
    }else{  
      best = apply(ColMat*SMat,1,max)
      ColMat[ColMat == 1 & ColMat*SMat != best] = 0
    }
    
    # Random selection among ties
    if(length(ColMat[rowSums(ColMat)>1,])>S) ColMat[rowSums(ColMat)>1,] <- t(apply(ColMat[rowSums(ColMat)>1,],1,pick))
    if(length(ColMat[rowSums(ColMat)>1,])==S) ColMat[rowSums(ColMat)>1,] <- pick(ColMat[rowSums(ColMat)>1,])
    
    
    ####################
    # Test if there is extinction
    ExtMat = matrix(0,nr=N,nc=S)
    randExt = matrix(runif(N*S,0,1),nr=N,nc=S)
    
    # Perform the test
    ExtMat[presMat == 1 & randExt < DMat] = -1
    
    ####################
    # Test if there is competitive exclusion
    CompDispMat = matrix(0,nr=N,nc=S)
    CompColMat = matrix(0,nr=N,nc=S)   		
    randComp1 = matrix(runif(N,0,1),nr=N,nc=S)
    randComp2 = matrix(runif(N*S,0,1),nr=N,nc=S)   	
    
    # Colonization by competitive exclusion    		
    CompColMat[presMat == 0 & rowSums(presMat) == 1 & rowSums(ExtMat) != -1 & randComp1 < pars$gamma & randComp2 < IMat] = 1
    best = apply(CompColMat*SMat,1,max)
    CompColMat[CompColMat == 1 & CompColMat*SMat != best] = 0		
    
    # Random selection among ties
    if(length(CompColMat[rowSums(CompColMat)>1,])>S) CompColMat[rowSums(CompColMat)>1,] <- t(apply(CompColMat[rowSums(CompColMat)>1,],1,pick))
    if(length(CompColMat[rowSums(CompColMat)>1,])==S) CompColMat[rowSums(CompColMat)>1,] <- pick(CompColMat[rowSums(CompColMat)>1,])
    
    # Displacement by competition
    CompDispMat[presMat == 1 & rowSums(CompColMat) == 1] = -1
    
    ######################
    # Apply changes
    presMat = presMat + ColMat + ExtMat + CompColMat + CompDispMat
    
    # Record regional occupancy
    Series[time,]=c(time,apply(presMat,2,sum)/N)
    
  }
  
  return(list(presMat,Series))
  
}

