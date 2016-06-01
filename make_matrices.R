# Matrix creation
makeMat <- function(pars,Env,model){
  EMat = matrix(rep(Env,each = S),nr=N,nc=S,byrow=TRUE)
  UMat = matrix(1,nr = N, nc = S)
  UMat = matrix(1,nr = N, nc = S)
  
  # Fecundity matrix
  uFMat = matrix(rep(pars$uF,each = N),nr=N,nc=S)
  sFMat = matrix(rep(pars$sF,each = N),nr=N,nc=S)  
  MaxFMat = matrix(rep(pars$maxF,each = N),nr=N,nc=S)  
  
  # Survival matrix
  uSMat = matrix(rep(pars$uS,each = N),nr=N,nc=S)
  sSMat = matrix(rep(pars$sS,each = N),nr=N,nc=S)  
  MaxSMat = matrix(rep(pars$maxS,each = N),nr=N,nc=S)  
  
  # Disturbance matrix
  uDMat = matrix(rep(pars$uD,each = N),nr=N,nc=S)
  sDMat = matrix(rep(pars$sD,each = N),nr=N,nc=S)  
  MaxDMat = matrix(rep(pars$maxD,each = N),nr=N,nc=S)  
  
  if(model == "Full") {
    FMat = gauss(uFMat,sFMat,MaxFMat,EMat)
    SMat = gauss(uSMat,sSMat,MaxSMat,EMat)
    DMat = 1 - gauss(uDMat,sDMat,MaxDMat,EMat)  
  }
  
  # Neutral model
  else if(model == "N") {
    FMat = matrix(rep(mean(MaxFMat),N*S),N,S)
    SMat = FMat
    DMat = 1 - matrix(rep(mean(MaxDMat),N*S),N,S)
  }
  
  # Species-sorting
  else if(model == "SS") {
    FMat = gauss(uFMat,sFMat,apply(MaxFMat,1,mean),EMat)
    SMat = gauss(uSMat,sSMat,apply(MaxSMat,1,mean),EMat)
    DMat = 1 - gauss(uDMat,sDMat,apply(MaxDMat,1,mean),EMat)
  }
  
  # Mass effect
  else if(model == "ME") {
    FMat = gauss(uFMat,sFMat,apply(MaxFMat,1,mean),EMat)
    SMat = gauss(uSMat,sSMat,apply(MaxSMat,1,mean),EMat)
    DMat = 1 - gauss(uDMat,sDMat,apply(MaxDMat,1,mean),EMat)
  }
  
  # Competition-colonization
  else if(model == "CC") {
    FMat = gauss(apply(uFMat,1,mean),apply(sFMat,1,mean),MaxFMat,EMat)
    SMat = gauss(apply(uSMat,1,mean),apply(sSMat,1,mean),MaxSMat,EMat)
    DMat = 1 - gauss(apply(uDMat,1,mean),apply(sDMat,1,mean),MaxDMat,EMat)
  }
  return(list(FMat,SMat,DMat))
}

colcompMat <- function(presMat,adjMat,SMat,FMat,DMat,model){
  IMat = matrix(nr=N,nc=S)
  if(model == "Full") {
    for(i in 1:S) {    
      matx = matrix(rep(SMat[,i], each=N),nr = N, nc = N, byrow=TRUE)  
      maty = matrix(rep(FMat[,i]*presMat[,i]*A,each=N),nr = N, nc = N, byrow=FALSE)  
      IMat[,i] = 1 - apply(1-adjMat*matx*maty,1,prod)  
    }  	
    nMat = pars$gamma*IMat/apply(IMat,1,sum)	
  }
  
  # Neutral model
  else if(model == "N") {
    for(i in 1:S) {  	
      matx = matrix(rep(SMat[,i],each=N),nr = N, nc = N, byrow=TRUE)	
      maty = matrix(rep(FMat[,i]*presMat[,i]*A,each=N),nr = N, nc = N, byrow=FALSE)	
      IMat[,i] = 1 - apply(1-adjMat*matx*maty,1,prod)	
    }		
    nMat = pars$gamma*IMat/apply(IMat,1,sum)
  }
  
  # Species-sorting
  else if(model == "SS") {
    for(i in 1:S) {  	
      matx = matrix(rep(SMat[,i],each=N),nr = N, nc = N, byrow=TRUE)	
      maty = matrix(rep(FMat[,i],each=N),nr = N, nc = N, byrow=FALSE)	
      IMat[,i] = 1 - apply(1-adjMat*matx*maty,1,prod)	
    }		
    nMat = 0*IMat/apply(IMat,1,sum)		
  }
  
  # Mass effect
  else if(model == "ME") {
    for(i in 1:S) {  	
      maty = matrix(rep(presMat[,i]*A,each=N),nr = N, nc = N, byrow=FALSE)  
      IMat[,i] = 1 - apply(1-adjMat*maty,1,prod)	
    }		
    nMat = 0*IMat/apply(IMat,1,sum)						
  }
  
  # Competition-colonization
  else if(model == "CC") {
    for(i in 1:S) {  	
      matx = matrix(rep(SMat[,i],each=N),nr = N, nc = N, byrow=TRUE)	
      maty = matrix(rep(FMat[,i]*presMat[,i]*A,each=N),nr = N, nc = N, byrow=FALSE)	
      IMat[,i] = 1 - apply(1-adjMat*matx*maty,1,prod)	
    }		
    nMat = pars$gamma*IMat/apply(IMat,1,sum)		
  }
  nMat[nMat=="NaN"]=0
  return(list(IMat,nMat))
}
