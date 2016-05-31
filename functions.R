### General functions
# gaussian
gauss = function(u,s,Max,E) Max*exp(-(u-E)^2/2/s^2)

# select among ties
pick <- function(x){
  x <- c(x)
  if(length(x[x!=0])>1) x[x!=0]=sample(c(1,rep(0,length(x[x!=0])-1)),length(x[x!=0]),replace=F)
  return(x)
} 

### Generate random spatial coordinates for each patch
land.aggr.rand <- function(aggr,N=N){
  
  # assess the aggregation of patches
  dis <- rep(aggr, N%/%aggr)
  if(sum(dis)!=N) dis <- c(dis, N-sum(dis)) 
  
  # select coordinates of patch aggregates
  XY <- data.frame(x=seq(50,950,10),y=seq(50,950,10))
  XY.pairs <- expand.grid(XY)
  XY.select <- XY.pairs[sample(1:nrow(XY.pairs),length(dis)),]
  
  # patch coordinates
  coord.mat <- 0
  dat <- seq(-50,50,1)
  x<-sample(dat,dis[1],replace=F) + XY.select[1,1]
  y<-sample(dat,dis[1],replace=F) + XY.select[1,2]
  coord.mat <- cbind(x,y)  
  for(i in 2:length(dis)){
    x<-sample(dat,dis[i],replace=F) + XY.select[i,1]
    y<-sample(dat,dis[i],replace=F) + XY.select[i,2]
    coord.mat <- rbind(coord.mat,cbind(x,y))
  }
  return(coord.mat)
}

### Generate random landscape characteristics
make.land <- function(N, aggr, r, model="Random") {
  if(model!="Random" & model!="Gradient" & model!="Patchy" & model!="Homogeneous") stop("model should be one of Random, Gradient, Patchy, Homogeneous")
  
  XY = land.aggr.rand(aggr=aggr,N=N) # spatial coordinates
  distMat = as.matrix(dist(XY, method = 'euclidean', upper = T, diag = T))
  if(model=="Random"){
    Env <- sample(1:N,N,replace=F)
  }else if(model=="Gradient"){
    Env <- c(1:N)[rank(distMat[,1])]
  }else if(model=="Patchy"){
    clus <- hclust(dist(XY, method = 'euclidean', upper = T, diag = T),"ward")
    clus <- cutree(clus, N%/%10)
    Env <- round(rank(clus,ties.method="average"),0)
  }else if(model=="Homogeneous"){
    Env <- rep(50,N)
  }
  diag(distMat) <- 10000
  adjMat = matrix(0, nr = N, nc = N)
  adjMat[distMat < r] = 1
  diag(adjMat) = 0
  
  MND.dis <- mean(apply(distMat,2,min))
  
  return(list(XY,Env,distMat,adjMat,MND.dis))
}

