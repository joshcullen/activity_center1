sample.coord=function(ac.coord,jump,dat,grid.coord,z,n.ac){
  ac.coord.orig=ac.coord.old=ac.coord
  for (i in 1:n.ac){
    for (j in 1:2){ #to sample each coordinate separately
      ac.coord.new=ac.coord.old
      ac.coord.new[i,j]=rnorm(1,mean=ac.coord.old[i,j],sd=jump[i,j])
      
      #get distances
      dist1.old=get.distance(ac.coord=ac.coord.old,grid.coord=grid.coord,
                             n.grid=n.grid,n.ac=n.ac)
      dist1.new=get.distance(ac.coord=ac.coord.new,grid.coord=grid.coord,
                             n.grid=n.grid,n.ac=n.ac)

      #get loglikel
      pold=get.loglikel(dist1=dist1.old,dat=dat,z=z,n.ac=n.ac,phi=phi)
      pnew=get.loglikel(dist1=dist1.new,dat=dat,z=z,n.ac=n.ac,phi=phi)
      
      #accept or reject MH
      k=acceptMH(p0=pold,p1=pnew,x0=ac.coord.old[i,j],x1=ac.coord.new[i,j],BLOCK=F)      
      ac.coord.old[i,j]=k$x
    }
  }
  list(ac.coord=ac.coord.old,accept=ac.coord.orig!=ac.coord.old)
}
#-----------------------------------
sample.phi=function(ac.coord,grid.coord,n.grid,n.ac,phi,jump,dat,z){
  dist1=get.distance(ac.coord=ac.coord,grid.coord=grid.coord,n.grid=n.grid,n.ac=n.ac)
  old=phi
  new=abs(rnorm(1,mean=old,sd=jump)) #reflection proposal around zero
  
  #get loglikel
  pold=get.loglikel(dist1=dist1,dat=dat,z=z,n.ac=n.ac,phi=old)
  pnew=get.loglikel(dist1=dist1,dat=dat,z=z,n.ac=n.ac,phi=new)
  
  #accept or reject MH
  k=acceptMH(p0=pold,p1=pnew,x0=old,x1=new,BLOCK=F)      
  list(phi=k$x,accept=k$accept)
}
#-----------------------------------
sample.z=function(ac.coord,grid.coord,n.grid,n.ac,n.tsegm,dat){
  #get distance
  dist1=get.distance(ac.coord=ac.coord,grid.coord=grid.coord,n.grid=n.grid,n.ac=n.ac)
  
  #get loglikel
  logl=matrix(NA,n.tsegm,n.ac)
  for (i in 1:n.ac){
    dist2=dist1[,i]
    prob=exp(-phi*dist2)
    prob1=matrix(prob/sum(prob),n.tsegm,n.grid,byrow=T)
    logl[,i]=rowSums(dat*log(prob1))    
  }
  maximo=matrix(apply(logl,1,max),n.tsegm,n.ac)
  logl=logl-maximo
  tmp=exp(logl)
  soma=matrix(rowSums(tmp),n.tsegm,n.ac)
  prob=tmp/soma
  
  z=rmultinom1(prob=prob,randu=runif(n.tsegm))
  z+1
}
#-----------------------------------
get.distance=function(ac.coord,grid.coord,n.grid,n.ac){
  res=matrix(NA,n.grid,n.ac)
  for (i in 1:n.ac){
    x2=(grid.coord$x-ac.coord$x[i])^2
    y2=(grid.coord$y-ac.coord$y[i])^2
    res[,i]=sqrt(x2+y2)
  }
  res
}
#-----------------------------------
get.loglikel=function(dist1,dat,z,n.ac,phi){
  res=rep(NA,n.ac)
  for (i in 1:n.ac){
    cond=z==i
    dat1=dat[cond,]
    dist2=dist1[,i]
    prob=exp(-phi*dist2)
    prob1=matrix(prob/sum(prob),sum(cond),n.grid,byrow=T)
    res[i]=sum(dat1*log(prob1))
  }
  sum(res)
}
#---------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#----------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}