sample.coord=function(ac.coord,jump,dat,grid.coord,z,n.ac,n.grid,phi){
  ac.coord.orig=ac.coord.old=ac.coord
  tmp=rnorm(n.ac*2,mean=ac.coord.old,sd=jump)
  ac.coord.prop=matrix(tmp,n.ac,2)
  
  for (i in 1:n.ac){
    cond=z==i
    n=sum(cond)
    dat1=dat[cond,]
    
    for (j in 1:2){ #to sample each coordinate separately
      ac.coord.new=ac.coord.old
      ac.coord.new[i,j]=ac.coord.prop[i,j]
        
      if (n>0){
        #get distances
        x2=(grid.coord[,'x']-ac.coord.old[i,'x'])^2
        y2=(grid.coord[,'y']-ac.coord.old[i,'y'])^2
        dist1.old=sqrt(x2+y2)
        
        x2=(grid.coord[,'x']-ac.coord.new[i,'x'])^2
        y2=(grid.coord[,'y']-ac.coord.new[i,'y'])^2
        dist1.new=sqrt(x2+y2)
          
        #get loglikel
        tmp=exp(-phi*dist1.old)
        tmp1=tmp/sum(tmp)
        lprob.old=matrix(log(tmp1),n,n.grid,byrow=T)
        pold=sum(dat1*lprob.old)
          
        tmp=exp(-phi*dist1.new)
        tmp1=tmp/sum(tmp)
        lprob.new=matrix(log(tmp1),n,n.grid,byrow=T)
        pnew=sum(dat1*lprob.new)
          
        #accept or reject MH
        k=acceptMH(p0=pold,p1=pnew,x0=ac.coord.old[i,j],x1=ac.coord.new[i,j],BLOCK=F)      
        ac.coord.old[i,j]=k$x
      }
    }
  }
  list(ac.coord=ac.coord.old,accept=ac.coord.orig!=ac.coord.old)
}
#-----------------------------------
sample.phi=function(ac.coord,grid.coord,n.grid,n.ac,phi,jump,dat,z){
  dist1=GetDistance(AcCoord=ac.coord,GridCoord=grid.coord, 
                    Ngrid=n.grid, Nac=n.ac)
  old=phi
  new=abs(rnorm(1,mean=old,sd=jump)) #reflection proposal around zero
  
  #get loglikel
  pold=get.loglikel(dist1=dist1,dat=dat,z=z,n.ac=n.ac,phi=old,n.grid=n.grid)
  pnew=get.loglikel(dist1=dist1,dat=dat,z=z,n.ac=n.ac,phi=new,n.grid=n.grid)
  
  #accept or reject MH
  k=acceptMH(p0=pold,p1=pnew,x0=old,x1=new,BLOCK=F)  
  logl=ifelse(k$accept==1,pnew,pold)
  list(phi=k$x,accept=k$accept,logl=logl)
}
#-----------------------------------
sample.z=function(ac.coord,grid.coord,n.grid,n.ac,n.tsegm,dat,phi,log.theta){
  #get distance
  dist1=GetDistance(AcCoord=ac.coord,GridCoord=grid.coord, 
                    Ngrid=n.grid, Nac=n.ac)
  
  #get loglikel
  logl=matrix(NA,n.tsegm,n.ac)
  for (i in 1:n.ac){
    dist2=dist1[,i]
    prob=exp(-phi*dist2)
    prob1=matrix(prob/sum(prob),n.tsegm,n.grid,byrow=T)
    logl[,i]=rowSums(dat*log(prob1))+log.theta[i]    
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
# get.distance=function(ac.coord,grid.coord,n.grid,n.ac){
#   res=matrix(NA,n.grid,n.ac)
#   for (i in 1:n.ac){
#     x2=(grid.coord$x-ac.coord$x[i])^2
#     y2=(grid.coord$y-ac.coord$y[i])^2
#     res[,i]=sqrt(x2+y2)
#   }
#   res
# }
#-----------------------------------
get.loglikel=function(dist1,dat,z,n.ac,phi,n.grid){
  res=rep(0,n.ac)
  for (i in 1:n.ac){
    cond=z==i
    n=sum(cond)
    if (n!=0){
      dat1=dat[cond,]
      dist2=dist1[,i]
      prob=exp(-phi*dist2)
      prob1=matrix(prob/sum(prob),n,n.grid,byrow=T)
      res[i]=sum(dat1*log(prob1))
    }
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
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------
sample.v=function(z,n.ac,gamma1){
  tmp=table(z)
  n=rep(0,n.ac)
  n[as.numeric(names(tmp))]=tmp
  
  seq1=n.ac:1
  tmp=cumsum(n[seq1])
  n.ge=tmp[seq1]
  n.ge1=n.ge[-1]
  v=rbeta(n.ac-1,n[-n.ac]+1,n.ge1+gamma1)
  c(v,1)
}