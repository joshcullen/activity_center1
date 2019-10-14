# rm(list=ls(all=TRUE))
set.seed(10)

setwd('U:\\GIT_models\\activity_center1')
library('Rcpp')
sourceCpp('aux1.cpp')
source('gibbs functions.R')

dat=read.csv('fake data.csv',as.is=T)
grid.coord=read.csv('fake data grid.csv',as.is=T)

#basic setup
n.tsegm=nrow(dat)
n.grid=nrow(grid.coord)
n.ac=3

#initial values
ac.coord=data.frame(x=c(5,50,90),
                    y=c(5,50,90))

#cluster membership
z=sample(1:n.ac,size=n.tsegm,replace=T)

#distance decay parameter
phi=0.0001

#matrices to store results
ngibbs=1000
nburn=ngibbs/2
store.coord=matrix(NA,ngibbs,n.ac*2)
store.z=matrix(NA,ngibbs,n.tsegm)
store.param=matrix(NA,ngibbs,1) #to store phi

#MH stuff
adaptMH=50
jump1=list(coord=matrix(10,n.ac,2),phi=0.2)
accept1=list(coord=matrix(0,n.ac,2),phi=0)

#gibbs sampler
for (i in 1:ngibbs){
  print(i)
  
  #sample coordinates
  tmp=sample.coord(ac.coord=ac.coord,jump=jump1$coord,dat=dat,
                   grid.coord=grid.coord,z=z,n.ac=n.ac)
  ac.coord=tmp$ac.coord
  accept1$coord=accept1$coord+tmp$accept
  # ac.coord=ac.coord.true
  
  #sample phi
  tmp=sample.phi(ac.coord=ac.coord,grid.coord=grid.coord,n.grid=n.grid,
                 n.ac=n.ac,phi=phi,jump=jump1$phi,dat=dat,z=z)
  phi=tmp$phi
  accept1$phi=accept1$phi+tmp$accept
  # phi=phi.true
  
  #sample z
  z=sample.z(ac.coord=ac.coord,grid.coord=grid.coord,
             n.grid=n.grid,n.ac=n.ac,n.tsegm=n.tsegm,dat=dat)
  # z=z.true
  
  #adapt MH
  if (i<nburn & i%%adaptMH==0){
    tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=adaptMH)
    jump1=tmp$jump1
    accept1=tmp$accept1
  }
    
  #store results
  store.coord[i,]=unlist(ac.coord)
  store.z[i,]=z
  store.param[i,]=phi
}
