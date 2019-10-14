rm(list=ls(all=TRUE))
set.seed(1)

#basic setup
n.tsegm=100
n.ac=3
n.grid=30
n=floor(runif(n.tsegm,min=0,max=900))

#spatial coordinates
grid.coord=data.frame(x=runif(n.grid,min=0,max=100),
                      y=runif(n.grid,min=0,max=100))
ac.coord.true=ac.coord=data.frame(x=runif(n.ac,min=0,max=100),
                                  y=runif(n.ac,min=0,max=100))

#cluster membership
z.true=z=sample(1:n.ac,size=n.tsegm,replace=T)

#distance decay parameter
phi.true=phi=0.1

#generate data
y=matrix(NA,n.tsegm,n.grid)
for (i in 1:n.tsegm){
  #calculate distance
  x2=(grid.coord$x-ac.coord$x[z[i]])^2
  y2=(grid.coord$y-ac.coord$y[z[i]])^2
  d=sqrt(x2+y2)
  
  #calculate probability
  tmp=exp(-phi*d)
  prob=tmp/sum(tmp)
  
  #draw results
  y[i,]=rmultinom(1,size=n[i],prob=prob)
}
image(y[z==3,])

setwd('U:\\GIT_models\\activity_center1')
write.csv(y,'fake data.csv',row.names=F)
write.csv(grid.coord,'fake data grid.csv',row.names=F)