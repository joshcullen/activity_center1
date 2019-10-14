gibbs.activity.center=function(dat,grid.coord,n.ac,ac.coord.init){
  #basic setup
  n.tsegm=nrow(dat)
  n.grid=nrow(grid.coord)

  #initial values
  ac.coord=ac.coord.init
  
  #cluster membership
  z=sample(1:n.ac,size=n.tsegm,replace=T)
  
  #distance decay parameter
  phi=0.0001
  
  #matrices to store results
  store.coord=matrix(NA,ngibbs,n.ac*2)
  store.z=matrix(NA,ngibbs,n.tsegm)
  store.param=matrix(NA,ngibbs,1) #to store phi
  store.logl=matrix(NA,ngibbs,1)
  
  #MH stuff
  adaptMH=50
  jump1=list(coord=matrix(10,n.ac,2),phi=0.2)
  accept1=list(coord=matrix(0,n.ac,2),phi=0)
  
  #gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    
    #sample coordinates
    tmp=sample.coord(ac.coord=ac.coord,jump=jump1$coord,dat=dat,
                     grid.coord=grid.coord,z=z,n.ac=n.ac,n.grid=n.grid,phi=phi)
    ac.coord=tmp$ac.coord
    accept1$coord=accept1$coord+tmp$accept
    # ac.coord=ac.coord.true
    
    #sample phi
    tmp=sample.phi(ac.coord=ac.coord,grid.coord=grid.coord,n.grid=n.grid,
                   n.ac=n.ac,phi=phi,jump=jump1$phi,dat=dat,z=z)
    phi=tmp$phi
    accept1$phi=accept1$phi+tmp$accept
    logl=tmp$logl
    # phi=phi.true
    
    #sample z
    z=sample.z(ac.coord=ac.coord,grid.coord=grid.coord,phi=phi,
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
    store.logl[i,]=logl
  }
  list(coord=store.coord,z=store.z,phi=store.param,logl=store.logl)  
}
