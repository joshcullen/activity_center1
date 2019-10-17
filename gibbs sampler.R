gibbs.activity.center=function(dat,grid.coord,n.ac,ac.coord.init,gamma1){
  #basic setup
  n.tsegm=nrow(dat)
  n.grid=nrow(grid.coord)
  grid.coord=data.matrix(grid.coord)
  
  #initial values
  ac.coord=data.matrix(ac.coord.init)
  z=sample(1:n.ac,size=n.tsegm,replace=T) #cluster membership
  phi=0.0001 #distance decay parameter
  theta=rep(1/n.ac,n.ac)
  
  #matrices to store results
  store.coord=matrix(NA,ngibbs,n.ac*2)
  store.z=matrix(NA,ngibbs,n.tsegm)
  store.param=matrix(NA,ngibbs,1) #to store phi
  store.logl=matrix(NA,ngibbs,1)
  store.theta=matrix(NA,ngibbs,n.ac)
  
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
               n.grid=n.grid,n.ac=n.ac,n.tsegm=n.tsegm,dat=dat,
               log.theta=log(theta))
    # z=z.true
    
    #sample theta
    v=sample.v(z=z,n.ac=n.ac,gamma1=gamma1)
    theta=rep(NA,n.ac)
    theta[1]=v[1]
    tmp=(1-v[1])
    for (j in 2:n.ac){
      theta[j]=v[j]*tmp
      tmp=tmp*(1-v[j])
    } 
      
    if (i<nburn & i%%adaptMH==0){
      #adapt MH
      tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=adaptMH)
      jump1=tmp$jump1
      if (jump1$phi<0.01) jump1$phi=0.01
      cond=jump1$coord<2; jump1$coord[cond]=2
      accept1=tmp$accept1
      
      #re-order data from time to time according to theta (largest to smallest)
      ordem=order(theta,decreasing=T)
      znew=rep(NA,n.tsegm)
      ac.coord=ac.coord[ordem,]
      for (j in 1:n.ac){
        cond=z==ordem[j]
        if (sum(cond)>0) znew[cond]=j
      }
      z=znew 
      theta=theta[ordem]
    }
    
    #store results
    store.coord[i,]=unlist(ac.coord)
    store.z[i,]=z
    store.param[i,]=phi
    store.logl[i,]=logl
    store.theta[i,]=theta
  }
  list(coord=store.coord,z=store.z,phi=store.param,logl=store.logl,theta=store.theta)  
}
