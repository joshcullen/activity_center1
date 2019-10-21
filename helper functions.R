assign.time.seg<- function(obs, breakpts, dat) {
  seg.n=apply(obs,1,sum)
  time.seg=rep(1:(length(breakpts)+1), seg.n)
  dat$time.seg<- time.seg
  
  dat
}
#------------------------------------------------
ind.func=function(x) {which(x==1)}
#------------------------------------------------
get.summary.stats=function(breakpt,dat,nloc){
  col.time1=which(colnames(dat)=='time1')
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  res=matrix(NA,n-1,nloc)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    tmp=dat[ind,-col.time1]
    
    ind.loc<- matrix(NA, nrow(tmp), 1)
    ind.loc<- apply(tmp, 1, FUN = ind.func)  #creates vector of locations occupied per observation
    
    tab<- tabulate(ind.loc)  #ensures that there are values for all nloc
    if (length(tab) < nloc) {
      tab<- c(tab,rep(0,nloc-length(tab)))
    }
    
    res[i-1,]=tab #takes count of the each location within given time segment
  }
  res
}
#------------------------------------------------
grid.summary.table=function(dat,crs){  #dat must already have time.seg assigned
  
  #create grid and extract coords per cell
  grid_5<- raster(extent(min(dat$utmlong), max(dat$utmlong),
                         min(dat$utmlat), max(dat$utmlat)) + 10000)
  res(grid_5)<- 5000
  proj4string(grid_5)<- crs
  grid_5[]<- 0
  grid.cell.locs<- coordinates(grid_5) %>% data.frame()
  names(grid.cell.locs)<- c("x", "y")
  grid.cell.locs$grid.cell<- 1:length(grid_5)
  grid.coord<- grid.cell.locs[grid.cell.locs$grid.cell %in% dat$grid.cell,]
  
  
  grid.coord
}
#------------------------------------------------
