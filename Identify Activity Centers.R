set.seed(1)

#load libraries and read important functions
library('Rcpp')
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)


sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')
source('helper functions.R')


#load data
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
#dat.list<- list(`1`=dat1, `12`=dat12, `19`=dat19, `27`=dat27)
obs<- read.csv("Occupancy Matrix for all Obs and Locs.csv", header = T, sep = ",")

obs1.seg<- read.csv("ID1 Seg x Loc.csv", header = T, sep = ',')
# obs12.seg<- read.csv("ID12 Seg x Loc.csv", header = T, sep = ',')
# obs19.seg<- read.csv("ID19 Seg x Loc.csv", header = T, sep = ',')
# obs27.seg<- read.csv("ID27 Seg x Loc.csv", header = T, sep = ',')


utm.crs<- CRS('+init=epsg:32617')
extent<- extent(min(dat$utmlong), max(dat$utmlong), min(dat$utmlat), max(dat$utmlat))
res<- 5000
buffer<- 10000

grid.coord<- grid.summary.table(dat=dat, crs=utm.crs, extent=extent, res=res, buffer=buffer)
dat<- left_join(dat, grid.coord, by="grid.cell") #add gridded locs to DF


#Define initial activity centers (obs > 100)
tmp<- which(apply(obs[,-1], 2, sum) > 100)
ac.coord.init<- grid.coord[tmp,-3]


#########################
### Run Gibbs sampler ###
#########################

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=nrow(ac.coord.init)
gamma1=0.1


#run gibbs sampler
options(warn=2)
res=gibbs.activity.center(dat=obs1.seg, grid.coord=grid.coord[,-3], n.ac=n.ac,
                          ac.coord.init=ac.coord.init, gamma1=gamma1)

#plot output and look at frequency of AC visitation
plot(res$logl,type='l')
plot(res$phi,type='l')

z=data.frame(z=res$z[ngibbs,], time.seg = 1:nrow(obs1.seg))
table(z$z)

## map

#Load world map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- sf::st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N

dat1<- dat %>% filter(id == 1)
dat1<- left_join(dat1, z, by = "time.seg")

ac.coords<- matrix(NA, 37, 2)
tmp<- res$coord[ngibbs,]

for (i in 1:37) {
  ac.coords[i,]<- c(tmp[i], tmp[i+37])
} 


ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_point(data = dat1, aes(x, y, color=z), size=1) +
  geom_point() +
  scale_color_viridis_c() +
  labs(x = "Easting", y = "Northing") +
  theme_bw()
