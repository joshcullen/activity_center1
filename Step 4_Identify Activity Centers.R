set.seed(1)

#load libraries and read important functions
library('Rcpp')
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(viridis)
library(progress)
library(rnaturalearth)
library(rnaturalearthdata)


sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')
source('helper functions.R')


#load data
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
obs<- get.summary.stats_obs(dat)

utm.crs<- CRS('+init=epsg:32617')
extent<- extent(min(dat$utmlong), max(dat$utmlong), min(dat$utmlat), max(dat$utmlat))
res<- 5000
buffer<- 10000

grid.coord<- grid.summary.table(dat=dat, crs=utm.crs, extent=extent, res=res, buffer=buffer)
dat<- left_join(dat, grid.coord, by="grid.cell") #add gridded locs to DF


#Define initial activity centers (obs > 100)
tmp<- which(apply(obs[,-1], 2, sum) > 100)
nobs<- colSums(obs[,tmp+1])
ac.coord.init<- grid.coord[tmp,-3]
ac.coord.init<- cbind(ac.coord.init, nobs)

#Lowest 20 of 37 initial ACs
ac.coord.init2<- arrange(ac.coord.init, nobs) %>% slice(n=1:20)


#########################
### Run Gibbs sampler ###
#########################

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=nrow(ac.coord.init2)
gamma1=0.1


#run gibbs sampler
pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)

res=gibbs.activity.center(dat=obs[,-1], grid.coord=grid.coord[,-3], n.ac=n.ac,
                          ac.coord.init=ac.coord.init2[,-3], gamma1=gamma1)

#plot output and look at frequency of AC visitation
plot(res$logl,type='l')
plot(res$phi,type='l')

ac<- data.frame(ac=res$z[ngibbs,])
table(ac$ac)
obs<- cbind(ac, obs)

############################
### Add ACs to Dataframe ###
############################

obs1<- obs %>% filter(id == 1) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)
obs12<- obs %>% filter(id == 12) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)
obs19<- obs %>% filter(id == 19) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)
obs27<- obs %>% filter(id == 27) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)

dat1<- dat %>% filter(id == 1) %>% left_join(obs1, by = "time.seg")
dat12<- dat %>% filter(id == 12) %>% left_join(obs12, by = "time.seg")
dat19<- dat %>% filter(id == 19) %>% left_join(obs19, by = "time.seg")
dat27<- dat %>% filter(id == 27) %>% left_join(obs27, by = "time.seg")

dat<- rbind(dat1, dat12, dat19, dat27)

#Calculate number of obs per AC
dat %>% filter(id==1) %>% dplyr::select(ac) %>% table()
dat %>% filter(id==12) %>% dplyr::select(ac) %>% table()
dat %>% filter(id==19) %>% dplyr::select(ac) %>% table()
dat %>% filter(id==27) %>% dplyr::select(ac) %>% table()



## Map

#Load world map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- sf::st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N

ac.coords<- matrix(NA, length(unique(ac$ac)), 2)
colnames(ac.coords)<- c("x","y")
tmp<- res$coord[ngibbs,]

for (i in 1:length(unique(ac$ac))) {
  ac.coords[i,]<- round(c(tmp[i], tmp[i+n.ac]), 0)
} 

ac.coords<- data.frame(ac.coords, ac=1:length(unique(ac$ac)))

# ACs and initial values
ggplot() +
  geom_point(data = ac.coord.init, aes(x, y, color = "A (initial: n=37)"), size = 3) +
  geom_point(data = ac.coord.init, aes(x, y, color = "B (initial: n=20)"), size = 3, pch = 1) +
  geom_point(data = ac.coords, aes(x, y, color = "C (model: n=20)"), size = 3, alpha = 0.5) +
  labs(x="Easting", y="Northing") +
  scale_color_manual("", values = c("grey40",viridis(n=5)[c(3,5)]))

# ACs and snail kite locs
ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_point(data = dat, aes(utmlong, utmlat, color=ac), size=0.5, alpha = 0.6) +
  geom_point(data = ac.coords, aes(x, y, color = ac), size = 4, pch = 1, stroke = 1) +
  scale_color_viridis_c() +
  labs(x = "Easting", y = "Northing") +
  theme_bw()


###################
### Save Output ###
###################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/cluster_tsegments_loc")

# write.csv(ac.coords, "Activity Center Coordinates.csv", row.names = F)
# write.csv(dat, "Snail Kite Gridded Data_AC.csv", row.names = F)

write.csv(ac.coords, "Activity Center Coordinates_reduced.csv", row.names = F)
