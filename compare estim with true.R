plot(res$logl,type='l')

plot(res$phi,type='l')
abline(h=phi.true,col='red')

k=data.frame(estim=res$z[ngibbs,],true1=z.true)
table(k)

ordem=c(6,2,4,3,5,1)
table(k)[ordem,]

n.ac=10
ac.coord=matrix(res$coord[ngibbs,],n.ac,2)
ac.coord[ordem,]
ac.coord.true

rango=range(c(ac.coord.true),ac.coord[ordem,])
plot(unlist(ac.coord.true),ac.coord[ordem,],xlim=rango,ylim=rango)
lines(rango,rango,col='red')

plot(ac.coord[1:6,1],ac.coord[1:6,2])