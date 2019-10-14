plot(res$logl,type='l')

plot(res$phi,type='l')
abline(h=phi.true,col='red')

k=data.frame(estim=res$z[ngibbs,],true1=z.true)
table(k)

ordem=c(2,3,1)
ac.coord=matrix(res$coord[ngibbs,],n.ac,2)
ac.coord[ordem,]
ac.coord.true

