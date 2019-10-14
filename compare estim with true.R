plot(store.param,type='l')
abline(h=phi.true,col='red')

k=data.frame(estim=z,true1=z.true)
table(k)

ordem=c(1,3,2)
ac.coord[ordem,]
ac.coord.true

