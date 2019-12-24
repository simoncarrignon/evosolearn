source("protomodels.R")

pop=cbind.data.frame(agent=1:100,parent=1:100,x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1),fitness=1:100)

env=c()
env$theta=environment(tstep,omega,delta)
env$b=b
env$K=K
a=1:n
pop=cbind.data.frame(agent=a,parent=a,x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1),fitness=rep(0,n))
parents=pop
meanf=c()
popsize=c()

socialLearning(childs,parents,type="average") - mean(parents$slp)

t=sample(tstep,1)
socialLearning(childs,parents,type="best",theta=env$theta[t]) - parents$ilp[which.min(abs(parents$ilp-env$theta[t]))] #should be 0

length(socialLearning(childs,parents,type="random"))-nrow(childs) #should be equal to zero



##Running test
epsilon=c(x=1,y=1,z=1)
sigma=c(s=1,y=1,z=1)
test=simpleEvoModel(100,500,omega = 2,delta = 4 ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma)

genes=c("x","y","z")

png("images/env_fit_pop.png",width=800,height=400)
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
plot(test$meanf,type="l",ylab="meanf")
plot(test$env,type="l",col="red",ylab="environment")
plot(test$popsize,type="l",col="blue",ylab="popsize")
dev.off()
dev.new()
png("images/allgenes.png",width=800,height=400)
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
sapply(genes,function(g)plot(sapply(test$allpop,function(i)mean(i[[g]])),ylab=paste("gene",g),type="l"))
dev.off()

dev.off()

cols=colorRampPalette(c("blue","yellow","red"))(500)[as.numeric(cut(test$meanf,breaks=500))]
plot3d(test$meanf,test$env,col=cols,pch=20)
plot(test$meanf,test$env,col=colorRampPalette(c("blue","yellow","red"))(5000),pch=20)

plot(test$meanf[1100:1300],type="l",ylab="meanf")
lapply(genes,function(g)sapply(test$allpop,function(i)c(sd=sd(i[[g]]),mean=mean(i[[g]]))))

pdf("allgenes.pdf" ,width=10,height=6)
genes=c("x","y","z")
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
sapply(genes,function(g)plot(sapply(test$allpop,function(i)mean(i[[g]])),ylab=paste("gene",g),type="l"))
dev.off()

pdf("env_fit_pop.pdf" ,width=10,height=6)
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
plot(test$meanf,type="l",ylab="meanf")
plot(test$env,type="l",col="red",ylab="environment")
plot(test$popsize,type="l",col="blue",ylab="popsize")
dev.off()


omegas=seq(0,3,.5)
allos_best=sapply(omegas,function(o)replicate(100,mean(simpleEvoModel(100,300,omega = o,delta = 2 ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,log=T)$pop$z)))
png("omegas_vs_z.png")
boxplot(allos_best,ylab="mean value of z",xlab=expression(omega),axes=F) 
axis(2)
axis(1,1:length(omegas),label = omegas)
box()
dev.off()

allos_rand=sapply(omegas,function(o)replicate(100,mean(simpleEvoModel(100,300,omega = o,delta = 2 ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,type="random",log=T)$pop$z)))
png("omegas_vs_z_random.png")
boxplot(allos_rand,ylab="mean value of z",xlab=expression(omega),axes=F) 
axis(2)
axis(1,1:length(omegas),label = omegas)
box()
dev.off()

library(rgl)
plot3d(test$meanf,test$env,col=cols,pch=20)

omegas=seq(-0.5,2.5,.5)
deltas=.0625*2^seq(0:6)

library(parallel)
names(omegas)=omegas
names(deltas)=deltas
cl <- makeForkCluster(4,outfile="")
osnds=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));mean(replicate(50,mean(simpleEvoModel(100,300,omega = o,delta = d ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,log=F)$pop$z,na.rm=T)))}))
osnds_rand=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));mean(replicate(50,mean(simpleEvoModel(100,300,type="random",omega = o,delta = d ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,log=F)$pop$z,na.rm=T)))}))
stopCluster(cl)

png("nonrand.png")
image(log(deltas),omegas,osnds,zlim=c(0,1))
dev.off()
png("rand.png")
image(log(deltas),omegas,osnds_rand,zlim=c(0,1))
dev.off()
