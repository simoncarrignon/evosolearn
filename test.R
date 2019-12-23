
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
epsilon=c(x=10,y=1,z=1)
sigma=c(s=2,y=2,z=2)
test=simpleEvoModel(100,500,omega = 2,delta = 4 ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma)


par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
plot(test$meanf,type="l",ylab="meanf")
plot(test$env,type="l",col="red",ylab="environment")
plot(test$popsize,type="l",col="blue",ylab="popsize")
dev.new()
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
sapply(genes,function(g)plot(sapply(test$allpop,function(i)mean(i[[g]])),ylab=paste("gene",g),type="l",ylim=c(0,1)))

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
sapply(genes,function(g)plot(sapply(test$allpop,function(i)mean(i[[g]])),ylab=paste("gene",g),type="l",ylim=c(0,1)))
dev.off()

pdf("env_fit_pop.pdf" ,width=10,height=6)
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
plot(test$meanf,type="l",ylab="meanf")
plot(test$env,type="l",col="red",ylab="environment")
plot(test$popsize,type="l",col="blue",ylab="popsize")
dev.off()
