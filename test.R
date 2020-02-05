source("protomodels.R")

pop=cbind.data.frame(agent=1:100,parent=1:100,x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1),w=1:100)

env=c()
env$theta=environment(tstep,omega,delta)
env$b=b
env$K=K
a=1:n
pop=cbind.data.frame(agent=a,parent=a,x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1),w=rep(0,n))
parents=pop
meanw=c()
popsize=c()

socialLearning(childs,parents,type="average") - mean(parents$slp)

t=sample(tstep,1)
socialLearning(childs,parents,type="best",theta=env$theta[t]) - parents$ilp[which.min(abs(parents$ilp-env$theta[t]))] #should be 0

length(socialLearning(childs,parents,type="random"))-nrow(childs) #should be equal to zero



##Running test
E=c(x=1,y=1,z=1)
sigma=c(s=1,y=1,z=1)
test=simpleEvoModel(100,500,omega = 2,delta = 4 ,b=2,K=200,mu=0.001,E=E,sigma=sigma)

genes=c("x","y","z")

png("images/env_fit_pop.png",width=800,height=400)
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
plot(test$meanw,type="l",ylab="meanw")
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

cols=colorRampPalette(c("blue","yellow","red"))(500)[as.numeric(cut(test$meanw,breaks=500))]
plot3d(test$meanw,test$env,col=cols,pch=20)
plot(test$meanw,test$env,col=colorRampPalette(c("blue","yellow","red"))(5000),pch=20)

plot(test$meanw[1100:1300],type="l",ylab="meanw")
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
plot(test$meanw,type="l",ylab="meanw")
plot(test$env,type="l",col="red",ylab="environment")
plot(test$popsize,type="l",col="blue",ylab="popsize")
dev.off()


omegas=seq(0,3,.5)
allos_best=sapply(omegas,function(o)replicate(100,mean(simpleEvoModel(100,300,omega = o,delta = 2 ,b=2,K=200,mu=0.001,E=E,sigma=sigma,log=T)$pop$z)))
png("omegas_vs_z.png")
boxplot(allos_best,ylab="mean value of z",xlab=expression(omega),axes=F) 
axis(2)
axis(1,1:length(omegas),label = omegas)
box()
dev.off()

allos_rand=sapply(omegas,function(o)replicate(100,mean(simpleEvoModel(100,300,omega = o,delta = 2 ,b=2,K=200,mu=0.001,E=E,sigma=sigma,type="random",log=T)$pop$z)))
png("omegas_vs_z_random.png")
boxplot(allos_rand,ylab="mean value of z",xlab=expression(omega),axes=F) 
axis(2)
axis(1,1:length(omegas),label = omegas)
box()
dev.off()

library(rgl)
plot3d(test$meanw,test$env,col=cols,pch=20)

omegas=seq(-0.5,2.5,.5)
deltas=.0625*2^seq(0:6)

library(parallel)
names(omegas)=omegas
names(deltas)=deltas

cl <- makeForkCluster(7,outfile="")
osnds_s=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));mean(replicate(50,mean(simpleEvoModel(100,10,omega = o,delta = d ,b=2,K=100,mu=0.001,E=E,sigma=sigma,log=F)$pop$z,na.rm=T)),na.rm=T)}))
osnds_rand=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));mean(replicate(50,mean(simpleEvoModel(200,1000,type="random",omega = o,delta = d ,b=2,K=200,mu=0.001,E=E,sigma=sigma,log=F)$pop$z,na.rm=T)),na.rm=T)}))
stopCluster(cl)


combined=expand.grid(omegas,deltas)
cl <- makeForkCluster(49,outfile="")

bigosnds_popsize=parApply(cl,combined,1,function(i){o=i[1];d=i[2];print(paste(o,d));mean(replicate(100,nrow(simpleEvoModel(1000,2000,omega = o,delta = d ,b=2,K=1000,mu=0.0001,E=E,sigma=sigma,log=F,type="best")$pop)),na.rm=T)})
bigosnds_rand_popsize=parApply(cl,combined,1,function(i){o=i[1];d=i[2];print(paste(o,d));mean(replicate(100,nrow(simpleEvoModel(1000,2000,omega = o,delta = d ,b=2,K=1000,mu=0.0001,E=E,sigma=sigma,log=F,type="random")$pop)),na.rm=T)})

tomat_rand_popsize=matrix(bigosnds_rand_popsize,nrow=7,ncol=7)
tomat_popsize=matrix(bigosnds_popsize,nrow=7,ncol=7)

matrixcol=hcl.colors(200, "YlOrRd", rev = TRUE)
png("nonrand.png")
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
image(log(deltas),omegas,osnds,zlim=c(0,1))
legend_image <- as.raster(matrix(rev(matrixcol), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'mean of the mean\n value of z')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

png("rand.png")
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
image(log(deltas),omegas,osnds_rand,zlim=c(0,1))
legend_image <- as.raster(matrix(rev(matrixcol), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'mean of the mean\n value of z')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

png("bignonrand.png")
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
image(log(deltas),omegas,t(tomat),zlim=c(0,1))
legend_image <- as.raster(matrix(rev(matrixcol), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'mean of the mean\n value of z')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

png("bigrand.png")
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
image(log(deltas),omegas,t(tomat_rand),zlim=c(0,1))
legend_image <- as.raster(matrix(rev(matrixcol), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'mean of the mean\n value of z')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

png("bignonrand_popsize.png")
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
image(log(deltas),omegas,t(tomat_popsize),zlim=c(0,1000))
legend_image <- as.raster(matrix(rev(matrixcol), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'mean of the mean\n population size\n (K=1000)')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1000,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

png("bigrand_popsize.png")
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
image(log(deltas),omegas,t(tomat_rand_popsize),zlim=c(0,1000))
legend_image <- as.raster(matrix(rev(matrixcol), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'mean of the mean\n population size\n (K=1000)')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1000,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()


##TESTING generate output:

n=100
tstep=200
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
t=simpleEvoModel(n,tstep,omega = 0,delta = 0 ,b=2,K=200,mu=c(x=0.01,y=0,z=0),E=c(x=.01,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,type="best",pop=pop)


statfun=c("mean","sd")
statvar=c("x","y","z","gp","ilp","p","w")
names(statvar)=statvar
names(statfun)=statfun

output=lapply(statfun,function(i)lapply(statvar,function(j)c()))
for( tt in 1:length(t$allpop))for(sf in statfun)for(sv in statvar)output[[sf]][[sv]]=c(output[[sf]][[sv]],match.fun(sf)(t$allpop[[tt]][,sv]))

output2=NULL
output2=updateOutput(output2,t$allpop[[tt]],statfun ,statvar)
for( tt in 1:length(t$allpop))output2=updateOutput(output2,t$allpop[[tt]],statfun ,statvar)

i=sample.int(tstep,1)
print(paste("testing if updateOutput work with timestep:",i))
test= ( output2$mean$p[i] == output$mean$p[i] &  output$mean$p[i] == mean(t$allpop[[i]]$p))
print(paste("test is:",test))
test= ( output2$sd$p[i] == output$sd$p[i] &  output$sd$p[i] == sd(t$allpop[[i]]$p))
print(paste("test is:",test))



##testing wrappers:
getVarXMeanW(n,tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=c(s=2^2,y=1,z=1),pop=pop)

#n,tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=c(s=10^e,y=1,z=1),pop=pop
replicateNTime(10,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=c(s=2^4,y=1,z=1),pop=pop)


#testing knockout genes:
n=100
tstep=200
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
t=simpleEvoModel(n,tstep,omega = 0,delta = 0 ,b=2,K=200,mu=c(x=0.01,y=0,z=0),E=c(x=.01,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,type="best",pop=pop,allpops=T)

#all following test should be true
testknockout <- function(res){
    tt=sample(tstep,1)
    if(!any(res$allpop[[tt]]$y == 0) ||!any(res$allpop[[tt]]$z == 0) ){print("genes are not knockout");return(FALSE);}
    any(( res$allpop[[tt]]$gp - res$allpop[[tt]]$ilp == res$allpop[[tt]]$gp - res$allpop[[tt]]$p) &&
    (res$allpop[[tt]]$gp - res$allpop[[tt]]$ilp == 0))
}




#testing getSummary and getTraj
testSummaries <- function(folder){
    onesubfolder= sample(list.dirs(folder),1)
    onesimu= sample(list.files(onesubfolder),1)
    tocheck=file.path(onesubfolder,onesimu)
    tocheck=gsub("/+","/",tocheck)
    print(paste("checking: ", tocheck))
    onetraj=getTraj(tocheck)
    laststep=sample.int(nrow(onetraj),1)
    summarytraj=getSummary(onetraj,nstep=laststep,vars=c("var_x","N","mean_w"))
    load(file.path(onesubfolder,"crossExplore.bin"))
    as.numeric(as.character(binded[binded[,"filename"] == tocheck,"mean_w"])) - summarytraj["mean_w"]
    as.numeric(as.character(binded[binded[,"filename"] == tocheck,"mean_w"])) - mean(onetraj[(nrow(onetraj)-laststep):(nrow(onetraj)),"mean_w"])
    summarytraj["var_x"] == mean(onetraj[(nrow(onetraj)-laststep):(nrow(onetraj)),"var_x"])
}
