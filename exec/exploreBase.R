source("protomodels.R")


n=200
omegas=seq(-0.5,2.5,.5)
deltas=.0625*2^seq(0:6)
names(omegas)=omegas
names(deltas)=deltas


##X and Z knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
n=200
omegas=seq(-0.5,2.5,.5)
deltas=.0625*2^seq(0:6)
names(omegas)=omegas
names(deltas)=deltas


##X and Z knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
t=evosolearn(n,200,omega = 0,delta = 0 ,b=2,K=1000,mu=c(x=0.01,y=0,z=0),E=c(x=1,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,sls="random",pop=pop)
plotAllVariable(t)



library(parallel)

tstep=100
cl <- makeForkCluster(4,outfile="")
extinctions_yzko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(evosolearn(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0,z=0),E=c(x=1,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,sls="random",pop=pop)$allpop)<tstep))}))

##Z knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=runif(n,0,1),z=rep(0,n)))
extinctions_zko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(evosolearn(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0.01,z=0),E=c(x=1,y=1,z=0),sigma=c(s=1,y=1,z=1),log=T,sls="random",pop=pop)$allpop)<tstep))}))


##no knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1)))
extinctions_nko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(evosolearn(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0.01,z=0.01),E=c(x=1,y=1,z=1),sigma=c(s=1,y=1,z=1),log=T,sls="random",pop=pop)$allpop)<tstep))}))


##y knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=runif(n,0,1)))
extinctions_yko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(evosolearn(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0,z=0.01),E=c(x=1,y=0,z=1),sigma=c(s=1,y=1,z=1),log=T,sls="random",pop=pop)$allpop)<tstep))}))

stopCluster(cl)

n=200
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=runif(n,0,1)))
t=evosolearn(n,1000,omega = 1,delta = 1 ,b=2,K=1000,mu=c(x=0.01,y=0,z=0.01),E=c(x=1,y=0,z=1),sigma=c(s=1,y=1,z=1),log=F,sls="random",pop=pop)


omegas=seq(1.5,2.5,.5)
deltas=.0625*2^seq(0:6)

names(omegas)=omegas
names(deltas)=deltas

extinctions=sapply(omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(50,length(evosolearn(n,200,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0,z=0),E=c(x=1,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,sls="random",pop=pop)$allpop)<500))}))


#GRAPH SERGEY

n=500
tstep=100

t=evosolearn(n,tstep,omega = 1.5,delta = 1 ,b=2,K=500,mu=c(x=0.01,y=0.01,z=0.01),E=c(x=1,y=1,z=1),sigma=c(s=2,y=2,z=2),log=T,sls="best")

png("allvariables.png",pointsize = 14,width=800,height=1200)
plotAllVariable(t)
dev.off()

png("allvariables_HDR.png",pointsize = 14,width=800,height=1200)
par(cex=2)
plotAllVariable(t,hdr=T)
dev.off()

png("comparisons.png",pointsize = 14,width=600,height=800)
par(mfrow=c(2,1),cex=1.2)
par(mar=c(2,4,2,4))
plot(sapply(1:tstep,function(it)mean(t$allpop[[it]]$p)),type="l",ylim=range(sapply(t$allpop,"[[","p")),ylab="p",main="Phenosls and theta",bty="n",)
lines(sapply(1:tstep,function(it)mean(t$allpop[[it]]$p)+sd(t$allpop[[it]]$p)),lty=3)
lines(sapply(1:tstep,function(it)mean(t$allpop[[it]]$p)-sd(t$allpop[[it]]$p)),lty=3)
par(new=T)
plot(t$env,type="l",col="blue",yaxt="n",xaxt="n",ylim=range(sapply(t$allpop,"[[","p")),ylab="",bty="n",)
axis(4,col="blue",col.axis="blue")
mtext(expression(theta),4,2,col="blue")

plot( sapply(1:tstep,function(it)mean(abs(t$allpop[[it]]$p-t$env[it]))),ylim=range(sapply(1:tstep,function(it)range(abs(t$allpop[[it]]$p-t$env[it])))),ylab=expression(group("|",theta - p,"|")),bty="n",type="l",main="Distance to theta and fitness")
lines(sapply(2:tstep,function(it)mean(abs(t$allpop[[it]]$p-t$env[it]))+sd(abs(t$allpop[[it]]$p-t$env[it]))),lty=3)
lines(sapply(1:tstep,function(it)mean(abs(t$allpop[[it]]$p-t$env[it]))-sd(abs(t$allpop[[it]]$p-t$env[it]))),lty=3)
par(new=T)
plot( sapply(1:tstep,function(it)mean(t$allpop[[it]]$w)),ylim=range(sapply(t$allpop,"[[","w")),type="l",col="red",yaxt="n",xaxt="n",bty="n",ylab="")
lines(sapply(1:tstep,function(it)mean(t$allpop[[it]]$w)+sd(t$allpop[[it]]$w)),col="red",lty=3)
lines(sapply(1:tstep,function(it)mean(t$allpop[[it]]$w)-sd(t$allpop[[it]]$w)),col="red",lty=3)
axis(4,col="red",col.axis="red")
mtext("w",4,2,col="red")
dev.off()


plot( sapply(1:tstep,function(it)mean(t$allpop[[it]]$p-t$env[it])),sapply(1:tstep,function(it)mean(t$allpop[[it]]$w)))

plot( sapply(1:tstep,function(it)mean(t$allpop[[it]]$p)),type="l")
points( sapply(1:tstep,function(it)mean(t$env[it])),type="l")



#### test 1/2/2020

#Fixed variable:

library(parallel)
cl <- makeForkCluster(7,outfile="")
testSelection=do.call("rbind",parLapply(cl,seqSigmas,function(e)
                                     {
                                         print(paste(repet,e))
                                         replicateNTime(repet,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=c(s=e,y=1,z=1),pop=pop)
                                     }))
stopCluster(cl)

final=as.data.frame(cbind(testSelection,sigma_s=rep(seqSigmas,each=repet)))

pdf("exploreSigmapdf",width=10)
par(mfrow=c(1,2))
plot(final[,3],final[,1],log="xy",xlab=colnames(final)[3],ylab=colnames(final)[1],pch=20,col=alpha(1,.2))
plot(final[,3],final[,2],log="x",xlab=colnames(final)[3],ylab=colnames(final)[2],pch=20,col=alpha(1,.2))
dev.off()

pdf("exploreMu.pdf",width=10)
par(mfrow=c(1,2))
plot(final[,3],final[,1],log="x",xlab=colnames(final)[3],ylab=colnames(final)[1],pch=20,col=alpha(1,.2))
plot(final[,3],final[,2],log="x",xlab=colnames(final)[3],ylab=colnames(final)[2],pch=20,col=alpha(1,.2))
dev.off()


###generalise exploration
explore=list()

omega=0
delta=0
n=1000
b=2
tstep=100
repet=10

allvarlists=list()
allvarlists[["mu"]]=10^seq(-4,-.1,length.out=50)
allvarlists[["K"]]=2^seq(8,15,length.out=50)
allvarlists[["E"]]=10^seq(-2,1.5,length.out=50)
allvarlists[["m"]]=10^seq(-3,1,length.out=50)
allvarlists[["sigma"]]=2^seq(-10,5,length.out=50)
#u=expand.grid(allvarlists)


cl <- makeForkCluster(50,outfile="")
for(gene in c("z")){
    ##Reset all variable and population
    pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)))
    for(var in names(allvarlists)){
        varlist=allvarlists[[var]]
        explore[[gene]][[var]]=do.call("rbind",
                                       parLapply(cl,varlist,function(v,var,gene)
                                                 {
                                                     sigma=c(s=1,y=1,z=1)
                                                     K=1000
                                                     mu=c(x=0,y=0,z=0)
                                                     E=c(x=0,y=0,z=0)
                                                     m=c(x=0,y=0,z=0)
                                                     #explore[[gene]]=list()
                                                     mu[gene]=1/(n)
                                                     m[gene]=0.1
                                                     E[gene]=0.1
                                                     if(var=="K")K=v
                                                     if(var=="sigma")sigma["s"]=v
                                                     if(var=="mu")mu[gene]=v
                                                     if(var=="m")m[gene]=v
                                                     if(var=="E")E[gene]=v
                                                     print(paste("gene",gene,repet,var,v))
                                                     do.call("rbind",lapply(1:repet,function(i)
                                                                            {
                                                                                if(gene == "x") pop[[gene]]=runif(n,-1,1)
                                                                                else pop[[gene]]=runif(n)
                                                                                getVarXMeanW(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m,gene)
                                                                            }))
                                                 },var=var,gene=gene))
                                       explore[[gene]][[var]]=as.data.frame(cbind(explore[[gene]][[var]],rep(varlist,each=repet)))
                                       colnames(explore[[gene]][[var]])[ncol(explore[[gene]][[var]])]=var
    }
stopCluster(cl)
save(file="firstExploration_uncrossed.bin",explore)

##X and Z knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
t=evosolearn(n,200,omega = 0,delta = 2 ,b=2,K=1000,mu=c(x=0.01,y=0,z=0),E=c(x=1,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,sls="random",pop=pop)
plotAllVariable(t)

##### testing simple predictions
n=5000
omegas=0
deltas=0
tstep=100
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)))
mus=c(0.001,.01,.1)
sigma=10
m=10

plot(mus,eq2830a(n,mus,sigma,10),ylim=c(0,10),type="l",col="red")
for(i in 1:50){

allvar=sapply(mus,function(mu){
t=evosolearn(n,tstep,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=m,y=0,z=0),mu=c(x=mu,y=0,z=0),E=c(x=0,y=0,z=0),sigma=c(s=sigma,y=1,z=1),log=T,sls="random",pop=pop)
t$sd$x[tstep]
                                       })

points(mus,allvar)
}


par(mforw=c(2,1))
sigma=10000
m=.1
check=binded[binded$sigma == sigma & binded$m == m & binded$E == .1 & binded$K == 2000,]
plot(mus,eq2833b(mean(check$N),mus,sqrt(sigma),m),col="red",pch=20,type="l")
lines(mus,check[,"var_x"])

for(sigma in c(0.1,1,10,100,1000))lines(mus,eq2833b(1000,mus,sigma,m),col="black",pch=20,xlab="my",ylab="mu",type="l")
lines(mus,check[,"var(x)"],col="green")
lines(mus,eq2833b(2000,mus,.6,m),col="red",pch=20)
lines(mus,eq2833b(2000,mus,sigma,m),col="red",pch=20)

mu=0.001
check=binded$x[binded$x$mu == mu & binded$x$m == m & binded$x$E == .1 & binded$x$K == 1000,]
sigmas=unique(na.exclude(binded$x$sigma))
varsigma=tapply(check[,"var(x)"],check$sigma,mean)
plot(sigmas,eq2833b(1000,mu,sigmas,m),col="red",pch=20,type="l",log="x")
lines(sigmas,varsigma)
points(check$sigma,check[,"var(x)"])
legend("topleft",c("eq28.33b","simulations"),col=c("red",1),lty=1)

mu=0.001
E=.1
K=1000
sigma=.1

check=binded[binded$mu == mu & binded$sigma == sigma & binded$m == m & binded$E == E & binded$K == K,]
plot(getTraj(check$filename,"N"))
mus=unique(na.exclude(binded$mu))
varmu=tapply(check[,"var_x"],check$mu,mean)
plot(mus,eq2833b(1000,mu=mus,sigma=500,m),col="red",pch=20,type="l",ylab="var(x)")
lines(mus,varmu)
points(check$mu,check[,"var_x"])

sigmas=10^(-5:5)
plot(sigmas,eq2829c(1000,mu,sigmas,m),col="red",pch=20,type="l",log="xy")
lines(sigmas,eq2833b(1000,mu,sigmas,m),col="green",pch=20,type="l")#xlim=c(0,1),ylim=c(0,.002))
lines(sigmas,eq2830a(1000,mu,sigmas,m),col="blue",pch=20,type="l")#xlim=c(0,1),ylim=c(0,.002))

##testing getTraj and stuff
folder="sixsixequalone/"
folder="unifPop50k/"
sumaries=getFullExperimentSummary(folder)
sumariesUnif=getFullExperimentSummary(folderUnif)
mus=unique(sumaries$mu)
ms=unique(sumaries$m)
m=0.1
E=0
K=1000
sigma=.1
nstep=3000


folder=c("sixsixequalone","unifPop50k","forpercpu")
allsummaries=lapply(folder,getFullExperimentSummary)
names(allsummaries)=folder

###======== environment changes:

folder=c("growingDelta","unifPop50k")
allsummaries=lapply(folder,getFullExperimentSummary)
names(allsummaries)=folder



    scale=list(var_x=list("0.1"=c(0,0.004),"0.2"=c(0,0.002),"0.4"=c(0,0.005),"10000"=c(0,0.15)),
               N=c(0,2100),
               mean_w=c(.2,1)
               )
    side=list(var_x="topleft",N="bottomright",mean_w="bottomright")


for(sigma in c(0.1,1,10,100,1000))lines(mus,eq2833b(1000,mus,sigma,m),col="black",pch=20,xlab="my",ylab="mu",type="l")
lines(mus,check[,"var(x)"],col="green")
plot(mus,eq2833b(.1,mus,K,m),col="red",pch=20,type="l")
lines(mus,eq2833b(1000,mus,1,m),col="red",pch=20,type="l")
lines(mus,eq2833b(1000,mus,10,m),col="red",pch=20,type="l")


par(mfrow=c(2,1),cex=1.2)
par(mar=c(2,4,2,4))
plot(signleExemple$mean$p,type="l",ylim=range(signleExemple$p),ylab="p",main="Phenosls and theta",bty="n",)
lines(signleExemple$mean$p-sqrt(signleExemple$mean$p),lty=3)
lines(signleExemple$mean$p+sqrt(signleExemple$mean$p),lty=3)
par(new=T)
plot(signleExemple$theta,type="l",col="blue",yaxt="n",xaxt="n",ylim=range(sapply(t$allpop,"[[","p")),ylab="",bty="n",)
axis(4,col="blue",col.axis="blue")
mtext(expression(theta),4,2,col="blue")

bigtest21=evosolearnM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)
system.time({bigtest22=evosolearnM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)})
system.time({bigtest23=evosolearnM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)})
system.time({bigtest24=evosolearnM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)})
system.time({bigtest11=evosolearn(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)})
system.time({bigtest12=evosolearn(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)})
system.time({bigtest13=evosolearn(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)})
system.time({bigtest14=evosolearn(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,sls="random",pop=pop)})
allvar=list(
            bigtest11$sd$x^2,
            bigtest12$sd$x^2,
            bigtest13$sd$x^2,
            bigtest14$sd$x^2,
            bigtest21$var$x,
            bigtest22$var$x,
            bigtest23$var$x,
            bigtest24$var$x
            )

par(mfrow=c(1,2))
cols=rainbow(8)
pdf("individualsrun.pdf")
plot(1,1,sls="n",xlim=c(1,10000),ylim=range(allvar),xlab="timestep",ylab="var(x)")
lapply(1:8,function(i)lines(allvar[[i]],col=cols[i]))
dev.off()

pdf("distributionvar.pdf")
plot(1,1,sls="n",ylim=c(0,500),xlim=c(0,.05),xlab="var(x)",ylab="density",main="distribution of variance for the 8000 last time step")
lapply(1:8,function(i)lines(density(allvar[[i]][2000:10000],from=0),col=cols[i],lwd=2))
dev.off()

#check different stats
boxplot(sapply(c(mean=mean,median=median,Mode=Mode,mode=hdrmode),function(f)sapply(allvar,function(i,b,e)f(i[b:e]),b=2000,e=10000)))

#check different windows
boxplot(sapply(seq(10,9000,100),function(b)sapply(allvar,function(i,b,e)Mode(i[b:e]),b=b,e=10000)))

mean(bigtest22$var$x[2000:10000])
mean(bigtest23$var$x[2000:10000])
mean(bigtest24$var$x[2000:10000])

hdr((bigtest11$sd$x[2000:10000])^2)$mode
hdr((bigtest12$sd$x[2000:10000])^2)$mode
hdr((bigtest13$sd$x[2000:10000])^2)$mode
hdr((bigtest14$sd$x[2000:10000])^2)$mode

dev.off()
pdf("distributionoflastvar.pdf")
plot(density((bigtest14$sd$x[2000:10000])^2))
lines(density(bigtest21$var$x[2000:10000]),ylim=c(0,400))
lines(density(bigtest22$var$x[2000:10000]),col="green")
lines(density(bigtest23$var$x[2000:10000]))
lines(density(bigtest24$var$x[2000:10000],from=0))

lines(density((bigtest11$sd$x[2000:10000])^2,from=0))
lines(density((bigtest12$sd$x[2000:10000])^2,from=0))
lines(density((bigtest13$sd$x[2000:10000])^2,from=0))
lines(density((bigtest14$sd$x[2000:10000])^2,from=0))

abline(v=Mode(bigtest21$var$x[2000:10000]))
abline(v=hdr(bigtest22$var$x[2000:10000])$mode,col="green")
abline(v=Mode(bigtest23$var$x[2000:10000]))
abline(v=Mode(bigtest24$var$x[2000:10000]))

abline(v=Mode((bigtest11$sd$x[2000:10000])^2))
abline(v=Mode((bigtest12$sd$x[2000:10000])^2))
abline(v=Mode((bigtest13$sd$x[2000:10000])^2))
abline(v=Mode((bigtest14$sd$x[2000:10000])^2))

test=evosolearn( 1000, 1000,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=100)


hdrmode <- function(d)hdr(d)$mode
