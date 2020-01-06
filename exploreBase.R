source("protomodels.R")


n=500
omegas=seq(-0.5,2.5,.5)
deltas=.0625*2^seq(0:6)
names(omegas)=omegas
names(deltas)=deltas


##X and Z knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
t=simpleEvoModel(n,200,omega = 0,delta = 2 ,b=2,K=1000,mu=c(x=0.01,y=0,z=0),E=c(x=1,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,type="random",pop=pop)
plotAllVariable(t)



library(parallel)

tstep=100
cl <- makeForkCluster(4,outfile="")
extinctions_yzko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(simpleEvoModel(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0,z=0),E=c(x=1,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,type="random",pop=pop)$allpop)<tstep))}))

##Z knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=runif(n,0,1),z=rep(0,n)))
extinctions_zko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(simpleEvoModel(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0.01,z=0),E=c(x=1,y=1,z=0),sigma=c(s=1,y=1,z=1),log=T,type="random",pop=pop)$allpop)<tstep))}))


##no knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1)))
extinctions_nko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(simpleEvoModel(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0.01,z=0.01),E=c(x=1,y=1,z=1),sigma=c(s=1,y=1,z=1),log=T,type="random",pop=pop)$allpop)<tstep))}))


##y knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=runif(n,0,1)))
extinctions_yko=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(10,length(simpleEvoModel(n,tstep,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0,z=0.01),E=c(x=1,y=0,z=1),sigma=c(s=1,y=1,z=1),log=T,type="random",pop=pop)$allpop)<tstep))}))

stopCluster(cl)

t=simpleEvoModel(n,tstep,omega = 1,delta = 1 ,b=2,K=1000,mu=c(x=0.01,y=0,z=0.01),E=c(x=1,y=0,z=1),sigma=c(s=1,y=1,z=1),log=T,type="random",pop=pop)


omegas=seq(-0.5,2.5,.5)
deltas=.0625*2^seq(0:6)

names(omegas)=omegas
names(deltas)=deltas

extinctions=sapply(omegas,function(o)sapply(deltas,function(d){print(paste(o,d));sum(replicate(50,length(simpleEvoModel(n,200,omega = o,delta = d ,b=2,K=1000,mu=c(x=0.01,y=0,z=0),E=c(x=1,y=0,z=0),sigma=c(s=1,y=1,z=1),log=T,type="random",pop=pop)$allpop)<500))}))


#GRAPH SERGEY

n=500
tstep=100

t=simpleEvoModel(n,tstep,omega = 1.5,delta = 1 ,b=2,K=500,mu=c(x=0.01,y=0.01,z=0.01),E=c(x=1,y=1,z=1),sigma=c(s=2,y=2,z=2),log=T,type="best")

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
plot(sapply(1:tstep,function(it)mean(t$allpop[[it]]$p)),type="l",ylim=range(sapply(t$allpop,"[[","p")),ylab="p",main="Phenotype and theta",bty="n",)
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
n=200
b=2
#pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)))
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
tstep=100
repet=100

K=500
sigma=c(s=1,y=1,z=1)
mu=c(x=0,y=0,z=0)
E=c(x=0,y=0,z=0)
m=c(x=0,y=0,z=0)

cl <- makeForkCluster(50,outfile="")
g="x"
    #pop[[g]]=runif(n)
    gene=g
    mu[gene]=1/(n)
    m[gene]=0.3
    E[gene]=0.01

    var="sigma"
    varlist=2^seq(-10,5,length.out=50)
    explore[[var]]=do.call("rbind",
                           parLapply(cl,varlist,function(v)
                                     {
                                         if(var=="K")K=v
                                         if(var=="sigma")sigma["s"]=v
                                         if(var=="mu")mu[gene]=v
                                         if(var=="m")m[gene]=v
                                         if(var=="E")E[gene]=v
                                         print(paste(repet,v))
                                         replicateNTime(repet,N=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
                                     }))
    explore[[var]]=as.data.frame(cbind(explore[[var]],rep(varlist,each=repet)))
    colnames(explore[[var]])[ncol(explore[[var]])]=var
    writeResults(explore[[var]],var,gene)



for(g in c("y","z")){
    pop[[g]]=runif(n)
    gene=g
    mu[gene]=1/(n)
    m[gene]=0.3
    E[gene]=0.1

    var="mu"
    varlist=10^seq(-10,-1,length.out=50)
    explore[[var]]=do.call("rbind",
                           parLapply(cl,varlist,function(v)
                                     {
                                         if(var=="K")K=v
                                         if(var=="sigma")sigma["s"]=v
                                         if(var=="mu")mu[gene]=v
                                         if(var=="m")m[gene]=v
                                         if(var=="E")E[gene]=v
                                         print(paste(repet,v))
                                         replicateNTime(repet,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
                                     }))
    explore[[var]]=as.data.frame(cbind(explore[[var]],rep(varlist,each=repet)))
    colnames(explore[[var]])[ncol(explore[[var]])]=var
    writeResults(explore[[var]],var,gene)

    var="m"
    varlist=10^seq(-3,1,length.out=50)
    explore[[var]]=do.call("rbind",
                           parLapply(cl,varlist,function(v)
                                     {
                                         if(var=="K")K=v
                                         if(var=="sigma")sigma["s"]=v
                                         if(var=="mu")mu[gene]=v
                                         if(var=="m")m[gene]=v
                                         if(var=="E")E[gene]=v
                                         print(paste(repet,v))
                                         replicateNTime(repet,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
                                     }))
    explore[[var]]=as.data.frame(cbind(explore[[var]],rep(varlist,each=repet)))
    colnames(explore[[var]])[ncol(explore[[var]])]=var
    writeResults(explore[[var]],var,gene)

    var="sigma"
    varlist=2^seq(-10,5,length.out=50)
    explore[[var]]=do.call("rbind",
                           parLapply(cl,varlist,function(v)
                                     {
                                         if(var=="K")K=v
                                         if(var=="sigma")sigma["s"]=v
                                         if(var=="mu")mu[gene]=v
                                         if(var=="m")m[gene]=v
                                         if(var=="E")E[gene]=v
                                         print(paste(repet,v))
                                         replicateNTime(repet,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
                                     }))
    explore[[var]]=as.data.frame(cbind(explore[[var]],rep(varlist,each=repet)))
    colnames(explore[[var]])[ncol(explore[[var]])]=var
    writeResults(explore[[var]],var,gene)

    var="E"
    varlist=10^seq(-2,1.5,length.out=50)
    explore[[var]]=do.call("rbind",
                           parLapply(cl,varlist,function(v)
                                     {
                                         if(var=="K")K=v
                                         if(var=="sigma")sigma["s"]=v
                                         if(var=="mu")mu[gene]=v
                                         if(var=="m")m[gene]=v
                                         if(var=="E")E[gene]=v
                                         print(paste(repet,v))
                                         replicateNTime(repet,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
                                     }))
    explore[[var]]=as.data.frame(cbind(explore[[var]],rep(varlist,each=repet)))
    colnames(explore[[var]])[ncol(explore[[var]])]=var
    writeResults(explore[[var]],var,gene)

    var="K"
    varlist=10^seq(1,10,length.out=50)
    explore[[var]]=do.call("rbind",
                           parLapply(cl,varlist,function(v)
                                     {
                                         if(var=="K")K=v
                                         if(var=="sigma")sigma["s"]=v
                                         if(var=="mu")mu[gene]=v
                                         if(var=="m")m[gene]=v
                                         if(var=="E")E[gene]=v
                                         print(paste(repet,v))
                                         replicateNTime(repet,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
                                     }))
    explore[[var]]=as.data.frame(cbind(explore[[var]],rep(varlist,each=repet)))
    colnames(explore[[var]])[ncol(explore[[var]])]=var
    writeResults(explore[[var]],var,gene)


}
    stopCluster(cl)


