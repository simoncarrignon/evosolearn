source("protomodels.R")


n=1000
omegas=seq(-0.5,2.5,.5)
deltas=.0625*2^seq(0:6)
names(omegas)=omegas
names(deltas)=deltas


##X and Z knockout
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
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
n=1000
b=2
#pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)))
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)))
tstep=1000
repet=10

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
                                         replicateNTime(repet,n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
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


     system.time({bigtest21=simpleEvoModelM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})
     system.time({bigtest22=simpleEvoModelM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})
     system.time({bigtest23=simpleEvoModelM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})
     system.time({bigtest24=simpleEvoModelM(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})
     system.time({bigtest11=simpleEvoModel(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})
     system.time({bigtest12=simpleEvoModel(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})
     system.time({bigtest13=simpleEvoModel(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})
     system.time({bigtest14=simpleEvoModel(n,10000,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=F,type="random",pop=pop)})

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
     plot(1,1,type="n",xlim=c(1,10000),ylim=range(allvar),xlab="timestep",ylab="var(x)")
     lapply(1:8,function(i)lines(allvar[[i]],col=cols[i]))
     dev.off()

     pdf("distributionvar.pdf")
     plot(1,1,type="n",ylim=c(0,500),xlim=c(0,.05),xlab="var(x)",ylab="density",main="distribution of variance for the 8000 last time step")
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

     
     hdrmode <- function(d)hdr(d)$mode
