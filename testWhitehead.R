source("protomodels.R")
getFirst <- function(allpopres){
    for(i in 2:length(allpopres$allpop)){
        ctsls=table(factor(allpopres$allpop[[i]][,"sls"],levels=c(-1,0,2,3)))
        if(length(ctsls)==1)break
    }
    return(ctsls)
}


q=0.8494 #sigma_s

delta=2^seq(-3,3)

omega=seq(-.5,2.5,.5)

sigma=c(0.08,0.25,.75) #sigma for whitead=>accuracy of behavior (E) for us
E_x=c(0.08,0.15,.25,.50,.75) #sigma for whitehead, E for us
C_v=1-c(0.01,0.03,0.05,0.07,0.1) #cost vertical learning
C_h=C_v-c(0.01,0.03,0.05,0.07,0.1)#cost horizont learning
C_i=C_h-c(0.01,0.03,0.05,0.07,0.1)#cost individual learning

u=expand.grid(sigma=E_x,C_v=C_v,C_i=C_i)

n=100
pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)))
pop=cbind(pop,sls=-1)  #sls -1 => genetic

pop[,"x"]=rnorm(n,0,sigma)
pop[26:50,"y"]=1
pop[26:50,"sls"]=0#sls 0 =>individual 
pop[51:100,"z"]=1
pop[51:75,"sls"]=2
pop[76:100,"sls"]=3


tstep=15000
K=1500

library(parallel)
cl <- makeForkCluster(3,outfile="")
allresults=lapply(1:7,function(d)
                  {
                      lapply(1:7,function(o)
                             {
                                 print(paste(o,d))
                                 env=environment(N=tstep,omega=omega[o],delta=delta[d])
                                 parLapply(cl,1:5,function(i)
                                           {
                                               print(i);
                                               p=unlist(u[sample(nrow(u),1),])
                                               result=getFirst(simpleEvoModel(n,tstep,omega = 0,delta = 0 ,b=2,K=K,mu=c(x=0,y=0,z=0),m=c(x=0,y=0,z=0),E=c(x=0,y=0,z=unname(p["sigma"]),sigma=c(s=exp(1),y=exp(unname(p["C_i"])),z=exp(unname(p["C_v"]))),log=F,sls="mixed",pop=pop,outputrate=1,allpops=T,theta=env))
                                               return(result)
                                           }
                                 )
                             }
                      )
                  }
)

stopCluster(cl)

summarized=lapply(rev(allresults),lapply,function(i)apply(sapply(i,function(i)(i>0)/(sum(i>0))),1,mean))
save(file=paste0("summarized.bin",summarized))


pdf("whitehead2007.pdf",width=6,height=6)
par(mfrow=c(7,7),mar=rep(.1,4),oma=c(6,6,1,1))
lapply(rev(summarized),lapply,function(i)pie(i,labels=NA))
par(new=T,mfrow=c(1,1),oma=rep(0,4),mar=c(4,4,1,1))
plot(1,1,ylim=c(0,1),xlim=c(0,1),type="n",xaxt="n",yaxt="n",xlab=expression(omega),ylab=expression(delta))
axis(1,at=seq(0,1,length.out=4),labels=omega[seq(1,7,length.out=4)])
axis(2,at=seq(0,1,length.out=4),labels=delta[seq(1,7,length.out=4)])
dev.off()



