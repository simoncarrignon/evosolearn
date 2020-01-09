allparameters=list()
lengthout=4
allparameters[["mu"]]=10^seq(-4,-.1,length.out=lengthout)
allparameters[["K"]]=2^seq(8,15,length.out=lengthout)
allparameters[["E"]]=2^seq(-10,2,length.out=lengthout)
allparameters[["m"]]=10^seq(-3,1,length.out=lengthout)
allparameters[["sigma"]]=2^seq(-10,5,length.out=lengthout)
parameters=as.data.frame(expand.grid(allparameters))

tablep=list()
lengthout=4
tablep[["mu"]]=paste0("$10^{",round(seq(-4,-.1,length.out=lengthout)),"}$")
tablep[["K"]]=paste0("$2^{",round(seq(8,15,length.out=lengthout)),"}$")
tablep[["E"]]=paste0("$2^{",round(seq(-10,2,length.out=lengthout)),"}$")
tablep[["m"]]=paste0("$10^{",round(seq(-3,1,length.out=lengthout)),"}$")
tablep[["sigma"]]=paste0("$2^{",round(seq(-10,5,length.out=lengthout)),"}$")
x=xtable(as.data.frame(tablep))
x

print(file="parameter_value.tex",x,type="latex",sanitize.text.function=function(x){x},include.rownames = FALSE,auto=T)
repet=4
parameters=parameters[rep(seq_len(nrow(parameters)),repet),]

omega=0
delta=0
n=1000
b=2
tstep=80
    mu=c(x=0,y=0,z=0)
    E=c(x=0,y=0,z=0)
    m=c(x=0,y=0,z=0)
    sigma=c(s=1,y=1,z=1)

cl <- makeForkCluster(50,outfile="")
for(gene in c("x","y","z")){
    ##Reset all variable and population
    pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)))
    mu=c(x=0,y=0,z=0)
    E=c(x=0,y=0,z=0)
    m=c(x=0,y=0,z=0)
    sigma=c(s=1,y=1,z=1)
    if(gene == "x") pop[[gene]]=runif(n,-1,1)
    else pop[[gene]]=runif(n)
    explore[[gene]]=do.call("rbind",
                            parLapply(cl,1:nrow(parameters),function(v,parameters,gene){
                                   print(paste0("g",gene,", sim #",v,"/",nrow(parameters)))
                                   K=parameters[v,"K"]
                                   mu[gene]=parameters[v,"mu"]
                                   m[gene]=parameters[v,"m"]
                                   E[gene]=parameters[v,"E"]
                                   sigma["s"]=parameters[v,"sigma"]
                                   getVarXMeanW(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m,gene=gene)
},parameters=parameters,gene=gene))
}
stopCluster(cl)
#save(file="crossExplore.bin",explore)
load(file="crossExplore.bin")
binded=lapply(explore,function(e)cbind.data.frame(e,parameters))




plot(cbind.data.frame(explore$x,parameters),log="xy")
dev.new()
plot(cbind.data.frame(explore$y,parameters),log="xy")
dev.new()
plot(cbind.data.frame(explore$z,parameters),log="xy")

tapply(binded$x[["var(x)"]],binded$x$sigma,mean,na.rm=T)
tapply(binded$x[["var(x)"]],binded$x$mu,mean,na.rm=T)


gene="z"
y="w"
alldata=binded

data=alldata[[gene]]
if(y=="var")
yl=paste0("var(",gene,")")
if(y=="w")
yl="mean(w)"

for(k in unique(data$K)){
    #dev.new()
    par(mfrow=c(4,4))
    par(mar=c(4,2,2,1))
    for(e in unique(data$E)){
        for(a in unique(data$sigma)){
            mt=parse(text=paste0("list(sigma==10^",round(log10(a)),",E==10^",round(log10(e)),",K==10^",round(log10(k)),")"))
            name=paste0("sigma=2exp",round(log2(a)),"_E=2exp",round(log2(e)),"_K=exp",round(log2(k)))
            print(name)
            png(paste0("images/explore_g",gene,"_",y,"_param_",name,".png"),pointsize=15)
            plot(1,1,type="n",xlim=range(data$mu),ylim=range(data[[yl]],na.rm=T),log="x",xlab=expression(mu),ylab=paste0(yl),main=mt)
            cols=rev(heat.colors(length(unique(data$m))))
            names(cols)=as.character(unique(data$m))
            for(b in unique(data$m)){
                subbdata=data[data$sigma == a & data$m ==b & data$E ==e & data$K == k,]
                means=tapply(subbdata[[yl]],subbdata$mu,mean,na.rm=T)
                lines(unique(data$mu),means,col=cols[as.character(b)],lwd=2)
                points(subbdata$mu,subbdata[[yl]],col=cols[as.character(b)],pch=20)
            }
            legend("topleft",legend=parse(text=paste("10^",round(log10(unique(data$m))))),lwd=2,col=cols,title=expression(m[x]),bty="n")
            dev.off()
        }
    }
}

par(mfrow=c(1,3))
for(g in c("x","y","z"))
            mt=parse(text=paste0("list(sigma==10^",round(log10(a)),",E==10^",round(log10(e)),",K==10^",round(log10(k)),")"))
            print(a)
            plot(1,1,type="n",xlim=range(data$mu),ylim=range(data[[yl]],na.rm=T),log="x",xlab=expression(mu),ylab=paste0("mean(",yl,")"),main=mt)
            cols=rev(heat.colors(length(unique(data$m))))
            names(cols)=as.character(unique(data$m))
            for(b in unique(data$m)){
                subbdata=data[data$sigma == a & data$m ==b & data$E ==e & data$K == k,]
                means=tapply(subbdata[[yl]],subbdata$mu,mean,na.rm=T)
                lines(unique(data$mu),means,col=cols[as.character(b)],lwd=2)
                points(subbdata$mu,subbdata[[yl]],col=cols[as.character(b)],pch=20)
            }

plot4dimensions <- function(alldata,y,dim1,dim2,dim3,dim4,gene){
    for(k in unique(paramane
    dev.new()
par(mfrow=c(4,4))
par(mar=c(4,2,2,1))
for(e in unique(alldata[[gene]][[dim1]])){
for(a in unique(parameters$sigma)){
    print(a)
plot(1,1,type="n",xlim=range(binded$x$mu),ylim=c(0,2),log= "x",xlab=expression(mu),ylab="mean(var(x))",main=parse(text=paste0("sigma==10^",round(log10(a)))))
cols=rev(heat.colors(length(unique(parameters$m))))

names(cols)=as.character(unique(parameters$m))
for(b in unique(parameters$m)){
    subbind=binded$x[binded$x$sigma == a & binded$x$m ==b & binded$x$E ==e,]
    means=tapply(subbind[["var(x)"]],subbind$mu,mean,na.rm=T)
    lines(unique(parameters$mu),means,col=cols[as.character(b)],lwd=2)
    points(subbind$mu,subbind[["var(x)"]],col=cols[as.character(b)],pch=20)
}
    legend("topleft",legend=parse(text=paste("10^",round(log10(unique(binded$x$m))))),lwd=2,col=cols,title=expression(m[x]))
}
}

}
