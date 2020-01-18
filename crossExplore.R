source("protomodels.R")
library(xtable)
library(parallel)
allparameters=list()
lengthout=10
allparameters[["mu"]]=c(0,0.001,0.01)
allparameters[["K"]]=c(500,1000,2000)
allparameters[["E"]]=c(0,.05,.1)
allparameters[["m"]]=c(.05,.1,.2)
allparameters[["sigma"]]=c(.1,.2,.4,10000)
allparameters[["delta"]]=0
parameters=as.data.frame(expand.grid(allparameters))

tablep=list()
tablep[["mu"]]=paste0("$10^{",round(seq(-4,-.1,length.out=lengthout)),"}$")
tablep[["K"]]=paste0("$2^{",round(seq(8,15,length.out=lengthout)),"}$")
tablep[["E"]]=paste0("$2^{",round(seq(-10,2,length.out=lengthout)),"}$")
tablep[["m"]]=paste0("$10^{",round(seq(-3,1,length.out=lengthout)),"}$")
tablep[["sigma"]]=paste0("$2^{",round(seq(-10,5,length.out=lengthout)),"}$")
tablep[["delta"]]=paste0("$10^{",round(seq(-1,1,length.out=lengthout)),"}$")
tp=as.data.frame(tablep)
colnames(tp)=c("$\\gm$","$K$","$E$","$m$","$\\gs_s$","$\\gd$")
x=xtable(tp)


print(file="tables/parameters_values.tex",x,type="latex",sanitize.text.function=identity,include.rownames = FALSE,auto=T)
repet=2
parameters=parameters[rep(seq_len(nrow(parameters)),repet),]

omega=0
n=1000
b=2
tstep=10000
mu=c(x=0,y=0,z=0)
E=c(x=0,y=0,z=0)
m=c(x=0,y=0,z=0)
sigma=c(s=1,y=1,z=1)
explore=list()

library(parallel)
genes=c("x","y","z")
cl <- makeForkCluster(5,outfile="")
for(gene in genes[1]){
    ##Reset all variable and population
    pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)),df=F)
    mu=c(x=0,y=0,z=0)
    E=c(x=0,y=0,z=0)
    m=c(x=0,y=0,z=0)
    sigma=c(s=1,y=1,z=1)
    if(gene == "x") pop[,gene]=runif(n,-1,1)
    else pop[,gene]=runif(n)
    explore[[gene]]=list()
    explore[[gene]]=do.call("rbind",
                            parLapply(cl,1:nrow(parameters),function(v,parameters){
                                      print(paste0("g",gene,", sim #",v,"/",nrow(parameters)))
                                      delta=parameters[v,"delta"]
                                      K=parameters[v,"K"]
                                      mu[gene]=parameters[v,"mu"]
                                      m[gene]=parameters[v,"m"]
                                      E[gene]=parameters[v,"E"]
                                      sigma["s"]=parameters[v,"sigma"]
                                      getVarXMeanW(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m,gene=gene,nstep=4000,sumstat=mean)
},parameters=parameters))
}
stopCluster(cl)
binded=lapply(explore,function(e)cbind.data.frame(e,parameters))
#save(file="crossExploreDeltas.bin",binded)
#load(file="crossExploreDeltas.bin")


for(gene in c("x","y","z")){
    for(var in c("w","var")){
        alltables=plotAlldimensions(alldata=binded,gene,var,write=T,dir=paste0("images/theta0/g",gene,"/"),pref="explore")

        for(i in 1:length(alltables)){
            alltables[[i]][]=vapply(alltables[[i]],function(m)paste0("\\includegraphics[width=3cm]{",m,"}"),character(1))
            rownames(alltables[[i]])=tablep[["E"]]
            colnames(alltables[[i]])=tablep[["sigma"]]

            addtorow <- list()
            addtorow$pos <- list(0, 0)
            if(var == "var")
                cap=paste("impact of the different variable on variance of gene ",gene)
            if(var == "w")
                cap=paste("impact of the different variable on fitness w")

            #header of the table
            addtorow$command <- c("& \\multicolumn{4}{c}{$\\gs_s$} \\\\\n",
                                  paste("E & ", paste(tablep[["sigma"]],collapse=" & "),"\\\\\n"))

            #print the table
            print(xtable(alltables[[i]],caption=cap),sanitize.text.function=identity,file=paste0("tables/","table_g",gene,"_K",i,"_",var,".tex"),include.colnames=F,add.to.row=addtorow)

        }

    }
}




                                  
