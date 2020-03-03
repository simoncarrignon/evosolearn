old <- Sys.time() 

args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 

ns=args[1]#first argument is the number of slave
nsm=args[2]#second argument the number of simulation 
mainfold=args[3] #third argument = name of the folder wher to store the results

if(is.na(mainfold) | mainfold=="") mainfold=Sys.info()['nodename']

fi=0
fold=paste0(mainfold,fi)
while(file.exists(fold)){
    fi=fi+1
    fold=paste0(mainfold,fi)
}

print(paste0("Abc will be stored in folder: ",fold))
dir.create(fold)


source("protomodels.R")
library(parallel)
allparameters=list()
allparameters[["mu"]]=c(0.001,0.01)
allparameters[["K"]]=c(500,1000,2000)
#allparameters[["E"]]=c(0,.05,.1)
allparameters[["E"]]=0
allparameters[["m"]]=c(.05,.1,.2)
#allparameters[["m"]]=.1
allparameters[["sigma"]]=c(1,2,4,10000)
allparameters[["delta"]]=c(0,.1,.2,.4,1)
#allparameters[["delta"]]=0
#allparameters[["vt"]]=c(0.001,0.02,.04,.08,.2)
allparameters[["vt"]]=0
allparameters[["omega"]]=c(0,2)
allparameters[["outputrate"]]=100
parameters=as.data.frame(expand.grid(allparameters))
repet=nsm
parameters=parameters[rep(seq_len(nrow(parameters)),repet),]


#omega=2
n=1000
b=2
tstep=10000
mu=c(x=0,y=0,z=0)
E=c(x=0,y=0,z=0)
m=c(x=0,y=0,z=0)
sigma=c(s=1,y=1,z=1)

genes=c("y")
pop=generatePop(n,distrib=list(x=rep(0,n),y=rep(0,n),z=rep(0,n)),df=F)
cl <- makeForkCluster(ns,outfile="")
for(gene in genes){
    ##Reset all variable and population
    mu=c(x=0,y=0,z=0)
    E=c(x=0,y=0,z=0)
    m=c(x=0,y=0,z=0)
    sigma=c(s=1,y=1,z=1)

    explore=do.call("rbind.data.frame",
                    parLapply(cl,1:nrow(parameters),function(v,parameters,gene,pop)
                              {
                                  print(paste0("g",gene,", sim #",v,"/",nrow(parameters)))
                                  delta=parameters[v,"delta"]
                                  omega=parameters[v,"omega"]
                                  K=parameters[v,"K"]
                                  outputrate=parameters[v,"outputrate"]
                                  vt=parameters[v,"vt"]
                                  mu[gene]=parameters[v,"mu"]
                                  m[gene]=parameters[v,"m"]
                                  E[gene]=parameters[v,"E"]
                                  sigma["s"]=parameters[v,"sigma"]
                                  if(gene == "x") pop[,gene]=runif(n,-1,1)
                                  else pop[,gene]=runif(n,0,1)
                                  fullmat=simpleEvoModel(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m,outputrate=outputrate,vt=vt)
                                  filename_mat=file.path(fold,paste0("fullmat",v,".bin"))
                                  save(file=filename_mat,fullmat)
                                  c(as.list(getSummary(fullmat,nstep=10,vars=c(paste0("var_",gene),"N","mean_w"))),filename=filename_mat)
                              },parameters=parameters,gene=gene,pop=pop)
                    )

}
stopCluster(cl)

binded=cbind(explore,parameters)
save(file=file.path(fold,"crossExplore.bin"),binded)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format

output=F
if(output){
    library(xtable)
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


    print(file="tables/parameters_values.tex",x,sls="latex",sanitize.text.function=identity,include.rownames = FALSE,auto=T)

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

}



                                  
