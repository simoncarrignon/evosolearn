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
allparameters[["K"]]=c(2000)
allparameters[["k_y"]]=c(1)
allparameters[["k_z"]]=c(2)
#allparameters[["K"]]=c(1000)
allparameters[["E"]]=c(0)#,.05,.1)
#allparameters[["E"]]=0
allparameters[["m"]]=c(.2)#,.5)
#allparameters[["m"]]=.1
allparameters[["sigma"]]=c(2,4)
allparameters[["delta"]]=c(1)
#allparameters[["delta"]]=0
allparameters[["vt"]]=c(0.002,.2)
#allparameters[["vt"]]=0
allparameters[["omega"]]=c(2)
#allparameters[["omega"]]=0
allparameters[["outputrate"]]=100
parameters=as.data.frame(expand.grid(allparameters))
repet=nsm
parameters=parameters[rep(seq_len(nrow(parameters)),repet),]


#omega=2
n=1000
b=2
tstep=40000
mu=c(x=0,y=0,z=0)
E=c(x=0,y=0,z=0)
m=c(x=0,y=0,z=0)
sigma=c(s=1,y=1,z=1)

genes=c("x","y","z")
cl <- makeForkCluster(ns,outfile="")

    ##Reset all variable and population
    mu=c(x=0,y=0,z=0)
    E=c(x=0,y=0,z=0)
    m=c(x=0,y=0,z=0)
    sigma=c(s=1,y=1,z=1)

    explore=do.call("rbind.data.frame",
                    parLapply(cl,1:nrow(parameters),function(v,parameters)
                              {
                                  print(paste0(paste0("g",genes,collapse=" "),", sim #",v,"/",nrow(parameters)))
                                  delta=parameters[v,"delta"]
                                  omega=parameters[v,"omega"]
                                  K=parameters[v,"K"]
                                  outputrate=parameters[v,"outputrate"]
                                  vt=parameters[v,"vt"]
                                  mu[genes]=parameters[v,"mu"]
                                  m[genes]=parameters[v,"m"]
                                  E[genes]=parameters[v,"E"]
                                  sigma=c(s=parameters[v,"sigma"],y=parameters[v,"sigma"]*parameters[v,"k_y"],z=parameters[v,"sigma"]*parameters[v,"k_y"]*parameters[v,"k_z"])
								  pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1)),df=F)
                                  fullmat=simpleEvoModel(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m,outputrate=outputrate,vt=vt,sls="random",allpop=T,prop=T)
                                  filename_mat=file.path(fold,paste0("fullmat",v,".bin"))
                                  summary=fullmat$summary
                                  save(file=filename_mat,summary)
                                  save(file=paste0(filename_mat,"pop"),fullmat)
                                  c(as.list(getSummary(summary,nstep=3,vars=c(paste0("mean_",genes),paste0("var_",genes),"N","mean_w","mean_p","theta","prop_y","prop_z"))),filename=filename_mat)
                              },parameters=parameters)
                    )

stopCluster(cl)

binded=cbind(explore,parameters)
save(file=file.path(fold,"crossExplore.bin"),binded)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format

                                  
