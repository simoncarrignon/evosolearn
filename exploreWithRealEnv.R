old <- Sys.time() 


args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 


ns=args[1]#first argument is the number of slave
nsm=args[2]#second argument the number of simulation 
mainfold=args[3] #third argument = name of the folder wher to store the results
env_i=args[4]
fun_i=args[5]
res=as.numeric(args[6])
lin=args[7]
sls=args[8]

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

print(paste("resolution should be",res))
##Manage Environment 

realdata=read.csv(paste0("data/",env_i,".csv"))
assign("f",get(fun_i))
if(fun_i != "interpolate"){
    new=f(realdata$dTsVscales,realdata$year,max(getDateResolution(realdata$year)))
    env=new$data
} else{
    realres=unique(getDateResolution(realdata$year)) #check if interpolation possible/usefull
    if(length(realres) == 1 && realres <= res){
        env=realdata$dTsVscales
    }else{
        if(lin){
            env=f(realdata$dTsVscales,realdata$year,res)$data
        }else{

            if(env_i == "vostok"){
                subss=getMean2(data=realdata$dTsVscales,year=realdata$year,by=max(getDateResolution(realdata$year)))
                env=interpolate(theta=subss$data,times=subss$year,finalres=20,delta=.8,omega=0.41)$data
            }

            else
                env=f(realdata$dTsVscales,realdata$year,res,omega=1.41,delta=0.8)$data
        }
    }
}

tstep=length(env)
allparameters=list()
allparameters[["mu"]]=10^(-5:-3)
allparameters[["K"]]=c(1000)
allparameters[["m"]]=(.2*(1:3))
allparameters[["E"]]=.2*(c(1,3,5))
#allparameters[["sigma"]]=(2^(0:4))[1]
allparameters[["sigma"]]=c(1,2)
#allparameters[["delta"]]=2^(0:4)
#allparameters[["vt"]]=(5^(0:4)*10^-3)[1:4]
#allparameters[["omega"]]=2^(-1:3)
allparameters[["outputrate"]]=ceiling(1/500*tstep)
allparameters[["k_z"]]=c(2,4,8)
allparameters[["k_y"]]=c(.3,.6,.9)
#allparameters[["sls"]]=c("random","best")
parameters=as.data.frame(expand.grid(allparameters))
repet=nsm
parameters=parameters[rep(seq_len(nrow(parameters)),repet),]
b=2
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
                parLapply(cl,1:nrow(parameters),function(v,parameters,env)
                          {
                              print(paste0(paste0("g",genes,collapse=" "),", sim #",v,"/",nrow(parameters)))
                              #delta=parameters[v,"delta"]
                              #omega=parameters[v,"omega"]
                              #vt=parameters[v,"vt"]
                              K=parameters[v,"K"]
                              n=parameters[v,"K"]
                              outputrate=parameters[v,"outputrate"]
                              mu[genes]=parameters[v,"mu"]
                              m[genes]=parameters[v,"m"]
                              E[genes]=parameters[v,"E"]
                              sigma["s"]=parameters[v,"sigma"]
			      #sls=parameters[v,"sls"]
                              pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
                              pop[,"x"]=rnorm(n,mean(env[1]),sd(env[1:5]))
                              sigma=c(s=parameters[v,"sigma"],y=parameters[v,"sigma"]*parameters[v,"k_y"],z=parameters[v,"sigma"]*parameters[v,"k_y"]*parameters[v,"k_z"])
                              fullmat=simpleEvoModel(n=n,tstep=tstep,omega = 0,delta = 0 ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m,outputrate=outputrate,vt=vt,sls=sls,allpop=F,repro="sex",prop=T,theta=env)
                              filename_mat=file.path(fold,paste0("fullmat",v,".rds"))
                              summary=fullmat
                              saveRDS(file=filename_mat,summary)
                              c(as.list(getSummary(summary,nstep=0,vars=c(paste0("mean_",genes),paste0("var_",genes),"N","mean_w","mean_p","theta"))),filename=filename_mat)
                          },parameters=parameters,env=env)
                )

stopCluster(cl)

binded=cbind(explore,parameters)
save(file=file.path(fold,"crossExplore.bin"),binded)
saveRDS(file=file.path(fold,"env.rds"),env)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


