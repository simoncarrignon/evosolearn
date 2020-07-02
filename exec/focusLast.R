old <- Sys.time() 

args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 

ns=args[1]#first argument is the number of slave
nsm=args[2]#second argument the number of simulation  #here will be 160 in any case
mainfold=args[3] #third argument = name of the folder wher to store the results
env_i=args[4] #for nwo will  be only NGRIP2 to avoid bad interpol
#fun_i=args[5]
#res=as.numeric(args[6])
#lin=args[7]
sls=args[5]

if(is.na(mainfold) | mainfold=="") mainfold=Sys.info()['nodename']

fi=0
fold=paste0(mainfold,fi)
while(file.exists(fold)){
    fi=fi+1
    fold=paste0(mainfold,fi)
}

print(paste0("Abc will be stored in folder: ",fold))
dir.create(fold)


source("R/corefunctions.R")
source("R/tools.R")
source("R/environment.R")
library(parallel)

##print(paste("resolution should be",res))
##Manage Environment 

#realdata=read.csv(paste0("data/",env_i,".csv"))
realdata=read.csv(paste0("report/data/ngrip2.csv"))
env=realdata$dTsVscales
#assign("f",get(fun_i))
#if(fun_i != "interpolate"){
#    new=f(realdata$dTsVscales,realdata$year,max(getDateResolution(realdata$year)))
#    env=new$data
#} else{
#    realres=unique(getDateResolution(realdata$year)) #check if interpolation possible/usefull
#    if(length(realres) == 1 && realres <= res){
#        env=realdata$dTsVscales
#    }else{
#        if(lin){
#            env=f(realdata$dTsVscales,realdata$year,res)
#        }else{
#
#            if(env_i == "vostok"){
#                ns=getMean2(data=realdata$dTsVscales,year=realdata$year,by=max(getDateResolution(realdata$year)))
#                env=interpolate(theta=ns$data,times=ns$year,finalres=20,delta=.8,omega=1)
#            }
#
#            else
#                env=f(realdata$dTsVscales,realdata$year,res,omega=1.41,delta=0.9894122)
#        }
#    }
#}

tstep=length(env)
allparameters=list()
allparameters[["mu_x"]]=10^runif(1,-5,-2)
allparameters[["mu_y"]]=10^runif(1,-5,-2)
allparameters[["mu_z"]]=10^runif(1,-5,-2)
allparameters[["K"]]=c(1000)
allparameters[["m_x"]]=runif(1)
allparameters[["m_y"]]=runif(1)
allparameters[["m_z"]]=runif(1)
allparameters[["E_x"]]=runif(1,.1,1.5)
allparameters[["E_y"]]=runif(1,.1,1.5)
#while(allparameters[["E_y"]]<allparameters[["E_x"]])allparameters[["E_y"]]=runif(1,.1,1.5)
allparameters[["E_z"]]=runif(1,.1,1.5)
#while(allparameters[["E_z"]]>allparameters[["E_y"]])allparameters[["E_z"]]=runif(1,.1,1.5)
allparameters[["sigma_s"]]=runif(1,1,10)
allparameters[["sigma_y"]]=runif(1,0,2)
while(allparameters[["sigma_y"]]>allparameters[["sigma_s"]])allparameters[["sigma_y"]]=runif(1,0,2)
allparameters[["sigma_z"]]=runif(1,1,100)
while(allparameters[["sigma_z"]]<allparameters[["sigma_y"]])allparameters[["sigma_z"]]=runif(1,1,100)

allparameters[["outputrate"]]=1

lengthtocheck=10000/20
parameters=as.data.frame(allparameters)
repet=160
parameters=parameters[rep(seq_len(nrow(parameters)),repet),]
genes=c("x","y","z")
cl <- makeForkCluster(ns,outfile="")

explore=do.call("rbind.data.frame",
                parLapply(cl,1:nrow(parameters),function(v,parameters,env)
                          {
                              print(paste0(paste0("g",genes,collapse=" "),", sim #",v,"/",nrow(parameters)))
                              K=parameters[v,"K"]
                              n=parameters[v,"K"]
                              outputrate=parameters[v,"outputrate"]
                              mu=parameters[v,paste("mu",genes,sep="_")]
                              names(mu)=genes
                              m=parameters[v,paste("m",genes,sep="_")]
                              names(m)=genes
                              E=parameters[v,paste("E",genes,sep="_")]
                              names(E)=genes
                              sigma=parameters[v,paste("sigma",c("s","y","z"),sep="_")]
                              names(sigma)=c("s","y","z")
                              pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
                              pop[,"x"]=rnorm(n,mean(env[1]),sd(env[1:5]))
                              fullmat=evosolearn(b=2,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m,outputrate=outputrate,sls="fitprop",allpop=F,repro="sex",theta=env)
                              filename_mat=file.path(fold,paste0("fullmat",v,".rds"))
                              saveRDS(file=filename_mat,fullmat[(nrow(fullmat)-lengthtocheck):nrow(fullmat),])
                              c(as.list(getSummary(fullmat,nstep=5,vars=c(paste0("mean_",genes),paste0("var_",genes),"N","mean_w","mean_p","theta"))),filename=filename_mat)
                          },parameters=parameters,env=env)
                )

stopCluster(cl)

binded=cbind(explore,parameters)
save(file=file.path(fold,"crossExplore.bin"),binded)
saveRDS(file=file.path(fold,"env.rds"),env)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


