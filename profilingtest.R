n=1000
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)

df <- "df.log"
Rprof(df, memory.profiling = TRUE)
testdf=testouille=simpleEvoModel(n,100,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=T,sls="random",pop=pop)
Rprof(NULL) ; 
sdf=summaryRprof(df)

ndf <- "ndf.log"
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
Rprof(ndf, memory.profiling = TRUE)
testM=simpleEvoModelM(n,100,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=T,sls="random",pop=pop)
Rprof(NULL) ; 
sndf=summaryRprof(ndf)

ndf <- "ndf.log"
Rprof(ndf, memory.profiling = TRUE)
testM=simpleEvoModel( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop)
Rprof(NULL) ; 
sondf=summaryRprof(ndf)

andf <- "andf.log"
Rprof(andf, memory.profiling = TRUE)
testB=simpleEvoModel( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=1)
Rprof(NULL) ; 
asondf=summaryRprof(andf)

oldfunction <- "oldfunction.log"
Rprof(oldfunction, memory.profiling = TRUE)
testB=simpleEvoModelOld( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=1)
Rprof(NULL) ; 
csondf=summaryRprof(oldfunction)

newfunction <- "newfunction.log"
Rprof(newfunction, memory.profiling = TRUE)
testB=simpleEvoModel( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=1)
Rprof(NULL) ; 
dsondf=summaryRprof(newfunction)

 microbenchmark(simpleEvoModel(100, 100,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=10),simpleEvoModel(100, 100,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=10))

