n=1000
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)

df <- "df.log"
Rprof(df, memory.profiling = TRUE)
evosolearn(n,100,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=T,sls="random",pop=pop)
Rprof(NULL) ; 
sdf=summaryRprof(df)

ndf <- "ndf.log"
pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
Rprof(ndf, memory.profiling = TRUE)
testM=evosolearnM(n,100,omega = 0,delta = 0 ,b=2,K=1000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=10,y=1,z=1),log=T,sls="random",pop=pop)
Rprof(NULL) ; 
sndf=summaryRprof(ndf)

ndf <- "ndf.log"
Rprof(ndf, memory.profiling = TRUE)
testM=evosolearn( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop)
Rprof(NULL) ; 
sondf=summaryRprof(ndf)

andf <- "andf.log"
Rprof(andf, memory.profiling = TRUE)
testB=evosolearn( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=1)
Rprof(NULL) ; 
asondf=summaryRprof(andf)

oldfunction <- "oldfunction.log"
Rprof(oldfunction, memory.profiling = TRUE)
testB=evosolearnOld( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=1)
Rprof(NULL) ; 
csondf=summaryRprof(oldfunction)

newfunction <- "newfunction.log"
Rprof(newfunction, memory.profiling = TRUE)
testB=evosolearn( 1000, 10001,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=1)
Rprof(NULL) ; 
dsondf=summaryRprof(newfunction)

 microbenchmark(evosolearn(100, 100,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=10),evosolearn(100, 100,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=.1,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=T,sls="random",pop=pop,outputrate=10))


n=1000
popsize=c(100,500,1000,2000,5000)

testN=lapply(popsize,function(n)
             {
                 pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=rep(0,n),z=rep(0,n)),df=F)
                 print(n)
                 microbenchmark(evosolearn(n=n,tstep=1000,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=0,z=0),mu=c(x=.001,y=0,z=0),E=c(x=0,y=0,z=0),sigma=c(s=1000,y=1,z=1),log=F,sls="random",pop=pop,outputrate=1),unit="s")
             }
)


testALLRandom=lapply(popsize,function(n)
               {
                   pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=runif(n),z=runif(n)),df=F)
                   print(n)
                   microbenchmark(evosolearn(n=n,tstep=1000,omega = 0,delta = 0 ,b=2,K=2000,m=c(x=.1,y=.1,z=.1),mu=c(x=.001,y=0.001,z=0.001),E=c(x=0,y=0,z=0),sigma=c(s=1000,y=10,z=10),log=F,sls="random",pop=pop,outputrate=1),unit="s")
               }
)



testALLBestK=lapply(popsize,function(k)
               {
                   N=100
                   pop=generatePop(N,distrib=list(x=runif(N,-1,1),y=runif(N),z=runif(N)),df=F)
                   aver=microbenchmark(evosolearn(n=N,tstep=200,omega = 0,delta = 0 ,b=2,K=k,m=c(x=.1,y=.1,z=.1),mu=c(x=.001,y=0.001,z=0.001),E=c(x=0,y=0,z=0),sigma=c(s=1000,y=10,z=10),log=F,sls="best",pop=pop,outputrate=1),unit="s",times=10)
               }
)

testALLBestKB=lapply(popsize,function(k)
               {
                   N=100
                   pop=generatePop(N,distrib=list(x=runif(N,-1,1),y=runif(N),z=runif(N)),df=F)
                   aver=microbenchmark(evosolearn(n=N,tstep=200,omega = 0,delta = 0 ,b=2,K=k,m=c(x=.1,y=.1,z=.1),mu=c(x=.001,y=0.001,z=0.001),E=c(x=0,y=0,z=0),sigma=c(s=1000,y=10,z=10),log=F,sls="best",pop=pop,outputrate=1),unit="s",times=10)
               }
)
