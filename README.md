# Evolution, social Learning and climate fluctuation

## Clone and run the model:

To clone this folder: 

```bash
git clone git@framagit.org:sc/pleistoclimate.git
```

then to test the model:

```bash
cd pleistoclimate/
```

then run `R` within this folder
and  from `R`

```R
source("protomodels.R")
```
to load the model and the functions it uses. To do a simple run of the model one can do:

```R
#setup the parameters
epsilon=c(x=1,y=1,z=1) #the standard deviation of the error associated with the expression of each phenotype (p' to p''')
sigma=c(s=2,y=2,z=2) #Selection strength
m=c(x=.3,y=.3,z=.3)
type="best" #type of copy for social learning
n=100
tstep=500

#run the simulation
test=simpleEvoModel(n = 100,tstep = 500,omega = 2,delta = 4 ,b = 2,K = 200,mu=0.001,epsilon = epsilon,sigma = sigma)
```

For know in the resulting list (here `test`), I return all the populations for all time steps. This is for now, to check if everything is good and how the different variables of interest (x,y,z,...) are distributed. Some summaries statistics are also available in `test$meanf`, `test$popsize`,... 

To plot the available summaries:

```R
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
plot(test$meanf,type="l",ylab="meanf")
plot(test$env,type="l",col="red",ylab="environment")
plot(test$popsize,type="l",col="blue",ylab="popsize")
```
![follow the link if not shown](images/env_fit_pop.png)


To get the mean and standard deviation of z for each time step:

```R
meansdz=sapply(test$allpop,function(i)c(sd(i[["z"]]),mean(i[["z"]])))
```
To check the evolution of the mean of each genes 

```R
genes=c("x","y","z")
par(mfrow=c(3,1))
par(mar=c(2,4,1,1))
sapply(genes,function(g)plot(sapply(test$allpop,function(i)mean(i[[g]])),ylab=paste("gene",g),type="l"))
```
![follow the link if image not shown](images/allgenes.png)

## Check environment 

Some quick exploration of the environment generation function

```R
par(mfrow=c(3,2))
t=10000 #number of timestep
par(mfrow=c(3,2),mar=c(2,2,1,1))
for(alpha in 0:2){ #alpha as used in YK95 is our $\omega$
    ts=TK95(t,alpha) #generate random noise with alpha =1
    plot(ts,type="l") #plot the environment
    y=getSpectrum(ts) #get spectrum of the environment generated
    x=1:length(y)
    plot(log(x),log(y))
    fit=lm(log(y)~log(x),cbind.data.frame(x=1:length(y),y=y)) #fit a linear model to check slope
    abline(fit,col="red") #visualise the fit and coefficient computed
    text(1,max(log(y)),paste("coef=",abs(round(fit$coefficients[2],2))),col="red")
}
```

![environment check](images/exploreEnv.png)


## Further exploration of the model

3d plot of the mean fitness given time and environment

```R
library(rgl)
plot3d(test$meanf,test$env,col=cols,pch=20)
```

In the next block we explore the impact of omega on the mean value of z at the end of the simulation. To do that we creat a list of omegas and for each of them we run 100 simulations: 
```R
omegas=seq(0,3,.5)
allos_best=sapply(omegas,function(o)replicate(100,mean(simpleEvoModel(100,300,omega = o,delta = 2 ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,log=T)$pop$z)))
boxplot(allos_best,ylab="mean value of z",xlab=expression(omega),axes=F) 
axis(2)
axis(1,1:length(omegas),label = omegas)
box()
```
![follow the link if image not shown](images/omegas_vs_z.png)

By default `simpleEvoModel` use the best mechanism to copy, we can compare this when using random by simply doing:

```R
allos_rand=sapply(omegas,function(o)replicate(100,mean(simpleEvoModel(100,300,omega = o,delta = 2 ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,type="random",log=T)$pop$z)))
boxplot(allos_rand,ylab="mean value of z",xlab=expression(omega),axes=F) 
axis(2)
axis(1,1:length(omegas),label = omegas)
box()
```
![follow the link if image not shown](images/omegas_vs_z_random.png)

Similarly we can explore the impact of omega and delta at the same time:
```R
osnds=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));mean(replicate(50,mean(simpleEvoModel(100,50,omega = o,delta = d ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,log=F)$pop$z)))}))
osnds_rand=parSapply(cl,omegas,function(o)sapply(deltas,function(d){print(paste(o,d));mean(replicate(50,mean(simpleEvoModel(100,50,type="random",omega = o,delta = d ,b=2,K=200,mu=0.001,epsilon=epsilon,sigma=sigma,log=F)$pop$z)))}))
```

![nonrand](images/nonrand.png)
![rand](images/rand.png)
