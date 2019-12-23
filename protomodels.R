source("environment.R")

#Return a list that stores different metrics of interest
simpleEvoModel <- function(n,tstep,epsilon=c(x=.01,y=.01,z=.01),sigma=c(s=1,y=1,z=1),omega,delta,b,K,mu,genes=c("x","y","z"),m=c(x=.3,y=.3,z=.3),type="best",log=T){
    env=c()
    env$theta=environment(tstep,omega,delta)
    env$b=b
    env$K=K
    a=1:n
    #Generate initial population (here all gene are randomly selected
    pop=cbind.data.frame(id=a,parent_id=a,x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1),fitness=rep(0,n))
    parents=NULL
    meanf=c()
    popsize=c()
    allpop=list()
    for( t in 1:tstep){
        if(log)print(paste(" timestep:",t))
        popsize=c(popsize,n)

        ##genetic phase
        e1=rnorm(n,0,epsilon['x'])
        pop$gp=pop$x+e1

        ##learning phase
        e2=rnorm(n,0,epsilon['y'])
        pop$ilp=pop$gp+pop$y*(env$theta[t]-pop$gp)+e2

        ##social learning phase
        P=socialLearning(pop,reference=parents,type=type,thetat=env$theta[t]) #get the list of which phenotype is socially copied by every agent
        e3=rnorm(n,0,sigma['z'])
        pop$slp=pop$ilp+pop$z*(P-pop$ilp)+e3

        ##phenotype check
        #pop$p=(1-pop$y)*(1-pop$z)*pop$x+(1-pop$z)*pop$y*env$theta[t]+pop$z*P+(e1*(1-pop$y)*(1-pop$z)+e2*(1-pop$z)+e3)
        #print(mean(pop$p - pop$slp))

        ##computation of the fitness
        pop$fitness = exp(-((pop$slp-env$theta[t])^2)/(2*sigma['s']^2)-((pop$y)^2)/(2*sigma['y']^2)-((pop$z)^2)/(2*sigma['z']^2))

        #selection
        selected=which(runif(n)<reproduction(pop$fitness,b,n,K))

        #reproduction
        nchilds=rpois(length(selected),b)
        childs=do.call("rbind",apply(cbind(selected,nchilds),1,function(s)do.call("rbind.data.frame",replicate(s[2],pop[s[1],],simplify=F))))
        if(is.null(childs))break
        childs$parent_id=childs$id
        childs$id=rownames(childs)

        #mutation

        newn=nrow(childs)
        for(g in genes){
            mutated=which(runif(newn)<mu)
            childs[mutated,g]=childs[mutated,g]+rnorm(length(mutated),0,m[g])
            if(g %in% c("y","z")){
                childs[mutated,g][childs[mutated,g]<0] = 0
                childs[mutated,g][childs[mutated,g]>1] = 1
            }
        }
        parents=pop[selected,] #we keep parents info (fitness,behavior, etc...) for next social learning
        oldpop=pop
        pop=childs
        n=newn
        meanf=c(meanf,mean(pop$fitness))
        allpop[[t]]=pop
    }
    return(list(meanf=meanf,env=env$theta,pop=pop,popsize=popsize,allpop=allpop))
}



#' @param w:vector with fitness
#' @param b:numeric value for the rate of birth
#' @param n:numeric value for the size of the pop (redundant/should be equal to length(w))
#' @param K:environment carrying capacity
#' @return: numeric vector of probabilities
reproduction <- function(w,b,n,K) 1/(1+(b-1)*(n/(K*w)))

#' @param newpop: a dataframe with fitness and agents ID
#' @param reference: a dataframe with phenotype and agents ID
#' @param type: a string define the type of copy to be done: \in {"parents","best","average","randon"}
#' @return: a unique numeric value or a vector of size nrow(newpop) with phenotypes to be copied 
socialLearning <- function(newpop,reference,thetat=NULL,type="random"){

    ##Checking for imature
    if(is.null(reference))reference=newpop #What happen for the first time step when the reference group doesn't have any final phenotype? should we choose phenotype before social learning? random social learning effect? 
    if(is.null(reference$slp))reference$slp=reference$ilp #What happen for the first time step when the reference group doesn't have any final phenotype? should we choose phenotype before social learning? random social learning effect? 
    if(anyNA(reference$slp))reference$slp[is.na(reference$slp)]=reference$ilp[is.na(reference$slp)] #if some of the reference group 

    if(type=="parents")
        return(reference$slp[match(newnewpop$parent_id,reference$id)])

    if(type=="best"){
        if(is.null(thetat))stop("when selecting best agents an environmental condition has to be given")
        best=which.min(abs(reference$slp-thetat))
        return(reference$slp[best]) #return the phenotype of the best individual in the reference group 
    }
    if(type=="average")
        return(mean(reference$slp))
    if(type=="random"){
        selected=sample(reference$id,nrow(newpop),replace=T) #we radomly assign a teacher for each individual of the new newpop
        return(reference$slp[match(selected,reference$id)])
    }
    stop("a type of copy should be chosen among parents,best,average,randon")
}



