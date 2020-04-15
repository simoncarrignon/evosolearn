source("environment.R")
source("tools.R")


#' @param w:vector with fitness
#' @param b:numeric value for the rate of birth
#' @param n:numeric value for the size of the pop (redundant/should be equal to length(w))
#' @param K:environment carrying capacity
#' @return: numeric vector of probabilities
selection <- function(w,b,n,K) 1/(1+(b-1)*(n/(K*w)))


#' @param ouptut a current output to be update or initialized if NULL
#' @param pop the current population upon wihch statistics have to be calculated
#' @param statfun a vector with the different function we apply on the population
#' @param statvar a vector with the different varaible we measure in the population
#' @param prop if the proporition of different strategies should be ouptut
updateOutputLine <- function(pop,statfun,statvar,getname=F,prop=T){
    res=c()
    if(getname){
        res=unlist(lapply(statfun,function(sf)lapply(statvar,function(sv)paste(sf,sv,sep="_"))))
        if(prop)res=c(res,c("prop_y","prop_z"))
    }
    else {
        res=unlist(lapply(statfun,function(sf)lapply(statvar,function(sv)match.fun(sf)(pop[,sv]))))
        if(prop)res=c(res,countStrategies(pop)[1:2])
    }
    return(res)
}


simpleEvoModel <- function(n,tstep,E=c(x=.01,y=.01,z=.01),sigma=c(s=1,y=1,z=1),omega,delta,b,K,mu=c(x=.3,y=.3,z=.3),genes=c("x","y","z"),m=c(x=.3,y=.3,z=.3),sls="best",log=F,pop=NULL,allpops=F,statfun=c("mean","var"),statvar=c("x","y","z","gp","ilp","p","w"),outputrate=1,vt=NULL,theta=NULL,prop=TRUE,reproduction="asex"){

    if(length(mu)==1)mu=c(x=mu,y=mu,z=mu)
	if(is.null(theta)){
	   theta=environment(tstep,omega,delta,vt)
	}
    tstep=min(length(theta),tstep)

    #Generate initial population (here all gene are randomly selected
    if(is.null(pop))pop=generatePop(n,distrib=list(x=runif(n,-1,1),y=runif(n,0,1),z=runif(n,0,1)),df=F)
    n=nrow(pop)

	#prepare outputs
	names(statfun)=statfun
	names(statvar)=statvar
    outputparam=c(E,sigma,omega,delta,b,K,mu,m)
    names(outputparam)=c(paste("E",names(E),sep="_"),paste("sigma",names(sigma),sep="_"),"omega","delta","b","K",paste("mu",names(mu),sep="_"),paste("m",names(m),sep="_"))
    outputparam=outputparam[c(1,4,8,9,10,11,14)]
    outputsnames=c("t",updateOutputLine(NULL,statfun,statvar,getname=T,prop=prop),"N","theta",names(outputparam))
    output=matrix(nrow=(tstep/outputrate),ncol=length(outputsnames))
    colnames(output)=outputsnames

    #output[1,]=c(1,updateOutputLine(pop,statfun,statvar,prop=prop),n,NA,outputparam)

    popsize=c()
    err1=E['x']>0 #precomputing error
    err2=E['y']>0
    err3=E['z']>0
    parents=NULL
    if(allpops)allpop=list()
    modt=1
    for( t in 1:tstep){
        if(log &&  ((t %% outputrate) == 0))print(paste(" timestep:",t))

        ##genetic phase
        if(err1)e1=rnorm(n,0,E['x'])else e1=0
        #e1=rnorm(n,0,E['x'])
        pop[,"gp"]=pop[,"x"]+e1

        ##learning phase
        if(err2)e2=rnorm(n,0,E['y'])else e2=0
        pop[,"ilp"]=pop[,"gp"]+pop[,"y"]*(theta[t]-pop[,"gp"])+e2

        ##social learning phase
        if(sum(pop[,"z"])>0)
            P=socialLearning(pop,reference=parents,sls=sls,thetat=theta[t]) #get the list of which phenotype is socially copied by every agent
        else
            P=0
        #e3=rnorm(n,0,E['z'])
        if(err3)e3=rnorm(n,0,E['z'])else e3=0
        pop[,"p"]=pop[,"ilp"]+pop[,"z"]*(P-pop[,"ilp"])+e3

        ##phenotype check
        #pop[,"p"]=(1-pop[,"y"])*(1-pop[,"z"])*pop[,"x"]+(1-pop[,"z"])*pop[,"y"]*theta[t]+pop[,"z"]*P+(e1*(1-pop[,"y"])*(1-pop[,"z"])+e2*(1-pop[,"z"])+e3)
        #print(mean(pop[,"p"] - pop[,"p"]))

        ##computation of the fitness
        pop[,"w"] = fitness(pop[,"p"],theta[t],pop[,"x"],pop[,"y"],pop[,"z"],sigma['s'],sigma['y'],sigma['z'])


        ##save the population state
        if((t %% outputrate) == 0 ){
            output[modt,]=c(t,updateOutputLine(pop,statfun,statvar,prop=prop),n,theta[t],outputparam)
            if(allpops)allpop[[modt+1]]=pop
            modt=modt+1 #maybe a way to calculate the indice of the output matrix witouth keeping this indice
        }


        #selection
        selected=which(runif(n)<selection(pop[,"w"],b,n,K))
        if(length(selected)<1)break

        #reproduction
        nchilds=rpois(length(selected),b)
        if(sum(nchilds)<1)break
        childs=matrix(nrow=sum(nchilds),ncol=ncol(pop))
        colnames(childs)=colnames(pop)
        c=1
        for( p in seq_along(selected)){
            if(nchilds[p]>0){
                for(i in c:(c+nchilds[p]-1)){
                    childs[i,]=pop[selected[p],]
                    if(reproduction=="sex"){
                        p2=sample(selected[-p],1)
                        cross=runif(3)<1/3
                        childs[i,c("x","y","z")][cross]=pop[p2,c("x","y","z")][cross]
                    }
                }
                c=c+nchilds[p]
            }
        }
        childs[,"parent_id"]=childs[,"id"]
        #if(sum(nchilds)<3)print((max(pop[,"id"])+1):(max(pop[,"id"])+1+sum(nchilds)))
        if(sum(nchilds)==1) childs[,"id"]=(max(pop[,"id"])+1)
        else childs[,"id"]=(max(pop[,"id"])+1):(max(pop[,"id"])+nrow(childs))


        #mutation

        newn=nrow(childs)
        for(g in genes){
            mutated=which(runif(newn)<mu[g])
            childs[mutated,g]=childs[mutated,g]+rnorm(length(mutated),0,m[g])
            if(g %in% c("y","z")){
                childs[mutated,g][childs[mutated,g]<0] = 0
                childs[mutated,g][childs[mutated,g]>1] = 1
            }
        }
        parents=pop[selected,,drop=F] #we keep parents info (fitness,behavior, etc...) for next social learning
        oldpop=pop
        pop=childs
        n=newn


    }
	if(allpops){output=list(summary=output,allpop=allpop)}
    return(output)
}

fitness <- function(p,theta,x,y,z,sigma_s,sigma_y,sigma_z){
        w = exp(-(p-theta)^2/(2*sigma_s^2)- y^2/(2*sigma_y^2)-z^2/(2*sigma_z^2))
        return(w)
}


#' @param newpop: a dataframe with fitness and agents ID
#' @param reference: a dataframe with phenotype and agents ID
#' @param sls: a string define the sls of copy to be done: \in {"parents","best","average","randon"}
#' @return: a unique numeric value or a vector of size nrow(newpop) with phenotypes to be copied 
socialLearning <- function(newpop,reference,thetat=NULL,sls="random"){

    ##Checking for imature
    if(is.null(reference))reference=newpop #What happen for the first time step when the reference group doesn't have any final phenotype? should we choose phenotype before social learning? random social learning effect? 
    if(is.null(reference[,"p"]))reference[,"p"]=reference[,"ilp"] #What happen for the first time step when the reference group doesn't have any final phenotype? should we choose phenotype before social learning? random social learning effect? 
    if(anyNA(reference[,"p"]))reference[,"p"][is.na(reference[,"p"])]=reference[,"ilp"][is.na(reference[,"p"])] #if some of the reference group 

    if(sls=="parents")
        return(unname(reference[,"p"][match(newpop[,"parent_id"],reference[,"id"])]))

    if(sls=="best"){
        if(is.null(thetat))stop("when selecting best agents an environmental condition has to be given")
        best=which.min(abs(reference[,"p"]-thetat))
        return(unname(unname(reference[,"p"][best]))) #return the phenotype of the best individual in the reference group 
    }
    if(sls=="average")
        return(mean(reference[,"p"]))
    if(sls=="random"){
        selected=sample(reference[,"id"],nrow(newpop),replace=T) #we radomly assign a teacher for each individual of the new newpop
        return(unname(reference[,"p"][match(selected,reference[,"id"])]))
    }
    if(sls=="mixed"){
        #parents Vertical
        selected_p = newpop[,"p"]

        #This could/should be done in a loop/apply for all sls

        best=which(newpop[,"sls"] == i_sls["best"])
        if(length(best)>0)
            selected_p[best]=reference[,"p"][which.min(abs(reference[,"p"]-thetat))]

        random=which(newpop[,"sls"] ==i_sls["random"])
        if(length(random)>0){
            randsel=sample(reference[,"id"],length(random),replace=T) #we radomly assign a teacher for each individual of the new newpop
            selected_p[random]= reference[,"p"][match(randsel,reference[,"id"])]
        }

        parentsls=which(newpop[,"sls"] == i_sls["parents"])
        if(length(parentsls)>0)
            selected_p[parentsls]=reference[,"p"][match(newpop[parentsls,"parent_id"],reference[,"id"])]

        average=which(newpop[,"sls"] == i_sls["average"])
        if(length(average)>0)selected_p[average]=mean(reference[,"p"])

        return(unname(selected_p))

    }
    stop("a sls should be chosen among parents,best,average,randon,mmixed")
}


#' @param distrib should be a list with 3 elements x yz)
generatePop <- function(n,distrib,df=F){
    if(length(distrib)!=3)stop("please give distribution for the 3 genes")
    if(sum(sapply(distrib,length))!= n * 3)stop("please give the full distribution (for each n individuals) for the 3 genes")
    a=1:n
    pop=c()
    if(df)pop=cbind.data.frame(id=a,parent_id=a,x=distrib$x,y=distrib$y,z=distrib$z,w=rep(0,n))
    else {
        pop=cbind(id=a,gp=rep(0,n),ilp=rep(0,n),p=rep(0,n),parent_id=a,x=distrib$x,y=distrib$y,z=distrib$z,w=rep(0,n))
    }
    return(pop)
}



#' @param ouptut a current output to be update or initialized if NULL
#' @param pop the current population upon wihch statistics have to be calculated
#' @param statfun a vector with the different function we apply on the population
#' @param statvar a vector with the different varaible we measure in the population
updateOutput <- function(output=NULL,pop,statfun,statvar){

    if(is.null(output) || sum(statvar %in% colnames(pop))<length(statvar))
        output=lapply(statfun,function(i)lapply(statvar,function(j)c()))
    else
        for(sf in statfun)
            for(sv in statvar)
                output[[sf]][[sv]]=c(output[[sf]][[sv]],match.fun(sf)(pop[,sv]))
    return(output)
}


#' @param pop the current population with characteristics of all individual
#' @param selected a vector with the id of the selected individuals
#' @param nchilds a vector with the number of children for each individual 
#' @return a matrix qiht the list of children
reproduction <- function(pop,selected,nchilds){
        childs=matrix(nrow=sum(nchilds),ncol=ncol(pop))
        colnames(childs)=colnames(pop)
        c=1
        if(sum(nchilds)==0)return(NULL)
        for( p in seq_along(selected)){
            if(nchilds[p]>0){
                childs[c:(c+nchilds[p]-1),]=t(replicate(nchilds[p],pop[selected[p],]))
                c=c+nchilds[p]
            }
        }
        childs[,"parent_id"]=childs[,"id"]
        #if(sum(nchilds)<3)print((max(pop[,"id"])+1):(max(pop[,"id"])+1+sum(nchilds)))
        if(sum(nchilds)==1) childs[,"id"]=(max(pop[,"id"])+1)
        else childs[,"id"]=(max(pop[,"id"])+1):(max(pop[,"id"])+nrow(childs))
        return(childs)
}
