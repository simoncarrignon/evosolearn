#' @export
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

i_sls <- 1:4 #a  global variable to store the indices of individidual socil laerning strategies
names(i_sls) <- c("best","parents","average","random")


#wrappers for exploration
#' @export
replicateNTime <- function(repet,n,tstep,omega,delta,b,K,mu,E,sigma,pop,m){

    do.call("rbind",replicate(repet,getVarXMeanW(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m),simplify=F))
}


#' get some useful summary stat from model run
#' @param nstep how many timestep should we use to compute the summary statitistic 
#' @param sumstat function used to compute the summary statistic
getVarXMeanW <- function(n,tstep,omega,delta,b,K,mu,E,sigma,pop,m,gene,nstep,sumstat){
    t=evosolearn(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
    res=c()
    pre=tstep-nstep
    if(is.null(t$var[[gene]][pre:tstep]))
        res=c(NA,NA)
    else
        res=c(sumstat(t$var[[gene]][pre:tstep]),sumstat(t$mean$w[pre:tstep]))

    names(res)=c(paste0("var(",gene,")"),"mean(w)")
    return(res)
}

#' get some useful summary stat from a model run
#' @param nstep how many timestep should we use to compute the summary statitistic 
#' @param sumstat function used to compute the summary statistic
#' @param fullmat a matrice with the result of a simulation
#' @param vars a vector with the variable we want to summarize using sumstat
#' @export
getSummary <- function(fullmat,nstep,sumstat=mean,vars){
    res=c()
    tstep=nrow(fullmat)
    pre=max(tstep-nstep,1)
    sapply(vars,function(v)sumstat(fullmat[pre:tstep,v]))
}

getTraj  <-  function(filename,var){
    rdsmat=NULL
    rdsmat=tryCatch(readRDS(filename),error=function(e)NULL)
    if(!is.null(rdsmat))
        fullmat=rdsmat
    else{
        fullmat=c()
        summary=c()
        rm(summary)
        rm(fullmat)
        rm("summary")
        rm("fullmat")
        load(file=as.character(filename))
    }
    if(!exists("fullmat")){
        fullmat=summary
    }
    if(length(fullmat)==2)fullmat=fullmat$summary
    if(var=="dist")
        res=abs(fullmat[,"mean_p"]-fullmat[,"theta"])
    else if(var=="distX")
        res=abs(fullmat[,"mean_x"]-fullmat[,"theta"])
    else{
        if(var %in% colnames(fullmat)){
            res=fullmat[,var]
        }
        else{
            res=rep(NA,nrow(fullmat))
        }
    }
    rm("fullmat")
    rm("summary")
    return(res)
}

writeResults <- function(final,var,gene){
    pdf(paste0("explore_",var,"_",gene,".pdf"),width=10)
    par(mfrow=c(1,2))
    plot(final[,3],final[,1],log="x",xlab=colnames(final)[3],ylab=colnames(final)[1],pch=20,col=alpha(1,.2))
    plot(final[,3],final[,2],log="x",xlab=colnames(final)[3],ylab=colnames(final)[2],pch=20,col=alpha(1,.2))
    dev.off()
}


#' @export
plotAlldimensions <- function(alldata,gene,y,dim1="K",dim2="E",dim3="sigma",dim4="m",dim5="mu",dir="images/",pref="explore",write=F,points=F){
    if(write)dir.create(dir,recursive = T)
    data=alldata[[gene]]
    if(y=="var")
        yl=paste0("var(",gene,")")
    if(y=="w")
        yl="mean(w)"
    storename=list()
    eldim1=unique(data[[dim1]])
    eldim2=unique(data[[dim2]])
    eldim3=unique(data[[dim3]])
    eldim4=unique(data[[dim4]])
    eldim5=unique(data[[dim5]])
    len1=length(eldim1)
    len2=length(eldim2)
    len3=length(eldim3)
    len4=length(eldim4)
    len5=length(eldim5)
    for(a in 1:len1){
        #dev.new()
        par(mfrow=c(4,4))
        par(mar=c(4,2,2,1))
        storename[[a]]=matrix(nrow=len2,ncol=len3)
        for(b in 1:len2){
            for(c in 1:len3){
                mt=parse(text=paste0("list(",dim1,"==2^",round(log2(eldim1[a])),",",dim2,"==2^",round(log2(eldim2[b])),",",dim3,"==2^",round(log2(eldim3[c])),")"))
                name=paste0(dim1,a,dim2,b,dim3,c)
                print(name)
                filename=paste0(dir,pref,"-g",gene,"-",y,"-param-",name,".png")
                print(filename)
                if(write)png(filename,pointsize=18)
                plot(1,1,sls="n",xlim=range(data[[dim5]]),ylim=range(data[[yl]],na.rm=T),xlab=parse(text=dim5),ylab=paste0(yl),main=mt)
                cols=rev(heat.colors(len4))
                for(d in 1:len4){
                    subbdata=data[data[[dim1]] == eldim1[a] & data[[dim2]] ==eldim2[b] & data[[dim3]] ==eldim3[c] & data[[dim4]] == eldim4[d],]
                    means=tapply(subbdata[[yl]],subbdata[[dim5]],mean,na.rm=T)
                    lines(eldim5,means,col=cols[d],lwd=3)
                    if(points)points(subbdata[[dim5]],subbdata[[yl]],col=cols[d],pch=20)
                }
                legend("topleft",legend=parse(text=paste("10^",round(log10(eldim4)))),lwd=3,col=cols,title=parse(text=paste0("m[",gene,"]")),bty="n")
                if(write)dev.off()
                storename[[a]][b,c]=filename
            }
        }
    }
    return(storename)
}

#' @export
printOne <- function(alldata,gene,y,K,E,sigma,m,mu,dir="images/",pref="one",write=F){
    data=alldata[[gene]]
    if(y=="var")
        yl=paste0("var(",gene,")")
    if(y=="w")
        yl="mean(w)"
    mt=parse(text=paste0("list(",dim1,"==2^",round(log2(eldim1[a])),",",dim2,"==2^",round(log2(eldim2[b])),",",dim3,"==2^",round(log2(eldim3[c])),")"))
    name=paste0(dim1,a,dim2,b,dim3,c)
    print(name)
    filename=paste0(pref,"_g",gene,"_",y,"_param_",name,".png")
    if(write)png(filename,pointsize=15)
    plot(1,1,sls="n",xlim=range(data[[dim5]]),ylim=range(data[[yl]],na.rm=T),log="x",xlab=expression(dim5),ylab=paste0(yl),main=mt)
    cols=rev(heat.colors(len4))
    for(d in 1:len4){
        subbdata=data[data[["K"]] == K & data[["E"]] == E & data[["sigma"]] ==sigma & data[["m"]] == m,]
        means=tapply(subbdata[[yl]],subbdata[["mu"]],mean,na.rm=T)
        lines(unique(subbdata[["mu"]]),means,col=cols[d],lwd=2)
        #points(subbdata[[dim5]],subbdata[[yl]],col=cols[d],pch=20)
    }
    #legend("topleft",legend=parse(text=paste("10^",round(log10(eldim4)))),lwd=2,col=cols,title=parse(text=paste0("m[",gene,"]")),bty="n")
    if(write)dev.off()
}

#' @export
plotTraj <- function(x,alltraj,col=1,ylim=NULL,xlim=NULL,add=F,lty=1,lwd=1.4,mean=F,...){
    #mn=apply(alltraj,1,mean)
    #sd=apply(alltraj,1,sd)
    qts=apply(alltraj,1,function(i)quantile(i,na.rm=T))
    if(is.null(xlim))xlim=c(0,nrow(alltraj))
    if(is.null(ylim))ylim=range(qts[c(2,4),])
    if(!add)plot(1,1,type="n",xlim=xlim,ylim=ylim,...)
    lines(x,qts[3,],col=col,lwd=lwd,lty=lty)
    if(!mean)lines(x,qts[4,],col=col,lwd=lwd/14,lty=lty)
    if(!mean)lines(x,qts[2,],col=col,lwd=lwd/14,lty=lty)
}

#' @export
addMeanSD <- function(x,y,col=1,plot=T){
    meany=tapply(y,x,mean)
    sdy=tapply(y,x,sd)
    if(plot)points(unique(x),meany,pch=20,col=col)
    if(plot)arrows(unique(x), meany+sdy, unique(x), meany-sdy,angle=90,code=3,length=.05,lwd=2,col=col)
    return(rbind(meany+sdy,meany,meany-sdy))
}

#' @export
getFullExperimentSummary <- function(fold,exclude=NULL){
    allfolds=list.dirs(fold,recursive=F)
    allsummary=c()
    for(f in allfolds){
        try({
        load(file.path(f,"crossExplore.bin"))
        #binded[,c("N","var_x","mean_w")]=apply(binded[,c("N","var_x","mean_w")],2,function(i)as.numeric(as.character(i)))
        if(!is.null(exclude))binded=binded[,which(colnames(binded)!="outputrate")]
        allsummary=rbind(allsummary,binded)
        })
    }
    return(allsummary)
}

#' @export
getSubsetWithTraj <- function(summarydataset,m,sigma,E,K,mu,delta=0,var="var_x",traj=T){
    res=list()
    res$summary=summarydataset[summarydataset$mu %in% mu & summarydataset$sigma %in% sigma & summarydataset$m %in% m & summarydataset$E %in% E & summarydataset$delta %in% delta & summarydataset$K %in% K,]
    if(traj)res$traj=sapply(res$summary$filename,getTraj,var=var)
    return(res)
}

#' Alpha
#'
#' A simple function to change the opacity of a color
#' @param  color the name or idea of a R color
#' @param  alpha a value in [0,1] defining the opacity wanted.
#' @export
alpha <- function(color,alpha) rgb(t(col2rgb(color)/255),alpha=alpha)

#' Shades
#'
#' A simple function to genarte shade of one color by changing its opacity
#' @param  color the name or idea of a R color
#' @param  n number or shades wanted
#' @export
shades<-function(color,n) sapply(seq(0,1,length.out=n+1),alpha,color=color)[2:(n+1)]

#' @export
eq2833a <- function(n,mu,sigma,m)return(mu*m^2*n)

#' @export
eq2830a <- function(n,mu,sigma,m)return((4*mu*sigma)/(1+sigma/(n*m^2)))

#' @export
eq2829c <- function(n,mu,sigma,m)return(m^2/4 * (sqrt(1+(16*mu*sigma)/(m^2))-1))

#' @export
eq2833b <- function(n,mu,sigma,m){
    gamma=m^2/(2*sigma^2)
    return(1*m^2*(gamma*n+1)/(4*gamma*n)*(sqrt(1+2*(gamma*n*4*mu*n)/(gamma*n+1)^2)-1))
}

#' countStrategies
#'
#' A function to count the number of pure social learning and pure individual learner
#' @param  upperlim the minimum value for the gene (y or z) to be considered as "pure" learning
#' @param  lowerlim the maximum value for the other gene (y or z) to be considered as "pure" learning
#' @export
countStrategies <- function(pop,lowerlim=.25,upperlim=.75){
    soc=sum(pop[,"z"] > upperlim & pop[,"y"]< lowerlim)/nrow(pop)
    il=sum(pop[,"y"] > upperlim & pop[,"z"]< lowerlim)/nrow(pop)
    mix=sum(pop[,"z"] > upperlim & pop[,"y"]> upperlim)/nrow(pop)
    return(c("prop_y"=il,"prop_z"=soc,"prop_yz"=mix,"prop_div"=(1-(il+soc+mix))))
}

#' plotMatrixStrateAndEn
#'
#' Function that allows to plot percentage of each different strategies below environment
#' @param sum_mat a dataframe with the percentage to plot, mean ,x
#' @param environment list  of values to be plot with sum_mat$mean_x  on top of the proportions
#'  The process that generate mat:
#' subb=droplevels(binded[binded$mu ==mu &binded$m ==m &binded$sigma ==s & binded$k_z ==kz & binded$k_y ==ky,])
#' nexpe=nrow(subb)
#' vars=c("prop_y","prop_z","mean_x","N","prop_yz")
#' names(vars)=vars
#'  mat_allexp= lapply(vars,function(i)matrix(NA,nexpe,length(theta)))
#'
#'  for(i in 1:50){
#'      load(as.character(subb$filename[i]))
#'      for(v  in vars) mat_allexp[[v]][i,]=summary[,v]
#'  }
#'  sum_mat=lapply(mat_allexp,function(m)apply(m,2,quantile,probs=c(.05,.5,.95),na.rm=T))
#' @export
plotMatrixStrateAndEn <- function(sum_mat,environment,epochs=NULL){
    nsteps=ncol(sum_mat[[1]])
    if(is.null(epochs))epochs=0:(nsteps-1)
    pureIl=rgb(1,.5,0)
    pureSl=rgb(0,.5,1)
    mixed=rgb(1,.5,1)
    par(fig=c(0,1,.8,1),mar=rep(0,4),oma=c(2,2,1,1))
    plot(environment,col="dark green",lwd=1,type="l",axes=F)
    par(new=T)
    plotTres(sum_mat$mean_x,ylim=range(sum_mat$mean_x,environment,na.rm=T,finite=T),xlim=c(0,nsteps))
    par(new=T)
    plotTres(sum_mat$mean_p,ylim=range(sum_mat$mean_p,environment,na.rm=T,finite=T),xlim=c(0,nsteps),col="red")
    mtext(expression(theta),4,0,col="dark green")
    mtext(expression(bar(x)),2,0,adj=.40,col="1")
    mtext(expression(bar(p)),2,0,adj=.70,col="red")
    par(fig=c(0,1,0,.8),new=T)
    plotTres(sum_mat$prop_y,col=pureIl,ylim=c(0,1),xlim=c(0,nsteps))
    par(new=T)
    plotTres(sum_mat$prop_z,col=pureSl,ylim=c(0,1),xlim=c(0,nsteps))
    par(new=T)
    plotTres(sum_mat$prop_yz,col=mixed,ylim=c(0,1),xlim=c(0,nsteps))
    ticsLab=seq(1,length(epochs),length.out=4)
    tics=seq(0,nsteps-1,length.out=4)
    axis(1,at=tics,labels=epochs[ticsLab]);axis(2);box()
}


#' @export
plotTres <- function(u,col=1,...){
    plot(u[2,],col=col,axes=F,type="l",...)
    lines(u[1,],col=col,lwd=.1)
    lines(u[3,],col=col,lwd=.1)
}


# trying to write a solid wrapper to take care of different file load/filename for uneique dataset stoage
#' @export
getUniqueExp <- function(filename){
    if(is.factor(filename))
        filename=as.character(filename)
    rdsmat=NULL
    rdsmat=tryCatch(readRDS(filename),error=function(e)NULL)
    if(!is.null(rdsmat))
        summary=rdsmat
    else{
        fullmat=c()
        summary=c()
        rm(summary)
        rm(fullmat)
        rm("summary")
        rm("fullmat")
        load(file=filename)
        if(!exists("summary")){
            summary=fullmat
        }
    }
    if(length(summary)==2)summary=summart$summary
    return(summary)
}

#' Calculate the CCFD following \emph{Vosoughi et al. (2018)}
#'@param size : sample of sizes
#'@return a order table with two column with first column : frequency and second column the probality of the frequency  
#' @export
ccfd <- function(size){
    total=length(size) 
    size=size[order(size)]
    counts=unique(size)
    id=unique(match(counts,size))
    p=sapply(id,function(i)length(size[i:total])) 
    p=p/total * 100
    x=counts
    y=p
    return(cbind(x,y))
}

#' randomCascades: A slightly (almost not) modified version of the original model by Alex & Damien
#' Original model is in tools/FNM.R
#' @return dataframe that allow to compare easily to the result of the model of the cascade
#' @export
randomCascades <- function(Nmin=100000,Nmax=50000,t_steps=50,y_max=1,runs=1,mu=0.00001,tau = 1,alpha=0,conformity=F,beta=-2,topfive=F,TF=5,C=.1,alberto=F){

    t<-seq(1,t_steps)
    N_vec=round(Nmin*exp(alpha*t))

    T <- as.integer(tau + t_steps) # total time_steps

    N=as.integer(N_vec[1])#inital N value

    old_pop <- seq(1,N) #inialise population
    base = N+1

    corpus <- list()
    vocab<- list()
    Z <-matrix(nrow=y_max,ncol=t_steps-1)

    #save all the values in one popualtion
    total_pop <- vector(mode='integer',length=N*t_steps)

    #beginning of t loop 
    for (t in 1:T){

        rnd = runif(N, min = 0, max = 1)

        ##### the assingment was set to pop, i've changed it to old_pop
        select_vec<-old_pop 

        pop<-old_pop # start new empty population


        if(alberto){
            pop=originalTopfive(pop,TF,C,rnd,mu,base)
        }
        else{ 
            if(conformity){ #introduction of a simple conformity bias 
                #freq=table(pop)
                #rnd=rnd/(freq[as.character(pop)])
                Bias <- conformityBias(old_pop,beta)
                pop[rnd>=mu] <- sample(old_pop, length(pop[rnd>=mu]), replace = TRUE,prob = Bias)
            }
            else{
                pop[rnd>=mu] <- sample(old_pop, length(pop[rnd>=mu]), replace = TRUE)
            }
            if(topfive){
                candidates=topfive(old_pop,TF,C)[runif(length(old_pop)*C)>=mu]
                pop[sample.int(length(old_pop),length(candidates))]= candidates
            }
        }


        pop[rnd<mu] <- seq(base,base+length(pop[rnd<mu])-1)


        base <- base + length(pop[rnd<mu]) # next trait to enter the population
        if (t > tau){
            N <- as.integer(N_vec[t-tau])
            s<-t-tau-1
            first<-(1+(s*N))
            last<-((s*N)+N)
            total_pop[first:last]<-pop
        }

        #rank_old <- rank_now
        old_pop <- pop # current pop becomes old_pop
        #end of t loop
    }

    MC <- sort(table(total_pop),decreasing = TRUE)
    df<-data.frame(rank=rank(MC,ties.method = "first"),MC)

    colnames(df)=c("rank","U","size")
    df$U=as.numeric(df$U)
    df
}


## Function to compute distances of all simulations  to an optimal distance
#' @export
distFivePercent <- function(data,clm){
    allexp=as.character(getUniqeSet(data,clm)$filename)

    y=sapply(allexp,function(f)readRDS(file=f)[,"mean_y"])
    z=sapply(allexp,function(f)readRDS(file=f)[,"mean_z"])
    x=sapply(allexp,function(f)readRDS(file=f)[,"mean_x"])

    nz=mean(apply(z,1,function(l)sum(l > .85,na.rm=T ))/ncol(z))
    s=0
    s=distT(nz,0.05)
    if(is.nan(s))s=1
    vary=sd(apply(y,1,var,na.rm=T))
    s=s+vary
    varx=sd(apply(x,1,var,na.rm=T))
    s=s+varx
    return(s)
}

#' @export
distT<-function(u,n) 1-(abs(u-n)/(sqrt(abs(u^2-n^2))))

#' @export
distU<-function(u,n) 1-(1/abs(u^2-n^2))


#' @export
getUniqeSet <- function(data,vect,out){
    sub=rep(TRUE,nrow(data))
    for(cn in names(vect)){
       sub=sub & (data[,cn] == vect[[cn]])
    }
    data[sub,]
}

#' read all experiments from a folder and return all concatened
#' prefix allows to putsome new folder that were add ed after simulation
#' it should be noted that in the output of a simulation the path of the full experiment is set when simulation are done. if this path change thus and if this function is call frome somewhere else thus pref is need
#' @export
getAllFromFolder <- function(foldername,pref=NULL,log=F){
    if(!is.null(pref))
        foldername=file.path(pref,foldername)
    allexpe=do.call("rbind",lapply(list.files(path=foldername,recursive=T,full.names=T,pattern="*cross*"),function(u){if(log)print(u);load(u);return(binded)}))
    if(!is.null(pref))
        allexpe$filename=file.path(pref,allexpe$filename)
    return(allexpe)
}
