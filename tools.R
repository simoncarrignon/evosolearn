

plotAllVariable <- function(results,hdr=F,vars=NULL,theta=NULL,t=NULL,...){
    varnames=NULL
    if(!is.null(vars))varnames=vars
    if(!is.null(results$allpop)){
        allsd=sapply(results$allpop,function(i)apply(i[,vars],2,sd))
        allmean=sapply(results$allpop,function(i)apply(i[,vars],2,mean))
        if(is.null(vars))varnames=rownames(allsd)
    }
    else{
        if(is.null(vars))varnames=names(results$allpop[[1]])
    }
    cols=c(rainbow(length(varnames)))
    names(cols)=varnames
    par(mfrow=c(length(vars)+1,1),mar=c(0,4,0,4))
    for(i in 1:(length(varnames))){
        varname=varnames[i]
        if(hdr){
            yvalues=lapply(results$allpop,"[[",varname)
            ylim=range(yvalues)
            if(varname %in% c("y","z"))ylim=c(0,1)
            hdr.boxplot(yvalues,border=NA,pch=".",outline=F,col=shades(cols[varname],3),prob=c(50,75,99),space=0,ylab=varname,ylim=ylim,h=.01,...)
            if(varname %in% c("x","ilp","gp","p") & !is.null(theta))
                lines(theta,col="red",lwd=1,lty=1)
        }
        else {
            if(!is.null(results$allpop)){
                meanvar=allmean[varname,]
                sdvar=allsd[varname,]
            }
            else{
                meanvar=results$mean[[varname]]
                sdvar=results$sd[[varname]]
            }
            lims=range(c(meanvar+sdvar,meanvar-sdvar))
            plot(meanvar,ylab=varname,sls="l",col=cols[varname],ylim=lims,xaxt='n',lwd=2,...)
            lines(meanvar+sdvar,ylab=varname,col=cols[varname],lwd=2,lty=3)
            lines(meanvar-sdvar,ylab=varname,col=cols[varname],lwd=2,lty=3)
        }
    }

    t=t[!is.na(t)]
    par(mar=c(2,4,0,4))
    varname=varnames[length(varnames)]
    plot(sapply(results$allpop,nrow),ylab="N",xaxt="n",type="l") 
    i=axis(1,labels=NA,col=NA)
    if(!is.null(t))axis(1,at=i,labels=t[seq(1,length(t),length.out=length(i))])
    par(new=T)
    plot(theta,xlab="",type="l",col="red",yaxt="n",xaxt="n",bty="n",ylab="")
    axis(4,col="red",col.axis="red")
    mtext(expression(theta),4,2,col="red")
}

#wrappers for exploration
replicateNTime <- function(repet,n,tstep,omega,delta,b,K,mu,E,sigma,pop,m){

    do.call("rbind",replicate(repet,getVarXMeanW(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m),simplify=F))
}


#' get some useful summary stat from model run
#' @param nstep how many timestep should we use to compute the summary statitistic 
#' @param sumstat function used to compute the summary statistic
getVarXMeanW <- function(n,tstep,omega,delta,b,K,mu,E,sigma,pop,m,gene,nstep,sumstat){
    t=simpleEvoModel(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
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
getSummary <- function(fullmat,gene,nstep,sumstat=mean,vars){
    res=c()
    tstep=nrow(fullmat)
    pre=max(tstep-nstep,1)
    sapply(vars,function(v)sumstat(fullmat[pre:tstep,v]))
}

getTraj  <-  function(filename,var){
    load(file=as.character(filename))
    if(var=="dist")
        abs(fullmat[,"mean_p"]-fullmat[,"theta"])
    else if(var=="distX")
        abs(fullmat[,"mean_x"]-fullmat[,"theta"])
    else
        fullmat[,var]
}

writeResults <- function(final,var,gene){
    pdf(paste0("explore_",var,"_",gene,".pdf"),width=10)
    par(mfrow=c(1,2))
    plot(final[,3],final[,1],log="x",xlab=colnames(final)[3],ylab=colnames(final)[1],pch=20,col=alpha(1,.2))
    plot(final[,3],final[,2],log="x",xlab=colnames(final)[3],ylab=colnames(final)[2],pch=20,col=alpha(1,.2))
    dev.off()
}


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

plotTraj <- function(x,alltraj,col=1,ylim=NULL,xlim=NULL,add=F,lty=1,lwd=1.4,mean=F,...){
    #mn=apply(alltraj,1,mean)
    #sd=apply(alltraj,1,sd)
    qts=apply(alltraj,1,function(i)quantile(i,na.rm=T))
    if(is.null(xlim))xlim=c(0,nrow(alltraj))
    if(is.null(ylim))ylim=range(qts[c(2,4),])
    if(!add)plot(1,1,type="n",xlim=xlim,ylim=ylim,...)
    lines(x,qts[3,],col=col,lwd=lwd,lty=lty)
    if(!mean)lines(x,qts[4,],lty=1,col=col,lwd=lwd/14,lty=lty)
    if(!mean)lines(x,qts[2,],lty=1,col=col,lwd=lwd/14,lty=lty)
}

addMeanSD <- function(x,y,col=1,plot=T){
    meany=tapply(y,x,mean)
    sdy=tapply(y,x,sd)
    if(plot)points(unique(x),meany,pch=20,col=col)
    if(plot)arrows(unique(x), meany+sdy, unique(x), meany-sdy,angle=90,code=3,length=.05,lwd=2,col=col)
    return(rbind(meany+sdy,meany,meany-sdy))
}

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
#' @param  alpha a value \in [0,1] defining the opacity wanted.
alpha <- function(color,alpha) rgb(t(col2rgb(color)/255),alpha=alpha)

#' Shades
#'
#' A simple function to genarte shade of one color by changing its opacity
#' @param  color the name or idea of a R color
#' @param  n number or shades wanted
shades<-function(color,n) sapply(seq(0,1,length.out=n+1),alpha,color=color)[2:(n+1)]

eq2833a <- function(n,mu,sigma,m)return(mu*m^2*n)

eq2830a <- function(n,mu,sigma,m)return((4*mu*sigma)/(1+sigma/(n*m^2)))

eq2829c <- function(n,mu,sigma,m)return(m^2/4 * (sqrt(1+(16*mu*sigma)/(m^2))-1))

eq2833b <- function(n,mu,sigma,m){
    gamma=m^2/(2*sigma^2)
    return(1*m^2*(gamma*n+1)/(4*gamma*n)*(sqrt(1+2*(gamma*n*4*mu*n)/(gamma*n+1)^2)-1))
}
