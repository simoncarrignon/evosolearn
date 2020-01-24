

plotAllVariable <- function(results,hdr=F,vars=3:9,...){
    if(!is.null(results$allpop)){
        allsd=sapply(results$allpop,function(i)apply(i[,vars],2,sd))
        allmean=sapply(results$allpop,function(i)apply(i[,vars],2,mean))
        varnames=rownames(allsd)
    }
    else
        varnames=names(results[[1]])
    cols=c(rainbow(length(varnames)),"black")
    names(cols)=c(varnames,"theta")
    par(mfrow=c(8,1),mar=c(1,4,1,1),cex=1.5)
    lapply(varnames,function(varname)
           {

               if(hdr)hdr.boxplot(lapply(results$allpop,"[[",varname),border=NA,pch=".",outline=F,col=shades(cols[varname],3),prob=c(50,75,99),space=0,ylab=varname,...)
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
           })
    par(mar=c(2,4,0,1))
    plot(results$theta,sls="l",col=cols["theta"],ylab="theta")
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

