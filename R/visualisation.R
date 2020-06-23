### Set some global variable
slpalette <- colorRampPalette(c(rgb(1,.5,0),rgb(0,.5,1)))
ilpalette <- colorRampPalette(c(rgb(0,.5,0),rgb(1,.5,0)))
survivalpallette <- colorRampPalette(c("grey","yellow","dark green"))
extinctpalette <- colorRampPalette(c("chartreuse4","white"))
pureSl <- rgb(0,.5,1)
pureIl <- rgb(1,.5,0)
pureGl <- rgb(0,.5,0)
###

#' Plot single run
#' 
#' Function that plot mean and standard for all values
#' 
#' @param restults the output of evosolearn()
#' @export
plotResults <- function(results,statfun=c("mean","sd"),statvar=c("w","p","x","y","z"),multi=F,N=T,theta=T,addrgb=T){
    allpop=F
    nlines=length(statvar)
    if(N)nlines=nlines+1
    if(theta)nlines=nlines+1
    if(addrgb)nlines=nlines+1
    vcol=rainbow(nlines)
    names(vcol)=statvar
    vcol["x"]=pureGl
    vcol["y"]=pureIl
    vcol["z"]=pureSl
    if(multi)nlines=nlines+2
    par(mfrow=c(nlines,1),mar=rep(.5,4),oma=c(5,5,2,1))

    if(typeof(results)=="list")
        allpop=T
    if(allpop){
        print("compute stats")
    }
    else{
        if(N){
            plot(results[,"N"],type="l",ylab="",main="",bty="n",xlab="t",axes=F,col="dark green")
            axis(2,xaxs="i")
            mtext(expression(N[e]),2,3)
        }
        for(var in statvar){
            sdvar=results[,paste0("sd_",var)]
            meanvar=results[,paste0("mean_",var)]
            if(var %in% c("y","z","w"))yrange=c(0,1)
            else yrange=range(meanvar-sdvar,meanvar+sdvar,na.rm=T)
            plot(meanvar,type="l",ylim=yrange,ylab="",main="",bty="n",xlab="t",axes=F,col=vcol[var])
            axis(2,xaxs="i")
            mtext(parse(text=paste0("bar(",var,")")),2,3,cex=.9)
            lines(meanvar+sdvar,lty=3,col=vcol[var])
            lines(meanvar-sdvar,lty=3,col=vcol[var])
        }
        if(addrgb){
            cls=matrix(nrow=1,ncol=nrow(results))
            cls[1,]=rgb(results[,"mean_y"],.5,results[,"mean_z"])
            plot.new()
            rasterImage(cls, 0, 0.2, 1, .8,interpolate=F,ylim=c(0,1))
            mtext("rgb(y,z)",2,0,cex=.8)
        }
        if(theta){
            plot(results[,"theta"],type="l",ylab="",main="",bty="n",xlab="t",axes=F,col="black")
            axis(2,xaxs="i")
            mtext(expression(theta),2,3)
        }
    }
    if(multi){
        plotDistThetaNW(results)
        plotThetaPhenotypes(results)
    }
    axis(1)
    mtext(expression(theta),2,3)
    mtext("Single run summary",3,0,outer=T)
}


#' Plot phenotype and  optimum
#' 
#' Function that plot mean phenotype a standard deviation
#' 
#' @param restults the output of evosolearn()
#' @export
plotThetaPhenotypes <- function(results,...){
    allpop=F
    par(mar=rep(4,4))
    if(typeof(results)=="list")
        allpop=T
    if(allpop){
        tstep=length(results$allpop)
        t=results
        plot(sapply(2:tstep,function(it)mean(t$allpop[[it]][,"p"])),type="l",ylim=range(sapply(t$allpop,"[[","p"),na.rm=T),ylab="p",main="Phenotype and theta",bty="n",...)
        lines(sapply(2:tstep,function(it)mean(t$allpop[[it]][,"p"])+sd(t$allpop[[it]][,"p"])),lty=3)
        lines(sapply(2:tstep,function(it)mean(t$allpop[[it]][,"p"])-sd(t$allpop[[it]][,"p"])),lty=3)
        par(new=T)
        plot(t$summary[,"theta"],type="l",col="blue",yaxt="n",xaxt="n",ylab="",bty="n",)
        axis(4,col="blue",col.axis="blue")
        mtext(expression(theta),4,2,col="blue")
    }
    else {
        sdr=results[,"sd_p"]
        meanr=results[,"mean_p"]
        plot(meanr,type="l",ylim=range(meanr-sdr,meanr+sdr,na.rm=T),ylab="p",main=expression(theta ~ "and p"),bty="n",xlab="t",...)
        lines(meanr+sdr,lty=3)
        lines(meanr-sdr,lty=3)
        par(new=T)
        plot(results[,"theta"],type="l",col="blue",yaxt="n",xaxt="n",ylim=range(meanr-sdr,meanr+sdr,na.rm=T),ylab="",bty="n",axes=F,xlab="")
        axis(4,col="blue",col.axis="blue")
        mtext(expression(theta),4,2,col="blue")
    }

}

#' Plot fitness and distance to optimum
#' 
#' Function that plot mean a standard deviation
#' 
#' @param restults the output of evosolearn()
#' @export
plotDistThetaNW <- function(results,...){
    allpop=F
    par(mar=rep(4,4))
    if(typeof(results)=="list")
        allpop=T
    if(allpop){
        plot( sapply(1:tstep,function(it)mean(abs(t$allpop[[it]]$p-t$env[it]))),ylim=range(sapply(1:tstep,function(it)range(abs(t$allpop[[it]]$p-t$env[it])))),ylab=expression(group("|",theta - p,"|")),bty="n",type="l",main="Distance to theta and fitness")
        lines(sapply(2:tstep,function(it)mean(abs(t$allpop[[it]]$p-t$env[it]))+sd(abs(t$allpop[[it]]$p-t$env[it]))),lty=3)
        lines(sapply(1:tstep,function(it)mean(abs(t$allpop[[it]]$p-t$env[it]))-sd(abs(t$allpop[[it]]$p-t$env[it]))),lty=3)
        par(new=T)
        plot( sapply(1:tstep,function(it)mean(t$allpop[[it]]$w)),ylim=range(sapply(t$allpop,"[[","w")),type="l",col="red",yaxt="n",xaxt="n",bty="n",ylab="")
        lines(sapply(1:tstep,function(it)mean(t$allpop[[it]]$w)+sd(t$allpop[[it]]$w)),col="red",lty=3)
        lines(sapply(1:tstep,function(it)mean(t$allpop[[it]]$w)-sd(t$allpop[[it]]$w)),col="red",lty=3)
        axis(4,col="red",col.axis="red")
        mtext("w",4,2,col="red")
    }
    else{
        sdr=results[,"sd_p"]
        meanr=results[,"mean_p"]
        sdw=results[,"sd_w"]
        meanw=results[,"mean_w"]
        theta=results[,"theta"]
        meanD=meanr-theta
        plot(meanD,type="l",ylim=range(meanD-sdr,meanD+sdr,na.rm=T),ylab=expression(group("|",theta - p,"|")),main=expression(theta - "p" ~ "and w"),bty="n",xlab="t",...)
        lines(meanD+sdr,lty=3)
        lines(meanD-sdr,lty=3)
        par(new=T)
        plot(meanw,type="l",col="red",yaxt="n",xaxt="n",bty="n",ylab="",axes=F,xlab="",ylim=c(0,1))
        lines(meanw+sdw,col="red",lty=3)
        lines(meanw-sdw,col="red",lty=3)
axis(4,col="red",col.axis="red")
mtext("w",4,2,col="red")
    }
}


#' 
plotAllVariable <- function(results,hdr=F,vars=NULL,theta=NULL,t=NULL,...){
    varnames=NULL
    allpop=F

    if(typeof(results)=="list")
        allpop=T
    if(!is.null(vars))varnames=vars
    if(allpop){
        allsd=sapply(results$allpop,function(i)apply(i[,vars],2,sd))
        allmean=sapply(results$allpop,function(i)apply(i[,vars],2,mean))
        if(is.null(vars))varnames=rownames(allsd)
    }
    else{
        if(is.null(vars))varnames=colnames(results)
    }
    cols=c(rainbow(length(varnames)))
    names(cols)=varnames
    par(mfrow=c(length(vars)+1,1),mar=c(0,4,2,4))
    for(i in 1:(length(varnames))){
        varname=varnames[i]
        if(hdr && allpop){
            yvalues=lapply(results$allpop,"[[",varname)
            ylim=range(yvalues)
            if(varname %in% c("y","z"))ylim=c(0,1)
            hdr.boxplot(yvalues,border=NA,pch=".",outline=F,col=cols[varname],prob=c(50,75,99),space=0,ylab=varname,ylim=ylim,h=.01,...)
            if(varname %in% c("x","ilp","gp","p") & !is.null(theta))
                lines(theta,col="red",lwd=1,lty=1)
        }
        else {
            if(allpop){
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
        par(mar=c(0,4,0,4))
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
