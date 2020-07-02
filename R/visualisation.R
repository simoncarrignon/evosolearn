#' Plot single run
#' 
#' RGB Scale
#'
#' Function that gives the rgb scale of a results
#' 
#' @param results the output of evosolearn()
#' @export
getRgbScale <- function(results){
    res=rep(NA,nrow(results))
    nas=is.na(results[,"mean_y"])
    res[!nas]=rgb(results[!nas,"mean_y"],.5,results[!nas,"mean_z"])
return(res)
}


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
#' @param statvar a vector with the different variable to plot 
#' @param multi if true, add some comparison at the end
#' @param N if true, add effective pop size
#' @param addrgb if true, add rgb scale
#' @param theta if true, add the optimum list
#' @param year if not null label theta axis using the range of years
#' @param results the output of evosolearn()
#' @export
plotResults <- function(results,statfun=c("mean","sd"),statvar=c("w","p","x","y","z"),multi=F,N=T,theta=T,addrgb=T,year=NULL){
	defpar=par()
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
    par(mfrow=c(nlines,1),mar=rep(.5,4),oma=c(5,5,2,1),xaxs="i")

    if(typeof(results)=="list")
        allpop=T
    if(allpop){
        print("compute stats")
    }
    else{
        if(N){
            plot(results[,"N"],type="l",ylab="",main="",bty="n",xlab="t",axes=F,col="dark green")
            axis(2)
            mtext(expression(N[e]),2,3)
        }
        for(var in statvar){
            sdvar=results[,paste0("sd_",var)]
            meanvar=results[,paste0("mean_",var)]
            if(var %in% c("y","z","w"))yrange=c(0,1)
            else yrange=range(meanvar-sdvar,meanvar+sdvar,na.rm=T)
            plot(meanvar,type="l",ylim=yrange,ylab="",main="",bty="n",xlab="t",axes=F,col=vcol[var])
            axis(2)
            mtext(parse(text=paste0("bar(",var,")")),2,3,cex=.9)
            lines(meanvar+sdvar,lty=3,col=vcol[var])
            lines(meanvar-sdvar,lty=3,col=vcol[var])
			box()
        }
        if(addrgb){
            cls=matrix(nrow=1,ncol=nrow(results))
            cls[1,]=getRgbScale(results)
            plot.new()
            rasterImage(cls, 0, 0.2, 1, .8,interpolate=F,ylim=c(0,1))
            mtext("rgb(y,z)",2,1,cex=.8)
        }
        if(theta){
            x=seq_along(results[,"theta"])
            if(!is.null(year))
                x=seq(min(year),max(year),length.out=length(x))
            plot(x,results[,"theta"],xlim=range(x,0),type="l",ylab="",main="",xlab="t",axes=F,col="black")
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
    par(mfrow=defpar$mfrow,mar=defpar$mar,oma=defpar$oma,xaxs=defpar$xaxs)
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
        plot(results[,"theta"],type="l",col="blue",yaxt="n",xaxt="n",ylab="",bty="n",axes=F,xlab="")
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
    plot(x,theta,xlab="",type="l",col="red",yaxt="n",xaxt="n",bty="n",xlab="y",ylab="")
    axis(4,col="red",col.axis="red")
    mtext(expression(theta),4,2,col="red")
}


#' @param limit number of step (ie generation) to visualise 
#' @param ind boolean, if TRUE, individual run are ploted
#' @param subnb if ind is TRUE, number of uindividual run to plot by experiments
#' @param limit number of step (ie generation) to visualise 
#' @param prop realtive width of the environment plot wrt the rgb matrix
#' @param varcol withc variable use forthe color
#' @param palette the color palette to use
printRGBpixels <- function(data,filename,ind=FALSE,tlimit=NULL,subnb=NULL,img.width=600,img.height=800,img.pointsize=14,env=NULL,prop=.1,ordered=NULL,varcol="rgb",palette=colorRampPalette(c("grey","yellow","dark green"))(1200),nexpe=20)
{

    mum=as.data.frame(expand.grid(m=unique(data$m),mu=unique(data$mu)))
    mum$prod=mum$m*mum$mu
    lkz=length(unique(data$k_z))
    lky=length(unique(data$k_y))
    le=length(unique(data$E))
    lm=length(unique(data$m))
    lmu=length(unique(data$mu))
    summary=getUniqueExp(data$filename[1])
    tsteps=nrow(summary)
    if(is.null(tlimit)){
        tlimit=tsteps
    }
    else{ 
        tlimit=ceiling(tlimit/unique(data$outputrate))
    }
    tlim=min(tsteps,tlimit)

    if(is.null(subnb))subnb=nexpe

    if(ind)
        bicpic=matrix(nrow=tsteps,ncol=lkz*lky*lm*lmu*subnb+1)
    else
        bicpic=matrix(nrow=tlim,ncol=lkz*lky*lm*lmu)
    bicpic[,]=NA
    p=1
    for(ky in unique(data$k_y)){
        for(kz in unique(data$k_z)){
            for(i in 1:nrow(mum)){
                mu=mum$mu[i]
                m=mum$m[i]
                subb=droplevels(data[data$mu ==mu &data$m ==m &data$E ==e & data$k_z ==kz & data$k_y ==ky,])
                nexpe=nrow(subb)
                if(ind){
                    subselect=sample.int(nexpe,subnb)
                    if(is.null(ordered))ordered="n"
                    if(ordered=="z"){
                    reorder=sapply(subselect,function(f)sum(getUniqueExp(subb$filename[f])[,"mean_z"],na.rm=T))
                    subselect=subselect[order(reorder)]
                    }
                    if(ordered=="y"){
                    reorder=sapply(subselect,function(f)sum(getUniqueExp(subb$filename[f])[,"mean_y"],na.rm=T))
                    subselect=subselect[order(reorder)]
                    }
                    for(f in subselect){
                        pxl=c()
                        summary=getUniqueExp(subb$filename[f])
                        na=which.max(is.na(summary[,1]))#find the first na ie when pop get extinct
                        if(na>1) summary=summary[1:(na-1),,drop=F]
                        if(varcol=="rgb"){
                            alphalvl=sapply(summary[,"N"]/1000,function(i)min(i,1))
                            pxl=rgb(summary[,"mean_y" ],.5,summary[,"mean_z"],alphalvl)
                        }
                        else{
                            pxl=palette[summary[,varcol]]
                        }

                        bicpic[1:length(pxl),p]=pxl
                        p=p+1
                        print(paste("individual run",p,"/",ncol(bicpic)))
                    }
                }
                else{ 
                    if(varcol=="rgb") vars=c("mean_y","mean_z")
                    else vars=varcol
                    names(vars)=vars
                    mat_allexp= lapply(vars,function(i)matrix(NA,nexpe,tlim))
                    for(i in 1:nexpe){
                        summary=getUniqueExp(subb$filename[i])
                        rng=(nrow(summary)-tlim+1):nrow(summary)
                        for(v  in vars){
                            if(length(mat_allexp[[v]][i,])!=length(summary[rng,v]))
                                print("expe is not the number of step expected")
                            else
                                mat_allexp[[v]][i,]=summary[rng,v]

                        }
                    }
                    sum_mat=lapply(mat_allexp,function(m)apply(m,2,mean,na.rm=T))
                    if(varcol=="rgb")
                        sum_mat$na=apply(mat_allexp$mean_y,2,function(i)sum(is.na(i)))
                    else
                        sum_mat$na=apply(mat_allexp[[varcol]],2,function(i)sum(is.na(i)))
                    sum_mat=lapply(sum_mat,function(m)m[!is.na(m)])
                    if(length(sum_mat$mean_y)>1){
                        if(is.na(sum_mat$na[1]))
                            pxl=NA
                        else{
                            #pxl=rgb(sum_mat$mean_y,.5,sum_mat$mean_z,alpha=1)
                        if(varcol=="rgb")    
                            pxl=rgb(sum_mat$mean_y,.5,sum_mat$mean_z,alpha=1-sum_mat$na/nexpe)
                        else    
                            pxl=alpha(palette[sum_mat[[varcol]]],alpha=1-sum_mat$na/nexpe)
                        }
                        bicpic[1:length(pxl),p]=pxl
                    }

                    p=p+1
                    print(paste("col",p,"/",ncol(bicpic)))
                }
            }
            scale=1 #to convert time scale in year
            png(filename,width=img.width,height=img.height,pointsize=img.pointsize)
            par(mar=rep(.5,4),oma=c(3,4,4,1))

            layout(matrix(c(1,2),ncol=2,nrow=1),width=c(prop,1-prop))

            plotVertEnv(env,scale)
            plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
            rasterImage(bicpic, 0, 0, 1, 1,interpolate=F)

            ## first level legend (from left to right)
            kyl=sapply(unique(data$k_y),function(u)as.expression(bquote(k[y] == .(u))))
            lkyl=length(kyl)
            mtext(kyl,side=3,line=1,at=seq((1/lkyl)/2,1-(1/lkyl)/2,length.out=length(kyl)),cex=.9)

            #mtext(bquote(E==.(e)),side=3,line=2)
            mtext(bquote((E[x]== .1 ~E[y]== 1 ~E[z]== .5  )),side=3,line=2,cex=.8)

            ## second level legend (from left to n*1/n)
            kzl=sapply(unique(data$k_z),function(u)as.expression(bquote(k[z] == .(u))))
            lkzl=length(kzl)
            mtext(kzl,side=3,line=0,at=seq(0+(1/(lkzl*lkyl))/2,lkzl/(lkzl*lkyl)-(1/(lkzl*lkyl))/2,length.out=lkzl),cex=.9)

            #mtext(mum$prod[seq(1,nrow(mum),length.out=3)],side=1,line=0,at=seq(0,1/9,length.out=3),cex=.7)
            #third level
            slvl3=1/ncol(bicpic) #scale of lvl3 legend
            if(ind)
                slvl3=slvl3*subnb
            nlvl3=length(mum$prod) #number of legend
            llvl3=mum$prod         #labels of legend
            tlvl3=seq(1/2*slvl3,slvl3*nlvl3-1/2*slvl3,length.out=nlvl3) #tick place of legend
            axis(1,at=tlvl3,label=llvl3,cex.axis=.6,las=2)
            mtext(bquote(mu %*% m),side=1,line=2,at=(slvl3*nlvl3-1/2*slvl3)/2,cex=.6)
            dev.off()
        }
    }

}

plotVertEnv <- function(realdata,scale=1,...){
    plot((realdata$dTsVscales)*scale,-realdata$year,type="l",yaxs="i",axes=F,xlab="T")
    axis(3,cex.axis=.6)
    axis(2,label=rev(round(sort(seq(max(realdata$year),min(realdata$year),length.out=5)*scale))),at=seq(min(-realdata$year),max(-realdata$year),length.out=5),cex.axis=.8)
    mtext("year BP",2,2,at=.5,outer=T,cex=.8)
    mtext("dt",3,1.5,at=.05,outer=T,cex=.7)
    #mtext(expression(delta^18*O),3,0,at=prop/2,outer=T)
}


setAxis <- function(data){
    mus=sapply(unique(data$mu),function(u)as.expression(bquote(mu == .(u))))
    tm=min(5,length(mus))
    nl=seq(1,length(mus),length.out=tm)
    axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

    ms=paste0("m=",unique(data$m))
    tm=min(4,length(ms))
    nl=seq(1,length(ms),length.out=tm)
    axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)

    e=sapply(unique(data$E),function(u)as.expression(bquote(E == .(u))))
    mtext(e,side=3,line=0,at=seq(.5/length(e),1-0.5*1/length(e),length.out=length(e)),cex=.9,outer=T)

    vsp=1/nlines

    ky=sapply(unique(data$k_y),function(u)as.expression(bquote(k[y] == .(u))))
    mtext(ky,side=2,line=2,at=seq(1.5*vsp,1-1.5*vsp,length.out=length(ky)),cex=.9,outer=T)

    vsp=1/nlines
    kz=sapply(unique(data$k_z),function(u)as.expression(bquote(k[z] == .(u))))
    mtext(kz,side=2,line=0,at=seq(vsp-.5*vsp,length(kz)*vsp-.5*vsp,length.out=length(kz)),cex=.9,outer=T)
}
