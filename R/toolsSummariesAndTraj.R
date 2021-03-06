
#' @export
plotAllVariableSummaries <- function(summaryresults,E,estimate=NULL,ylim=NULL,var="var_x"){
    Ks=sort(unique(summaryresults$K))
    mus=sort(unique(summaryresults$mu))
    ms=sort(unique(summaryresults$m))
    sigmas=sort(unique(summaryresults$sigma))
    deltas=sort(unique(summaryresults$delta))
    omegas=sort(unique(summaryresults$omega))
    if(is.null(omegas))omega=0
    colsK=rev(colorRampPalette(c("red","blue"))(length(Ks)))

    par(mfrow=c(length(sigmas),length(mus)),mar=rep(1,4),oma=c(4,4,1,4))

    for(sigma in sigmas){
        for(mu in mus){
            subset=getSubsetWithTraj(summaryresults,sigma=sigma,mu=mu,m=ms,E=E,K=Ks,delta=deltas,traj=F)$summary
            subdeltas=sort(unique(subset$delta))
            getallmetrics=lapply(deltas,function(d)
                                 {
                                     singlesetup=getSubsetWithTraj(subset,sigma=sigma,mu=mu,m=ms,E=E,K=Ks,delta=d,traj=F)$summary
                                     if(sum(lengths(singlesetup))==0)return(NULL)
                                     allmeans=tapply(singlesetup[,var],singlesetup[,c("K","m")],mean,na.rm=T)
                                     allsds=tapply(singlesetup[,var],singlesetup[,c("K","m")],sd,na.rm=T)
                                     allestimates=NULL
                                     noselection=NULL
                                     if(!is.null(estimate)){
                                         allpopMean=tapply(singlesetup$N,singlesetup[,c("K","m")],mean,na.rm=T)
                                         allestimates=allpopMean
                                         if(sigma>1000)
                                             noselection=allpopMean
                                         for(k in as.character(Ks))
                                             for(im in as.character(ms)){
                                                 allestimates[k,im]=estimate(allpopMean[k,im],mu,sigma,as.numeric(im))
                                                 if(sigma>1000)
                                                     noselection[k,im]=eq2833a(allpopMean[k,im],mu,sigma,as.numeric(im))
                                             }
                                     }
                                     allranges=range(allmeans+allsds,allmeans-allsds,allestimates)
                                     return(list(mean=allmeans,sd=allsds,range=allranges,estimate=allestimates,noselection=noselection)) 
                                 }
            )
            if(sum(lengths(getallmetrics))==0)break
            getallmetrics=getallmetrics[lengths(getallmetrics)>0]
            res=barplot(getallmetrics[[1]]$mean,beside=T,legend=F,col="white",plot=F)
            if(is.null(ylim)){
            ranges=range(lapply(getallmetrics,"[[","range"),na.rm=T)
            if(any(is.infinite(ranges)))ranges=c(-1,1)
            }
            else{ ranges=ylim}
            plot(0:max(res),ylim=ranges,type="n",axes=F,xlab="",ylab="")
            axis(1,at=res[2,],labels=sapply(ms,function(i)as.expression(bquote(m[x] == .(i)))),lwd=0)

            axis(2)
            mtext(var,2,2,cex=1)
            mtext(bquote(mu == .(mu) ~ sigma[s] == .(sigma) ),3,0)
            box()
            for(id in 1:length(getallmetrics)){
                d=getallmetrics[[id]]
                nres=res+id/10
                if(!is.null(estimate)){
                    points(as.vector(nres),as.vector(d$estimate),pch=21,bg="white",cex=1)
                    if(sigma>1000) points(as.vector(nres),as.vector(d$noselection),pch=22,bg="white",cex=1)
                }
                arrows(as.vector(nres), d$mean+d$sd, as.vector(nres), d$mean-d$sd,angle=90,code=3,length=.01,lwd=1.5,col=colsK,lty=id)
            }
            if(!is.null(estimate)){
                if(sigma>1000){
                    legend("topleft",
                           legend=c(paste0("K=",Ks),"Hermisson","no selection",sapply(subdeltas,function(d)as.expression(bquote(delta==.(d))))),
                           col=c(colsK,1,1,rep(1,length(subdeltas))),
                           lty=c(rep(1,length(Ks)),NA,NA,seq_along(subdeltas)),
                           pch=c(rep(NA,length(Ks)),21,22,rep(NA,length(subdeltas))),
                           )
                }
                else{
                    legend("topleft",
                           legend=c(paste0("K=",Ks),"Hermisson",sapply(subdeltas,function(d)as.expression(bquote(delta==.(d))))),
                           col=c(colsK,1,rep(1,length(subdeltas))),
                           lty=c(rep(1,length(Ks)),NA,seq_along(subdeltas)),
                           pch=c(rep(NA,length(Ks)),21,rep(NA,length(subdeltas))),
                           )
                }
            }
            else{
                legend("topleft",
                       legend=c(paste0("K=",Ks),sapply(subdeltas,function(d)as.expression(bquote(delta==.(d))))),
                       col=c(colsK,rep(1,length(subdeltas))),
                       lty=c(rep(1,length(Ks)),seq_along(subdeltas)),
                       pch=c(rep(NA,length(Ks)),rep(NA,length(subdeltas))),
                       )
            }
        }
    }

}


#this function use the filename stored in an summarized ouput  to get the individual trajectory for each experimetn an trace them.
#' @export
plotAllTrajVar <- function(summaryresults,obs="N",m,E,ylim=NULL,legside="topleft"){
    Ks=sort(unique(summaryresults$K))
    mus=sort(unique(summaryresults$mu))
    ms=sort(unique(summaryresults$m))
    sigmas=sort(unique(summaryresults$sigma))
    deltas=sort(unique(summaryresults$delta))
    omega=sort(unique(summaryresults$omega))
    if(is.null(omega))omega=0

    par(mfrow=c(length(sigmas),length(mus)))

    #different colos will represent different Ks
    cols=colorRampPalette(c("red","blue"))(length(Ks))
    names(cols)=Ks

    for(sigma in sigmas){
        for(mu in mus){

            subsetraj=getSubsetWithTraj(summaryresults,sigma=sigma,mu=mu,m=m,E=E,K=Ks,delta=deltas,var=obs,traj=T)
            param=bquote(sigma[s] == .(sigma) ~ E == .(E) ~ m == .(m) ~ mu == .(mu) ~ omega == .(omega) )

            if(is.null(dim(subsetraj$traj))){#some simulation have differen lenght
                minlen=min(lengths(subsetraj$traj))
                subsetraj$traj=sapply(subsetraj$traj,function(i){
                            if(length(i)>minlen){
                                rescale=seq(1,length(i),length.out=minlen)
                                return(i[rescale])
                            }
                            return(i)
               })
            }

            a=apply(subsetraj$traj,1,function(l)quantile(l[(length(l)/2):length(l)],na.rm=T))
            maxts_n=max(subsetraj$summary$extinction) # get real scale of observation

            if(is.null(ylim))yl=range(a,na.rm=T)
            else yl =ylim
            if(any(is.infinite(yl)))yl=c(-1,1)

            #param=bquote("full traj" ~ .(obs) ~ "with" ~ sigma[s] == .(sigma) ~ E == .(E) ~ m == .(m) ~ mu == .(mu) )

            #prepare main plot area

            plot(1,1,xlab="time",ylab=obs,ylim=yl,xlim=c(1,maxts_n),type="n",main=param)
            for(id in 1:length(deltas)){
                d=deltas[id]

                singlesetup=lapply(Ks,function(k)getSubsetWithTraj(subsetraj$summary,sigma=sigma,mu=mu,m=m,E=E,K=k,delta=d,var=obs,traj=T)$traj)

                if(sum(lengths(singlesetup))==0)break
                #adjust length of single setup to ouptut
                maxts=nrow(singlesetup[[1]])
                rangesummaries=seq(1,maxts_n,length.out=maxts)

                lapply(1:length(Ks),function(sub)
                       {
                           plotTraj(rangesummaries,singlesetup[[sub]],col=cols[sub],xlim=range(rangesummaries),add=T,mean=T,lty=id,lwd=.8)
                       })
            }
            legend(legside,
                   legend=c(paste0("K=",Ks),sapply(deltas,function(d)as.expression(bquote(delta==.(d))))),
                   col=c(cols,rep(1,length(deltas))),
                   lty=c(rep(1,length(Ks)),seq_along(deltas))
                   )
        }

    }

}


## this funciton take all experiment from one or more folders and return it as a unique dataframe
#' @export
getAlllSummaries <- function(folder,exclude=NULL){
    allsummaries=lapply(folder,getFullExperimentSummary,exclude=exclude)
    names(allsummaries)=folder
    allexperiments=lapply(allsummaries,function(curr)
                          {
                              allTl=lapply(curr$filename,getTraj,var="N")
                              allLen=sapply(allTl,function(i)length(i[!is.na(i)]))
                              cbind(curr,extinction=allLen)
                          }
    )
    df.allexp=do.call("rbind",allexperiments)
    return(df.allexp)
}

#' @export
updateScale <- function(df.exp){
    df.exp$extinction=(df.exp$extinction-1)*unique(df.exp$outputrate)+1
    return(df.exp)
}


#plot tow lines: distance to optimum and fintes
#' @export
plotDistVsFitness <- function(summaryresults,ylim=NULL,legside="topleft",...){
    Ks=sort(unique(summaryresults$K))
    mus=sort(unique(summaryresults$mu))
    ms=sort(unique(summaryresults$m))
    sigmas=sort(unique(summaryresults$sigma))
    deltas=sort(unique(summaryresults$delta))
    omega=sort(unique(summaryresults$omega))

    subsetraj=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="dist",traj=T)
    subsetrajFitness=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="mean_w",traj=T)

    if(is.null(dim(subsetraj$traj))){#some simulation have differen lenght
        minlen=min(lengths(subsetraj$traj))
        subsetraj$traj=sapply(subsetraj$traj,function(i){
                              if(length(i)>minlen){
                                  rescale=seq(1,length(i),length.out=minlen)
                                  return(i[rescale])
                              }
                              return(i)
    })
    }

    a=apply(subsetraj$traj,1,function(l)quantile(l,na.rm=T))
    b=apply(subsetrajFitness$traj,1,function(l)quantile(l,na.rm=T))
    maxts_n=max(subsetraj$summary$extinction) # get real scale of observation

                maxts=nrow(subsetraj$traj)
                rangesummaries=seq(1,maxts_n,length.out=maxts)

    if(is.null(ylim))yl=range(a[c(2,4),],na.rm=T)
    else yl =ylim
    if(any(is.infinite(yl)))yl=c(-1,1)

    #param=bquote("full traj" ~ .(obs) ~ "with" ~ sigma[s] == .(sigma) ~ E == .(E) ~ m == .(m) ~ mu == .(mu) )

    #prepare main plot area
    par(mar=c(2,4,2,4))
    plot(1,1,xlab="time",ylim=yl,xlim=c(1,maxts_n),type="n",ylab=expression(group("|",theta - p,"|")))
    qts=a
    lines(rangesummaries,qts[3,],col="black",lwd=1,lty=1)
    lines(rangesummaries,qts[4,],lty=2,col="black",lwd=.5)
    lines(rangesummaries,qts[2,],lty=2,col="black",lwd=.5)
    par(new=T)
    if(is.null(ylim))yl=range(b,na.rm=T)
    else yl =ylim
    if(any(is.infinite(yl)))yl=c(-1,1)
    plot(1,1,xlab="",type="l",ylim=yl,xlim=c(1,maxts_n),col="red",yaxt="n",xaxt="n",bty="n",ylab="",...)
    qts=b
    axis(4,col="red",col.axis="red")
    mtext("w",4,2,col="red")
    lines(rangesummaries,qts[3,],col="red",lwd=1,lty=1)
    lines(rangesummaries,qts[4,],lty=2,col="red",lwd=.5)
    lines(rangesummaries,qts[2,],lty=2,col="red",lwd=.5)


}



#plot tow lines: distance to optimum and fintes
#' @export
plot3Genes <- function(summaryresults,ylim=NULL,side="topright",...){
    Ks=sort(unique(summaryresults$K))
    mus=sort(unique(summaryresults$mu))
    ms=sort(unique(summaryresults$m))
    sigmas=sort(unique(summaryresults$sigma))
    deltas=sort(unique(summaryresults$delta))
    omega=sort(unique(summaryresults$omega))

    subsetrajX=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="mean_x",traj=T)
    subsetrajY=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="mean_y",traj=T)
    subsetrajZ=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="mean_z",traj=T)
    subsetrajT=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="distX",traj=T)

        subsetraj=subsetrajX
    if(is.null(dim(subsetraj$traj))){#some simulation have differen lenght
        minlen=min(lengths(subsetraj$traj))
        subsetraj$traj=sapply(subsetraj$traj,function(i){
                              if(length(i)>minlen){
                                  rescale=seq(1,length(i),length.out=minlen)
                                  return(i[rescale])
                              }
                              return(i)
    })
    }

    x=apply(subsetrajX$traj,1,function(l)quantile(l,na.rm=T))
    y=apply(subsetrajY$traj,1,function(l)quantile(l,na.rm=T))
    z=apply(subsetrajZ$traj,1,function(l)quantile(l,na.rm=T))
    t=apply(subsetrajT$traj,1,function(l)quantile(l,na.rm=T))
    maxts_n=max(subsetrajX$summary$extinction) # get real scale of observation

    maxts=nrow(subsetrajX$traj)
    rangesummaries=seq(1,maxts_n,length.out=maxts)

    if(is.null(ylim))yl=range(z[c(2,4),],y[c(2,4),],na.rm=T)
    else yl =ylim
    if(any(is.infinite(yl)))yl=c(-1,1)

    #prepare main plot area
    par(mar=c(2,4,2,4))
    plot(1,1,xlab="time",ylim=yl,xlim=c(1,maxts_n),type="n",ylab="learning",...)
    lines(rangesummaries,y[3,],col="black",lwd=1,lty=1)
    lines(rangesummaries,y[4,],lty=2,col="black",lwd=.5)
    lines(rangesummaries,y[2,],lty=2,col="black",lwd=.5)
    lines(rangesummaries,z[3,],col="green",lwd=1,lty=1)
    lines(rangesummaries,z[4,],lty=2,col="green",lwd=.5)
    lines(rangesummaries,z[2,],lty=2,col="green",lwd=.5)
    legend(side,legend=c("individual","social","genetic"),lty=1,col=c("black","green","red"))
    par(new=T)
    if(is.null(ylim))yl=range(t[c(2,4),],na.rm=T)
    else yl =ylim
    if(any(is.infinite(yl)))yl=c(-1,1)
    plot(1,1,xlab="",type="l",ylim=yl,xlim=c(1,maxts_n),col="red",yaxt="n",xaxt="n",bty="n",ylab="")
    axis(4,col="red",col.axis="red")
    mtext(expression(group("|",theta - x,"|")),4,3,col="red")
    lines(rangesummaries,t[3,],col="red",lwd=1,lty=1)
    lines(rangesummaries,t[4,],lty=2,col="red",lwd=.5)
    lines(rangesummaries,t[2,],lty=2,col="red",lwd=.5)


}


#plot tow lines: distance to optimum and fintes
#' @export
plotProp <- function(summaryresults,ylim=NULL,side="topright",...){
    Ks=sort(unique(summaryresults$K))
    mus=sort(unique(summaryresults$mu))
    ms=sort(unique(summaryresults$m))
    sigmas=sort(unique(summaryresults$sigma))
    deltas=sort(unique(summaryresults$delta))
    omega=sort(unique(summaryresults$omega))

    subsetrajY=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="prop_y",traj=T)
    subsetrajZ=getSubsetWithTraj(summaryresults,sigma=sigmas,mu=mus,m=ms,E=0,K=Ks,delta=deltas,var="prop_z",traj=T)

        subsetraj=subsetrajX
    if(is.null(dim(subsetraj$traj))){#some simulation have differen lenght
        minlen=min(lengths(subsetraj$traj))
        subsetraj$traj=sapply(subsetraj$traj,function(i){
                              if(length(i)>minlen){
                                  rescale=seq(1,length(i),length.out=minlen)
                                  return(i[rescale])
                              }
                              return(i)
    })
    }

    y=apply(subsetrajY$traj,1,function(l)quantile(l,na.rm=T))
    z=apply(subsetrajZ$traj,1,function(l)quantile(l,na.rm=T))
    maxts_n=max(subsetrajY$summary$extinction) # get real scale of observation

    maxts=nrow(subsetrajY$traj)
    rangesummaries=seq(1,maxts_n,length.out=maxts)

    if(is.null(ylim))yl=range(z[c(2,4),],y[c(2,4),],na.rm=T)
    else yl =ylim
    if(any(is.infinite(yl)))yl=c(-1,1)

    #prepare main plot area
    par(mar=c(2,4,2,4))
    plot(1,1,xlab="time",ylim=yl,xlim=c(1,maxts_n),type="n",ylab="learning",...)
    lines(rangesummaries,y[3,],col="black",lwd=1,lty=1)
    lines(rangesummaries,y[4,],lty=2,col="black",lwd=.5)
    lines(rangesummaries,y[2,],lty=2,col="black",lwd=.5)
    lines(rangesummaries,z[3,],col="green",lwd=1,lty=1)
    lines(rangesummaries,z[4,],lty=2,col="green",lwd=.5)
    lines(rangesummaries,z[2,],lty=2,col="green",lwd=.5)
    legend(side,legend=c("individual","social"),lty=1,col=c("black","green"))
}



#this function use the filename stored in an summarized ouput  to get the individual trajectory for each experimetn an trace them.
#' @export
plotAllProp <- function(summaryresults,obs="prop",m,E,ylim=NULL,legside="topleft"){
    Ks=sort(unique(summaryresults$K))
    mus=sort(unique(summaryresults$mu))
    ms=sort(unique(summaryresults$m))
    sigmas=sort(unique(summaryresults$sigma))
    deltas=sort(unique(summaryresults$delta))
    omega=sort(unique(summaryresults$omega))
    if(is.null(omega))omega=0

    par(mfrow=c(length(sigmas),length(mus)))

    #different colos will represent different Ks
    cols=colorRampPalette(c("red","blue"))(length(Ks))
    names(cols)=Ks

    for(sigma in sigmas){
        for(mu in mus){

            subsetrajY=getSubsetWithTraj(summaryresults,sigma=sigma,mu=mu,m=m,E=E,K=Ks,delta=deltas,var="prop_y",traj=T)
            param=bquote(sigma[s] == .(sigma) ~ E == .(E) ~ m == .(m) ~ mu == .(mu) ~ omega == .(omega) )

            maxts_n=max(subsetrajY$summary$extinction) # get real scale of observation

            #prepare main plot area

            plot(1,1,xlab="time",ylab=obs,ylim=c(0,1),xlim=c(1,maxts_n),type="n",main=param)
            for(id in 1:length(deltas)){
                d=deltas[id]

                singlesetupY=lapply(Ks,function(k)getSubsetWithTraj(subsetrajY$summary,sigma=sigma,mu=mu,m=m,E=E,K=k,delta=d,var="prop_y",traj=T)$traj)
                singlesetupZ=lapply(Ks,function(k)getSubsetWithTraj(subsetrajY$summary,sigma=sigma,mu=mu,m=m,E=E,K=k,delta=d,var="prop_z",traj=T)$traj)
                singlesetupYZ=lapply(Ks,function(k)getSubsetWithTraj(subsetrajY$summary,sigma=sigma,mu=mu,m=m,E=E,K=k,delta=d,var="prop_yz",traj=T)$traj)

                if(sum(lengths(singlesetupY))==0)break
                #adjust length of single setup to ouptut
                maxts=nrow(singlesetupY[[1]])
                rangesummaries=seq(1,maxts_n,length.out=maxts)

                plotTraj(rangesummaries,singlesetupY[[1]],col="red",xlim=range(rangesummaries),add=T,lty=1,lwd=.8)
                plotTraj(rangesummaries,singlesetupYZ[[1]],col="orange",xlim=range(rangesummaries),add=T,lty=1,lwd=.8)
                plotTraj(rangesummaries,singlesetupZ[[1]],col="green",xlim=range(rangesummaries),add=T,lty=1,lwd=.8)
            }
            legend(legside,
                   legend=c("individual","social","mixed",sapply(deltas,function(d)as.expression(bquote(delta==.(d))))),
                   col=c("green","red","orange",rep(1,length(deltas))),
                   lty=c(rep(1,3),seq_along(deltas))
                   )
        }

    }

}
