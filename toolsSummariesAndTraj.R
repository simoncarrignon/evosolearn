### barplot version:


colsK=rev(colorRampPalette(c("red","blue"))(length(Ks)))

#for(i in folder){
i="deltas"
pdf(paste0(i,"_allvariance.pdf"),title=i,width=10,heigh=10,pointsize=12)
par(mfrow=c(4,2),mar=c(4,4,1,1))
#sigmascale=c(0.01,0.01,0.01,1.7)
#names(sigmascale)=sigmas
allscale=list()
for(sigma in unique(sumaries$sigma)){
    allscale[[as.character(sigma)]]=list()
    for(mu in mus){
        allsubset=getSubsetWithTraj(allsummaries[["unifPop50k"]],sigma=sigma,mu=mu,m=ms,E=E,K=Ks,traj=F)
        allmeans=tapply(allsubset$summary$var_x,allsubset$summary[,c("K","m")],mean)
        allpopMean=tapply(allsubset$summary$N,allsubset$summary[,c("K","m")],mean)
        allEstimates=allpopMean
        for(k in as.character(Ks))
            for(im in as.character(ms))
                allEstimates[k,im]=eq2833b(allpopMean[k,im],mu,sigma,as.numeric(im))
        allsds=tapply(allsubset$summary$var_x,allsubset$summary[,c("K","m")],sd)
        res=barplot(allmeans,beside=T,legend=F,col="white",plot=F)
        allscale[[as.character(sigma)]][[as.character(mu)]]=c(0,max(allmeans+allsds,allEstimates))
        ulul=plot(0:max(res),ylim=allscale[[as.character(sigma)]][[as.character(mu)]],type="n",axes=F,xlab="",ylab="")
        arrows(as.vector(res), allmeans+allsds, as.vector(res), allmeans-allsds,angle=90,code=3,length=.05,lwd=1.5,col=colsK)
        points(as.vector(res),as.vector(allEstimates),pch=21,bg="white",cex=1)
        axis(1,at=res[2,],labels=sapply(ms,function(i)as.expression(bquote(m[x] == .(i)))),lwd=0)

        axis(2)
        mtext(bquote(var[x]),2,2,cex=1)
        mtext(bquote(mu == .(mu) ~ sigma[s] == .(sigma) ),1,3)
        box()
        for(id in 1:length(deltas)){
            d=deltas[id]
            allsubset=getSubsetWithTraj(sumaries,sigma=sigma,mu=mu,m=ms,E=E,K=Ks,delta=d,traj=F)
            allmeans=tapply(allsubset$summary$var_x,allsubset$summary[,c("K","m")],mean)
            allpopMean=tapply(allsubset$summary$N,allsubset$summary[,c("K","m")],mean)
            allEstimates=allpopMean
            for(k in as.character(Ks))
                for(im in as.character(ms))
                    allEstimates[k,im]=eq2833b(allpopMean[k,im],mu,sigma,as.numeric(im))
            allsds=tapply(allsubset$summary$var_x,allsubset$summary[,c("K","m")],sd)
            nres=res+id/10
            arrows(as.vector(nres), allmeans+allsds, as.vector(nres), allmeans-allsds,angle=90,code=3,length=.05,lwd=1.5,col=colsK,lty=(id+1))
            points(as.vector(nres),as.vector(allEstimates),pch=21,bg="white",cex=1)
        }
    }
    legend("topleft",
           legend=c(paste0("K=",Ks),"Hermisson",sapply(c(0,deltas),function(d)as.expression(bquote(delta==.(d))))),
           col=c(colsK,1,rep(1,length(deltas)+1)),
           lty=c(rep(1,length(Ks)),NA,1,seq_along(deltas)+1),
           pch=c(rep(NA,length(Ks)),21,rep(NA,length(deltas)+1)),
           )
}
dev.off()



#this function use the filename stored in an summarized ouput  to get the individual trajectory for each experimetn an trace them.
plotAllTrajVar <- function(summaryresults,obs="N",m,E,ylim=NULL,legside="topleft"){

    Ks=sort(unique(summaryresults$K))
    mus=sort(unique(summaryresults$mu))
    ms=sort(unique(summaryresults$m))
    sigmas=sort(unique(summaryresults$sigma))
    deltas=sort(unique(summaryresults$delta))

    par(mfrow=c(length(sigmas),length(mus)))
    for(sigma in sigmas){
        for(mu in mus){

            subsetraj=getSubsetWithTraj(summaryresults,sigma=sigma,mu=mu,m=ms,E=E,K=Ks,delta=deltas,var=obs,traj=T)
            maxts_n=max(subsetraj$summary$extinction) # get real scale of observation

            if(is.null(ylim))ylim=range(subsetraj$traj,na.rm=T)
            param=bquote("full traj" ~ .(obs) ~ "with" ~ sigma[s] == .(sigma) ~ E == .(E) ~ m == .(m) ~ mu == .(mu) )

            #prepare main plot area
            plot(1,1,xlab="time",ylab=obs,ylim=ylim,xlim=c(1,maxts_n),type="n",main=param)
            legend(legside,
                   legend=c(paste0("K=",Ks),sapply(deltas,function(d)as.expression(bquote(delta==.(d))))),
                   col=c(cols,rep(length(deltas))),
                   lty=c(rep(1,length(Ks)),seq_along(deltas))
                   )

            for(id in 1:length(deltas)){
                d=deltas[id]

                singlesetup=lapply(Ks,function(k)getSubsetWithTraj(subsetraj$summary,sigma=sigma,mu=mu,m=m,E=E,K=k,delta=d,var=obs,traj=T)$traj)

                if(sum(lengths(singlesetup))==0)break
                #adjust length of single setup to ouptu
                maxts=nrow(singlesetup[[1]])
                rangesummaries=seq(1,maxts_n,length.out=maxts)

                a=sapply(singlesetup,quantile,na.rm=T)[c(2,4),]
                cols=colorRampPalette(c("red","blue"))(length(singlesetup))
                names(cols)=allsubmu
                lapply(1:length(allsubmu),function(sub)
                       {
                           plotTraj(rangesummaries,singlesetup[[sub]],col=cols[sub],xlim=range(rangesummaries),add=T,mean=T,lty=id,lwd=.8)
                       })
            }
        }

    }

}


## this funciton take all experiment from one or more folders and return it as a unique dataframe
getAlllSummaries <- function(folder){
    folder=c("omega1growingDelta")
    allsummaries=lapply(folder,getFullExperimentSummary)
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

updateScale <- function(df.exp){
    df.exp$extinction=(df.exp$extinction-1)*unique(df.exp$outputrate)+1
    return(df.exp)
}

