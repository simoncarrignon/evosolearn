    folder=c("growingDelta","unifPop50k")
    allsummaries=lapply(folder,getFullExperimentSummary)
    names(allsummaries)=folder

    Ks=unique(allsummaries[[1]]$K)
    mus=unique(allsummaries[[1]]$mu)
    ms=unique(allsummaries[[1]]$m)
    sigmas=unique(allsummaries[[1]]$sigma)
    deltas=unique(allsummaries[[1]]$delta)
    E=0


    scale=list(var_x=list("0.1"=c(0,0.004),"0.2"=c(0,0.002),"0.4"=c(0,0.005),"10000"=c(0,0.15)),
               N=c(200,2100),
               mean_w=c(.2,1)
               )
    side=list(var_x="topleft",N="bottomright",mean_w="bottomright")
#for(exp in folder[1]){
    sumaries=allsummaries[["growingDelta"]]
    for(obs in c("var_x","mean_w","N")){
        tmp=exp
        if(exp=="sixsixequalone")tmp="startzero"
        tmp="delta"
        pdf(paste0(tmp,"_traj",obs,".pdf"),width=10,height=12)
        par(mfrow=c(4,2))
        for(sigma in unique(sumaries$sigma)){
            for(mu in unique(sumaries$mu)){
                ### plot traj for delta=0
                allsubset=getSubsetWithTraj(allsummaries[["unifPop50k"]],sigma=sigma,mu=mu,m=ms,E=E,K=Ks,traj=T)
                maxts=nrow(allsubset$traj)
                rangesummaries=seq(1,maxts,length.out=maxts_n)
                options(scipen=999)
                allsubmu=lapply(Ks,function(k)getSubsetWithTraj(allsummaries[["unifPop50k"]],sigma=sigma,mu=mu,m=m,E=E,K=k,delta=0,var=obs)$traj[rangesummaries,])
                a=sapply(allsubmu,quantile)[c(2,4),]
                param=bquote("full traj" ~ .(obs) ~ "with" ~ sigma[s] == .(sigma) ~ E == .(E) ~ m == .(m) ~ mu == .(mu) )
                if(is.null(scale[obs]))
                    ylim=range(a[1,],a[2,])
                else{
                    if(is.list(scale[[obs]]))
                       ylim=scale[[obs]][[as.character(sigma)]]
                   else
                       ylim=scale[[obs]]
                }
                plot(1,1,type="n",xlab="time step",ylab=obs,ylim=ylim,xlim=range(rangesummaries),main=param)
                cols=colorRampPalette(c("red","blue"))(length(allsubmu))
                names(cols)=allsubmu
                lapply(1:length(allsubmu),function(sub)
                       {
                           plotTraj(rangesummaries,allsubmu[[sub]],col=cols[sub],xlim=range(rangesummaries),add=T,mean=T,lwd=.8)
                       })
                legend(side[[obs]],
                       legend=c(paste0("K=",Ks),sapply(c(0,deltas),function(d)as.expression(bquote(delta==.(d))))),
                       col=c(cols,rep(1,length(deltas)+1)),
                       lty=c(rep(1,length(Ks)),1,seq_along(deltas)+1)
                       )
                ### plot traj for all other deltas
                for(id in 1:length(deltas)){
                    d=deltas[id]
                    allsubset=getSubsetWithTraj(allsummaries[["growingDelta"]],sigma=sigma,mu=mu,m=ms,E=E,K=Ks,delta=d,traj=T)
                    maxts_n=nrow(allsubset$traj)
                    rangesummaries=1:maxts_n
                    rangesummaries_scale=seq(1,maxts,length.out=maxts_n)
                    allsubmu=lapply(Ks,function(k)getSubsetWithTraj(allsummaries[["growingDelta"]],sigma=sigma,mu=mu,m=m,E=E,K=k,delta=d,var=obs,traj=T)$traj[rangesummaries,])
                    a=sapply(allsubmu,quantile,na.rm=T)[c(2,4),]
                    cols=colorRampPalette(c("red","blue"))(length(allsubmu))
                    names(cols)=allsubmu
                    lapply(1:length(allsubmu),function(sub)
                           {
                               plotTraj(rangesummaries_scale,allsubmu[[sub]],col=cols[sub],xlim=range(rangesummaries),add=T,mean=T,lty=id,lwd=.8)
                           })
                }
            }
        }
        dev.off()
    }
#}


    ### barplot version:
    

    colsK=rev(colorRampPalette(c("red","blue"))(length(Ks)))

    #for(i in folder){
    i="deltas"
        pdf(paste0(i,"_allvariance.pdf"),title=i,width=10,heigh=10,pointsize=12)
        par(mfrow=c(4,2),mar=c(4,4,1,1))
        #sigmascale=c(0.01,0.01,0.01,1.7)
        sumaries=allsummaries[["growingDelta"]]
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
    #}



