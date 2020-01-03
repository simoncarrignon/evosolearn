

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
                   plot(meanvar,ylab=varname,type="l",col=cols[varname],ylim=lims,xaxt='n',lwd=2,...)
                   lines(meanvar+sdvar,ylab=varname,col=cols[varname],lwd=2,lty=3)
                   lines(meanvar-sdvar,ylab=varname,col=cols[varname],lwd=2,lty=3)
               }
           })
    par(mar=c(2,4,0,1))
    plot(results$theta,type="l",col=cols["theta"],ylab="theta")
}
}

