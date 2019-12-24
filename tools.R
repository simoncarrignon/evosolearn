

plotAllVariable <- function(results){
    allsd=sapply(results$allpop,function(i)apply(i[,3:9],2,sd))
    allmean=sapply(results$allpop,function(i)apply(i[,3:9],2,mean))
    cols=rainbow(nrow(allsd)+1)
    names(cols)=c(rownames(allsd),"theta")
    par(mfrow=c(8,1),mar=c(1,4,1,1))
    lapply(rownames(allsd),function(varname)
           {
               lims=range(c(allmean[varname,]+allsd[varname,],allmean[varname,]-allsd[varname,]))
               plot(allmean[varname,],ylab=varname,type="l",col=cols[varname],ylim=lims,xaxt='n')
               lines(allmean[varname,]+allsd[varname,],ylab=varname,col=cols[varname],lty=3)
               lines(allmean[varname,]-allsd[varname,],ylab=varname,col=cols[varname],lty=3)
           })
    par(mar=c(2,4,0,1))
    plot(results$env,type="l",col=cols["theta"],ylab="theta")
}



plotAllHDR <- function(results,var="w"){
    allsd=sapply(results$allpop,function(i)apply(i[,3:9],2,sd))
    cols=rainbow(nrow(allsd)+1)
    par(mfrow=c(8,1),mar=c(1,4,1,1))
    lapply(rownames(allsd),function(varname)
           {
               hdr.boxplot(lapply(results$allpop,"[[",varname),border=NA,pch=".",outline=F,col=c(cols[varname],"darkgreen"),space=0,ylab=varname)
           })
    par(mar=c(2,4,0,1))
    plot(results$env,type="l",col=cols["theta"],ylab="theta")
}

