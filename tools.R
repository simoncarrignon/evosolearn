

plotAllVariable <- function(results,hdr=F,vars=3:9){
    allsd=sapply(results$allpop,function(i)apply(i[,vars],2,sd))
    allmean=sapply(results$allpop,function(i)apply(i[,vars],2,mean))
    cols=c(rainbow(nrow(allsd)),"black")
    names(cols)=c(rownames(allsd),"theta")
    par(mfrow=c(8,1),mar=c(1,4,1,1))
    lapply(rownames(allsd),function(varname)
           {
	
			   if(hdr)hdr.boxplot(lapply(results$allpop,"[[",varname),border=NA,pch=".",outline=F,col=shades(cols[varname],3),prob=c(50,75,99),space=0,ylab=varname)
			   else {
				   lims=range(c(allmean[varname,]+allsd[varname,],allmean[varname,]-allsd[varname,]))
				   plot(allmean[varname,],ylab=varname,type="l",col=cols[varname],ylim=lims,xaxt='n')
				   lines(allmean[varname,]+allsd[varname,],ylab=varname,col=cols[varname],lty=3)
				   lines(allmean[varname,]-allsd[varname,],ylab=varname,col=cols[varname],lty=3)
			   }
		   })
    par(mar=c(2,4,0,1))
    plot(results$env,type="l",col=cols["theta"],ylab="theta")
}

