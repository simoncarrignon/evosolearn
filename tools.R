

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

#wrappers for exploration
replicateNTime <- function(repet,n,tstep,omega,delta,b,K,mu,E,sigma,pop,m){

    do.call("rbind",replicate(repet,getVarXMeanW(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m),simplify=F))
}


getVarXMeanW <- function(n,tstep,omega,delta,b,K,mu,E,sigma,pop,m){
    t=simpleEvoModel(n=n,tstep=tstep,omega = omega,delta = delta ,b=b,K=K,mu=mu,E=E,sigma=sigma,pop=pop,m=m)
    if(is.null(t$sd$x[tstep]))
        return(c("var(x)"=NA,"mean(w)"=NA))
    else
        return(c("var(x)"=t$sd$x[tstep],"mean(w)"=t$mean$w[tstep]))
}

writeResults <- function(final,var,gene){
    pdf(paste0("explore_",var,"_",gene,".pdf"),width=10)
    par(mfrow=c(1,2))
    plot(explore[[var]][,3],final[,1],log="x",xlab=colnames(final)[3],ylab=colnames(final)[1],pch=20,col=alpha(1,.2))
    plot(final[,3],final[,2],log="x",xlab=colnames(final)[3],ylab=colnames(final)[2],pch=20,col=alpha(1,.2))
    dev.off()
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

