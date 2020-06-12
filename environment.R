
#Function from tuneR package
#Based on Timmer and Koening 1995: https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract
#' @param alpha: the slope of the power distribution (called omega in whitehead)
#' @param N: the length of the timeseri to 
TK95 <- function(N, alpha = 1){ 
    f <- seq(from=0, to=pi, length.out=(N/2+1))[-c(1,(N/2+1))] # Fourier frequencies
    f_ <- 1 / f^alpha # Power law
    RW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the real part
    IW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the imaginary part
    fR <- complex(real = c(rnorm(1), RW, rnorm(1), RW[(N/2-1):1]), 
                  imaginary = c(0, IW, 0, -IW[(N/2-1):1]), length.out=N)
     # Those complex numbers that are to be back transformed for Fourier Frequencies 0, 2pi/N, 2*2pi/N, ..., pi, ..., 2pi-1/N 
     # Choose in a way that frequencies are complex-conjugated and symmetric around pi 
     # 0 and pi do not need an imaginary part
    reihe <- fft(fR, inverse=TRUE) # go back into time domain
    return(Re(reihe)) # imaginary part is 0
}

#return spectrum of a timeserie
getSpectrum <- function(x){
	N <- length(x)
    pw=8
	M <- 2^pw # zero-pad total length
    while(M<N){
        pw=pw+1
        M <- 2^pw
    }
	freq <- seq(0, 0.5, by = 1/M)
	x.zp <- c(x, rep(0, M-N))
	S.pgram <- (1/N)*abs(fft(x.zp)[1:(M/2+1)])^2
    S.pgram
}

#return a 1/f environmen with curve = omega and sd = delta and increase rate vt
environment <- function(N,omega,delta,vt=NULL){
    ts=TK95(N,omega)
    ts=delta*ts/sd(ts)
    if(!is.null(vt))ts = ts + vt*1:N
    return(ts)
}


gauss <- function(N,mean,delta,v){
    sapply(1:N,function(t)rnorm(1,mean*v,delta))
}


exploreEnvironmentsProperties <- function(){

    #get mean correlation
    mean(replicate(10,{
                   t=environment(10000,2,1)
                   vtm=t[1:(length(t)-1)]
                   vt=t[2:length(t)]
                   cov(vt,vtm)}))


    pdf("env_corel.pdf",width=8,height=4.5)
    par(mfrow=c(1,3))

    for(omega in 0:2){
        t=environment(50000,omega,.1)
        vtm=t[1:(length(t)-1)]
        vt=t[2:length(t)]
        print(cov(vt,vtm))
        cv=round(cov(vt,vtm),digit=3)
        plot((1:50000)[seq(1,50000,100)],t[seq(1,50000,100)],main=bquote(omega == .(omega)~","~delta==.1~", cov("*theta[t]*","*theta[t-1]*")" == .(cv)),ylab=bquote(theta[t]),,xlab="t",type="l",lwd=1)
    }
    dev.off()
}



#Function interpolate
#To add intermediate value between two theta (using straight linear interpolation)
#' @param theta: the original theta vector to complete
#' @param times: the time associated to the thetas
#' @param finalres: the resolution, in the unit used in the vector times, at which new data had to be generate. if times is in year and finalres is 1, the vector thetat will correspond to data ofr each eyar, if time is in thousand of years (as in LR04 stack) and finalres is 0.5 thus the output will correspond to one point every 500 years.
#' @param omega: autocorrelation of the noise use to generate interpolate datapoints
#' @param delta: variance of the noise use to generate interpolate datapoints
#' @return a vector of dimension: seq(1,length(times),by=finalres)
interpolate <- function(theta,times,finalres,delta=0,omega=0){
    newtheta=c()
    for(i in 1:(length(times)-1)){
        rangesyears=seq(times[i],(times[i+1]),by=finalres)
        newt=seq(theta[i],theta[i+1],length.out=length(rangesyears))
        names(newt)=rangesyears
        if(omega > 0 && delta >0){ #add generated noise between unknown data points
            noise=environment(length(newt),omega=omega,delta=delta) # we generate noise between the two datapoints, including them
            newt[2:(length(newt-1))]=newt[2:(length(newt-1))]+noise #we add the generated noise only to the new datapoints
        }
        newtheta=c(newtheta,newt[-length(newt)])
    }
    newtheta=c(newtheta,theta[length(theta)])
    names(newtheta)[length(newtheta)]=times[length(times)]
    year=as.numeric(names(newtheta))
    names(newtheta)=NULL
    return(cbind.data.frame(data=newtheta,year=year))
}

fillgape <- function(){

    realdata=read.csv("data/theta_real.csv")
    theta=realdata$permille
	### calculate mean for smaller time windows
	#the order and signe of the timeseries dosent matter: omega(a) == omega(-rev(a))

	delta_w=sapply(1:(length(theta)-w),function(i)sd(theta[i:(i+w)]))
	plot(delta,xlab="time",main=paste0("windows size=",w*500," years (",w,"x 500yr)"),ylab=expression(delta),type="l")


	omega_w=sapply(1:(length(theta)-w),function(i)
				   {
					   y=getSpectrum(theta[i:(i+w)]) #get spectrum of the environment generated
					   x=1:length(y)
					   fit=lm(log(y)~log(x),cbind.data.frame(x=1:length(y),y=y)) #fit a linear model to check slope
					   return(abs(fit$coefficients[2]))
				   }
	)

    newt=interpolate(-realdata$permille,-realdata$years.BP.2000,finalres=-.25)
	#linear interpolation
    newy=seq(min(-realdata$years.BP.2000),max(-realdata$years.BP.2000),.25)
	#linear + noise
    newtN=interpolate(-realdata$permille,-realdata$years.BP.2000,finalres=-.25,omega=omw,delta=dtw)
    plot(-realdata$years.BP.2000[1:100],realdata$permille[1:100],type="l",lwd=10)
    plot(-realdata$years.BP.2000,-realdata$permille,type="l")
    plot(newy,newt,type="l",lwd=2,col="red")
    par(mfrow=c(1,2))
    plot(realdata$years.BP.2000,realdata$permille,type="l")
    sub=seq(1,length(newy),length.out=length(newy)/5)

    plot(-realdata$years.BP.2000,-realdata$permille,type="l",lwd=10,xlim=c(-300,0))
    lines(rev(newy),newtN,type="l",col="red",lwd=3)
    lines(rev(newy),newt,type="l",col="green")

	w=10

	omega_w=sapply(1:(length(theta)-w),function(i)getOmega(theta[i:(i+w)]) )

}

getOmega <- function(t){
	y=getSpectrum(t) #get spectrum of the environment generated
	x=1:length(y)
	fit=lm(log(y)~log(x),cbind.data.frame(x=1:length(y),y=y)) #fit a linear model to check slope
	return(abs(fit$coefficients[2])[[1]])

}

#Convert suite of delta 18 O and convert it to temperature
convertEpstein <- function(d) return(16.5-4.5*d+.14*d^2)

#get the lowest resolution of the
getDateResolution <- function(years) years[2:length(years)]-years[1:(length(years)-1)]

getLast <- function(x)return(x[1])
getFirst <- function(x)return(x[length(x)])
#getMean <- function(x)mean(x)

applySampling <- function(x,y,fun)tapply(y,x,fun)

#return resempled vector of date and data
# data a vector of measure
# year a vector of date
# by the ne sampling rate in year
getMean2 <- function(data,year,by){
    newyears=rev(seq(max(year),min(year),-by))
    newdata=sapply(1:length(newyears),function(y,data,oldyear)
                   {
                       if(y==1)
                           slice= data[ oldyear <=   newyears[y] ]
                       else
                           slice= data[ oldyear <=  newyears[y] & oldyear >  newyears[y-1]]
                       mean(slice)
                   },data=data,oldyear=year)
    newyearsmean=seq(min(newyears)+.5*by,max(newyears),by)
    return(cbind.data.frame(data=newdata[-1],year=newyearsmean))
}

getClosest <- function(data,year,by){
    newyears=rev(seq(max(year),min(year),-by))
    newdata=sapply(newyears,function(y,data,oldyear)
                   {
                       if(length(data[oldyear == y])==1)
                           return(data[oldyear == y])
                       upper=min(oldyear[oldyear > y])
                       lower=max(oldyear[oldyear < y])
                       if(abs(upper-y)<abs(lower-y))
                           dy=upper
                       else
                           dy=lower
                       data[oldyear==dy]
                   },data=data,oldyear=year)
    return(cbind.data.frame(data=newdata,year=newyears))
}

multitaper  <- function(){
    library(multitaper)
    x <- -ngrip$year
    y <- ngrip$dTsVscales
    vars <- multitaperTrend(xd=y, B=0.05, deltat=50.0, t.in=x)
    plot(x,y,type="l")
    lines(x,vars[[1]]+vars[[2]]*vars[[3]],type="l",col="red")

}
  
getSpectrumMTM <- function(data,freq,nw=2.0,k=3,...){
    require(multitaper)
res=unique(getDateResolution(freq))
if(length(res)>1)res=mean(res)
    data.ts <- ts(data=data, start=min(freq), freq=1/res) 
    data.spec <- spec.mtm(data.ts, nw=nw, k=k,plot=F) 
    return(data.spec)
}

#' Get the fit between two vector, possibility to limit the fit to a region of the data (using u and sp)
getOmega2 <- function(x,y,u=NULL,sp=.5){
    w=T
    if(!is.null(u)){
        w=x > u-sp & x < u+sp #window to comput beta (omega for us)
    }
    fitw=lm(y[w] ~ x[w])
    return(abs(fitw$coefficients[[2]]))
}

plotSpectrum <- function(data,year,...){
    logyx=getLogSpec(data,year)
    plot(logyx,axes=F,type="l",xlab="Cycle Length (years)",ylab=expression(log[10]*P(f)),...)
    axis(1,at=seq(-6,-2,by=1),label=round(1/10^seq(-6,-2,by=1)))
    axis(2)
}

#' Get the log10 of the spectrum density and frequencies of a vector
getLogSpec <- function(data,year){
    m.spec=getSpectrumMTM(data,year)
    x=log10(m.spec$freq)
    y=log10(m.spec$spec)
    return(cbind(freq=x,spec=y))
}

#' apply the function f on subsets of size w of the vector data
getWindows<-function(data,w,f)sapply(1:(length(data)-w),function(i,w)f(data[i:(i+w)]),w=w) 

#' getOmegaWrap return the slope of the spectrum of a vector
#' @param data a vector 
#' @return a double 
getOmegaWrap <- function(data){
    sp=getSpectrumMTM(data,1:length(data))
    return(getOmega2(log10(sp$freq[2:length(sp$freq)]),log10(sp$spec[2:length(sp$freq)])))
}

