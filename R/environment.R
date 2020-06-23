
#Function from tuneR package
#Based on Timmer and Koening 1995: https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract
#' @param alpha: the slope of the power distribution (called omega in whitehead)
#' @param N: the length of the timeserie to generate
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

#return a 1/f environment of N steps with sd = delta and the slope of it's spectrum decompoistion = -omega + possibility to increase the mean of the environment at a rate vt (in wich case the slop may not be -omega)
environment <- function(N,omega,delta,vt=NULL){
    ts=TK95(N,omega)
    ts=delta*ts/sd(ts)
    if(!is.null(vt))ts = ts + vt*1:N
    return(ts)
}

# white noise
gauss <- function(N,mean,delta,v){
    sapply(1:N,function(t)rnorm(1,mean*v,delta))
}



#Function interpolate
#To add intermediate value between two theta (using straight linear interpolation)
#' @param theta: the original theta vector to complete
#' @param times: the time associated to the thetas
#' @param finalres: the resolution, in the unit used in the vector times, at which new data had to be generate. if times is in year and finalres is 1, the vector thetat will correspond to data ofr each eyar, if time is in thousand of years (as in LR04 stack) and finalres is 0.5 thus the output will correspond to one point every 500 years.
#' @param omega: autocorrelation of the noise use to generate interpolate datapoints
#' @param delta: variance of the noise used for bigger gap
#' @param delta2: variance of the noise used for smaller gaps
#' @return a two columns data frame. Each column is a of dimension: seq(1,length(times),by=finalres), one with the years and other with associated theta
interpolate <- function(theta,times,finalres,delta=0,omega=0,delta2=NULL){
    if(is.null(delta2))delta2=delta
    newtheta=c()
    for(i in 1:(length(times)-1)){
        if(abs(times[i]-(times[i+1])) > finalres){ #new points needed
            rangesyears=seq(times[i],(times[i+1]),by=finalres) #generate new time points between the two original ones given our new resolution
            newt=seq(theta[i],theta[i+1],length.out=length(rangesyears)) #generate new , regularly spaced ( ie linear) data points between the original one
            names(newt)=rangesyears
            if(5*finalres < abs(times[i]-(times[i+1]))){#original resolution above threshold: 1/f noise
                n=0
                if(length(newt) %%2 != 0 ) n =1
                noise=environment(length(newt)+n,omega=omega,delta=delta) # we generate the noise needed between the two datapoints included + 1 if the length is uneven
            }
            else{ #original resolution below threshold: white noise
                noise=rnorm(length(newt),0,delta2)
            }
            newt[2:(length(newt-1))]=newt[2:(length(newt-1))]+noise[2:(length(newt-1))] #we add the generated noise only to the new datapoints, excluding the original ones
        }
        else{ #original resolution high enough, we keep original datapoints
            newt=theta[i:(i+1)]
            names(newt)=times[i:(i+1)]
        }
        newtheta=c(newtheta,newt[-length(newt)]) #exluce the last original point (will be include next step)
    }
    newtheta=c(newtheta,theta[length(theta)]) #include the last original point
    names(newtheta)[length(newtheta)]=times[length(times)]
    year=as.numeric(names(newtheta))
    names(newtheta)=NULL
    return(cbind.data.frame(data=newtheta,year=year))
}


#return the coefficient of a linear fit to the spectrum decomposition of a time serie t
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

##wrapers to functionalise simple subsetting call
getLast <- function(x)return(x[1])
getFirst <- function(x)return(x[length(x)])
applySampling <- function(x,y,fun)tapply(y,x,fun)
###

#' create an evenly spaced vector of date and associate to it the mean of the datapoints falling within the new intervals from the original dataset 
#' @param data a vector of measure
#' @param year a vector of date
#' @param by the ne sampling rate in year
#' @return a two columns data frame. Each column is of dimension: seq(min(year),max(year),by=by), one with the years and other with associated theta
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

#' create an evenly spaced vector of date and associate to it the closest datapoint from the original dataset 
#' @param data a vector of measure
#' @param year a vector of date
#' @param by the ne sampling rate in year
#' @return a two columns data frame. Each column is of dimension: seq(min(year),max(year),by=by), one with the years and other with associated theta
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

##wrapper to calculate the spectrum decomposition of a timeserie using Multitaper Methods
getSpectrumMTM <- function(data,freq,nw=2.0,k=3,...){
    require(multitaper)
    res=unique(getDateResolution(freq))
    if(length(res)>1)res=mean(res)
    data.ts <- ts(data=data, start=min(freq), freq=1/res) 
    data.spec <- spec.mtm(data.ts, nw=nw, k=k,plot=F) 
    return(data.spec)
}

#' Get the absolute slope of a linear fit between two vectors, possibility to limit the fit to a region of the data (using u and sp)
getOmega2 <- function(x,y,u=NULL,sp=.5){
    w=T
    if(!is.null(u)){
        w=x > u-sp & x < u+sp #window to comput beta (omega for us)
    }
    fitw=lm(y[w] ~ x[w])
    return(abs(fitw$coefficients[[2]]))
}

plotSpectrum <- function(data,year,add=F,...){
    logyx=getLogSpec(data,year)
    if(add)
        lines(logyx,...)
    else{ 
        plot(logyx,axes=F,type="l",xlab="Cycle Length (years)",ylab=expression(log[10]*P(f)),...)
        axis(1,at=seq(-6,-2,by=1),label=round(1/10^seq(-6,-2,by=1)))
        axis(2)
    }
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

#' getOmegaWrap return the slope of the spectrum decomposition of a vector
#' @param data a vector 
#' @return a double 
getOmegaWrap <- function(data){
    sp=getSpectrumMTM(data,1:length(data))
    return(getOmega2(log10(sp$freq[2:length(sp$freq)]),log10(sp$spec[2:length(sp$freq)])))
}

