
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
}

getSpectrum <- function(x){
	S.pgram <- (1/length(x))*abs(fft(x))^2
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
