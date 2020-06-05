source("environment.R")
setAxis <- function(data){
    mus=sapply(unique(data$mu),function(u)as.expression(bquote(mu == .(u))))
    tm=min(5,length(mus))
    nl=seq(1,length(mus),length.out=tm)
    axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

    ms=paste0("m=",unique(data$m))
    tm=min(4,length(ms))
    nl=seq(1,length(ms),length.out=tm)
    axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)

    e=sapply(unique(data$E),function(u)as.expression(bquote(E == .(u))))
    mtext(e,side=3,line=0,at=seq(.5/length(e),1-1/length(e),length.out=length(e)),cex=.9,outer=T)

    vsp=1/nlines

    ky=sapply(unique(data$k_y),function(u)as.expression(bquote(k[y] == .(u))))
    mtext(ky,side=2,line=2,at=seq(1.5*vsp,1-1.5*vsp,length.out=length(ky)),cex=.9,outer=T)

    vsp=1/nlines
    kz=sapply(unique(data$k_z),function(u)as.expression(bquote(k[z] == .(u))))
    mtext(kz,side=2,line=0,at=seq(vsp-.5*vsp,length(kz)*vsp-.5*vsp,length.out=length(kz)),cex=.9,outer=T)
}

plotVertEnv <- function(realdata,scale,...){
    plot((realdata$dTsVscales)*scale,-realdata$year,type="l",yaxs="i",axes=F,...)
    plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
    axis(2,label=rev(round(sort(seq(0,min(realdata$year),length.out=5)*scale))),at=seq(0,1,length.out=5),outer=T)
}

slpalette=colorRampPalette(c(rgb(1,.5,0),rgb(0,.5,1)))
ilpalette=colorRampPalette(c(rgb(0,.5,0),rgb(1,.5,0)))
survivalpallette=colorRampPalette(c("grey","yellow","dark green"))
extinctpalette=colorRampPalette(c("chartreuse4","white"))
pureSl=slpalette(2)[2]
pureIl=ilpalette(2)[2]


namesRealEnv=c("epica","lr04","ls16","vostok","ngrip","martrat")
namesFun=c("getMean","getLast","getFirst")
allexp=list(
#list(folder="exploreRealEnvRANDOMNoisy",
#idexpe="stackLR04InterpolLinear")
#,
#list(folder="exploreRealEnvRANDOMNoisyHigherRes",
#idexpe="stackLR04InterpolNoisy")
#,
list(folder="epicagetLast",
idexpe="epicagetLast")
)

combiname=expand.grid(namesRealEnv,namesFun)
allexp=apply(combiname,1,function(i)list(folder=paste0(i,collapse=""),idexpe=paste0(i,collapse=""),env=i[1],fun=i[2]))

for( exp in allexp){

folder=exp$folder
idexpe=exp$idexpe

ns=length(list.files(path=folder,recursive=T,full.names=T,pattern="*cross*"))

binded=do.call("rbind",lapply(list.files(path=folder,recursive=T,full.names=T,pattern="*cross*"),function(u){print(u);load(u);return(binded)}))
load(file=as.character(binded$filename[1]))
nsteps=nrow(summary)

nlines=length(unique(binded$k_z))*length(unique(binded$k_y))
ncols=length(unique(binded$E))

binded2=binded
for(s in unique(binded2$sigma)){

    binded=binded2[binded2$sigma==s,]

    png(paste0("images/",idexpe,"3Genes_sigma",s,"_E.png"),width=500,height=600,pointsize=16)

    par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
    for(ky in rev(unique(binded$k_y))){
        for(kz in rev(sort(unique(binded$k_z)))){
            for(e in unique(binded$E)){
                subb=binded[binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,]
                par(mar=rep(0,4))
                resy=tapply(subb$mean_y,subb[,c("mu","m")],mean,na.rm=T)
                resz=tapply(subb$mean_z,subb[,c("mu","m")],mean,na.rm=T)

                pxl=resz
                pxl[,]=NA
                for(i in 1:nrow(resy)){
                    for(j in 1:ncol(resy)){
                        if(is.nan(resy[i,j]))pxl[i,j]=NA
                        else
                            pxl[i,j]=rgb(red=resy[i,j],green=.5,blue=resz[i,j])
                    }
                }
                pxl=pxl[nrow(pxl):1,]
                plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
                rasterImage(pxl, 0, 0, 1, 1,interpolate=F)
                box()
            }
        }
    }
    setAxis(binded)
    mtext("mean genotype",1,2,at=0.1,outer=T,adj=0,cex=1.2)

    dev.off()
    ################ pop size


    png(paste0("images/",idexpe,"N_sigma",s,"_E.png"),width=500,height=600,pointsize=16)

    par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
    for(ky in rev(unique(binded$k_y))){
        for(kz in rev(sort(unique(binded$k_z)))){
            for(e in unique(binded$E)){
                subb=binded[binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,]
                par(mar=rep(0,4))
                res=tapply(subb$N,subb[,c("mu","m")],mean,na.rm=T)
                image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$N,na.rm=T),col=survivalpallette(4000))
            }
        }
    }
    setAxis(binded)
    par(xpd=NA)
    top=10
    points(rep(1.3,8),seq(top,top-.5,length.out=8),cex=2,pch=22,col=NA,bg=rev(survivalpallette(8)))
    text(rep(1.5,3),seq(top,top-.5,length.out=3),paste0(round(seq(max(binded$N,na.rm=T),0,length.out=3))),cex=1.2,adj=0)
    text(1.3,top+.25,expression(bar(N)[e]),cex=1.2,adj=0)
    points(1.5,top-.5-.5,cex=2,pch=22,col=1,bg="white")
    text(1.5,top-.5-.5,"NA",cex=1.2,adj=0)

    mtext("Mean population size",1,2,at=0.1,outer=T,adj=0,cex=1.2)
    dev.off()

    png(paste0("images/",idexpe,"NA_sigma",s,"_E.png"),width=500,height=600,pointsize=16)

    nexp=c()
    par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
    for(ky in rev(unique(binded$k_y))){
        for(kz in rev(sort(unique(binded$k_z)))){
            for(e in unique(binded$E)){
                subb=binded[binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,]
                par(mar=rep(0,4))
                res=tapply(subb$mean_y,subb[,c("mu","m")],function(i)sum(is.na(i)))
                nexp=min(tapply(subb$mean_y,subb[,c("mu","m")],length))
                image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(0,nexp),col=extinctpalette(nexp))
            }
        }
    }
    setAxis(binded)
    par(xpd=NA)
    points(rep(1.3,8),seq(top,top-.5,length.out=8),cex=2,pch=22,col=NA,bg=extinctpalette(8))
    text(rep(1.5,3),seq(top,top-.5,length.out=3),paste0(round(seq(0,nexp,length.out=3))),cex=1.2,adj=0)
    text(1.5,top+.25,"#ext",cex=1.2,adj=0)
    mtext("Number of Extinctions",1,2,at=0.1,outer=T,adj=0,cex=1.2)

    dev.off()


}

for(s in unique(binded2$sigma)){

    binded=binded2[binded2$sigma==s,]
    mum=as.data.frame(expand.grid(m=unique(binded$m),mu=unique(binded$mu)))
    mum$prod=mum$m*mum$mu
    lkz=length(unique(binded$k_z))
    lky=length(unique(binded$k_y))
    le=length(unique(binded$E))
    lm=length(unique(binded$m))
    lmu=length(unique(binded$mu))
    tsteps=nsteps
    for(e in unique(binded$E)){
        bicpic=matrix(nrow=tsteps,ncol=lkz*lky*lm*lmu)
        bicpic[,]=NA
        p=0
        for(ky in unique(binded$k_y)){
            for(kz in unique(binded$k_z)){
                for(i in 1:nrow(mum)){
                    mu=mum$mu[i]
                    m=mum$m[i]
                    subb=droplevels(binded[binded$mu ==mu &binded$m ==m &binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,])
                    nexpe=nrow(subb)
                    vars=c("mean_y","mean_z")
                    names(vars)=vars
                    mat_allexp= lapply(vars,function(i)matrix(NA,nexpe,tsteps))
                    for(i in 1:nexpe){
                        load(as.character(subb$filename[i]))
                        for(v  in vars){
                            if(length(mat_allexp[[v]][i,])!=length(summary[,v]))
                                print(i)
                            else
                                mat_allexp[[v]][i,]=summary[,v]

                        }
                    }
                    sum_mat=lapply(mat_allexp,function(m)apply(m,2,mean,na.rm=T))
                    sum_mat$na=apply(mat_allexp$mean_y,2,function(i)sum(is.na(i)))
                    sum_mat=lapply(sum_mat,function(m)m[!is.na(m)])
                    if(length(sum_mat$mean_y)>1){
                        if(is.na(sum_mat$na[1]))
                            pxl=NA
                        else
                            #pxl=rgb(sum_mat$mean_y,.5,sum_mat$mean_z,alpha=1)
                            pxl=rgb(sum_mat$mean_y,.5,sum_mat$mean_z,alpha=1-sum_mat$na/nexpe)
                        bicpic[1:length(pxl),p]=pxl
                    }

                    p=p+1
                    print(paste("col",p,"/",ncol(bicpic)))
                }
                scale=1 #to convert time scale in year
                png(paste0("images/vertical_",idexpe,"_sigma",s,"_E",e,".png"),width=600)
                par(mar=rep(0,4),oma=c(1,4,4,1))
                prop=.1
                layout(matrix(c(1,2),ncol=2,nrow=1),width=c(prop,1-prop))
                env=exp$env
                env=read.csv(paste0("data/",exp$env,".csv"))
                plotVertEnv(env,scale)
                rasterImage(bicpic, 0, 0, 1, 1,interpolate=F)
                mtext("year BP",3,2,at=0,outer=T)
                mtext(expression(delta[18]*O),3,0,at=prop/2,outer=T)

                kyl=sapply(unique(binded$k_y),function(u)as.expression(bquote(k[y] == .(u))))
                mtext(kyl,side=3,line=2,at=seq(0+(1/3)/2,1-(1/3)/2,length.out=length(kyl)),cex=.9)

                kzl=sapply(unique(binded$k_z),function(u)as.expression(bquote(k[z] == .(u))))
                mtext(kzl,side=3,line=0,at=seq(0+(1/9)/2,3/9-(1/9)/2,length.out=length(kzl)),cex=.9)
                mtext(bquote(E==.(e)),side=3,line=3)
                dev.off()
            }
        }
    }

}

for(s in unique(binded2$sigma)){
    subnum=20

    binded=binded2[binded2$sigma==s,]
    mum=as.data.frame(expand.grid(m=unique(binded$m),mu=unique(binded$mu)))
    mum$prod=mum$m*mum$mu
    lkz=length(unique(binded$k_z))
    lky=length(unique(binded$k_y))
    le=length(unique(binded$E))
    lm=length(unique(binded$m))
    lmu=length(unique(binded$mu))
    tsteps=nsteps
    for(e in unique(binded$E)){
        bicpic=matrix(nrow=tsteps,ncol=lkz*lky*lm*lmu*subnum+1)
        bicpic[,]=NA
        p=1
        for(ky in unique(binded$k_y)){
            for(kz in unique(binded$k_z)){
                for(i in 1:nrow(mum)){
                    mu=mum$mu[i]
                    m=mum$m[i]
                    subb=droplevels(binded[binded$mu ==mu &binded$m ==m &binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,])
                    nexpe=nrow(subb)
                    for(f in sample.int(nexpe,subnum)){
                        pxl=c()
                        load(as.character(subb$filename[f]))
                        na=which.max(is.na(summary[,1]))#find the first na ie when pop get extinct
                        if(na>1) summary=summary[1:(na-1),,drop=F]
                        if(na==1)pxl=NA
                        else
                            pxl=rgb(summary[,"mean_y" ],.5,summary[,"mean_z"],alpha=1)
                            #pxl=rgb(summary[,"mean_y" ],.5,summary[,"mean_z"],alpha=.5+1/2*min(1,summary[,"N"]/1000))
                        bicpic[1:length(pxl),p]=pxl
                        p=p+1
                        print(paste("col",p,"/",ncol(bicpic)))
                    }
                }
                scale=1000 #to convert time scale in year
                png(paste0("images/verticalM_",idexpe,"_sigma",s,"_E",e,".png"),width=1800,height=1400,pointsize=42)
                par(mar=rep(0,4),oma=c(1,4,4,1))
                prop=.1
                layout(matrix(c(1,2),ncol=2,nrow=1),width=c(prop,1-prop))
                env=read.csv(paste0("data/",exp$env,".csv"))
                plotVertEnv(env,scale,lwd=4)
                rasterImage(bicpic, 0, 0, 1, 1,interpolate=F)
                mtext("year BP",3,2,at=0,outer=T)
                mtext(expression(delta[18]*O),3,0,at=prop/2,outer=T)

                kyl=sapply(unique(binded$k_y),function(u)as.expression(bquote(k[y] == .(u))))
                mtext(kyl,side=3,line=2,at=seq(0+(1/3)/2,1-(1/3)/2,length.out=length(kyl)),cex=.9)

                kzl=sapply(unique(binded$k_z),function(u)as.expression(bquote(k[z] == .(u))))
                mtext(kzl,side=3,line=0,at=seq(0+(1/9)/2,3/9-(1/9)/2,length.out=length(kzl)),cex=.9)
                mtext(bquote(E==.(e)),side=3,line=3)
                dev.off()
            }
        }
    }

}
}
