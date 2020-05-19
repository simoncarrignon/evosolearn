slpalette=colorRampPalette(c(rgb(0,.5,0),rgb(0,.5,1)))
ilpalette=colorRampPalette(c(rgb(0,.5,0),rgb(1,.5,0)))
survivalpallette=colorRampPalette(c("grey","yellow","dark green"))
extinctpalette=colorRampPalette(c("chartreuse4","white"))
pureSl=slpalette(2)[2]
pureIl=ilpalette(2)[2]

realdata=read.csv("data/41586_2004_BFnature02805_MOESM1_ESM.csv")
theta=rev(tapply(realdata$permille,realdata$years.BP.2000,mean))
plot(theta)
binded=do.call("rbind",lapply(list.files(path="exploreRealEnvBESTsigMas",recursive=T,full.names=T,pattern="*cross*"),function(u){print(u);load(u);return(binded)}))
binded=binded[binded$sigma==2,]
nonoise=binded
binded=do.call("rbind",lapply(list.files(path="exploreRealEnvBESTNoisy",recursive=T,full.names=T,pattern="*cross*"),function(u){print(u);load(u);return(binded)}))
binded=rbind(nonoise,binded)

nlines=3*3
ncols=length(unique(binded$E))

setAxis <- function(){
    mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
    tm=5
    nl=seq(1,length(mus),length.out=tm)
    axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

    ms=paste0("m=",unique(binded$m))
    tm=4
    nl=seq(1,length(ms),length.out=tm)
    axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)

    e=sapply(unique(binded$E),function(u)as.expression(bquote(E == .(u))))
    mtext(e,side=3,line=0,at=seq(.5/length(e),3.5/length(e),length.out=length(e)),cex=.9,outer=T)

    vsp=1/nlines

    ky=sapply(unique(binded$k_y),function(u)as.expression(bquote(k[y] == .(u))))
    mtext(ky,side=2,line=2,at=seq(1.5*vsp,1-1.5*vsp,length.out=length(ky)),cex=.9,outer=T)

    vsp=1/nlines
    kz=sapply(unique(binded$k_z),function(u)as.expression(bquote(k[z] == .(u))))
    mtext(kz,side=2,line=0,at=seq(vsp-.5*vsp,length(kz)*vsp-.5*vsp,length.out=length(kz)),cex=.9,outer=T)
}


png(paste0("exploringRealEnv3Genes_E.png"),width=500,height=600,pointsize=16)

par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
for(ky in rev(unique(binded$k_y))){
    for(kz in rev(sort(unique(binded$k_z)))){
        for(s in unique(binded$E)){
            subb=binded[binded$E ==s & binded$k_z ==kz & binded$k_y ==ky,]
            par(mar=rep(0,4))
            resy=tapply(subb$mean_y,subb[,c("mu","m")],mean,na.rm=T)
            resz=tapply(subb$mean_z,subb[,c("mu","m")],mean,na.rm=T)

            pxl=resz
            pxl[,]=NA
            for(i in 1:nrow(resy)){
                for(j in 1:ncol(resy)){
                    print(paste(i,j))
                    if(is.nan(resy[i,j]))pxl[i,j]=NA
                    else
                        pxl[i,j]=rgb(red=resy[i,j],green=.5,blue=resz[i,j])
                }
            }
            pxl=pxl[nrow(pxl):1,]
            plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
            rasterImage(pxl, 0, 0, 1, 1,interpolate=F)
        }
    }
}
setAxis()
mtext("mean genotype",1,2,at=0.1,outer=T,adj=0,cex=1.2)

dev.off()
################ pop size


png(paste0("exploringRealEnvN_E.png"),width=500,height=600,pointsize=16)

par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
for(ky in rev(unique(binded$k_y))){
    for(kz in rev(sort(unique(binded$k_z)))){
        for(s in unique(binded$E)){
            subb=binded[binded$E ==s & binded$k_z ==kz & binded$k_y ==sg,]
            par(mar=rep(0,4))
            res=tapply(subb$N,subb[,c("mu","m")],mean,na.rm=T)
            image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$N,na.rm=T),col=survivalpallette(4000))
        }
    }
}
setAxis()
par(xpd=NA)
top=10
points(rep(1.3,8),seq(top,top-.5,length.out=8),cex=2,pch=22,col=NA,bg=rev(survivalpallette(8)))
text(rep(1.5,3),seq(top,top-.5,length.out=3),paste0(round(seq(max(binded$N,na.rm=T),0,length.out=3))),cex=1.2,adj=0)
text(1.3,top+.25,expression(bar(N)[e]),cex=1.2,adj=0)
points(1.5,top-.5-.5,cex=2,pch=22,col=1,bg="white")
text(1.5,top-.5-.5,"NA",cex=1.2,adj=0)

mtext("Mean population size",1,2,at=0.1,outer=T,adj=0,cex=1.2)
dev.off()

png(paste0("exploringRealEnv_E.png"),width=500,height=600,pointsize=16)

nexp=c()
par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
for(ky in rev(unique(binded$k_y))){
    for(kz in rev(sort(unique(binded$k_z)))){
        for(s in unique(binded$E)){
            subb=binded[binded$E ==s & binded$k_z ==kz & binded$k_y ==ky,]
            par(mar=rep(0,4))
            res=tapply(subb$mean_y,subb[,c("mu","m")],function(i)sum(is.na(i)))
            nexp=min(tapply(subb$mean_y,subb[,c("mu","m")],length))
            image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(0,nexp),col=extinctpalette(nexp))
        }
    }
}
setAxis()
par(xpd=NA)
points(rep(1.3,8),seq(top,top-.5,length.out=8),cex=2,pch=22,col=NA,bg=extinctpalette(8))
text(rep(1.5,3),seq(top,top-.5,length.out=3),paste0(round(seq(0,nexp,length.out=3))),cex=1.2,adj=0)
text(1.5,top+.25,"#ext",cex=1.2,adj=0)
mtext("Number of Extinctions",1,2,at=0.1,outer=T,adj=0,cex=1.2)

dev.off()

