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


for(s in unique(binded$sigma)){
    for(kz in unique(binded$k_z)){
        for(ky in unique(binded$k_y)){
            #ky=1
            #kz=1
            #par(mfrow=c(4,5),mar=rep(1,4),oma=c(4,3,0,0))
            #par(mfrow=c(2,1))
            for(mu in unique(binded$mu)){
                for(m in unique(binded$m)){
                    subb=droplevels(binded[binded$mu ==mu &binded$m ==m &binded$sigma ==s & binded$k_z ==kz & binded$k_y ==ky,])
                    nexpe=nrow(subb)
                    vars=c("prop_y","prop_z","mean_x","N","prop_yz")
                    names(vars)=vars
                    mat_allexp= lapply(vars,function(i)matrix(NA,nexpe,length(theta)))

                    for(i in 1:50){
                        load(as.character(subb$filename[i]))
                        for(v  in vars) mat_allexp[[v]][i,]=summary[,v]
                    }
                    sum_mat=lapply(mat_allexp,function(m)apply(m,2,quantile,probs=c(.05,.5,.95),na.rm=T))
                    fname=paste0("traj_s",s,"_m",m,"_mu",mu,"_ky",ky,"_kz",kz,".png")
                    #png(fname,width=450,height=450,pointsize=28)
                    plotMatrixStrateAndEn(sum_mat,theta)
                    #dev.off()
                }
            }
        }
    }
}




################ les deux genes

nlines=length(unique(binded$delta))*length(unique(binded$omega))
ncols=length(unique(binded$vt))*length(unique(binded$k_y))


png(paste0("exploringRealEnv3Genes.png"),width=450,height=500,pointsize=16)
nlines=3
ncols=3
par(mfrow=c(nlines,ncols),oma=c(6,5,2,5))
for(s in unique(binded$sigma)){
    for(kz in rev(sort(unique(binded$k_z)))){
        for(ky in unique(binded$k_y)){
            subb=binded[binded$sigma ==s & binded$k_z ==kz & binded$k_y ==ky,]
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
            #print(res)
            #kimage(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$mean_z,na.rm=T),col=mypalette(4000))
            #text(allval[,1],allval[,2],allval[,3],cex=.3) 
        }
    }
}
mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
tm=5
nl=seq(1,length(mus),length.out=tm)
axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

ms=paste0("m=",unique(binded$m))
tm=4
nl=seq(1,length(ms),length.out=tm)
axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
ky=sapply(unique(binded$k_y),function(u)as.expression(bquote(k[y] == .(u))))
mtext(ky,side=3,line=0,at=seq(.5/length(ky),2.5/length(ky),length.out=length(ky)),cex=.9,outer=T)
kz=sapply(unique(binded$k_z),function(u)as.expression(bquote(k[z] == .(u))))
mtext(kz,side=2,line=2,at=seq(.5/length(kz),2.5/length(kz),length.out=length(kz)),cex=.9,outer=T)

dev.off()
################ pop size


png(paste0("exploringRealEnvN.png"),width=450,height=500,pointsize=16)

par(mfrow=c(nlines,ncols),oma=c(6,5,2,5))
for(s in unique(binded$sigma)){
    for(kz in rev(sort(unique(binded$k_z)))){
        for(sg in unique(binded$k_y)){
            subb=binded[binded$sigma ==s & binded$k_z ==kz & binded$k_y ==sg,]
            par(mar=rep(0,4))
            res=tapply(subb$N,subb[,c("mu","m")],mean,na.rm=T)
            image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$N,na.rm=T),col=survivalpallette(4000))
        }
    }
}
mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
tm=5
nl=seq(1,length(mus),length.out=tm)
axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)
ms=paste0("m=",unique(binded$m))
tm=4
nl=seq(1,length(ms),length.out=tm)
axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
ky=sapply(unique(binded$k_y),function(u)as.expression(bquote(k[y] == .(u))))
mtext(ky,side=3,line=0,at=seq(.5/length(ky),2.5/length(ky),length.out=length(ky)),cex=.9,outer=T)
kz=sapply(unique(binded$k_z),function(u)as.expression(bquote(k[z] == .(u))))
mtext(kz,side=2,line=2,at=seq(.5/length(kz),2.5/length(kz),length.out=length(kz)),cex=.9,outer=T)

par(xpd=NA)
points(rep(1.5,8),seq(3,2.8,length.out=8),cex=2,pch=22,col=NA,bg=rev(survivalpallette(8)))
text(rep(1.8,3),seq(3,2.8,length.out=3),paste0(round(seq(max(binded$N,na.rm=T),0,length.out=3),digits=2)),cex=1.2,adj=0)
text(1.8,3.11,expression(bar(N)[e]),cex=1.2,adj=0)
points(2,2.6,cex=2,pch=22,col=1,bg="white")
text(1.3,2.6,"NA",cex=1.2,adj=0)

dev.off()

png(paste0("exploringRealEnv.png"),width=450,height=500,pointsize=16)

par(mfrow=c(nlines,ncols),oma=c(6,5,2,5))
nexp=c()
for(s in unique(binded$sigma)){
    for(kz in rev(sort(unique(binded$k_z)))){
        for(ky in unique(binded$k_y)){
            subb=binded[binded$sigma ==s & binded$k_z ==kz & binded$k_y ==ky,]
            par(mar=rep(0,4))
            res=tapply(subb$mean_y,subb[,c("mu","m")],function(i)sum(is.na(i)))
            nexp=min(tapply(subb$mean_y,subb[,c("mu","m")],length))
            image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(0,nexp),col=extinctpalette(nexp))
        }
    }
}
mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
tm=5
nl=seq(1,length(mus),length.out=tm)
axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)
ms=paste0("m=",unique(binded$m))
tm=4
nl=seq(1,length(ms),length.out=tm)
axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
ky=sapply(unique(binded$k_y),function(u)as.expression(bquote(k[y] == .(u))))
mtext(ky,side=3,line=0,at=seq(.5/length(ky),2.5/length(ky),length.out=length(ky)),cex=.9,outer=T)
kz=sapply(unique(binded$k_z),function(u)as.expression(bquote(k[z] == .(u))))
mtext(kz,side=2,line=2,at=seq(.5/length(kz),2.5/length(kz),length.out=length(kz)),cex=.9,outer=T)

par(xpd=NA)
points(rep(1.5,8),seq(3,2.8,length.out=8),cex=2,pch=22,col=NA,bg=rev(extinctpalette(nexp)))
text(rep(1.8,3),seq(3,2.8,length.out=3),paste0(round(seq(0,nexp,length.out=3),digits=2)),cex=1.2,adj=0)
text(1.8,3.11,expression(bar(N)[e]),cex=1.2,adj=0)

dev.off()
