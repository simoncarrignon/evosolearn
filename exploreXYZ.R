#slpalette=colorRampPalette(c("dark green","yellow","dark blue"))
#ilpalette=colorRampPalette(c("dark green","yellow", "red"))
slpalette=colorRampPalette(c(rgb(0,.5,0),rgb(0,.4,1)))
ilpalette=colorRampPalette(c(rgb(0,.5,0),rgb(1,.4,0)))

bindedAll=do.call("rbind",lapply(list.files(path="exploreX1000T1SEX",recursive=T,full.names=T),function(u){load(u);return(binded)}))




sls=c("BEST","RANDOM")

for(s in sls){
binded=do.call("rbind",lapply(list.files(path=paste0("exploreXYZ1000T1",s),recursive=T,full.names=T),function(u){load(u);return(binded)}))


    ################ les deux genes
    for(kz in unique(binded$k_z)){

        png(paste0("exploringXYadaptation_K1000_t3000_",s,"_yz_k",kz,".png"),width=850,height=1200,pointsize=28)

        nlines=length(unique(binded$delta))*length(unique(binded$omega))
        ncols=length(unique(binded$vt))*length(unique(binded$k_y))

        par(mfrow=c(nlines,ncols),oma=c(5,5,2,5))
        for(d in rev(unique(binded$delta))){
            for(o in rev(unique(binded$omega))){
                for(v in unique(binded$vt)){
                    for(sg in unique(binded$k_y)){
                        subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$k_z ==kz & binded$k_y ==sg,]
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
        }
        mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
        tm=3
        nl=seq(1,length(mus),length.out=tm)
        axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

        ms=paste0("m=",unique(binded$m))
        tm=3
        nl=seq(1,length(ms),length.out=tm)
        axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
        #mtext(unique(binded$delta),side=2,line=4*10+1,at=(0:3)+1.5)

        delt=sapply(unique(binded$delta),function(u)as.expression(bquote(delta == .(u))))
        mtext(delt,side=2,,line=2,at=seq(.5/length(delt),4.5/length(delt),length.out=length(delt)),cex=.9,outer=T)

        omeg=sapply(unique(binded$omega),function(d)as.expression(bquote(.(d))))
        mtext(omeg,side=2,line=0,at=seq(.5/nlines,4.5/nlines,length.out=length(omeg)),cex=.9,outer=T)
        mtext(expression(omega),side=2,line=1,at=0.02,cex=1,outer=T)

        sig=sapply(unique(binded$k_y),function(d)as.expression(bquote(.(d))))
        mtext(sig,side=1,line=1.5,at=seq(.5/ncols,2.5/ncols,length.out=length(sig)),cex=1,outer=T)
        mtext(expression(k[y]*": "),side=1,line=1.5,at=-.5/ncols,cex=1,outer=T)

        vs=paste0("v=",unique(binded$vt))
        mtext(vs,side=1,line=3,at=seq(.1,.9,length.out=length(vs)),cex=.85,outer=1)
        par(xpd=NA)
        points(rep(2,8),seq(30,28,length.out=8),cex=2,pch=22,col=NA,bg=rev(mypalette(8)))
        text(rep(2.3,3),seq(30,28,length.out=3),paste0(round(seq(max(binded$mean_y,na.rm=T),0,length.out=3),digits=2)),cex=1.2,adj=0)
        text(2.3,31,expression(bar(y)),cex=1.2,adj=0)
        points(2,26,cex=2,pch=22,col=1,bg="white")
        text(2.3,26,"NA",cex=1.2,adj=0)
        dev.off()
    }

    ################ pop size

    for(kz in unique(binded$k_z)){
        png(paste0("exploringXYadaptation_K1000_t3000_",s,"_N_k",kz,".png"),width=850,height=1200,pointsize=28)

        nlines=length(unique(binded$delta))*length(unique(binded$omega))
        ncols=length(unique(binded$vt))*(length(unique(binded$k_y))+1)
        par(mfrow=c(nlines,ncols),oma=c(5,5,2,5))
        for(d in rev(unique(binded$delta))){
            for(o in rev(unique(binded$omega))){
                for(v in unique(binded$vt)){
                    subb=bindedAll[bindedAll$omega == o & bindedAll$delta == d & bindedAll$vt == v & bindedAll$sigma ==2,]
                    par(mar=rep(0,4))
                    res=tapply(subb$N,subb[,c("mu","m")],mean,na.rm=T)
                    #print(res)
                    image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(bindedAll$N,na.rm=T),col=survivalpallette(4000))
                    box(col="violet",lwd=6)
                    for(sg in unique(binded$k_y)){
                        subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$k_y ==sg & binded$k_z ==kz ,]
                        par(mar=rep(0,4))
                        res=tapply(subb$N,subb[,c("mu","m")],mean,na.rm=T)
                        #print(res)
                        image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$N,na.rm=T),col=survivalpallette(4000))
                        #text(allval[,1],allval[,2],allval[,3],cex=.3) 
                    }
                }
            }
        }
        mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
        tm=3
        nl=seq(1,length(mus),length.out=tm)
        axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

        ms=paste0("m=",unique(binded$m))
        tm=3
        nl=seq(1,length(ms),length.out=tm)
        axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
        #mtext(unique(binded$delta),side=2,line=4*10+1,at=(0:3)+1.5)

        delt=sapply(unique(binded$delta),function(u)as.expression(bquote(delta == .(u))))
        mtext(delt,side=2,,line=2,at=seq(.5/length(delt),4.5/length(delt),length.out=length(delt)),cex=.9,outer=T)

        omeg=sapply(unique(binded$omega),function(d)as.expression(bquote(.(d))))
        mtext(omeg,side=2,line=0,at=seq(.5/nlines,4.5/nlines,length.out=length(omeg)),cex=.9,outer=T)
        mtext(expression(omega),side=2,line=1,at=0.02,cex=1,outer=T)

        sig=sapply(unique(binded$k_y),function(d)as.expression(bquote(.(d))))
        mtext(sig,side=1,line=1.5,at=1/ncols+seq(.5/ncols,2.5/ncols,length.out=length(sig)),cex=1,outer=T)
        mtext(expression(k[y]*": "),side=1,line=1.5,at=0.5/ncols,cex=1,outer=T)
        mtext("no Y",side=3,line=.5,at=seq(0.5/ncols,1,4/ncols),cex=.7,outer=T)

        vs=paste0("v=",unique(binded$vt))
        mtext(vs,side=1,line=3,at=seq(.1,.9,length.out=length(vs)),cex=.85,outer=1)

        par(xpd=NA)
        points(rep(2,8),seq(30,28,length.out=8),cex=2,pch=22,col=NA,bg=rev(survivalpallette(8)))
        text(rep(2.3,3),seq(30,28,length.out=3),paste0(round(seq(max(binded$N,na.rm=T),0,length.out=3),digits=2)),cex=1.2,adj=0)
        text(2.3,31,expression(bar(N)[e]),cex=1.2,adj=0)
        points(2,26,cex=2,pch=22,col=1,bg="white")
        text(2.3,26,"NA",cex=1.2,adj=0)
        points(2,25,cex=2,pch=22,col="violet",bg="white",lwd=3)
        text(2.3,25,"no Y",cex=1.2,adj=0)

        dev.off()
    }

    ################ num extinctions
    for(kz in unique(binded$k_z)){
        png(paste0("exploringXYadaptation_K1000_t3000_NA_",s,"_k",kz,".png"),width=850,height=1200,pointsize=28)

        nlines=length(unique(binded$delta))*length(unique(binded$omega))
        ncols=length(unique(binded$vt))*length(unique(binded$k_y))
        par(mfrow=c(nlines,ncols),oma=c(5,5,2,5))
        for(d in rev(unique(binded$delta))){
            for(o in rev(unique(binded$omega))){
                for(v in unique(binded$vt)){
                    for(sg in unique(binded$k_y)){
                        subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$k_y ==sg & binded$k_z ==kz,]
                        par(mar=rep(0,4))
                        res=tapply(subb$mean_y,subb[,c("mu","m")],function(i)sum(is.na(i)))
                        #print(res)
                        image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(0,20),col=extinctpalette(20))
                        #text(allval[,1],allval[,2],allval[,3],cex=.3) 
                    }
                }
            }
        }
        mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
        tm=3
        nl=seq(1,length(mus),length.out=tm)
        axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

        ms=paste0("m=",unique(binded$m))
        tm=3
        nl=seq(1,length(ms),length.out=tm)
        axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
        #mtext(unique(binded$delta),side=2,line=4*10+1,at=(0:3)+1.5)

        delt=sapply(unique(binded$delta),function(u)as.expression(bquote(delta == .(u))))
        mtext(delt,side=2,,line=2,at=seq(.5/length(delt),4.5/length(delt),length.out=length(delt)),cex=.9,outer=T)

        omeg=sapply(unique(binded$omega),function(d)as.expression(bquote(.(d))))
        mtext(omeg,side=2,line=0,at=seq(.5/nlines,4.5/nlines,length.out=length(omeg)),cex=.9,outer=T)
        mtext(expression(omega),side=2,line=1,at=0.02,cex=1,outer=T)

        sig=sapply(unique(binded$k_y),function(d)as.expression(bquote(.(d))))
        mtext(sig,side=1,line=1.5,at=seq(.5/ncols,2.5/ncols,length.out=length(sig)),cex=1,outer=T)
        mtext(expression(k[y]*": "),side=1,line=1.5,at=-.5/ncols,cex=1,outer=T)

        vs=paste0("v=",unique(binded$vt))
        mtext(vs,side=1,line=3,at=seq(.1,.9,length.out=length(vs)),cex=.85,outer=1)
        par(xpd=NA)
        points(rep(2,8),seq(30,28,length.out=8),cex=2,pch=22,col=NA,bg=rev(extinctpalette(8)))
        text(rep(2.3,3),seq(30,28,length.out=3),paste0(round(seq(20,0,length.out=3),digits=2)),cex=1.2,adj=0)
        text(2.3,31,"$extinct",cex=1.2,adj=0)
        #points(2,26,cex=2,pch=22,col=1,bg="white")
        #text(2.3,26,"NA",cex=1.2,adj=0)
        dev.off()
    }


    ################ y 
    for(kz in unique(binded$k_z)){
        png(paste0("exploringXYadaptation_K1000_t3000_",s,"_y_k",kz,".png"),width=850,height=1200,pointsize=28)

        nlines=length(unique(binded$delta))*length(unique(binded$omega))
        ncols=length(unique(binded$vt))*(length(unique(binded$k_y)))
        par(mfrow=c(nlines,ncols),oma=c(5,5,2,5))
        for(d in rev(unique(binded$delta))){
            for(o in rev(unique(binded$omega))){
                for(v in unique(binded$vt)){
                    for(sg in unique(binded$k_y)){
                        subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$k_y ==sg & binded$k_z ==kz ,]
                        par(mar=rep(0,4))
                        res=tapply(subb$mean_y,subb[,c("mu","m")],mean,na.rm=T)
                        image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$mean_y,na.rm=T),col=ilpalette(4000))
                        #text(allval[,1],allval[,2],allval[,3],cex=.3) 
                    }
                }
            }
        }
        mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
        tm=3
        nl=seq(1,length(mus),length.out=tm)
        axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

        ms=paste0("m=",unique(binded$m))
        tm=3
        nl=seq(1,length(ms),length.out=tm)
        axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
        #mtext(unique(binded$delta),side=2,line=4*10+1,at=(0:3)+1.5)

        delt=sapply(unique(binded$delta),function(u)as.expression(bquote(delta == .(u))))
        mtext(delt,side=2,,line=2,at=seq(.5/length(delt),4.5/length(delt),length.out=length(delt)),cex=.9,outer=T)

        omeg=sapply(unique(binded$omega),function(d)as.expression(bquote(.(d))))
        mtext(omeg,side=2,line=0,at=seq(.5/nlines,4.5/nlines,length.out=length(omeg)),cex=.9,outer=T)
        mtext(expression(omega),side=2,line=1,at=0.02,cex=1,outer=T)

        sig=sapply(unique(binded$k_y),function(d)as.expression(bquote(.(d))))
        mtext(sig,side=1,line=1.5,at=1/ncols+seq(.5/ncols,2.5/ncols,length.out=length(sig)),cex=1,outer=T)
        mtext(expression(k[y]*": "),side=1,line=1.5,at=0.5/ncols,cex=1,outer=T)
        mtext("no Y",side=3,line=.5,at=seq(0.5/ncols,1,4/ncols),cex=.7,outer=T)

        vs=paste0("v=",unique(binded$vt))
        mtext(vs,side=1,line=3,at=seq(.1,.9,length.out=length(vs)),cex=.85,outer=1)

        par(xpd=NA)
        points(rep(2,8),seq(30,28,length.out=8),cex=2,pch=22,col=NA,bg=rev(ilpalette(8)))
        text(rep(2.3,3),seq(30,28,length.out=3),paste0(round(seq(max(binded$mean_y,na.rm=T),0,length.out=3),digits=2)),cex=1.2,adj=0)
        text(2.3,31,expression(bar(y)),cex=1.2,adj=0)
        points(2,26,cex=2,pch=22,col=1,bg="white")
        text(2.3,26,"NA",cex=1.2,adj=0)
        #points(2,25,cex=2,pch=22,col="violet",bg="white",lwd=3)
        #text(2.3,25,"no Y",cex=1.2,adj=0)

        dev.off()
    }

    ################ z 
    for(kz in unique(binded$k_z)){
        png(paste0("exploringXYadaptation_K1000_t3000_",s,"_z_k",kz,".png"),width=850,height=1200,pointsize=28)

        nlines=length(unique(binded$delta))*length(unique(binded$omega))
        ncols=length(unique(binded$vt))*(length(unique(binded$k_y)))
        par(mfrow=c(nlines,ncols),oma=c(5,5,2,5))
        for(d in rev(unique(binded$delta))){
            for(o in rev(unique(binded$omega))){
                for(v in unique(binded$vt)){
                    for(sg in unique(binded$k_y)){
                        subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$k_y ==sg & binded$k_z ==kz ,]
                        par(mar=rep(0,4))
                        res=tapply(subb$mean_z,subb[,c("mu","m")],mean,na.rm=T)
                        image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$mean_z,na.rm=T),col=slpalette(4000))
                        #text(allval[,1],allval[,2],allval[,3],cex=.3) 
                    }
                }
            }
        }
        mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
        tm=3
        nl=seq(1,length(mus),length.out=tm)
        axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

        ms=paste0("m=",unique(binded$m))
        tm=3
        nl=seq(1,length(ms),length.out=tm)
        axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)
        #mtext(unique(binded$delta),side=2,line=4*10+1,at=(0:3)+1.5)

        delt=sapply(unique(binded$delta),function(u)as.expression(bquote(delta == .(u))))
        mtext(delt,side=2,,line=2,at=seq(.5/length(delt),4.5/length(delt),length.out=length(delt)),cex=.9,outer=T)

        omeg=sapply(unique(binded$omega),function(d)as.expression(bquote(.(d))))
        mtext(omeg,side=2,line=0,at=seq(.5/nlines,4.5/nlines,length.out=length(omeg)),cex=.9,outer=T)
        mtext(expression(omega),side=2,line=1,at=0.02,cex=1,outer=T)

        sig=sapply(unique(binded$k_y),function(d)as.expression(bquote(.(d))))
        mtext(sig,side=1,line=1.5,at=1/ncols+seq(.5/ncols,2.5/ncols,length.out=length(sig)),cex=1,outer=T)
        mtext(expression(k[y]*": "),side=1,line=1.5,at=0.5/ncols,cex=1,outer=T)
        mtext("no Y",side=3,line=.5,at=seq(0.5/ncols,1,4/ncols),cex=.7,outer=T)

        vs=paste0("v=",unique(binded$vt))
        mtext(vs,side=1,line=3,at=seq(.1,.9,length.out=length(vs)),cex=.85,outer=1)

        par(xpd=NA)
        points(rep(2,8),seq(30,28,length.out=8),cex=2,pch=22,col=NA,bg=rev(slpalette(8)))
        text(rep(2.3,3),seq(30,28,length.out=3),paste0(round(seq(max(binded$mean_z,na.rm=T),0,length.out=3),digits=2)),cex=1.2,adj=0)
        text(2.3,31,expression(bar(z)),cex=1.2,adj=0)
        points(2,26,cex=2,pch=22,col=1,bg="white")
        text(2.3,26,"NA",cex=1.2,adj=0)
        #points(2,25,cex=2,pch=22,col="violet",bg="white",lwd=3)
        #text(2.3,25,"no Y",cex=1.2,adj=0)

        dev.off()
    }
}
