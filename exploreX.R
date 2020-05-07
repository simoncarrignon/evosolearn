survivalpallette=colorRampPalette(c("grey","yellow","dark green"))
extinctpalette=colorRampPalette(c("chartreuse","white"))

N=c(2000,1000)
sls=c("BEST","RANDOM")
for(n in N){
binded=do.call("rbind",lapply(list.files(path=paste0("exploreX",n,"T1",sls),recursive=T,full.names=T),function(u){load(u);return(binded)}))
#binded[is.na(binded)]=0

png(paste0("exploringXadaptation_K",n,"_t3000_big.png"),width=1200,height=1200,pointsize=28)


nlines=length(unique(binded$delta))*length(unique(binded$omega))
ncols=length(unique(binded$vt))*length(unique(binded$sigma))

par(mfrow=c(nlines,ncols),oma=c(5,5,2,5))

for(d in rev(unique(binded$delta))){
    for(o in rev(unique(binded$omega))){
        for(v in unique(binded$vt)){
            for(sg in unique(binded$sigma)){
                subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$sigma ==sg,]
                #print(subb)
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
mtext(delt,side=2,,line=2,at=seq(0.1,.9,length.out=length(delt)),cex=.9,outer=T)

omeg=sapply(unique(binded$omega),function(d)as.expression(bquote(.(d))))
mtext(omeg,side=2,line=0,at=seq(0.02,.18,length.out=length(omeg)),cex=.9,outer=T)
mtext(expression(omega),side=2,line=1,at=0.02,cex=1,outer=T)

sig=sapply(unique(binded$sigma),function(d)as.expression(bquote(.(d))))
mtext(sig,side=1,line=1.5,at=seq(0.02,.18,length.out=length(sig)),cex=1,outer=T)
mtext(expression(sigma*": "),side=1,line=1.5,at=0.0,cex=1,outer=T)

vs=paste0("v=",unique(binded$vt))
mtext(vs,side=1,line=3,at=seq(.1,.9,length.out=length(vs)),cex=.85,outer=1)
par(xpd=NA)
points(rep(2,8),seq(30,28,length.out=8),cex=2,pch=22,col=NA,bg=rev(survivalpallette(8)))
text(rep(2.3,3),seq(30,28,length.out=3),paste0(round(seq(max(binded$N,na.rm=T),0,length.out=3),digits=-2)),cex=1.2,adj=0)
text(2.3,31,"N",cex=1.2,adj=0)
points(2,26,cex=2,pch=22,col=1,bg="white")
text(2.3,26,"NA",cex=1.2,adj=0)

dev.off()

png(paste0("exploringXadaptation_K",n,"_t3000_NA.png"),width=1200,height=1200,pointsize=28)

nlines=length(unique(binded$delta))*length(unique(binded$omega))
ncols=length(unique(binded$vt))*length(unique(binded$sigma))
par(mfrow=c(nlines,ncols),oma=c(5,5,2,5))
for(d in rev(unique(binded$delta))){
    for(o in rev(unique(binded$omega))){
        for(v in unique(binded$vt)){
            for(sg in unique(binded$sigma)){
                subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$sigma ==sg,]
                #print(subb)
                par(mar=rep(0,4))
                res=tapply(subb$mean_x,subb[,c("mu","m")],function(i)sum(is.na(i)))
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

sig=sapply(unique(binded$sigma),function(d)as.expression(bquote(.(d))))
mtext(sig,side=1,line=1.5,at=seq(.5/ncols,2.5/ncols,length.out=length(sig)),cex=1,outer=T)
mtext(expression(sigma*": "),side=1,line=1.5,at=-.5/ncols,cex=1,outer=T)

vs=paste0("v=",unique(binded$vt))
mtext(vs,side=1,line=3,at=seq(.1,.9,length.out=length(vs)),cex=.85,outer=1)
par(xpd=NA)
points(rep(2,8),seq(30,28,length.out=8),cex=2,pch=22,col=NA,bg=rev(extinctpalette(8)))
text(rep(2.3,3),seq(30,28,length.out=3),paste0(round(seq(20,0,length.out=3),digits=2)),cex=1.2,adj=0)
text(2.3,31,"#extinct",cex=1.2,adj=0)
points(2,26,cex=2,pch=22,col=1,bg="white")
text(2.3,26,"NA",cex=1.2,adj=0)

dev.off()
}

