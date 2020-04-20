
binded=do.call("rbind",lapply(list.files(path="exploreX1000T4",recursive=T,full.names=T),function(u){load(u);return(binded)}))
binded[is.na(binded)]=0

allval=cbind(expand.grid(seq(0,1,length=4),seq(0,1,length=4)), unlist(lapply(unique(binded$m),function(i)i*unique(binded$mu))))

png("exploringXadaptation_K1000_t12000.png",width=1200,height=1200,pointsize=24)
par(mfrow=c(16,16),oma=c(5,4,2,5))
for(d in rev(unique(binded$delta))){
for(o in rev(unique(binded$omega))){
for(v in unique(binded$vt)){
for(sg in unique(binded$sigma)){
    subb=binded[binded$omega == o & binded$delta == d & binded$vt == v & binded$sigma ==sg,]
    par(mar=rep(0,4))
res=tapply(subb$N,subb[,c("mu","m")],mean)
#print(res)
image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$N),col=colorRampPalette(c("red","yellow","dark green"))(4000))
#text(allval[,1],allval[,2],allval[,3],cex=.3) 
}
}
}
}
mus=sapply(unique(binded$mu),function(u)as.expression(bquote(mu == .(u))))
axis(1,label=mus,at=seq(0,1,length.out=length(mus)),cex=.6,las=3)
ms=paste0("m=",unique(binded$m))
axis(4,label=ms,at=seq(0,1,length.out=length(ms)),las=1,cex=.6)
#mtext(unique(binded$delta),side=2,line=4*10+1,at=(0:3)+1.5)
omeg=sapply(unique(binded$omega),function(u)as.expression(bquote(omega == .(u))))
mtext(omeg,side=2,line=5*10+2,at=seq(0.5,4.5,length.out=length(omeg)),cex=.7)
delt=sapply(unique(binded$delta),function(d)as.expression(bquote(delta == .(d))))
mtext(delt,side=2,line=5*10+2+1,at=seq(2.5,18.5,length.out=length(delt)),cex=.85)
sig=sapply(unique(binded$sigma),function(d)as.expression(bquote(sigma == .(d))))
mtext(sig,side=1,line=1,at=seq(-19.5,-15.5,length.out=4),cex=.7)
vs=paste0("v=",unique(binded$vt))
mtext(vs,side=1,line=2,at=seq(-17,-1.5,length.out=length(vs)),cex=.85)
par(xpd=NA)
points(rep(2,8),seq(20,18,length.out=8),cex=2,pch=22,col=NA,bg=rev(colorRampPalette(c("red","yellow","dark green"))(8)))
text(rep(2.3,3),seq(20,18,length.out=3),paste0("N=",seq(max(binded$N),0,length.out=3)),cex=.7,adj=0)

dev.off()
