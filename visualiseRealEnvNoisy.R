slpalette=colorRampPalette(c(rgb(0,.5,0),rgb(0,.5,1)))
ilpalette=colorRampPalette(c(rgb(0,.5,0),rgb(1,.5,0)))
survivalpallette=colorRampPalette(c("grey","yellow","dark green"))
extinctpalette=colorRampPalette(c("chartreuse4","white"))
pureSl=slpalette(2)[2]
pureIl=ilpalette(2)[2]

realdata=read.csv("data/theta_real.csv")
theta=rev(realdata$permille)
plot(theta)
binded=do.call("rbind",lapply(list.files(path="exploreRealEnvBESTsigMas",recursive=T,full.names=T,pattern="*cross*"),function(u){print(u);load(u);return(binded)}))
binded=binded[binded$sigma==2,]
nonoise=binded
binded=do.call("rbind",lapply(list.files(path="exploreRealEnvBESTNoisy/",recursive=T,full.names=T,pattern="*cross*"),function(u){print(u);load(u);return(binded)}))
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



mum=as.data.frame(expand.grid(m=unique(binded$m),mu=unique(binded$mu)))
mum$prod=mum$m*mum$mu
lkz=length(unique(binded$k_z))
lky=length(unique(binded$k_y))
le=length(unique(binded$E))
lm=length(unique(binded$m))
lmu=length(unique(binded$mu))
for(e in unique(binded$E)){
	bicpic=matrix(nrow=length(theta),ncol=lkz*lky*lm*lmu)
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
				mat_allexp= lapply(vars,function(i)matrix(NA,nexpe,length(theta)))
				for(i in 1:nexpe){
					load(as.character(subb$filename[i]))
					for(v  in vars) mat_allexp[[v]][i,]=summary[,v]
				}
				sum_mat=lapply(mat_allexp,function(m)apply(m,2,mean,na.rm=T))
				sum_mat$na=apply(mat_allexp$mean_y,2,function(i)sum(is.na(i)))
				sum_mat=lapply(sum_mat,function(m)m[!is.na(m)])
				if(length(sum_mat$mean_y)>1){
					pxl=rgb(sum_mat$mean_y,.5,sum_mat$mean_z,alpha=1-sum_mat$na/nexpe)
					bicpic[1:length(pxl),p]=pxl
				}

				p=p+1
				print(paste("col",p,"/",ncol(bicpic)))
			}
			png(paste0("images/vertical_E",e,".png"),width=600)
			par(mar=rep(0,4),oma=c(1,4,4,1))
			prop=.1
			layout(matrix(c(1,2),ncol=2,nrow=1),width=c(prop,1-prop))
			plot(rev(realdata$permille),rev(realdata$years.BP.2000),type="l",yaxs="i",axes=F)
			plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
			axis(2,label=realdata$years.BP.2000[seq(1,length(realdata$years.BP.2000),length.out=5)],at=seq(0,1,length.out=5),outer=T)
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


