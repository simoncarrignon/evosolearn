source("environment.R")
source("tools.R")
setAxis <- function(data){
    mus=sapply(sort(unique(data$mu)),function(u)as.expression(bquote(mu == .(u))))
    tm=min(5,length(mus))
    nl=seq(1,length(mus),length.out=tm)
    axis(1,label=mus[nl],at=seq(0,1,length.out=length(nl)),cex.axis=1.1,las=3)

    ms=paste0("m=",sort(unique(data$m)))
    tm=min(4,length(ms))
    nl=seq(1,length(ms),length.out=tm)
    axis(4,label=ms[nl],at=seq(0,1,length.out=length(nl)),las=1,cex.axis=1.1)

    e=sapply(sort(unique(data$E)),function(u)as.expression(bquote(E == .(u))))
    mtext(e,side=3,line=0,at=seq(.5/length(e),1-0.5*1/length(e),length.out=length(e)),cex=.9,outer=T)

    vsp=1/nlines

    ky=sapply(sort(unique(data$k_y)),function(u)as.expression(bquote(k[y] == .(u))))
    mtext(ky,side=2,line=2,at=seq(1.5*vsp,1-1.5*vsp,length.out=length(ky)),cex=.9,outer=T)

    vsp=1/nlines
    kz=sapply(sort(unique(data$k_z)),function(u)as.expression(bquote(k[z] == .(u))))
    mtext(kz,side=2,line=0,at=seq(vsp-.5*vsp,length(kz)*vsp-.5*vsp,length.out=length(kz)),cex=.9,outer=T)
}

plotVertEnv <- function(realdata,scale=1,...){
    plot((realdata$dTsVscales)*scale,-realdata$year,type="l",yaxs="i",axes=F)
    plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
    axis(2,label=rev(round(sort(seq(max(realdata$year),min(realdata$year),length.out=5)*scale))),at=seq(0,1,length.out=5),outer=T)
                mtext("year BP",3,2,at=0,outer=T)
                mtext(expression(delta^18*O),3,0,at=prop/2,outer=T)
}

slpalette=colorRampPalette(c(rgb(1,.5,0),rgb(0,.5,1)))
ilpalette=colorRampPalette(c(rgb(0,.5,0),rgb(1,.5,0)))
scprepalette=colorRampPalette(c("white","dark red"))
survivalpallette=colorRampPalette(c("grey","yellow","dark green"))
extinctpalette=colorRampPalette(c("chartreuse4","white"))
pureSl=slpalette(2)[2]
pureIl=ilpalette(2)[2]


combiname=expand.grid(namesRealEnv,namesFun)
allexp=apply(combiname,1,function(i)list(folder=paste0(i,collapse=""),idexpe=paste0(i,collapse=""),env=i[1],fun=i[2]))

for( exp in allexp){
    print(exp)

folder="ngrip2_Slsrandom_DIST/"
idexpe="dist"

ns=length(list.files(path=folder,recursive=T,full.names=T,pattern="*cross*"))

distances=getAllFromFolder(folder)
onesum=readRDS(file=as.character(binded$filename[1]))
nsteps=nrow(onesum)


columnparam=c("mu","K","m","E","sigma","k_z","k_y")

simpleDS=unique(distances[,columnparam])
allscore=apply(simpleDS,1,distFivePercent,data=distances)
distances=cbind(simpleDS,allscore)

clm=unique(binded2[2050,columnparam])

binded2=distances
for(s in unique(binded2$sigma)){
    png(paste0("images/",idexpe,"Score_sigma",s,"_E.png"),width=500,height=600,pointsize=16)

    binded=binded2[binded2$sigma==s,]
nlines=length(unique(binded$k_z))*length(unique(binded$k_y))
ncols=length(unique(binded$E))

    par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
    for(ky in rev(sort(unique(binded$k_y)))){
        for(kz in rev(sort(unique(binded$k_z)))){
            for(e in unique(binded$E)){
                subb=binded[binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,]
                par(mar=rep(0,4))
                res=tapply(subb$allscore,subb[,c("mu","m")],unique,na.rm=T)
                image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$allscore,na.rm=T),col=scprepalette(100))
            }
        }
    }
    setAxis(binded)
    #par(xpd=NA)
	#usr=par('usr')
    #top=usr[4]+usr[4]*nlines
    #points(rep(1.3,8),seq(top,top-.5,length.out=8),cex=2,pch=22,col=NA,bg=rev(survivalpallette(8)))
    #text(rep(1.5,3),seq(top,top-.5,length.out=3),paste0(round(seq(max(binded$N,na.rm=T),0,length.out=3))),cex=1.2,adj=0)
    #text(1.3,top+.25,expression(bar(N)[e]),cex=1.2,adj=0)
    #points(1.5,top-.5-.5,cex=2,pch=22,col=1,bg="white")
    #text(1.5,top-.5-.5,"NA",cex=1.2,adj=0)

    #mtext("Mean population size",1,2,at=0.1,outer=T,adj=0,cex=1.2)
    dev.off()
}


png("posteriorKY.png")
boxplot(distances$allscore ~ distances$k_y)
dev.off()
png("posteriorKZ.png")
boxplot(distances$allscore ~ distances$k_z)
dev.off()
png("posteriorE.png")
boxplot(distances$allscore ~ distances$E)
dev.off()
