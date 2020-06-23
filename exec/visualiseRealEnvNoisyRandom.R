source("environment.R")
source("tools.R")
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
    mtext(e,side=3,line=0,at=seq(.5/length(e),1-0.5*1/length(e),length.out=length(e)),cex=.9,outer=T)

    vsp=1/nlines

    ky=sapply(unique(data$k_y),function(u)as.expression(bquote(k[y] == .(u))))
    mtext(ky,side=2,line=2,at=seq(1.5*vsp,1-1.5*vsp,length.out=length(ky)),cex=.9,outer=T)

    vsp=1/nlines
    kz=sapply(unique(data$k_z),function(u)as.expression(bquote(k[z] == .(u))))
    mtext(kz,side=2,line=0,at=seq(vsp-.5*vsp,length(kz)*vsp-.5*vsp,length.out=length(kz)),cex=.9,outer=T)
}



#' @param limit number of step (ie generation) to visualise 
#' @param ind boolean, if TRUE, individual run are ploted
#' @param subnb if ind is TRUE, number of uindividual run to plot by experiments
#' @param limit number of step (ie generation) to visualise 
#' @param prop realtive width of the environment plot wrt the rgb matrix
#' @param varcol withc variable use forthe color
#' @param palette the color palette to use
printRGBpixels <- function(data,filename,ind=FALSE,tlimit=NULL,subnb=NULL,img.width=600,img.height=800,img.pointsize=14,env=NULL,prop=.1,ordered=NULL,varcol="rgb",palette=colorRampPalette(c("grey","yellow","dark green"))(1200))
{

    mum=as.data.frame(expand.grid(m=unique(data$m),mu=unique(data$mu)))
    mum$prod=mum$m*mum$mu
    lkz=length(unique(data$k_z))
    lky=length(unique(data$k_y))
    le=length(unique(data$E))
    lm=length(unique(data$m))
    lmu=length(unique(data$mu))
    summary=getUniqueExp(data$filename[1])
    tsteps=nrow(summary)
    if(is.null(tlimit)){
        tlimit=tsteps
    }
    else{ 
        tlimit=ceiling(tlimit/unique(data$outputrate))
    }
    tlim=min(tsteps,tlimit)

    if(is.null(subnb))subnb=nexpe

    if(ind)
        bicpic=matrix(nrow=tsteps,ncol=lkz*lky*lm*lmu*subnb+1)
    else
        bicpic=matrix(nrow=tlim,ncol=lkz*lky*lm*lmu)
    bicpic[,]=NA
    p=1
    for(ky in unique(data$k_y)){
        for(kz in unique(data$k_z)){
            for(i in 1:nrow(mum)){
                mu=mum$mu[i]
                m=mum$m[i]
                subb=droplevels(data[data$mu ==mu &data$m ==m &data$E ==e & data$k_z ==kz & data$k_y ==ky,])
                nexpe=nrow(subb)
                if(ind){
                    subselect=sample.int(nexpe,subnb)
                    if(is.null(ordered))ordered="n"
                    if(ordered=="z"){
                    reorder=sapply(subselect,function(f)sum(getUniqueExp(subb$filename[f])[,"mean_z"],na.rm=T))
                    subselect=subselect[order(reorder)]
                    }
                    if(ordered=="y"){
                    reorder=sapply(subselect,function(f)sum(getUniqueExp(subb$filename[f])[,"mean_y"],na.rm=T))
                    subselect=subselect[order(reorder)]
                    }
                    for(f in subselect){
                        pxl=c()
                        summary=getUniqueExp(subb$filename[f])
                        na=which.max(is.na(summary[,1]))#find the first na ie when pop get extinct
                        if(na>1) summary=summary[1:(na-1),,drop=F]
                        if(varcol=="rgb"){
                            alphalvl=sapply(summary[,"N"]/1000,function(i)min(i,1))
                            pxl=rgb(summary[,"mean_y" ],.5,summary[,"mean_z"],alphalvl)
                        }
                        else{
                            pxl=palette[summary[,varcol]]
                        }

                        bicpic[1:length(pxl),p]=pxl
                        p=p+1
                        print(paste("individual run",p,"/",ncol(bicpic)))
                    }
                }
                else{ 
                    if(varcol=="rgb") vars=c("mean_y","mean_z")
                    else vars=varcol
                    names(vars)=vars
                    mat_allexp= lapply(vars,function(i)matrix(NA,nexpe,tlim))
                    for(i in 1:nexpe){
                        summary=getUniqueExp(subb$filename[i])
                        rng=(nrow(summary)-tlim+1):nrow(summary)
                        for(v  in vars){
                            if(length(mat_allexp[[v]][i,])!=length(summary[rng,v]))
                                print("expe is not the number of step expected")
                            else
                                mat_allexp[[v]][i,]=summary[rng,v]

                        }
                    }
                    sum_mat=lapply(mat_allexp,function(m)apply(m,2,mean,na.rm=T))
                    if(varcol=="rgb")
                        sum_mat$na=apply(mat_allexp$mean_y,2,function(i)sum(is.na(i)))
                    else
                        sum_mat$na=apply(mat_allexp[[varcol]],2,function(i)sum(is.na(i)))
                    sum_mat=lapply(sum_mat,function(m)m[!is.na(m)])
                    if(length(sum_mat$mean_y)>1){
                        if(is.na(sum_mat$na[1]))
                            pxl=NA
                        else{
                            #pxl=rgb(sum_mat$mean_y,.5,sum_mat$mean_z,alpha=1)
                        if(varcol=="rgb")    
                            pxl=rgb(sum_mat$mean_y,.5,sum_mat$mean_z,alpha=1-sum_mat$na/nexpe)
                        else    
                            pxl=alpha(palette[sum_mat[[varcol]]],alpha=1-sum_mat$na/nexpe)
                        }
                        bicpic[1:length(pxl),p]=pxl
                    }

                    p=p+1
                    print(paste("col",p,"/",ncol(bicpic)))
                }
            }
            scale=1 #to convert time scale in year
            png(filename,width=img.width,height=img.height,pointsize=img.pointsize)
            par(mar=rep(.5,4),oma=c(3,4,4,1))

            layout(matrix(c(1,2),ncol=2,nrow=1),width=c(prop,1-prop))

            plotVertEnv(env,scale)
            rasterImage(bicpic, 0, 0, 1, 1,interpolate=F)
            mtext("year BP",3,2,at=0,outer=T)
            mtext(expression(delta^18*O),3,0,at=prop/2,outer=T)

            ## first level legend (from left to right)
            kyl=sapply(unique(data$k_y),function(u)as.expression(bquote(k[y] == .(u))))
            lkyl=length(kyl)
            mtext(kyl,side=3,line=1,at=seq((1/lkyl)/2,1-(1/lkyl)/2,length.out=length(kyl)),cex=.9)

            mtext(bquote(E==.(e)),side=3,line=2)

            ## second level legend (from left to n*1/n)
            kzl=sapply(unique(data$k_z),function(u)as.expression(bquote(k[z] == .(u))))
            lkzl=length(kzl)
            mtext(kzl,side=3,line=0,at=seq(0+(1/(lkzl*lkyl))/2,lkzl/(lkzl*lkyl)-(1/(lkzl*lkyl))/2,length.out=lkzl),cex=.9)

            #mtext(mum$prod[seq(1,nrow(mum),length.out=3)],side=1,line=0,at=seq(0,1/9,length.out=3),cex=.7)
            #third level
            slvl3=1/ncol(bicpic) #scale of lvl3 legend
            if(ind)
                slvl3=slvl3*subnb
            nlvl3=length(mum$prod) #number of legend
            llvl3=mum$prod         #labels of legend
            tlvl3=seq(1/2*slvl3,slvl3*nlvl3-1/2*slvl3,length.out=nlvl3) #tick place of legend
            axis(1,at=tlvl3,label=llvl3,cex.axis=.6,las=2)
            mtext(bquote(mu %*% m),side=1,line=2,at=-1/9)
            dev.off()
        }
    }

}

plotVertEnv <- function(realdata,scale=1,...){
    plot((realdata$dTsVscales)*scale,-realdata$year,type="l",yaxs="i",axes=F)
    plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
    axis(2,label=rev(round(sort(seq(max(realdata$year),min(realdata$year),length.out=5)*scale))),at=seq(0,1,length.out=5),outer=T)
}

slpalette=colorRampPalette(c(rgb(1,.5,0),rgb(0,.5,1)))
ilpalette=colorRampPalette(c(rgb(0,.5,0),rgb(1,.5,0)))
survivalpallette=colorRampPalette(c("grey","yellow","dark green"))
extinctpalette=colorRampPalette(c("chartreuse4","white"))
pureSl=slpalette(2)[2]
pureIl=ilpalette(2)[2]


namesRealEnv=c("ls16","lr04","epica","vostok","ngrip","martrat")
#namesRealEnv=c("vostok")
namesFun=c("getMean2","getClosest")
namesFun=c("interpolate_LinTRUERes20RANDOM","interpolate_LinFALSERes20BEST")
#namesFun=c("interpolate_LinTRUERes20SlsrandomTRIS","interpolate_LinTRUERes20SlsbestTRIS")
#namesFun=c(
#           "interpolate_LinTRUERes20SlsrandomBIS","interpolate_LinTRUERes20SlsbestBIS",
#           "interpolate_LinTRUERes20SlsrandomTRIS","interpolate_LinTRUERes20SlsbestTRIS",
#           "interpolate_LinTRUERes20SlsrandomQUIS","interpolate_LinTRUERes20SlsbestQUIS"
#           )

combiname=expand.grid(namesRealEnv,namesFun)
allexp=apply(combiname,1,function(i)list(folder=paste0(i,collapse=""),idexpe=paste0(i,collapse=""),env=i[1],fun=i[2]))

for( exp in allexp){
    print(exp)

folder=exp$folder
idexpe=exp$idexpe

ns=length(list.files(path=folder,recursive=T,full.names=T,pattern="*cross*"))

binded=do.call("rbind",lapply(list.files(path=folder,recursive=T,full.names=T,pattern="*cross*"),function(u){print(u);load(u);return(binded)}))
summary=getUniqueExp(binded$filename[1])
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
	usr=par('usr')
    top=usr[4]+usr[4]*nlines
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


for( exp in allexp){
    folder=exp$folder
    idexpe=exp$idexpe
    try({
    binded=do.call("rbind",lapply(list.files(path=folder,recursive=T,full.names=T,pattern="*cross*"),function(u){print(u);load(u);return(binded)}))
    summary=getUniqueExp(binded$filename[1])
    nsteps=nrow(summary)

    nlines=length(unique(binded$k_z))*length(unique(binded$k_y))
    ncols=length(unique(binded$E))

    binded2=binded
    for(s in unique(binded2$sigma)){

        binded=binded2[binded2$sigma==s,]
        env=read.csv(paste0("data/",exp$env,".csv"))
        for(e in unique(binded$E)){
            printRGBpixels(data=binded,filename=paste0("images/vertical_",idexpe,"_sigma",s,"_E",e,".png"),ind=F,env=env,img.width=800)
            printRGBpixels(data=binded,filename=paste0("images/verticalM_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=1800,img.height=1400,img.pointsize=42)
            printRGBpixels(data=binded,filename=paste0("images/verticalMy_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=1800,img.height=1400,img.pointsize=42,ordered="y")
            printRGBpixels(data=binded,filename=paste0("images/verticalMz_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=1800,img.height=1400,img.pointsize=42,ordered="z")
        }
    }
    })
    print(exp$idexpe)

}
                    



