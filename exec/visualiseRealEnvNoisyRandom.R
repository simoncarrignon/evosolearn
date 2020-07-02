source("../R/environment.R")
source("../R/tools.R")



namesRealEnv=c("lr04","epica","vostok","ngrip","martrat")
namesRealEnv=c("ngrip2")
namesRealEnv=c("ngrip2","vostok","martrat","ls16","ngrip","epica")

namesFun=c("interpolate_LinTRUERes20RANDOM","interpolate_LinFALSERes20RANDOM","interpolate_LinTRUERes20BEST","interpolate_LinFALSERes20BEST","interpolate_LinTRUERes20SlsrandomBIS","interpolate_LinTRUERes20SlsbestBIS")

#namesFun=c("interpolate_LinTRUERes20SlsrandomTRIS","interpolate_LinTRUERes20SlsbestTRIS")
namesFun=c(
           #"interpolate_LinTRUERes20SlsrandomBIS","interpolate_LinTRUERes20SlsbestBIS",
		   "interpolate_SlsfitpropLARGERPARAM"
           #"interpolate_SlsrandomNEWF"#,"interpolate_LinTRUERes20SlsbestTRIS",
           #"interpolate_LinFALSERes20SlsrandomLOW",
           #"interpolate_LinFALSERes20SlsrandomBIS2"
           )

combiname=expand.grid(namesRealEnv,namesFun)
allexp=apply(combiname,1,function(i)list(folder=paste0(i,collapse=""),idexpe=paste0(i,collapse=""),env=i[1],fun=i[2]))

for( exp in allexp){
    print(exp)

    idexpe=exp$idexpe
    binded=getAllFromFolder(exp$folder,pref="output")
    if(!is.null(binded)){

        ns=length(list.files(path=folder,recursive=T,full.names=T,pattern="*cross*"))

        summary=getUniqueExp(binded$filename[1])
        nsteps=nrow(summary)

        nlines=length(unique(binded$k_z))*length(unique(binded$k_y))
        ncols=length(unique(binded$E))

        binded2=binded
        for(s in unique(binded2$sigma)){

            binded=binded2[binded2$sigma==s,]

            png(paste0("images/",idexpe,"3Genes_sigma",s,"_E.png"),width=5000,height=600,pointsize=16)

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
                        pxl=pxl[nrow(pxl):0,]
                        plot(c(0, 1), c(0, 1), type = "n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
                        rasterImage(pxl, -1, 0, 1, 1,interpolate=F)
                        box()
                    }
                }
            }
            setAxis(binded)
            mtext("mean genotype",1,2,at=0.1,outer=T,adj=0,cex=1.2)

            dev.off()
            ################ pop size


            png(paste0("images/",idexpe,"N_sigma",s,"_E.png"),width=5000,height=600,pointsize=16)

            par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
            for(ky in rev(unique(binded$k_y))){
                for(kz in rev(sort(unique(binded$k_z)))){
                    for(e in unique(binded$E)){
                        subb=binded[binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,]
                        par(mar=rep(0,4))
                        res=tapply(subb$N,subb[,c("mu","m")],mean,na.rm=T)
                        image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(binded$N,na.rm=T),col=survivalpallette(3999))
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

            png(paste0("images/",idexpe,"NA_sigma",s,"_E.png"),width=5000,height=600,pointsize=16)

            nexp=c()
            par(mfrow=c(nlines,ncols),oma=c(6,6,2,5))
            for(ky in rev(unique(binded$k_y))){
                for(kz in rev(sort(unique(binded$k_z)))){
                    for(e in unique(binded$E)){
                        subb=binded[binded$E ==e & binded$k_z ==kz & binded$k_y ==ky,]
                        par(mar=rep(0,4))
                        res=tapply(subb$mean_y,subb[,c("mu","m")],function(i)sum(is.na(i)))
                        nexp=min(tapply(subb$mean_y,subb[,c("mu","m")],length))
                        image(t(res),xaxt="n",yaxt="n",ylab="",xlab="",zlim=range(-1,nexp),col=extinctpalette(nexp))
                    }
                }
            }
            setAxis(binded)
            par(xpd=NA)
            points(rep(1.3,8),seq(top,top-.5,length.out=8),cex=2,pch=22,col=NA,bg=extinctpalette(8))
            text(rep(1.5,3),seq(top,top-.5,length.out=3),paste0(round(seq(0,nexp,length.out=3))),cex=1.2,adj=0)
            text(1.5,top+.25,"#ext",cex=1.2,adj=0)
            mtext("Number of Extinctions",0,2,at=0.1,outer=T,adj=0,cex=1.2)

            dev.off()


        }
    }
}


for( exp in allexp){
    idexpe=exp$idexpe
    folder=exp$folder
    try({
    binded=getAllFromFolder(folder,pref="output")
    summary=getUniqueExp(binded$filename[1])
    nsteps=nrow(summary)

    nlines=length(unique(binded$k_z))*length(unique(binded$k_y))
    ncols=length(unique(binded$E))

    binded2=binded
    for(s in unique(binded2$sigma)){

        binded=binded2[binded2$sigma==s,]
        env=read.csv(paste0("../report/data/",exp$env,".csv"))
        for(e in unique(binded$E)){
            printRGBpixels(data=binded,filename=paste0("images/vertical_",idexpe,"_sigma",s,"_E",e,".png"),ind=F,env=env,img.width=800)
            printRGBpixels(data=binded,filename=paste0("images/verticalM_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=18000,img.height=1400,img.pointsize=42)
            printRGBpixels(data=binded,filename=paste0("images/verticalMy_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=18000,img.height=1400,img.pointsize=42,ordered="y")
            printRGBpixels(data=binded,filename=paste0("images/verticalMz_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=18000,img.height=1400,img.pointsize=42,ordered="z")
            printRGBpixels(data=binded,filename=paste0("../report/images/vertical_N_",idexpe,"_sigma",s,"_E",e,".png"),ind=F,env=env,img.width=800,varcol="N")
            printRGBpixels(data=binded,filename=paste0("../report/images/vertical_N_M_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=18000,img.height=1400,img.pointsize=42,varcol="N")
            printRGBpixels(data=binded,filename=paste0("../report/images/vertical_N_My_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=18000,img.height=1400,img.pointsize=42,ordered="y",varcol="N")
            printRGBpixels(data=binded,filename=paste0("../report/images/vertical_N_Mz_",idexpe,"_sigma",s,"_E",e,".png"),ind=T,env=env,img.width=18000,img.height=1400,img.pointsize=42,ordered="z",varcol="N")
        }
    }
    })
    print(exp$idexpe)

}
                    



