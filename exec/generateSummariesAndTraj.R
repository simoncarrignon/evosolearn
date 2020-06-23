source("toolsSummariesAndTraj.R")




nullOmeg=getAlllSummaries(c("growingDelta","unifPop50k"))

    allsummaries=lapply(folder,getFullExperimentSummary)
    names(allsummaries)=folder
    allexperiments=lapply(allsummaries,function(curr)
                          {
                              allTl=lapply(curr$filename,getTraj,var="N")
                              allLen=sapply(allTl,function(i)length(i[!is.na(i)]))
                              cbind(curr,extinction=allLen)
                          }
    )
    nullOmeg=do.call("rbind",allexperiments)
nullOmeg=updateScale(nullOmeg)
plotAllTrajVar(nullOmeg,obs="N",m=.01,E=0)
plotAllVariableSummaries(nullOmeg,E=0,estimate=eq2833b)

allexperiments$growingDelta$extinction=(allexperiments$growingDelta$extinction-1)*100+1
allexperiments$unifPop50k$extinction=(allexperiments$unifPop50k$extinction-1)*10+1

omega=getAlllSummaries(c("omega1growingDelta"))
omega=updateScale(omega)

omega2=getAlllSummaries(c("omega2growingDelta"))
omega2=updateScale(omega2)

plotAllTrajVar(omega,obs="N",m=.05,E=0)

pdf("allvarianceOmega1.pdf",width=11,height=12)
plotAllVariableSummaries(omega,E=0,estimate=eq2833b)
dev.off()
pdf("allvarianceOmega2.pdf",width=11,height=12)
plotAllVariableSummaries(omega2,E=0,estimate=eq2833b)
dev.off()

pdf("traj_distoptimumOmega2.pdf",width=11,height=12)
plotAllTrajVar(omega2,obs="dist",m=.05,E=0)
dev.off()

pdf("traj_distoptimumOmega1.pdf",width=11,height=12)
plotAllTrajVar(omega,obs="dist",m=.05,E=0)
dev.off()

plotAllTrajVar(omega2,obs="dist",m=.05,E=0)
plotAllVariableSummaries(allcur,E=0,estimate=eq2833b)
## need to adjust time step as the resolutions is not the same for all setups

allcur$growingDelta$extinction=(allcur$growingDelta$extinction-1)*100+1
allcur$unifPop50k$extinction=(allcur$unifPop50k$extinction-1)*10+1
allcur=do.call("rbind",allcur)

Ks=unique(current$K)
mus=unique(current$mu)
ms=unique(current$m)
sigmas=unique(current$sigma)
deltas=unique(current$delta)
E=0

options(scipen=999)

scale=list(var_x=list("0.1"=c(0,0.004),"0.2"=c(0,0.002),"0.4"=c(0,0.005),"10000"=c(0,0.15)),
           N=c(0,2100),
           mean_w=c(.2,1)
           )
side=list(var_x="topleft",N="bottomright",mean_w="bottomright")

mean(replicate(10,{
t=environment(10000,2,1)
vtm=t[1:(length(t)-1)]
vt=t[2:length(t)]
cov(vt,vtm)
           }))


