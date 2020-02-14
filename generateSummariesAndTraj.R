source("toolsSummariesAndTraj.R")

allcur$growingDelta$extinction=(allcur$growingDelta$extinction-1)*100+1
allcur$unifPop50k$extinction=(allcur$unifPop50k$extinction-1)*10+1



folder=c("growingDelta","unifPop50k")
omega=getAlllSummaries(c("omega1growingDelta"))
omega=updateScale(omega)
plotAllTrajVar(omega,obs="N",m=.05,E=0)


plotAllTrajVar(allcur,obs="N",m=.01,E=0)
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

tmp="delta"
for(obs in c("var_x","mean_w","N")){
}
dev.off()
}


