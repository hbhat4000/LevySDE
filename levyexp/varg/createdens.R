rm(list=ls(all=TRUE))

load('x.RData')

thresh=10

numnor=sum(abs(x)<thresh)
dennor=length(x)
nor=c(numnor,dennor)
save(nor,file='nor10.RData')

denBR=density(x[abs(x)<thresh])
save(denBR,file='denBR10.RData')

thresh=40

numnor=sum(abs(x)<thresh)
nor=c(numnor,dennor)
save(nor,file='nor40.RData')

denBR=density(x[abs(x)<thresh])
save(denBR,file='denBR40.RData')

