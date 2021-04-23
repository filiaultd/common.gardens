##bbrachi
##4/10/2018

## ---- fitlat

library(ggmap)
library(grid)

####local adaptation on the mean variation of fitness in the N and in the South.

##read in data
phen=read.table("./res/means_NvsS.txt", sep="\t", h=T, stringsAsFactors = T)
clim=read.table("./data/worldclim_swedish_acc.txt", sep="\t", h=T, stringsAsFactor=T )
clim_exp=read.table("./data/worldclim_swedish_exp.txt", , sep="\t", h=T)
colnames(clim)[1]="id"

##merge
phen=merge(phen, clim[, c(1, 3, 5,6)])

##plot and estimate correlations
bootCor=function(x, N=1000,method="spearman"){
    res=c()
    for(i in 1:N){
        s=sample(1:nrow(x), size=nrow(x), replace=T)
        res=c(res, cor.test(x[s,1], x[s,2], method=method)$estimate)
    }
    return(res)
}

rN=bootCor(phen[,c("lat", "N_fitness")])
rS=bootCor(phen[,c("lat", "S_fitness")])

##reorder factor levels for the region

phen$region = factor(phen$region,levels(phen$region)[c(2, 1, 3)])

pdf("./figures/fitness_N_S_lat.pdf", paper="special", height=6, width=4)
m=matrix(c(1,1, 2,2,3,3, 3, 3,4,4,4,5), ncol=4, byrow=T)
layout(m)
par(mar=c(3, 3, 1, 1), mgp=c(1.5, 0.5, 0))
plot(phen$lat, phen$N_fitness, cex=1, col="dodgerblue3", pch=16, xlab="latitude", ylab="fitness in the North", ylim=c(0, 0.012))
legend("topright", "A", bty="n", cex=1.5)
plot(phen$lat, phen$S_fitness, cex=1, col="firebrick3", pch=16, xlab="latitude", ylab="fitness in the South", ylim=c(0, 0.012))
legend("topright", "B", bty="n", cex=1.5)
hist(rN, xlim=range(c(rN, rS)), ylim=c(0, 100), col="dodgerblue3", breaks=50, main="", xlab="rho")
q=quantile(rN, c(0.025, 0.957))
b=mean(rN)
segments(q, 0,q, 110, col="dodgerblue3", lwd=2, lty=2)
segments(b, 0,b, 110, col="dodgerblue3", lwd=2, lty=1)
hist(rS, xlim=range(c(rN, rS)), col="firebrick3", add=T, breaks=50)
q=quantile(rS, c(0.025, 0.957))
b=mean(rS)
segments(b, 0,b, 110, col="firebrick3", lwd=2, lty=1)
segments(q, 0,q, 110, col="firebrick3", lwd=2, lty=2)
box()
legend("topright", "C", bty="n", cex=1.5)
mypal=colorRampPalette(c("firebrick3", "dodgerblue3"))
cols=mypal(200)
r=range(phen$lat, na.rm=T)
plot(phen$N_fitness, phen$S_fitness, col=cols[findInterval(phen$lat, seq(r[1],r[2],length=200))], pch=16, cex=1, xlab="fitness in the North", ylab="fitness in the South")
abline(0,1, lwd=1, col=1, lty=2)
x=cor.test(phen$N_fitness, phen$S_fitness, method="spearman")
legend("topright", "D", cex=1.5,bty="n")
legend("topleft",legend=bquote(paste(rho, "=", .(round(x$estimate, 3)), "; ", italic(p),"-value = ", .(format(x$p.value, scientific=T, digits=3)))), cex=0.8, bty="n")
par(mar=c(3, 0.2 ,1, 1))
plot(0, 0, xlim=c(0,2), ylim=c(0, 1), type="n", bty="n", axes=F, xlab="", ylab="")
for(i in 1:200){
    rect(0,(i-1)/200,0.3,i/200, col=cols[i], border=FALSE, lwd=0)
}
a=0
for(i in seq(r[1], r[2], length.out=6)){
    text(0.5, a, round(i, 2), cex=1, adj=-0.5)
    a=a+1/5
}
dev.off()

## ---- end-of-fitlat

