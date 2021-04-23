## ---- ff1

load("./data/accs.per.plot.Rdata")
er=accs.per.plot

fit=read.table("./res/means_NvsS.txt", h=T, sep="\t")

sites=names(er[,1,1])
plots=names(er[1,,1])

acc=read.table("./data/acc_list.txt", h=T, sep="\t")
acc=droplevels(acc[acc$tubes<=200,])

##are there accessions that are sample more than 2 in each plot in more than one experiments

L=list()
for(s in sites){
    sub=er[s,,]
    rs=rowSums(sub)
    sub=sub[rs>0,]
    x=apply(sub, 2, function(x){sum(x>0)>=(length(x)-1)})
    #print(s)
    print(acc[match(names(x[x==T]), acc$lines),1:2])
    L[[s]]=names(x[x==T])
}
dom=unlist(L)
t1=table(dom)

fittest=names(t1[t1>=3])

for(s in sites){
    sub=er[s,,]
    print(100-100*sum(sub)/350)
}


pdf("./figures/fit_vs_overall_freq.pdf", paper="special", width=5, height=10, pointsize=8)
layout(matrix(1:8, ncol=2, byrow=T), width=c(2, 1))
par(mar=c(4, 4, 1, 1))
for(s in sites){
    fa=data.frame(id=colnames(er[s,,]), freq=colSums(er[s,,]))
    fa=fa[fa[,2]>0,]
    res=merge(fa, fit, by="id", all=F)
    col=rep("Dodgerblue", nrow(res))
    col[res$id%in%fittest]="firebrick"
    plot(res$N_fitness, res$S_fitness, cex=sqrt(res$freq), pch=16, col=col, xlab="fitness in the North", ylab="fitness in the South")
    abline(0, 1)
    subf=fit[fit$id%in%res$id==F,]
    points(subf$N_fitness, subf$S_fitness, pch=16, col="gold3")
    legend("topright", s, bty="n")
    bp=acc
    bp$neversampled=acc$lines%in%subf$id
    t1=table(bp$region, bp$neversampled)
    S=colSums(t1)
    t1[,1]=t1[,1]/S[1]
    t1[,2]=t1[,2]/S[2]
    barplot(t1, legend=T, names=c(paste("Sampled\n(N=", S[1], ")", sep=""), paste("Not Sampled\n(N=", S[2], ")", sep="")))
}
dev.off()

## ---- end-of-ff1

## sub=er[1,,]
## acc$neversampled=acc$lines%in%colnames(sub)==F
## t1=table(acc$region, acc$neversampled)
## barplot(t1, legend=T)



