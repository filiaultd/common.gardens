##for each region, find the allele which is beneficial, for each SNP, compute it's average latitude and determine it's cost in the other region.

##select the top 1% SNPs

N=round(0.0005*nrow(p)) ##1%
top=p[rank(p$p_wald)<=N+1,]


##for each SNP, add the average latitute of the alternate allele
coord=read.table("./data/worldclim_swedish_acc.txt", h=T, sep="\t", stringsAsFactor=F)[,c(1, 5:6)]
##order like the SNP file to save time
row.names(coord)=coord[,1]
coord=coord[paste(fam[,1]),]

benefN=data.frame()
benefS=data.frame()

deletN=data.frame()
deletS=data.frame()


##for the North

for(i in 1:nrow(top)){
    r=geno[paste(top$rs[i]),]
    g=r[-1:-3]
    a=r[2:3]
    eN=top$beta_1[i]
    eS=top$beta_2[i]
    benefN[i, "rs"]=paste(top$rs[i])
    if(eN>0){benefN[i, "benefN"]="ALT"
        sign=1; all=2}else{benefN[i, "benefN"]="REF"
                  sign=-1; all=0
              }
    benefN[i, "avLatbN"]=mean(coord[names(g)[g==all],"lat"], na.rm=T)
    benefN[i, "effbN"]=top$beta_1[i]*sign
    benefN[i, "Southeff"]=top$beta_2[i]*sign
    ###
    benefS[i, "rs"]=paste(top$rs[i])
    if(eS>0){benefS[i, "benefS"]="ALT"
        sign=1; all=2}else{benefS[i, "benefS"]="REF"
                         sign=-1; all=0
                     }
    benefS[i, "avLatbS"]=mean(coord[names(g)[g==all],"lat"], na.rm=T)
    benefS[i, "effbS"]=top$beta_2[i]*sign
    benefS[i, "Northeff"]=top$beta_1[i]*sign

    ##deletrious

    deletN[i, "rs"]=paste(top$rs[i])
    if(eN<0){deletN[i, "deletN"]="ALT"
        sign=1; all=2}else{deletN[i, "deletN"]="REF"
                         sign=-1; all=0
                     }
    deletN[i, "avLatdN"]=mean(coord[names(g)[g==all],"lat"], na.rm=T)
    deletN[i, "effdN"]=top$beta_1[i]*sign
    deletN[i, "Southeff"]=top$beta_2[i]*sign
 
    deletS[i, "rs"]=paste(top$rs[i])
    if(eS<0){deletS[i, "deletS"]="ALT"
        sign=1; all=2}else{deletS[i, "deletS"]="REF"
                         sign=-1; all=0
                     }
    deletS[i, "avLatdS"]=mean(coord[names(g)[g==all],"lat"], na.rm=T)
    deletS[i, "effdS"]=top$beta_1[i]*sign
    deletS[i, "Northeff"]=top$beta_2[i]*sign
 
 

}


saveRDS(benefN, "./res/benefN.rds")
saveRDS(benefS, "./res/benefS.rds")

pdf("./figures/half_donut.pdf", paper="special", width=8, height=5)

layout(matrix(c(1:6), ncol=3, byrow=2))
par(mar=c(4, 4,1, 1))
plot(benefN$effbN, benefN$Southeff, pch=16, col="Dodgerblue", cex=0.8, xlab="effect of beneficial allele in the North", ylab="effect of allele in the South", xlim=c(0, 0.0012), ylim=c(-1.2e-3, 1.2e-3))
abline(0, 0)
legend("bottomright", legend=sum(benefN$Southeff<0), bty="n")
legend("topright", legend=sum(benefN$Southeff>0), bty="n")
plot(benefS$effbS, benefS$Northeff, pch=16, col="Firebrick", cex=0.8,  xlab="effect of beneficial allele in the South", ylab="effect of allele in the North", xlim=c(0, 0.0012), ylim=c(-1.2e-3, 1.2e-3))
abline(0, 0)
legend("bottomright", legend=sum(benefS$Northeff<0), bty="n")
legend("topright", legend=sum(benefS$Northeff>0), bty="n")
plot(density(benefN$avLatbN), col="Dodgerblue", xlab="average latitude of beneficial alleles", main="")
points(density(benefS$avLatbS), type="l", col="firebrick")
legend("topright", lty=c(1, 1), col=c("Dodgerblue", "firebrick"), legend=c("North", "South"), bty="n")
plot(deletN$effdN, deletN$Southeff, pch=16, col="Dodgerblue", cex=0.8, xlab="effect of deleticial allele in the North", ylab="effect of allele in the South", ylim=c(-1.2e-3, 1.2e-3))
abline(0, 0)
legend("bottomright", legend=sum(deletN$Southeff<0), bty="n")
legend("topright", legend=sum(deletN$Southeff>0), bty="n")
plot(deletS$effdS, deletS$Northeff, pch=16, col="Firebrick", cex=0.8,  xlab="effect of deleticial allele in the South", ylab="effect of allele in the North", ylim=c(-1.2e-3, 1.2e-3))
abline(0, 0)
legend("bottomright", legend=sum(deletS$Northeff<0), bty="n")
legend("topright", legend=sum(deletS$Northeff>0), bty="n")
plot(density(deletN$avLatdN), col="Dodgerblue", xlab="average latitude of deleticial alleles", main="")
points(density(deletS$avLatdS), type="l", col="firebrick")
legend("topright", lty=c(1, 1), col=c("Dodgerblue", "firebrick"), legend=c("North", "South"), bty="n")
dev.off()


##second option: polarize by strongest least negative effect

b=top[, c("beta_1", "beta_2", "af")]

x=t(apply(b,1, function(x){
    y=as.numeric(x)
    if(all(y[1:2]>0)){return(y)}
    if(all(y[1:2]<0)){return(c(-1*y[1:2], 1-y[3]))}
    if(sum(y[1:2]>0)==1){
        w=which(abs(y[1:2])==max(abs(y[1:2])))
        if(y[w]<0){return(c(-1*y[1:2], 1-y[3]))}else{return(y)}
    }
}))

segments(0, -100, 0, 100)
segments(-100, 0, 100, 0)

N=round(0.005*nrow(p)) ##1%
top=p[rank(p$p_wald)<=N+1,]

col=rep("dodgerblue", nrow(top))
col[top$af>0.5]="firebrick"
plot(top$beta_1, top$beta_2, col=col, pch=16, cex=1)
segments(0, -100, 0, 100)
segments(-100, 0, 100, 0)

col=rep("dodgerblue", nrow(top))
col[x[,3]>0.5]="firebrick"
plot(x[,1:2], col=col, pch=16, cex=1)
segments(0, -100, 0, 100)
segments(-100, 0, 100, 0)



##





##     }else{benef[i, "benefN"]="REF"
##         benef[i, "avLatbN"]=mean(coord[names(g)[g==0],"lat"], na.rm=T)
##         benef[i, "effbN"]=top$beta_1[i]*(-1)
##     }
##     if(eS>0){benef[i, "benefS"]="ALT"
##         benef[i, "avLatbS"]=mean(coord[names(g)[g==2],"lat"], na.rm=T)
##         benef[i, "effbS"]=top$beta_2[i]
##     }else{benef[i, "benefS"]="REF"
##         benef[i, "avLatbS"]=mean(coord[names(g)[g==0],"lat"], na.rm=T)
##         benef[i, "effbS"]=top$beta_2[i]*(-1)
##     }
##     delet[i, "rs"]=paste(top$rs[i])
##     if(eN<0){delet[i, "deletN"]="ALT"
##         delet[i, "avLatdN"]=mean(coord[names(g)[g==2],"lat"], na.rm=T)
##         delet[i, "effdN"]=top$beta_1[i]
##     }else{delet[i, "deletN"]="REF"
##         delet[i, "avLatdN"]=mean(coord[names(g)[g==0],"lat"], na.rm=T)
##         delet[i, "effdN"]=top$beta_1[i]*(-1)
##     }
##     if(eS<0){delet[i, "deletS"]="ALT"
##         delet[i, "avLatdS"]=mean(coord[names(g)[g==2],"lat"], na.rm=T)
##         delet[i, "effdS"]=top$beta_2[i]
##     }else{delet[i, "deletS"]="REF"
##         delet[i, "avLatdS"]=mean(coord[names(g)[g==0],"lat"], na.rm=T)
##         delet[i, "effdS"]=top$beta_2[i]*(-1)
##     }
## }

## saveRDS(benef, "./res/beneficial_allele_latitute.rds")
## saveRDS(delet, "./res/deleterious_allele_latitute.rds")
  

##






Northern=acclist$lines[acclist$region=="N Sweden"]

refacc=sample(Northern, 1)

##make vector to repolarize SNP effects

snp=geno[,refacc]
v=(geno[,refacc]-1)*-1

vf=rep(1, nrow(top))
vf[top$af>=0.5]=-1

##repolarized SNP effects

p$beta_1V=p$beta_1*v
p$beta_2V=p$beta_2*v

    
##select the top 1% SNPs

N=round(0.01*nrow(p))
top=p[rank(p$p_wald)<=N+1,]
##for each SNP, add the average latitute of the alternate allele
coord=read.table("./data/worldclim_swedish_acc.txt", h=T, sep="\t", stringsAsFactor=F)[,c(1, 5:6)]
##order like the SNP file to save time
row.names(coord)=coord[,1]
coord=coord[paste(fam[,1]),]

top$av_lat=apply(geno[paste(top$rs),(-1:-3)], 1, function(snp, coord){return(mean(coord[snp!="0", "lat"], na.rm=T))}, coord=coord)

x=colorRampPalette(c("firebrick","grey80", "Dodgerblue"))
cols <- x(10)[as.numeric(cut(top$av_lat,breaks = 10))]

layout(matrix(c(2, 1, 1, 3, 1, 1, 1, 1, 4, 1, 1, 5), ncol=4, byrow=T))
plot(top$beta_1, top$beta_2, pch=16, cex=1.2, xlab="SNP fitness effect in the North", ylab="SNP fitness effect in the South", col=cols)
segments(0, -100, 0, 100, lty=2)
segments(-100, 0, 100, 0, lty=2)
plot(density(top$av_lat))
points(density(na.omit(top$av_lat[top$beta_1<0 & top$beta_2>0])), type="l", col=2)
plot(density(top$av_lat))
points(density(na.omit(top$av_lat[top$beta_1>0 & top$beta_2>0])), type="l", col=2)
plot(density(top$av_lat))
points(density(na.omit(top$av_lat[top$beta_1>0 & top$beta_2<0])), type="l", col=2)
plot(density(top$av_lat))
points(density(na.omit(top$av_lat[top$beta_1<0 & top$beta_2<0])), type="l", col=2)
hist(top$av_lat)




## compute the average latitute of the allele that increase fitness in each region.
## once beneficial in each region

##for allele frequency. 




