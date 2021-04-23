## Benjamin Brachi
## 07/10/2018

library(readr)
library(snpStats)

## read in the mlmm results.

## ---- readGWA

pheno="N_fitness_S_fitness"
x=unlist(strsplit(scan(paste("./GWA/output_gemma_loco/out_loco_", pheno,sep=""), what="character")[1], ","))[7]
f=substring(gsub("]]", "", x), 2, (nchar(x)-3))
p=read.delim(f, header=T, sep="\t")

p$score=-log10(p$p_wald)
p=p[p$af>=0.1,]
p$fdr_pval=p.adjust(p$p_wald, "fdr")

##top SNPs

## ---- end-of-readGWA

##read the SNP matrix

pref="./GWA/snps/sweden_200_MAF10"

##read the genotypes
gz=gzfile(paste(pref, ".bimbam.geno.gz", sep=""))
geno=read.delim(gz,header=F, sep=",") 
fam=read.delim(paste(pref, ".fam", sep=""), h=F, sep=" ")
map=read.delim(paste(pref, ".bim", sep=""), h=F, sep="\t")

row.names(geno)=map[,2]
colnames(geno)[-1:-3]=fam[,1]
colnames(geno)[1:3]=c("snp", "ALT", "REF")

##
acclist=read.table("./data/acc_list.txt", sep="\t", h=T, stringsAsFactor=F)


##polarize SNP by effect of most common allele in nonSWE. 

top=p[rank(p$p_wald)<1001,]

##cluster SNPs in top to keep one per peak

map=cbind(top[,c("rs", "chr", "ps")], nb)
map=map[order(map$ps),]
map=map[order(map$chr),]

res=data.frame(matrix(ncol=ncol(map)+3))
colnames(res)=c(colnames(map), c("n_snps", "start", "end"))
fill=1
chrs=unique(map$chr)

for(chr in chrs){
    sub=map[map$chr==chr,]
    if(nrow(sub)==1){sub$clust=1; 
        res[fill,1:ncol(map) ]=sub[sub$rank==min(sub$rank),1:4]
        res[fill,(ncol(map)+1):(ncol(map)+3)] =c(nrow(sub), range(sub$pos))
        fill=fill+1
    }else{
        h=hclust(dist(sub$ps))
        sub$clust=cutree(h, h=40000)
        for(u in unique(sub$clust)){
            subsub=sub[sub$clust==u,]
            res[fill,1:ncol(map) ]=subsub[abs(subsub$beta_1)==max(abs(subsub$beta_1)) | abs(subsub$beta_2)==max(abs(subsub$beta_2)), 1:4][1,]
            res[fill,(ncol(map)+1):(ncol(map)+3)] =c(nrow(subsub), range(subsub$ps))
            fill=fill+1
        }
    }
}



afww=readRDS("./res/freq_SNPs_nonSWE_1001g.rds")

afww$rs=paste(afww$chr,"_", afww$pos,sep="")
afww=afww[,c("rs", "af_ALT_ww")]

top=merge(top, afww, by="rs", all.x=T, all.y=F)
          
##how many missing SNPs
sum(is.na(top$af_ALT_ww))
##11 missing top SNP, is acceptable

##if the alternate allele has an allele frequency below 0.5, outside of sweden, switch the effects signs
sw=rep(1, nrow(top))
sw[top$af_ALT_ww<0.5]=-1



##for each SNP, add the average latitute of the alternate allele
coord=read.table("./data/worldclim_swedish_acc.txt", h=T, sep="\t", stringsAsFactor=F)[,c(1, 5:6)]
##order like the SNP file to save time
row.names(coord)=coord[,1]
coord=coord[paste(fam[,1]),]

nb=top[, c("beta_1", "beta_2")]*sw
row.names(nb)=top$rs
nb$afSW=top$af
nb$afSW[sw==-1]=1-nb$afSW[sw==-1]
nb$afGL=top$af_ALT_ww
nb$afGL[sw==-1]=1-top$af_ALT_ww[sw==-1]

##subset the swedish SNPs
x=geno[row.names(nb), -1:-3]
acclist=acclist[acclist$lines!=".",]
row.names(acclist)=acclist$lines
regions=acclist[colnames(x), "region"]

for(s in 1:nrow(top)){
    if(sw[s]==1){subacc=na.omit(colnames(x)[x[s,]==2])
        #nb$afCA[s]=top$af[s]
        #nb$afGL[s]=top$af_ALT_ww[s]
    }else{
        subacc=na.omit(colnames(x)[x[s,]==0])
        #nb$afCA[s]=1-top$af[s]
        #nb$afGL[s]=top$af_ALT_ww[s]
    }
    nb$avlatCA[s]=mean(coord$lat[coord$lines%in%subacc], na.rm=T)
    nb$avlatRA[s]=mean(coord$lat[coord$lines%in%subacc==F], na.rm=T)
    if(sw[s]==1){all=2}else{all=0}
    nb$dfreqNorth[s]=sum(x[s, regions=="N Sweden"]==all, na.rm=T)/length(na.omit(x[s, regions=="N Sweden"]))-nb$afGL[s]
    nb$dfreqSouth[s]=sum(x[s, regions=="S Sweden"]==all, na.rm=T)/length(na.omit(x[s, regions=="S Sweden"]))-nb$afGL[s]
    nb$dfreqCentral[s]=sum(x[s, regions=="C Sweden"]==all, na.rm=T)/length(na.omit(x[s, regions=="C Sweden"]))-nb$afGL[s]
}


wh=list()
nb$quad[nb[,1]>0 & nb[,2]>0]=1; wh[[1]]=c(1, 1)
nb$quad[nb[,1]>0 & nb[,2]<0]=3; wh[[3]]=c(1, -1)
nb$quad[nb[,1]<0 & nb[,2]<0]=2; wh[[2]]=c(-1, -1)
nb$quad[nb[,1]<0 & nb[,2]>0]=4; wh[[4]]=c(-1, 1)
tb=100*table(nb$quad)/nrow(nb)
tc=c(0.00018)

cols=c("gold", "grey20", "forestgreen","blue")
nb$col=cols[nb$quad]

pdf("./figures/donut_af1001g.pdf", paper="special", width=12, height=5)
m1=matrix(c(0, 0, 1, 0, 0, 0, 3, 1, 2,0, 1, 1, 1, 1,1, 0, 4, 1, 5, 0,0, 0, 1, 0, 0) , ncol=5, byrow=T)
m2=m1+5
m2[m2==5]=0
m=cbind(m1, m2)
layout(m, width=c(0.1, 1, 1, 1, 0.1, 0.1,1, 1, 1, 0.1), heights=c(0.1, 1, 1, 1, 0.1)) 
par(cex.lab=1.2, mar=c(5, 5, 1, 1))
plot(nb[,1:2], col=nb$col, pch=16, cex=1, xlab="effect of globally common allele on fitness in the North"
   , ylab="effect of globally common allele on fitness in the South", xlim=c(-2e-3, 2e-3), ylim=c(-2e-3, 2e-3))
segments(0, -100, 0, 100)
segments(-100, 0, 100, 0)
for(q in 1:4){
    c1=(wh[[q]]*tc)[1]; c2=(wh[[q]]*tc)[2];
    text(c1, c2, labels=paste(round(tb[q], 2),"%", sep=""))
}
freqCols=c("freqSouth","freqNorth")
namesCols=c("S", "N")
##add the distributions of 
sub=nb[nb[,1]>0 & nb[,2]>0,]
boxplot(sub[, c("afSW", "afGL")], outline=F)
stripchart(sub[, c("afSW", "afGL")], col="Dodgerblue", pch=16, cex=0.5, method="jitter", jitter=0.1, vertical=T, add=T)
sub=nb[nb[,1]<0 & nb[,2]>0,]
plot(density(na.omit(sub$dfreqSouth)), col="firebrick", xlim=c(-1, 1), main="", xlab="")
points(density(na.omit(sub$dfreqNorth)), type="l", col="Dodgerblue")
segments(0, 100, 0,0) 
sub=nb[nb[,1]<0 & nb[,2]<0,]
boxplot(sub[, c("afSW", "afGL")], outline=F)
stripchart(sub[, c("afSW", "afGL")], col="Dodgerblue", pch=16, cex=0.5, method="jitter", jitter=0.1, vertical=T, add=T)
sub=nb[nb[,1]>0 & nb[,2]<0,]
plot(density(na.omit(sub$dfreqSouth)), col="firebrick", xlim=c(-1, 1), main="", xlab="")
points(density(na.omit(sub$dfreqNorth)), type="l", col="Dodgerblue")
segments(0, 100, 0,0)
##
nb=na.omit(nb)
plot(nb$dfreqSouth, nb$dfreqNorth, col=nb$col, pch=16, xlab="Delta freq: South-Global", ylab="Delta freq: North-Global")
segments(0, -100, 0, 100, lty=2)
abline(0, 1, , lty=3)
segments(-100, 0, 100, 0, lty=2)
nb$fq[nb$dfreqSouth>0 & nb$dfreqNorth>0]=1
nb$fq[nb$dfreqSouth<0 & nb$dfreqNorth>0]=2
nb$fq[nb$dfreqSouth>0 & nb$dfreqNorth<0]=4
nb$fq[nb$dfreqSouth<0 & nb$dfreqNorth<0]=3
for(q in 1:4){
    sub=nb[nb$fq==q,]
    t=table(sub$quad)
    pie(t, col=cols)
}
dev.off()


library(yarrr)




