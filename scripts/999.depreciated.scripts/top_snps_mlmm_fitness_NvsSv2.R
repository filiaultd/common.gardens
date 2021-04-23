
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
p$fdr_pval=p.adjust(p$p_wald,"fdr")

##top SNPs

## ---- end-of-readGWA

##read the SNP matrix

pref="./GWA/snps/sweden_200_MAF10"
##read the genotypes
gz=gzfile(paste(pref, ".bimbam.geno.gz", sep=""))
geno=read.delim(gz,header=F, sep=",")
fam=read.delim(paste(pref,".fam", sep=""), h=F, sep=" ")
map=read.delim(paste(pref, ".bim",sep=""), h=F, sep="\t")

row.names(geno)=map[,2] colnames(geno)[-1:-3]=fam[,1]
colnames(geno)[1:3]=c("snp", "ALT", "REF")

## acclist=read.table("./data/acc_list.txt", sep="\t", h=T,
stringsAsFactor=F)

##polarize SNP by effect of most common allele in nonSWE.

top=p[rank(p$p_wald)<1001,]

wh=list()
top$quad[top[,"beta_1"]>0 & top[,"beta_2"]>0]=1; wh[[1]]=c(1, 1)
top$quad[top[,"beta_1"]>0 & top[,"beta_2"]<0]=3; wh[[3]]=c(1, -1)
top$quad[top[,"beta_1"]<0 & top[,"beta_2"]<0]=2; wh[[2]]=c(-1, -1)
top$quad[top[,"beta_1"]<0 & top[,"beta_2"]>0]=4; wh[[4]]=c(-1, 1)
tb=100*table(top$quad)/nrow(top)
tc=c(0.00018)


##
CE=top[(top$beta_1*top$beta_2)>0,]
DE=top[(top$beta_1*top$beta_2)<0,]

##for SNPs with different effects in North and South
##find the allele with the positive effect and see if it's more frequent in the north or the South.

##for each SNP, add the average latitute of the alternate allele
coord=read.table("./data/worldclim_swedish_acc.txt", h=T, sep="\t", stringsAsFactor=F)[,c(1, 5:6)]
##order like the SNP file to save time
row.names(coord)=coord[,1]
coord=coord[paste(fam[,1]),]

##read acclist and sort it by fam
acclist=read.table("./data/acc_list.txt", h=T, sep="\t", stringsAsFactors=F)
acclist=droplevels(acclist[acclist$lines!=".",])
row.names(acclist)=acclist$lines
acclist=acclist[paste(fam[,1]),]

##is the allele with positive effect in the South (beta_2>0) 
for(s in 1:nrow(DE)){
    nb=DE[s, c("beta_1", "beta_2")]
    w=which(abs(nb[1,])==max(abs(nb[1,])))
    if(nb[w]>0){all=2}else{all=0; nb=-1*nb}
    if(w==1){home="N Sweden";away="S Sweden"}else{home="S Sweden"; away="N Sweden"}
    ghome=na.omit(as.numeric(geno[paste(DE$rs[s]),paste(acclist$lines[acclist$region==home])]))
    gaway=na.omit(as.numeric(geno[paste(DE$rs[s]),paste(acclist$lines[acclist$region==away])]))       
    DE$fhome[s]=sum(ghome==all)/length(ghome)
    DE$faway[s]=sum(gaway==all)/length(gaway)
}


##for CE 

fWW=readRDS("./res/freq_SNPs_nonSWE_1001g.rds")
row.names(fWW)=paste(fWW[,1], "_" ,fWW[,2], sep="")

for(s in 1:nrow(CE)){
    nb=CE[s, c("beta_1", "beta_2")]
    w=which(abs(nb[1,])==max(abs(nb[1,])))
    if(all(nb<0)){all=0; CE$fSWE[s]=1-CE$af[s]; CE$fWW[s]=1-fWW[paste(CE$rs[s]), 3]}else{all=2; CE$fSWE[s]=CE$af[s]; CE$fWW[s]=fWW[paste(CE$rs[s]), 3]}
}


##make a plot 

#cols=c("gold", "grey20", "forestgreen","blue")
CE$col="grey30"##cols[CE$quad]
DE$col="grey30"##cols[DE$quad]

pdf("./figures/home_vs_away_freq.pdf", paper="special", width=8, height=5)
par(mfrow=c(1, 2))
plot(DE$fhome, DE$faway, col=DE$col, pch=16, xlab="Home frequency", ylab="Away frequency", main="Frequency of the allele\nwith the highest positive effect")
abline(0, 1, lty=2)
plot(CE$fSWE, CE$fWW, col=DE$col, pch=16, xlab="Frequency in Sweden", ylab="Frequency in everywhere else", main="Frequency of allele with consistent\npositive effect across region")
abline(0, 1, lty=2)
dev.off()
