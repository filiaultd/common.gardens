##bbrachi
##04/10/2018

## ---- readaf

load("./006.10k.permutation.parsed.output/all.dat.Rdata")
#load("./007.calculate.afd.from.start/pol.AFD.out.Rdata")
af=all.dat##pol.AFD.out
rm(all.dat)##pol.AFD.out)

##subset the extrem afs

ext=af[af$North==0 | af$South==0,]

x=round(table(N=(af$North==0), S=(af$South==0))*100/nrow(af), 2)

kable(x, caption="Contengency table with percentage of SNPs for which no more extrem allele frequency changes were observed in 10,000 permutations.")

##ext2=af[(abs(af$N.mean)-4*af$N.sd>0 | abs(af$S.mean)-4*af$S.sd>0) & (af$North==0 | af$South==0),]

## ---- end-of-readaf

## ---- readGWA
## read GWA res from the mlmm.

pheno="N_fitness_S_fitness"
x=unlist(strsplit(scan(paste("./GWA/output_gemma_loco/out_loco_", pheno,sep=""), what="character")[1], ","))[7]
f=substring(gsub("]]", "", x), 2, (nchar(x)-3))
p=read.delim(f, header=T, sep="\t")
p$score=-log10(p$p_wald)
p=p[p$af>=0.1,]
p$fdr_pval=p.adjust(p$p_wald, "fdr")

##top SNPs
top=p[p$fdr_pval<0.1,]
saveRDS(top,"./res/mlmm_top_SNPs_fdr10.rds")

## ---- end-of-readGWA

## ---- mlmmvsafchange

top=readRDS("./res/mlmm_top_SNPs_fdr10.rds")
af$rs=paste(af$Chromosome, "_", af$Position, sep="")

ov=af[af$rs%in%top$rs,c("rs", "N.mean", "N.sd", "S.mean", "S.sd", "North", "South")]
kable(ov, caption="mlmm top SNPs in the SNE.")

## ---- end-of-mlmmvsafchange

## ## ---- correffaf

##add the SNP effects inthe GWA to the af

##for the top 100/200 SNPs, do we see correlations between SNP effects and allele frequency changes

gaf=merge(af, p, by="rs")
gaf$GWAsignif=(gaf$fdr_pval<=0.1)

top=gaf[rank(gaf$p_wald)<=201,]


pdf("./figures/top_snps_effects_vs_afchange.pdf", paper="special", width=8, height=6)
par(mfrow=c(2,2), mar=c(4, 4, 1, 1))
for(a in c("N.mean", "S.mean")){
    AF=top[,a]
    if(a=="N.mean"){xname="AF North"}else{xname="AF South"}
    for(b in c("beta_1", "beta_2")){
        if(b=="beta_1"){yname="GWA North"}else{yname="GWA South"}
        GWA=top[,b]
        plot(GWA,AF, pch=16, col="Dodgerblue", xlab=xname, ylab=yname)
        C=cor.test(GWA,AF)
        m=lm(AF~GWA)
        sm=summary(m)
        if(sm$coefficients[2, 4]<0.05){abline(m)}
        legend("bottomleft", legend=paste("Pearson\'s corr:", round(C$estimate, 2), "\npval=", round(C$p.value, 3)), bty="n")
    }
}
dev.off()




## ---- end-of-correffaf
