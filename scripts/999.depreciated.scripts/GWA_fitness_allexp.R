## Benjamin Brachi
## 02/10/2018

##the idea is to map traits measured in the common gardens

## ---- prepmap
library(gplots)
library(plyr)
library(mixOmics)
library(lme4)
library(igraph)
library(Hmisc)
library(vegan)
library(readr)
library(mashr)
library(Rmosek)
library(REBayes)
require(GenomicRanges)
require(TxDb.Athaliana.BioMart.plantsmart22)
##require(VariantAnnotation)
require(snpStats)
require(dbscan)
require(reshape2)
txdb=TxDb.Athaliana.BioMart.plantsmart22
genes=genes(txdb)
source("./scripts/GWA_functions.R")
source("./scripts/modif_heatmap.R")

set.seed(123)

##read the blups for heritable hubs
blups=read.table("./res/blups.txt", h=T, sep="\t")
row.names(blups)=blups[,1]
blups=blups[,-1]

##select the phenotypes of interests
trait="fitness"

blups=blups[,grep(trait, colnames(blups))]
##remove columns that have only zeros (which result from computing blups for non heritable traits

blups=blups[,apply(blups, 2, function(x){sum(x,na.rm=T)!=0 & sum(is.na(x))<(length(x)/2)})]

##regress out area from fecundity to make new late life fitness trait. 

pref="./GWA/snps/sweden_200_MAF10.bimbam"

##make a phenotype file
acclist=scan(paste(pref, ".acclist", sep=""), sep=",")



##make a table of blups in the same order as the genotype
phen=blups[paste(acclist), ]
phenfile=paste(pref, "_blups_", trait, ".phen.txt", sep="")
write.table(phen, phenfile, sep="\t", row.names=F, col.names=F)

## ---- end-of-prepmap

## ---- mapfec

##run gemma using gemma-wrapper
system("mkdir ./GWA/output_gemma_loco/")
system(paste("gemma-wrapper --loco --json --cache-dir ./GWA/.gemma-cache -- -gk -g ", pref, ".geno.gz -a ", pref, ".map -p ", phenfile, " > ./GWA/output_gemma_loco/K.json", sep=""))

for(n in 1:ncol(phen)){
    pheno=colnames(phen)[n]
    system(paste("gemma-wrapper --loco --json --cache-dir ./GWA/.gemma-cache --input ./GWA/output_gemma_loco/K.json -- -lmm 4 -g ", pref, ".geno.gz -a ", pref, ".map -p ", phenfile, " -n ", n , " > ./GWA/output_gemma_loco/out_loco_", pheno, sep=""))
}

## ---- end-of-maphh


## ---- mlmm8fit

##run gemma using gemma-wrapper
n=1:8
sel=paste(n, collapse=" ")
pheno="fit8exp"
system(paste("gemma-wrapper --loco --json --cache-dir ./GWA/.gemma-cache --input ./GWA/output_gemma_loco/K.json -- -lmm 4 -g ", pref, ".geno.gz -a ", pref, ".map -p ", phenfile, " -n ", sel , " > ./GWA/output_gemma_loco/out_loco_", pheno, sep=""))


## ---- end-of-mlmm8fit


## ---- geteff

ph2=data.frame(phen=colnames(phen), ph2=NA)
gwsignif=list()
for(n in 1:ncol(phen)){
    pheno=colnames(phen)[n]
    x=unlist(strsplit(scan(paste("./GWA/output_gemma_loco/out_loco_", pheno,sep=""), what="character")[1], ","))[7]
    f=substring(gsub("]]", "", x), 2, (nchar(x)-3))
    f2=gsub("assoc", "log", f)
    ph2$ph2[n]=ph2gemma(f2)
    p=read.delim(f, header=T, sep="\t")
    ##make a manhattan plot
    colnames(p)[3]="pos"
    colnames(p)=gsub("p_lrt", "pval", colnames(p))
    p$score=-log10(p$pval)
    p=p[p$af>=0.1,]
    p$fdr_pval=p.adjust(p$pval, "fdr")
    gwsignif[[pheno]]=p[p$fdr_pval<=0.05,]
    jpeg(paste("./GWA/manhattans/lmm_gwa_", pheno, ".jpeg", sep=""),res=600, unit="in", width=8, height=4, pointsize=12)
    manhattan(p)
    dev.off()
    if(n==1){
        beta=p[, c("rs","beta")]
        se=p[, c("rs", "se")]
    }else{
        beta=merge(p[, c("rs","beta")], beta, by="rs", all=T)
        se=merge(p[, c("rs","se")], se, by="rs", all=T)
    }
}


saveRDS(ph2, "./res/ph2_fitness.rds")
saveRDS(gwsignif, "./res/list_signif_SNPs_fitness.rds")

colnames(beta)=c("rs", colnames(phen))
colnames(se)=c("rs", colnames(phen))
row.names(beta)=beta$rs
row.names(se)=se$rs
beta=beta[, -1]
se=se[,-1]

saveRDS(beta, "./GWA/mashr/beta_fitness.rds")
saveRDS(se, "./GWA/mashr/se_fitness.rds")

## ---- end-of-geteff

## ---- mashrfitNvsS

beta=readRDS("./GWA/mashr/beta_fitness.rds")
se=readRDS("./GWA/mashr/se_fitness.rds")

##only keep a subset of the traits, as many are highly correlated
##remove traits that are not particularly interesting.


##try mashr
beta=as.matrix(na.omit(beta))
se=as.matrix(na.omit(se))
##
beta=beta[row.names(beta)%in%row.names(se),]
se=se[row.names(se)%in%row.names(beta),]
se=se[row.names(beta),]
data = mash_set_data(beta, se)
## identify a set of strong tests
out1by1=paste("./GWA/mashr/mashr_1by1_fitness_per_exp.rds", sep="")
                                        #if(file.exists(out1by1)==F){
m.1by1 = mash_1by1(mash_set_data(data$Bhat,data$Shat))
    saveRDS(m.1by1, out1by1)
#}else{
#m.1by1=readRDS(out1by1)
#}
strong.set = get_significant_results(m.1by1,0.1)
random.set=sample(1:nrow(data$Bhat), size=50000, replace=F)
data.temp = mash_set_data(data$Bhat[random.set,],data$Shat[random.set,])
Vhat = estimate_null_correlation(data.temp, z_thresh=2, apply_lower_bound=F)
rm(data.temp)
data.random = mash_set_data(data$Bhat[random.set,],data$Shat[random.set,],V=Vhat)
data.strong = mash_set_data(data$Bhat[strong.set,],data$Shat[strong.set,], V=Vhat)
U.pca = cov_pca(data.strong, ncol(beta))
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
saveRDS(m, paste("./GWA/mashr/mashr_50000_random_snps_fitness.rds", sep=""))
m=readRDS(paste("./GWA/mashr/mashr_50000_random_snps_fitness.rds", sep=""))
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
saveRDS(m2, paste("./GWA/mashr/mashr_strongset_fitness.rds", sep=""))
m2=readRDS(paste("./GWA/mashr/mashr_strongset_fitness.rds", sep=""))


## ---- end-of-runmashr

snps=names(get_significant_results(m2, thresh = 0.05, conditions = NULL,sig_fn = get_lfsr))

##cluster the SNPs by regions and plot the effects for each region

map=data.frame(rs=snps, stringsAsFactors=F)
map$chr=as.integer(substring(map$rs, 1, 1))
map$pos=as.integer(substring(map$rs, 3, ))

##sort the SNPs and cluster within each chromosomes
map=map[order(map$pos, decreasing=F), ]
map=map[order(map$chr, decreasing=F), ]
##add the significance rank
map$rank=match(map$rs, snps)


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
        h=hclust(dist(sub$pos))
        sub$clust=cutree(h, h=40000)
        for(u in unique(sub$clust)){
            subsub=sub[sub$clust==u,]
            res[fill,1:ncol(map) ]=subsub[subsub$rank==min(subsub$rank),1:4]
            res[fill,(ncol(map)+1):(ncol(map)+3)] =c(nrow(subsub), range(subsub$pos))
            fill=fill+1
        }
    }
}


## only keep cluster that include the best 10

#res=res[res$rank<=50,]
saveRDS(res, paste("./res/associated_loci_", paste(r, collapse="_"), ".rds", sep=""))

top=res$rs

b=m2$result$PosteriorMean[paste(top),]
flsr=m2$result$lfsr[paste(top),]
fldr=m2$result$lfdr[paste(top),]

b[flsr>=0.01]=0
b[fldr>=0.01]=0

##fix the column names to match the names of sites in the paper

change_exp_names=function(x){
    exps1=c("ULL", "RAT", "RAM", "ADA")
    exps2=c("S1", "S2", "N1", "N2")
    for(i in 1:4){
        x=gsub(exps1[i], exps2[i], x)
    }
    return(x)
}

colnames(b)=change_exp_names(colnames(b))
colnames(b)=gsub(paste("_",r, sep=""), "", colnames(b))
 
x=colorRampPalette(c("Dodgerblue", "green", "gold", "Firebrick"))
cr=x(200)
pdf(paste("./figures/heatmap_mashr_", paste(r, collapse="_"),".pdf", sep=""), paper="special", width=7, height=7,pointsize=8)
heatmap.2(b, scale="none", col=cr, margins=c(8, 8))
dev.off()

