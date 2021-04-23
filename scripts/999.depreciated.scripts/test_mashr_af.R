## ---- mashraf
library(mashr)

beta=af[,c("N.mean", "S.mean")]
se=af[,c("N.sd", "S.sd")]*2 ## two standard deviation... 

row.names(beta)=af$rs
row.names(se)=af$rs

beta=as.matrix(na.omit(beta))
se=as.matrix(na.omit(se))
##
beta=beta[row.names(beta)%in%row.names(se),]
se=se[row.names(se)%in%row.names(beta),]
se=se[row.names(beta),]
data = mash_set_data(beta, se)
## identify a set of strong tests
out1by1=paste("./GWA/mashr/mashr_1by1_af.rds", sep="")
#if(file.exists(out1by1)==F){
    m.1by1 = mash_1by1(mash_set_data(data$Bhat,data$Shat))
    saveRDS(m.1by1, out1by1)
#}else{
#    m.1by1=readRDS(out1by1)
#}
strong.set = get_significant_results(m.1by1,0.0001)
random.set=sample(1:nrow(data$Bhat), size=100000, replace=F)
data.temp = mash_set_data(data$Bhat[random.set,],data$Shat[random.set,])
Vhat = estimate_null_correlation(data.temp, z_thresh=2, apply_lower_bound=F)
rm(data.temp)
data.random = mash_set_data(data$Bhat[random.set,],data$Shat[random.set,],V=Vhat)
data.strong = mash_set_data(data$Bhat[strong.set,],data$Shat[strong.set,], V=Vhat)
U.pca = cov_pca(data.strong, ncol(beta))
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
saveRDS(m, paste("./GWA/mashr/mashr_50000_random_snps_af.rds", sep=""))
m=readRDS(paste("./GWA/mashr/mashr_50000_random_snps_af.rds", sep=""))
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
saveRDS(m2, paste("./GWA/mashr/mashr_strongset_af.rds", sep=""))
m2=readRDS(paste("./GWA/mashr/mashr_strongset_af.rds", sep=""))


snps=names(get_significant_results(m2, thresh = 0.001, conditions = NULL,sig_fn = get_lfsr))

##cluster the SNPs by regions and plot the effects for each region

map=data.frame(rs=snps, stringsAsFactors=F)
map$chr=as.integer(substring(map$rs, 1, 1))
map$pos=as.integer(substring(map$rs, 3, ))

##sort the SNPs and cluster within each chromosomes
map=map[order(map$pos, decreasing=F), ]
map=map[order(map$chr, decreasing=F), ]
##add the significance rank
map$rank=match(map$rs, snps)

map=map[map$rank<=1000,]

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

## ---- end-of-mashraf


## ---- clusteraf

res=data.frame(matrix(ncol=ncol(map)+3))
colnames(res)=c(colnames(map), c("n_snps", "start", "end"))
fill=1
chrs=unique(map$chr)
for(chr in chrs){
    sub=ext[ext$Chromosome==chr,]
    if(nrow(sub)==1){sub$clust=1; 
        res[fill,1:ncol(ext) ]=sub[sub$rank==min(sub$rank),1:4]
        res[fill,(ncol(ext)+1):(ncol(ext)+3)] =c(nrow(sub), range(sub$Position))
        fill=fill+1
    }else{
        h=hclust(dist(sub$Position))
        sub$clust=cutree(h, h=40000)
        for(u in unique(sub$clust)){
            subsub=sub[sub$clust==u,]
            res[fill,1:ncol(ext) ]=subsub[subsub$rank==min(subsub$rank),1:4]
            res[fill,(ncol(ext)+1):(ncol(ext)+3)] =c(nrow(subsub), range(subsub$Position))
            fill=fill+1
        }
    }
}



## ---- end-of-clusteraf

