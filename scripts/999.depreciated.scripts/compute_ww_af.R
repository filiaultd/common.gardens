##bbrachi
##8/10/2018
library(readr)

##compute allele frequency of SNPs outside of sweden

gf=gzfile("./GWA/snps/1001g/filtered_snps.csv.gz")
snps=read_delim(gf, delim=",", col_names=T)


##get the geographical origin of accessions from somewhere.

acc=read.table("./GWA/snps/1001g/accessions1001g.csv", sep=",", h=T)

##remove the swedish accessions from the dataset

nonSWE=paste(acc[acc$country!="SWE","tg_ecotypeid"])

##subset the snps

snps=snps[,c("chr", "pos", nonSWE)]

##compute allele frequency of ALT alleles

rs=rowSums(snps[,-1:-2]) 

af=rs/length(nonSWE)

##write it for use

map=data.frame(snps[,1:2], af_ALT_ww=af)

saveRDS(map, "./res/freq_SNPs_nonSWE_1001g.rds")
