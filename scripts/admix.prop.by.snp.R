## summarize the admixture groups that each SNP occurs in in the 1001g data
## DLF 20SEPT19

setwd("/lustre/scratch/projects/field_experiments/adaptation.sweden/common.gardens")
library(data.table)

###############################
### prep SNP and admixture data
###############################


### snp data
snp.dat <- fread("/lustre/scratch/projects/field_experiments/001.common.reference.files/010.1001g.paper.snps/1001genome.sweden.pos.gt.table",stringsAsFactors=FALSE)

### sed 's/0|0/0/g; s/1|1/1/g; s/.\/./2/g; s/0|1/0.5/g; s/1|0/0.5/g' to generate this data
### so 0=ref, 1=alt, 2=NA, 0.5=het

### need column names
vcf.head <- read.table("/lustre/scratch/projects/field_experiments/001.common.reference.files/010.1001g.paper.snps/1001g.vcf.header.txt", comment.char="")
vcf.head <- vcf.head[1,c(1,2,4,5,10:ncol(vcf.head))]
colnames(snp.dat) <- as.character(unlist(vcf.head[1,]))
colnames(snp.dat)[1:4] <- c("chr", "pos", "ref", "alt")
#save(snp.dat, file="1001g.snp.dat.Rdat")

### prep admixture data
setwd("/lustre/scratch/projects/field_experiments/adaptation.sweden/common.gardens")
### admixture groups for 1001g data - from the paper:
#"We defined nine groups based on these clusters and assigned each individual to a group if more than 60% of its genome
#derived from the corresponding cluster. The 135 individuals not
#matching this criterion were labeled ‘‘Admixed.’’"
admix <- read.csv(file="./data/1001genomes-accessionsand1001genomes-admixture.csv", stringsAsFactors=FALSE)
# remove USA accessions
admix <- admix[admix$country!="USA",]
admix.m <- admix[,grep("K9", colnames(admix))]
rownames(admix.m) <- admix$id

max.a <- apply(admix.m,1, function(x)(max(x, na.rm=TRUE)))
a.group <- rep(NA, nrow(admix.m))
for(up in 1:nrow(admix.m)){
  up.dat <- admix.m[up,]
  a.group[up] <- which(up.dat==max.a[up])-1
}
names(a.group)<- rownames(admix.m)
a.group <- cbind(admix[,c(1:4)], max.a,a.group, admix.m)
a.group <- a.group[which(a.group$max.a>0.6),] ##998, which is not what is reported in 1001g paper???
ag <- a.group$a.group
names(ag) <- a.group$id

### for one snp, get proportion each admixture group
ad.levs <- unique(ag)
ad.levs <- ad.levs[order(ad.levs)]
up.snp <- unlist(snp.dat[1,])

ad.count <- function(up.snp){
  up.snp <- up.snp[5:ncol(snp.dat)]
  up.r <- up.snp[up.snp==0]
  r.admix <- ag[names(ag)%in%names(up.r)]
  r.admix <- table(factor(r.admix,levels=ad.levs))
  up.a <- up.snp[up.snp==1]
  a.admix <- ag[names(ag)%in%names(up.a)]
  a.admix <- table(factor(a.admix,levels=ad.levs))
  up.h <- up.snp[up.snp==0.5]
  h.admix <- ag[names(ag)%in%names(up.h)]
  h.admix <- table(factor(h.admix,levels=ad.levs))

  r.admix <- r.admix+h.admix
  a.admix <- a.admix+h.admix

  out.dat <- c(r.admix,a.admix)
  return(out.dat)
}

admix.count <- apply(snp.dat, 1, ad.count)

#save(admix.count, file="admix.count.Rdat")

### also need chr and position of SNPs

admix.pos <- paste(unlist(snp.dat[,1]),unlist(snp.dat[,2]), sep="_")

admix.count <- t(admix.count)
rownames(admix.count) <- admix.pos
colnames(admix.count)[1:9] <- paste("R",colnames(admix.count)[1:9], sep="")
colnames(admix.count)[10:18] <- paste("A",colnames(admix.count)[10:18], sep="")


#####################################################
###### polarize this by ancestral/derived############
#####################################################

### ancestral and derived alleles were determined in /projects/field_experiments/001.common.reference.files/004.genus.vcfs

p.snps <- read.table("./data/Polarized.snps.txt",colClasses="character")
rownames(p.snps) <- gsub("Chr","",rownames(p.snps))

### 000 is cases without calls in either ref or alt species 299151
### 100 is reference ancestral  1018851 (ie alt is derived allele)
### 010 is alt ancestral  111968 
### 001 is both ref and alt alleles occur in species 351834

p.snps$rs <- rownames(p.snps)
colnames(p.snps)[1] <- "history"

ac <- merge(admix.count,p.snps, by="row.names") ###1755305 SNPs


### need 100 and 010 cases
### 100 can stay the way they are
### (going to plot the effect of the derived allele)
### for 010, need to flip both the betas and the lat degrees between centroids.

ac <- ac[ac$history %in% c("100", "010"),]  ### only 1129233 snps left
### make this a bit less unwieldy
acs <- split(ac, ac$history)
acs.alt <- acs[[1]]
### reverse columns for these guys
acs.alt <- acs.alt[,c(1,11:19,2:10,20:21)]
colnames(acs.alt)[11:19] <- paste("D", 0:8, sep="")

acs.ref <- acs[[2]]
colnames(acs.ref)[11:19] <- paste("D", 0:8, sep="")
colnames(acs.ref)[2:10] <- paste("A", 0:8, sep="")

## put everthing back together
admix.prop.pol <- rbind(acs.alt, acs.ref)
admix.prop.pol <- admix.prop.pol[,-1]
save(admix.prop.pol, file="./data/admix.prop.pol.Rdat")









