## summarize the admixture groups that each SNP occurs in in the 1001g data
## for new cluster cbe
## DLF 04Nov19

##srun --qos=short --partition=c --cpus-per-task=1 --mem=100gb --pty bash
##ml r/3.5.1-foss-2018b
##ml r-bundle-bioconductor/3.8-foss-2018b-r-3.5.1

setwd("/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/")
library(data.table)

###############################
### prep SNP and admixture data
###############################


### snp data
snp.dat <- fread("/groups/nordborg/projects/field_experiments/001.common.reference.files/010.1001g.paper.snps/1001genome.sweden.pos.gt.table",stringsAsFactors=FALSE)

### sed 's/0|0/0/g; s/1|1/1/g; s/.\/./2/g; s/0|1/0.5/g; s/1|0/0.5/g' to generate this data
### so 0=ref, 1=alt, 2=NA, 0.5=het

### need column names
vcf.head <- read.table("/groups/nordborg/projects/field_experiments/001.common.reference.files/010.1001g.paper.snps/1001g.vcf.header.txt", comment.char="")
vcf.head <- vcf.head[1,c(1,2,4,5,10:ncol(vcf.head))]
colnames(snp.dat) <- as.character(unlist(vcf.head[1,]))
colnames(snp.dat)[1:4] <- c("chr", "pos", "ref", "alt")
#save(snp.dat, file="1001g.snp.dat.Rdat")

### prep admixture data
#setwd("/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/")
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
#admix.count <- apply(snp.dat[1:10000,], 1, ad.count)

#save(admix.count, file="admix.count.Rdat")

### also need chr and position of SNPs

admix.pos <- paste(unlist(snp.dat[,1]),unlist(snp.dat[,2]), sep="_")

admix.count <- t(admix.count)
#rownames(admix.count) <- admix.pos[1:10000]
rownames(admix.count) <- admix.pos
colnames(admix.count)[1:9] <- paste("R",colnames(admix.count)[1:9], sep="")
colnames(admix.count)[10:18] <- paste("A",colnames(admix.count)[10:18], sep="")

save(admix.count, file="./data/admix.count.Rdat")