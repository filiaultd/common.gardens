## summarize the allele counts per K group using Fernando's Swedish data
## 25AUG20

##srun --qos=short --partition=c --cpus-per-task=1 --mem=100gb --pty bash
##ml r/3.5.1-foss-2018b
##ml r-bundle-bioconductor/3.8-foss-2018b-r-3.5.1

setwd("/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/")
library(data.table)

all.gts.file <- "/groups/nordborg/projects/field_experiments/001.common.reference.files/006.field.SNPs.Fernando/02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.table"
header.file <- "/groups/nordborg/projects/field_experiments/001.common.reference.files/006.field.SNPs.Fernando/02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.csv"
### this is actually still a csv
### ran:
### sed 's/-1/NA/g; s/2/0.5/g' to change encoding of 02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.csv to .table
### so 0=ref, 1=alt, 2=NA, 0.5=het
### original encoding is -1=NA, 0=ref, 1=alt, 2=het
### unfortunately, this also impacts column names and position information, so need to recover those from .csv file
kg.file <- "./data/K.groups.txt"

###############################
### prep SNP and admixture data
###############################

### snp data
all.gts <- fread(all.gts.file,stringsAsFactors=FALSE,colClasses="numeric",header=TRUE)
hs <- read.csv(header.file, nrows=1, header=TRUE)
colnames(all.gts) <- gsub("X","",colnames(hs))
pos <- fread(header.file,stringsAsFactors=FALSE,colClasses="numeric",header=TRUE)
pos <- pos[,.(Chromosome, Position)]

### K.groups data
K.groups <- read.table(kg.file,stringsAsFactors=FALSE)
colnames(K.groups) <- "K.group"
K.groups$id <- rownames(K.groups)

kcol <- brewer.pal(6, "Paired")[c(6,1,2,5,3,4)]
ki <- data.frame(matrix(c("1","S1","2","S2","3","C","4","N1","5","N2","6","B"),ncol=2,byrow=TRUE))
colnames(ki) <- c("K.group","K.name")
K.groups <- merge(K.groups, ki, by="K.group")
K.groups$K.name <- factor(K.groups$K.name, levels=c("B","S1","S2","C","N1","N2"))

### subset snp data to lines in K.group data
all.gts <- all.gts[,K.groups$id,with=FALSE]

#####################################
### get counts of alt and ref for each K group at each SNP
#####################################

### for one snp, get proportion each admixture group
ad.levs <- levels(K.groups$K.name)
#up.snp <- unlist(snp.dat[1,])

ad.count <- function(up.snp){
  up.r <- up.snp[up.snp==0]
  r.admix <- K.groups[K.groups$id%in%names(up.r),]
  r.admix <- table(factor(r.admix$K.name,levels=ad.levs))
  up.a <- up.snp[up.snp==1]
  a.admix <- K.groups[K.groups$id%in%names(up.a),]
  a.admix <- table(factor(a.admix$K.name,levels=ad.levs))
  up.h <- up.snp[up.snp==0.5]
  h.admix <- K.groups[K.groups$id%in%names(up.h),]
  h.admix <- table(factor(h.admix$K.name,levels=ad.levs))
  
  r.admix <- r.admix+h.admix
  a.admix <- a.admix+h.admix
  
  out.dat <- c(r.admix,a.admix)
  return(out.dat)
}

admix.count <- apply(snp.dat, 1, ad.count)
#admix.count <- apply(snp.dat[1:10000,], 1, ad.count)

#####################################
### reformat, add position dat, and output
#####################################

admix.count <- t(admix.count)
#rownames(admix.count) <- admix.pos[1:10000]
admix.pos <- paste(pos[,Chromosome], pos[,Position], sep="_")
rownames(admix.count) <- admix.pos
colnames(admix.count)[1:6] <- paste("R",colnames(admix.count)[1:6], sep="")
colnames(admix.count)[7:12] <- paste("A",colnames(admix.count)[7:12], sep="")
allele.count <- admix.count
save(allele.count, file="./data/K.group.allele.count.Rdat")


#####################################
### get allele frequencies per K group
#####################################

rc <- allele.count[,1:6]
ac <- allele.count[,7:12]
tc <- rc + ac

rf <- rc/tc
af <- ac/tc