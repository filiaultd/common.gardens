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
### unfortunately, this also impacts column names, so need to recover those from .csv file

###############################
### prep SNP and admixture data
###############################

### snp data

all.gts <- fread(all.gts.file,stringsAsFactors=FALSE,colClasses="numeric",header=TRUE)
hs <- read.csv(header.file, nrows=1)