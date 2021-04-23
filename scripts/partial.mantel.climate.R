
### script to run Hancock-style partial Mantel tests for v2.climate variables and PCs of climate variables
### with 200 swedish lines in experiment
### method developed in 37.climate.associations.mantel.Rmd

args = commandArgs(trailingOnly=TRUE)
up.col <- as.numeric(args[[1]])+3
print(up.col)

library(ecodist)
source("/groups/nordborg/projects/field_experiments/adaptation_sweden/select.and.resequence/002.scripts/00.allele.freq.change.fxns.R")
setwd("/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/")

## 1. load data
# genotypes
# loading SNPs using a script from 00.allele.freq.change.fxns.R
all.gts.file <- "/groups/nordborg/projects/field_experiments/001.common.reference.files/006.field.SNPs.Fernando/02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.csv"
## this is the data set used for genotyping these samples (from Fernando)
min.MAF <- 0.1
max.NAF <- 0.05 ## remove if NA rate is higher than this
### set these so agree with Ben's GWAS filtering so our SNPs overlap
gts <- prep.SAF(all.gts.file=all.gts.file, min.MAF=min.MAF, max.NAF=max.NAF)
# reads in file, recodes to NA, 0 (ref), 0,5 (het), 1 (alt)
# Angela's original method didn't do minor allele frequency filtering, but I am still including it

#experimental lines
lines <- read.table("./data/acc_list.txt",stringsAsFactors=FALSE, fill=TRUE,sep="\t", header=TRUE, nrows=200)
exp.lines <- lines$lines
#exp.lines%in%colnames(gts) looks OK

# K.matrix
# as calculated in plink from these data
K <- read.table("./data/K.matrix.200Swedes.txt")
## need accession names
pref="./GWA/snps/sweden_200_MAF10.bimbam"
acclist=scan(paste(pref, ".acclist", sep=""), sep=",")
rownames(K) <- acclist
colnames(K) <- acclist
K <- as.matrix(K)
K.dist <- dist(K)
## distance matrix from K

# climate variables
# bioclim v2, BIO1 to BIO19
bioc <- read.table("./data/bioclim.v2.200.experimental.lines.txt", header=TRUE, stringsAsFactors=FALSE)

# PCs of this data - use first 4 in analysis
b.pcs <- read.table("./data/experimental.pcas.txt", header=TRUE, stringsAsFactors=FALSE)
b.pcs <- b.pcs[rownames(b.pcs)%in%exp.lines,]

bioc <- merge(bioc, b.pcs, by.x="lines", by.y="id", all=TRUE)

print("files loaded")

## 2. partial mantel function
### snp.no is row number of up.gts (genotypes)
### bc.no is the column number of the bioclim variable

pm.test <- function(snp.no, bc.no, up.gts){
  # get distance matrices
  gt.up <- t(up.gts[snp.no,])
  up.pos <- gt.up[1:2,]
  gt.up <- gt.up[-(1:2),]
  gt.d <- dist(gt.up)
  
  bc.up <- bioc[,bc.no]
  #names(bc.up)==names(gt.up)  looks ok.
  names(bc.up) <- bioc$lines
  bc.d <- dist(bc.up)
  
  #rownames(K)==names(bc.up) also looks OK
  
  # run test
  test.mt <- mantel(gt.d ~ bc.d + K.dist, mrank=TRUE)
  return(test.mt)
}

## 3. run for BIOCLIM data
## columns is up.col (command line arguement [[1]])
## run by chromosome in case walltime is insufficient

for(up.chr in 1:5){
  up.gts <- gts[gts$Chromosome==up.chr,]
  up.var <- colnames(bioc)[up.col]
  out.file <- paste("./data/001.partial.mantel.climate/",up.var,".CHR",up.chr,".partial.mantel.output.Rdat", sep="")
  test.m <- matrix(NA, nrow=nrow(up.gts), ncol=6)
  #for(up.snp in 1:100){
  for(up.snp in 1:nrow(up.gts)){
    up.mt <- pm.test(snp.no=up.snp, bc.no=up.col, up.gts=up.gts)
    test.m[up.snp,] <- up.mt[1:6]
  }
  save(test.m, file=out.file)
}



