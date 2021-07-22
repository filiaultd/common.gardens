## ----setup, include=FALSE------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
up.rot <- as.numeric(args[[1]])
print(up.rot)

knitr::opts_chunk$set(echo = TRUE)
#library(dgof)
library("dgof", lib.loc="/users/daniele.filiault/Rpackages")
#library(ggplot2)
#library(gridExtra)
#library(ggpmisc)
#library(tidyr)
#library(ggpubr)
#library(cowplot)
#library(dplyr)
#library(viridis)
source("./50.compare.scan.nonscan.subsets.functions.R")

outfile <- paste("./data/53.data/s.scan.rot",up.rot,"Rdat", sep=".")

## ----get chromosome lengths----------------------------------------------

# get chromosome lengths
len <- read.table("../../001.common.reference.files/001.TAIR10.genome/TAIR10_all.fa.fai",stringsAsFactors=FALSE, nrows=7)
len <- len[1:5,]
len.cs <- cumsum(len[,2])
len.cs <- c(0,len.cs)
len.max <- max(len.cs)


## ----load AF and AFD data------------------------------------------------
### allele frequency data from script 46.allele.freq.differences.Kgroups
### this is needed for add.gwas() function in script 50.  It should have been explicitly coded as an input variable in the function, but that would now require changing a bunch of scripts that work.  If I need to revisit these, I will rewrite the function.
load("./data/pop.af.dat.Rdat")
load("./data/allele.freq.GWAS.snps.Rdata") #afg


## ----prep AFD data-------------------------------------------------------
### combine AFD datasets
colnames(pop.af.dat)[1:2] <- c("chrom", "pos")
pos <- do.call(rbind,strsplit(rownames(afg),"_"))
colnames(pos) <- c("chrom", "pos")
pos <- as.data.frame(pos)
afg <- cbind(afg, pos)
af.dat <- merge(pop.af.dat, afg, by=c("chrom", "pos"), all=TRUE)


## ----fxns for rotation and KS tests--------------------------------------

###################################
### add relative position (relpos) to any dataframe with "chrom" and "pos" columns
###################################
#up.dat is dataframe to use, len.cs is length to add to each chromosome
relpos <- function(up.dat, len.cs){
  ud.s <- split(up.dat, up.dat$chrom)
  for(chr in names(ud.s)){
    up.s <- ud.s[[chr]]
    up.s$rel.pos <- up.s$pos + len.cs[as.numeric(chr)]
    ud.s[[chr]] <- up.s
  }
  ud.s <- do.call(rbind, ud.s)
  return(ud.s)
}

###############################################
### rotates a vector by a certain number of bp
################################################
### pos is a vector of positions to rotate
### bp.slide is the number of bases to rotate
### max.bp is the maximum positions in the genome

genome.rotate <- function(pos,bp.slide, max.bp){
	pos.r <- pos+bp.slide
	pos.rr <- sapply(pos.r, function(x){
		if(x>max.bp){x <- x-max.bp}
		return(x)
	})
	return(pos.rr)
}

#####################################################
### Kolmogorov-Smirnov test between two distributions
#####################################################

 ks.test.column <- function(a.dat, b.dat, var.name){
  aval <- a.dat[,colnames(a.dat)==var.name]
  bval <- b.dat[,colnames(b.dat)==var.name]
  out.test <- suppressWarnings(ks.test(aval, bval, alternative = "greater"))
  ### use caution with suppressing warnings!  I did it here to keep my logfile on the cluster from slowing everything down. 
}

################################################
### do KS tests of observed data AF, AFD, home.beta for a set of experiments
############################################
#ss.dat <- up.ss.dat
#scan.name <- up.scan.name
ks.all <- function(ss.dat, scan.name){
  ks.out <- matrix(NA, ncol=6, nrow=4)
  exps <- unique(ss.dat$exp)

  for(up in 1:length(exps)){
    up.dat <- ss.dat[ss.dat$exp==exps[up],]
    up.af.test <- ks.test.column(a.dat=up.dat, b.dat=up.dat[up.dat$scan==TRUE,], var.name="ahome")
    up.af.sum <- unlist(up.af.test[c("statistic","p.value")])
    names(up.af.sum) <- c("statistic", "p.value")
    names(up.af.sum) <- paste("af", names(up.af.sum), sep=".")
    
    up.afd.test <- ks.test.column(a.dat=up.dat, b.dat=up.dat[up.dat$scan==TRUE,], var.name="ha.afd")
    up.afd.sum <- unlist(up.afd.test[c("statistic","p.value")])
    names(up.afd.sum) <- c("statistic", "p.value")
    names(up.afd.sum) <- paste("afd", names(up.afd.sum), sep=".")
    
    up.beta.test <- ks.test.column(a.dat=up.dat, b.dat=up.dat[up.dat$scan==TRUE,], var.name="home.beta")
    up.beta.sum <- unlist(up.beta.test[c("statistic","p.value")])
    names(up.beta.sum) <- c("statistic", "p.value")
    names(up.beta.sum) <- paste("beta", names(up.beta.sum), sep=".")
    
    out.sum <- c(up.af.sum, up.afd.sum,up.beta.sum)
    ks.out[up,] <- out.sum
    colnames(ks.out) <- names(out.sum)
  }
  ks.out <- as.data.frame(ks.out)
  ks.out$exp <- exps
  ks.out$scan <- scan.name
  return(ks.out)
}
 
#####################################
### do genome rotation for KS stats
#####################################
 
ks.genome.rotation <- function(nrotations, win.size, ss.dat, len.cs, len.max, home.beta){
  # get rotation values, set up relative positions
  rot.bp <- sample(1:len.max-1, nrotations)
  colnames(ss.dat)[1] <- "chrom"
  ss.dat.rel <- relpos(up.dat=ss.dat, len.cs=len.cs)
  
  #constuct output file
  rot.ks <- as.list(rep(NA,nrotations))
  #do rotations
  for(up in 1:nrotations){
    print(up)
    up.r <- rot.bp[[up]]
    #rotate selection scan positions
    ss.dat.up <- ss.dat.rel
    ss.dat.up$rel.pos <- genome.rotate(pos=ss.dat.rel$rel.pos, bp.slide=up.r, max.bp=len.max)
    up.rot.dat <- scan.snps.rel.pos(ss.dat=ss.dat.up, win.size=win.size, home.beta=home.beta)
    #do rotated KS test
    up.scan.ks <- ks.all(ss.dat=up.rot.dat, scan.name=up.scan.name)
    up.scan.ks$rotation <- as.character(up)
    rot.ks[[up]] <- up.scan.ks
  }
  rot.ks <- do.call(rbind, rot.ks)
  return(rot.ks)
} 
 

## ----polarize betas south and get AF and AFD bins------------------------
### GWAS run in gemma in /groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/28.BLUP.GWAS.Rmd

gwas.res.files.s <- c("./res/gemma_marginal/gemma_lmm_blup_RAT_2011.rds", "./res/gemma_marginal/gemma_lmm_blup_RAT_2012.rds", "./res/gemma_marginal/gemma_lmm_blup_ULL_2011.rds","./res/gemma_marginal/gemma_lmm_blup_ULL_2012.rds")

home.allele="ASOUTH"
away.allele="ANORTH"

s.genome.dat <- as.list(1:4)

for(up.fn in 1:4){
  up.f <- gwas.res.files.s[up.fn]
  print(up.f)    
  up.dat <- add.gwas(up.gwa.file=up.f, home.allele=home.allele, away.allele=away.allele)
  up.pheno <- get.p.name(up.gwa.file=up.f)
  up.short.name <- short.p.name(p.name=up.pheno)
  # get AFD bins
  breaks <- seq(0.5,1,0.1)
  # specify interval/bin labels
  tags <- c("[.5-.6)","[.6-.7)","[.7-.8)","[.8-.9)","[.9-1)")  # bucketing values into bins north/south
  up.dat$ahome.bins <- cut(up.dat$ahome, breaks=breaks, include.lowest=TRUE, right=FALSE, labels=tags)
  up.dat$ahome.bins <- factor(up.dat$ahome.bins, levels = c(tags, "[1]"),ordered = TRUE)
  # add fixed bins
  up.dat[up.dat$ahome==1, 35] <- "[1]"
  up.dat$exp <- up.pheno
  up.dat <- up.dat[,c("chrom", "pos","exp","home.beta","ahome", "ahome.bins", "ha.afd", "ha.bins")]
  s.genome.dat[[up.fn]] <- up.dat
}
s.genome.dat <- do.call(rbind, s.genome.dat)
s.genome.dat <- relpos(up.dat=s.genome.dat, len.cs=len.cs)



## ----ks tests selection scan data - South experiments--------------------
gwas.files <- gwas.res.files.s
win.size=10000
nrotations <- 40

scannames <- c("swedishgenomes","huber","hortonFST","hortonCLR","hortonPHS","priceBayescan","priceAFDLD")

s.sweep.rot <- as.list(1:length(scannames))

### 1. swedish genomes
up.scan <- 1
up.scan.name <- scannames[up.scan]
up.scan.file <- "./data/003.selection.scans/Swedish.genomes.paper.sweeps.simple.csv"
sg.dat <- read.csv(up.scan.file, stringsAsFactors=FALSE)
sg.dat <- sg.dat[,1:2]
colnames(sg.dat) <- c("chr", "pos")
ss.dat <- sg.dat
# get observed data
up.ss.dat <- scan.snps(ss.dat=ss.dat, win.size=win.size, home.beta=s.genome.dat)
colnames(up.ss.dat)[1] <- "chr"
obs.scan.ks <- ks.all(ss.dat=up.ss.dat, scan.name=up.scan.name)
obs.scan.ks$rotation <- "observed"
# get rotated data
rot.scan.ks <- ks.genome.rotation(nrotations=nrotations, win.size=win.size, ss.dat=ss.dat, len.cs=len.cs, len.max=len.max, home.beta=s.genome.dat)
# output data
s.sweep.rot[[up.scan]] <- rbind(obs.scan.ks, rot.scan.ks)
save(s.sweep.rot, file=outfile)

### 2. Huber scans
up.scan <- 2
up.scan.name <- scannames[up.scan]
up.scan.file <- "./data/003.selection.scans/christian.intervals.csv"
ns.dat <- read.csv(up.scan.file, stringsAsFactors=FALSE)
ss.dat <- ns.dat[,1:2]
# get observed data
up.ss.dat <- scan.snps(ss.dat=ss.dat, win.size=win.size, home.beta=s.genome.dat)
colnames(up.ss.dat)[1] <- "chr"
obs.scan.ks <- ks.all(ss.dat=up.ss.dat, scan.name=up.scan.name)
obs.scan.ks$rotation <- "observed"
# get rotated data
rot.scan.ks <- ks.genome.rotation(nrotations=nrotations, win.size=win.size, ss.dat=ss.dat, len.cs=len.cs, len.max=len.max, home.beta=s.genome.dat)
# output data
s.sweep.rot[[up.scan]] <- rbind(obs.scan.ks, rot.scan.ks)
save(s.sweep.rot, file=outfile)

### 3. Horton FST
up.scan <- 3
up.scan.name <- scannames[up.scan]
up.scan.file <- "./data/003.selection.scans/selection_scans_horton/global.all.fsts"
fst.dat <- read.table(up.scan.file, stringsAsFactors=FALSE, header=TRUE)
fst.dat <- fst.dat[order(fst.dat$fst, decreasing=TRUE),]
#take top 1%
fst.dat <- fst.dat[1:(nrow(fst.dat)*0.01),]
fst.dat <- fst.dat[,2:3]
colnames(fst.dat) <- c("chr", "pos")
ss.dat <- fst.dat
# get observed data
up.ss.dat <- scan.snps(ss.dat=ss.dat, win.size=win.size, home.beta=s.genome.dat)
colnames(up.ss.dat)[1] <- "chr"
obs.scan.ks <- ks.all(ss.dat=up.ss.dat, scan.name=up.scan.name)
obs.scan.ks$rotation <- "observed"
# get rotated data
rot.scan.ks <- ks.genome.rotation(nrotations=nrotations, win.size=win.size, ss.dat=ss.dat, len.cs=len.cs, len.max=len.max, home.beta=s.genome.dat)
# output data
s.sweep.rot[[up.scan]] <- rbind(obs.scan.ks, rot.scan.ks)
save(s.sweep.rot, file=outfile)

### 4. Horton CLR
up.scan <- 4
up.scan.name <- scannames[up.scan]
up.scan.file <- "./data/003.selection.scans/selection_scans_horton/allCLR.txt"
clr.dat <- read.table(up.scan.file, stringsAsFactors=FALSE, header=TRUE)
clr.dat <- clr.dat[order(clr.dat$meanCLR, decreasing=TRUE),]
#take top 1%
clr.dat <- clr.dat[1:(nrow(clr.dat)*0.01),]
clr.dat <- clr.dat[,1:2]
colnames(clr.dat) <- c("chr", "pos")
ss.dat <- clr.dat
# get observed data
up.ss.dat <- scan.snps(ss.dat=ss.dat, win.size=win.size, home.beta=s.genome.dat)
colnames(up.ss.dat)[1] <- "chr"
obs.scan.ks <- ks.all(ss.dat=up.ss.dat, scan.name=up.scan.name)
obs.scan.ks$rotation <- "observed"
# get rotated data
rot.scan.ks <- ks.genome.rotation(nrotations=nrotations, win.size=win.size, ss.dat=ss.dat, len.cs=len.cs, len.max=len.max, home.beta=s.genome.dat)
# output data
s.sweep.rot[[up.scan]] <- rbind(obs.scan.ks, rot.scan.ks)
save(s.sweep.rot, file=outfile)

### 5. Horton PHS
up.scan <- 5
up.scan.name <- scannames[up.scan]
up.scan.file <- "./data/003.selection.scans/selection_scans_horton/phsscores.pos.txt"
phs.dat <- read.table(up.scan.file, stringsAsFactors=FALSE, header=TRUE)
phs.dat <- phs.dat[order(phs.dat$PHS, decreasing=TRUE),]
#take top 1%
phs.dat <- phs.dat[1:(nrow(phs.dat)*0.01),]
phs.dat <- phs.dat[,1:2]
colnames(phs.dat) <- c("chr", "pos")
ss.dat <- phs.dat
# get observed data
up.ss.dat <- scan.snps(ss.dat=ss.dat, win.size=win.size, home.beta=s.genome.dat)
colnames(up.ss.dat)[1] <- "chr"
obs.scan.ks <- ks.all(ss.dat=up.ss.dat, scan.name=up.scan.name)
obs.scan.ks$rotation <- "observed"
# get rotated data
rot.scan.ks <- ks.genome.rotation(nrotations=nrotations, win.size=win.size, ss.dat=ss.dat, len.cs=len.cs, len.max=len.max, home.beta=s.genome.dat)
# output data
s.sweep.rot[[up.scan]] <- rbind(obs.scan.ks, rot.scan.ks)
save(s.sweep.rot, file=outfile)


### 6. Price Bayescan
up.scan <- 6
up.scan.name <- scannames[up.scan]
up.scan.file <- "./data/003.selection.scans/price_2020/BAYESCAN_FDR_FINAL"
pbayes.dat <- read.table(up.scan.file, stringsAsFactors=FALSE, header=TRUE)
pbayes.dat <- pbayes.dat[pbayes.dat$FDR<0.1,]
pbayes.dat <- pbayes.dat[,1:2]
colnames(pbayes.dat) <- c("chr", "pos")
ss.dat <- pbayes.dat
# get observed data
up.ss.dat <- scan.snps(ss.dat=ss.dat, win.size=win.size, home.beta=s.genome.dat)
colnames(up.ss.dat)[1] <- "chr"
obs.scan.ks <- ks.all(ss.dat=up.ss.dat, scan.name=up.scan.name)
obs.scan.ks$rotation <- "observed"
# get rotated data
rot.scan.ks <- ks.genome.rotation(nrotations=nrotations, win.size=win.size, ss.dat=ss.dat, len.cs=len.cs, len.max=len.max, home.beta=s.genome.dat)
# output data
s.sweep.rot[[up.scan]] <- rbind(obs.scan.ks, rot.scan.ks)
save(s.sweep.rot, file=outfile)

### 7. Price AFD.LD
up.scan <- 7
up.scan.name <- scannames[up.scan]
up.scan.file <- "./data/003.selection.scans/price_2020/AFD_LD_Italy_Sweden_pops"
pal.dat <- read.table(up.scan.file, stringsAsFactors=FALSE, header=TRUE)
pal.dat <- pal.dat[pal.dat$AFD>0.7 & pal.dat$LD>0.19,] #cutoff -> AFD>0.7 and LD>0.19 (95th percentiles)
pal.dat <- pal.dat[,1:2]
colnames(pal.dat) <- c("chr", "pos")
ss.dat <- pal.dat
# get observed data
up.ss.dat <- scan.snps(ss.dat=ss.dat, win.size=win.size, home.beta=s.genome.dat)
colnames(up.ss.dat)[1] <- "chr"
obs.scan.ks <- ks.all(ss.dat=up.ss.dat, scan.name=up.scan.name)
obs.scan.ks$rotation <- "observed"
# get rotated data
rot.scan.ks <- ks.genome.rotation(nrotations=nrotations, win.size=win.size, ss.dat=ss.dat, len.cs=len.cs, len.max=len.max, home.beta=s.genome.dat)
# output data
s.sweep.rot[[up.scan]] <- rbind(obs.scan.ks, rot.scan.ks)
save(s.sweep.rot, file=outfile)

### 8. output data
#save(s.sweep.rot, file="./data/53.data/s.sweep.rot.Rdat")




