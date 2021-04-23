library(plyr)
library(mixOmics)
library(lme4)
library(igraph)
library(Hmisc)
library(vegan)



wd="/group/bergelson-lab/bbrachi_data/sweden/adaptation/GWA"
prefix="/group/bergelson-lab/bbrachi_data/sweden/microbiota/GWA/snps/sweden"


name_file= paste(getwd(),"/GWA/phen/phen_CG.tsv", sep="")
write.table(colnames(phen)[-1],name_file, col.names=F, row.names=F, sep="\t", quote=F)
pheno_file= paste(getwd(), "/GWA/phen/blups_CG.tsv", sep="")
write.table(phen[,-1],pheno_file, col.names=F, row.names=F, sep="\t", quote=F)

names=colnames(phen[,-1])


cat("", file=paste(wd, "/submit/submit_all_hh_GWAs.sh",sep=""), append=F)
for(p in names){
    phen=p
    n=match(p, names)
    cat(paste("
#!/bin/bash
#PBS -N lmm_GWA_", phen,"
#PBS -d ", wd,"
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=12gb
#PBS -q batch
#PBS -l walltime=3:00:00
#PBS -o ", wd, "/logs/GWA_lr_hh_", phen, ".log
#PBS -e ", wd, "/logs/GWA_lr_hh_", phen, ".err
module add gemma
cd ", wd, "
gemma -gk 1 -g ", prefix, ".bimbam.geno -p ", phen_file, " -n ", n, " -o K_", phen, "
gemma -lmm 2 -g ", prefix, ".bimbam.geno -a ", prefix, ".bimbam.map -p ", phen_file, " -n ", n ," -k ./output/K_", phen, ".cXX.txt -o out_", phen, "
",sep=""), file=paste(wd, "/submit/submit_GWA_", phen, ".sh", sep=""))
    cat(paste("qsub ", wd, "/submit/submit_GWA_", phen, ".sh\n", sep=""), file=paste(wd, "/submit/submit_all_lr_hh_GWAs.sh", sep=""), append=T)
}

