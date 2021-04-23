
setwd("/group/bergelson-lab/project_sweden/sweden/adaptation")

prefix="/group/bergelson-lab/project_sweden/sweden/adaptation/GWA/snps/"
GWAfolder="/group/bergelson-lab/project_sweden/sweden/adaptation/GWA/"

phen=read.table("./res/means_NvsS.txt", sep="\t",h=T)
##set the data type (blups, means...)
dt="means_NvsS"
system(paste("mkdir ./GWA/", dt, sep=""))
system(paste("mkdir ./GWA/", dt, "/phen",sep=""))
system(paste("mkdir ./GWA/", dt, "/submit",sep=""))
system(paste("mkdir ./GWA/", dt, "/logs",sep=""))

##order the phenotypes to match genotypes

acc=read.table(paste(prefix, "/sweden_200_MAF10.bimbam.acclist", sep=""), h=F)
colnames(acc)="id"
p=merge(acc, phen, by="id", sort=F)

##restrict to traits to interests

write.table(data.frame(p)[,-1], paste("./GWA/", dt,"/phen/GWA_phen_mlmm.txt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")

system(paste("rm ./GWA/", dt,"/submit/submit_all_mlmm_fit_NvsS.sh", sep=""))

ns=" 1 2"
cat(paste("#!/bin/bash
#PBS -N GWAmlmm_", dt,"
#PBS -d ",GWAfolder,"/", dt,"
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=6gb
#PBS -l walltime=2:00:00
#PBS -o ",GWAfolder,"/", dt, "/log.txt
#PBS -e ",GWAfolder,"/", dt, "/error.txt
module load gcc/6.2.0
module add gemma
cd ",GWAfolder,"/",dt,"
gemma -g ", prefix, "sweden_200_MAF10.bimbam.geno -gk 1 -p ",GWAfolder, "/", dt,"/phen/GWA_phen_mlmm.txt -o K_sweden_fit_NvsS","
gemma -lmm -n ", ns, " -g ", prefix, "sweden_200_MAF10.bimbam.geno -k ",GWAfolder,"/", dt, "/output/K_sweden_fit_NvsS.cXX.txt -p ",GWAfolder,"/", dt,"/phen/GWA_phen_mlmm.txt -a ", prefix, "sweden_200_MAF10.bimbam.map -o mlmm_fit_NvsS","
gzip  ./output/mlmm_fit_NvsS_.*.txt
", sep="", collapse=" "),file=paste("./GWA/", dt,"/submit/mlmm_fit_NvsS.qsub", sep=""))

