#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=00:15:00
#SBATCH --output=%A_%a.blup.limix.stdout
#SBATCH --array=1-8

# set environment #
ml anaconda3/2019.03
source /software/2020/software/anaconda3/2019.03/etc/profile.d/conda.sh
conda activate ~/.conda/envs/limix

# define variables #
export i=$SLURM_ARRAY_TASK_ID
N=`expr $i - 1`
PHENO_FILE="/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/fitness.BLUPS.with.Kgroup.txt"
GENO_FILE="/groups/nordborg/projects/field_experiments/001.common.reference.files/006.field.SNPs.Fernando/imputed/02_2.3M_200Swedes.biallelic.imputed.filtered.bed"
K_FILE="/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/K.matrix.200Swedes.labels.txt"
ALLPHENOS=(RAM_2011 ADA_2011 RAT_2011 ULL_2011 RAM_2012 ADA_2012 RAT_2012 ULL_2012)
UP_PHENO=${ALLPHENOS[$N]}
MAR_PY="/users/daniele.filiault/limix/003.Daniele.scripts/limix.marginal.py"

echo $N
echo $UP_PHENO


# run marginal GWAS #

python $MAR_PY -pf $PHENO_FILE -gf $GENO_FILE -kf $K_FILE -pn $UP_PHENO
