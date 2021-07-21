#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=long
#SBATCH --time=2-00:00:00
#SBATCH --output=%A_%a.gea.south.out
#SBATCH --array=1-25

export UPVAR=$SLURM_ARRAY_TASK_ID

cd /groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens

ml r/3.5.1-foss-2018b
ml r-bundle-bioconductor/3.8-foss-2018b-r-3.5.1

Rscript --vanilla 53.rotate.gea.south.R  $UPVAR

