#!/usr/bin/env bash

# Lines beginning with #SBATCH are instructions for SLURM, this line is a comment :)
# For example:

# This job requests from SLURM to allocate 1 node.
#SBATCH --nodes=1
# On that node, it will run 4 tasks, each with 1 core and 1 GB of memory.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=long
#SBATCH --time=4-00:00:00
#SBATCH --output=%A_%a.stdout
#SBATCH --array=1-23 #this is how to make an array should be 1-23, try with 4 to start

export UPVAR=$SLURM_ARRAY_TASK_ID

echo $UPVAR

cd /groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/scripts

ml r/3.5.1-foss-2018b
ml r-bundle-bioconductor/3.8-foss-2018b-r-3.5.1

Rscript --vanilla partial.mantel.climate.R $UPVAR

# after job is finished, see file my.stdout for results
