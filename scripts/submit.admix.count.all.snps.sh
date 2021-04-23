#!/usr/bin/env bash

# Lines beginning with #SBATCH are instructions for SLURM, this line is a comment :)
# For example:

# This job requests from SLURM to allocate 1 node.
#SBATCH --nodes=1
# On that node, it will run 4 tasks, each with 1 core and 1 GB of memory.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=0-08:00:00
# And it will place the output of the commands into my.stdout file
#SBATCH --output=admix.count.output

cd /groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/scripts

ml r/3.5.1-foss-2018b
ml r-bundle-bioconductor/3.8-foss-2018b-r-3.5.1

Rscript admix.count.all.snps.R

# after job is finished, see file my.stdout for results
