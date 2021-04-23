#!/bin/bash

#PBS -P aquilegia
#PBS -N admix.prop
#PBS -j oe
#PBS -o admix.prop.log
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=100gb

module load R

Rscript /lustre/scratch/projects/field_experiments/adaptation.sweden/common.gardens/scripts/admix.prop.by.snp.R

