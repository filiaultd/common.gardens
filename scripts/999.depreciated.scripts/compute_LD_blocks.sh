#!/bin/bash

#PBS -d /group/bergelson-lab/project_sweden/sweden/adaptation/
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=12gb
#PBS -l walltime=96:00:00
#PBS -o /group/bergelson-lab/project_sweden/sweden/adaptation/logs/LD_blocks.log
#PBS -e /group/bergelson-lab/project_sweden/sweden/adaptation/logs/LD_blocks.err
#PBS -M benjamin.brachi@inra.fr
#PBS -m e

module load gcc/6.2.0
module add plink
module add gemma

wd="/group/bergelson-lab/project_sweden/sweden/genomes/002.Swedes220.SNPs"
vcf="Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.vcf"
adaptation="/group/bergelson-lab/project_sweden/sweden/adaptation/GWA/snps"

##compute the haplotype blocks

plink --bfile $wd/sweden_200_MAF10_filt2 --allow-no-sex -maf 0.10 --missing-code N,-9,0,NA --no-pheno -blocks 'no-pheno-req' --blocks-max-kb 50 --blocks-min-maf 0.10 --out $adaptation/sweden_200_MAF10

