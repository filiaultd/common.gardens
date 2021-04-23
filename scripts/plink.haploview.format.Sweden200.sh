## run interactively on cbe
## uses files from Ben's prep_SNPs.sh file that was meant to prep SNPs for GWAS
## convert format so can be used in Haploview

srun --qos=short --partition=c --cpus-per-task=1 --mem=10gb --pty bash

ml plink/1.9b_6.10-x86_64

plink --bfile ./sweden_200_MAF10_filt2 --recode --out sweden_200_MAF10_filt2_HVformat


