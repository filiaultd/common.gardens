## run interactively on cbe
## uses files from Ben's prep_SNPs.sh file that was meant to prep SNPs for GWAS


srun --qos=short --partition=c --cpus-per-task=1 --mem=10gb --pty bash

ml plink/1.9b_6.10-x86_64

#plink --bfile ./sweden_200_MAF10_filt2 --freq

#plink --bfile $wd/sweden_200_MAF10_filt2 --allow-no-sex -maf 0.10 --missing-code N,-9,0,NA --no-pheno -recode bimbam --keep-allele-order --out $bimbam2

#plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 200 --out 
#plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 400 --out plink.300kb
## for some reason???

#plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 600 --out plink.600kb

## DLF 04Nov20
## running with some additional parameter combinations to try to improve block definition
## first is exact settings Ben uses in microbiome paper

plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 200 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --out plink.ben
## this reduces the high end of the "strong LD" 
## "Two variants are normally considered by this procedure to be in "strong LD" if the bottom of the 90% D-prime confidence interval is greater than 0.70, and the top of the confidence interval is at least 0.98."

## I would like to try a few more tweeks, as I think this isn't enough to define good blocks here.

## 1. increase blocks max to 400kb (this seemed to help a bit before)
plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 400 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --out plink.ben.400kb

## 2. "By default, this procedure treats confidence interval tops smaller than 0.90 as strong evidence for historical recombination; use --blocks-recomb-highci to adjust this threshold." So if I am understanding correctly, we want to decrease this a bit to reduce "recombination" detection (our blocks are apparently too short).
plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 400 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --blocks-recomb-highci 0.8 --out plink.ben.400kb.recomb8

## "Normally, the number of "strong LD" pairs within a haploblock must be more than 0.95 times the total number of informative pairs (i.e. either "strong LD" or 'recombination'). This threshold can be adjusted with --blocks-inform-frac." For our case, I think that this should be lowered appreciably.  Let's try a few iterations

plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 200 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --blocks-inform-frac 0.9 --out plink.ben.inform9

plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 200 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --blocks-inform-frac 0.8 --out plink.ben.inform8

plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 200 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --blocks-inform-frac 0.5 --out plink.ben.inform5

## From playing with these results in 24.haplotype.block.thin.Rmd, it seems like max-kb should be increased

plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 400 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --blocks-inform-frac 0.9 --out plink.ben400.inform9

plink --bfile ./sweden_200_MAF10_filt2 --blocks 'no-pheno-req' --blocks-max-kb 400 --blocks-strong-lowci 0.7005 --blocks-strong-highci 0.9005 --blocks-inform-frac 0.8 --out plink.ben400.inform8

