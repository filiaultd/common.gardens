# module load gcc/6.2.0
# module add plink
# module add gemma

module load GCC/6.3.0-2.27
module load PLINK/1.9b_4.1-x86_64
module load GEMMA/0.97-foss-2017a

wd="/lustre/scratch/projects/field_experiments/001.common.reference.files/006.field.SNPs.Fernando"
vcf="02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.vcf"
adapt=$wd

plink --vcf $wd/$vcf --make-bed --out $wd/sweden --missing-genotype N --chr 1-5  --allow-extra-chr --double-id --keep-allele-order

##the keep-allele-order is important because it allows to maintain the genotypes in ref/alt format instead of minor/major. This will allow to interprete coefficient signs from GWAs along the the change of allele frequency of the

head -n 200 $wd/sweden.fam > $wd/keep.txt
head $wd/keep.txt
tail $wd/keep.txt

mkdir $adapt
recodeAT=$wd/recodeAT/sweden_200_MAF10
outbed=$wd/bed/sweden_200_MAF10
bimbam=$wd/bimbam/sweden_200_MAF10
mkdir $wd/recodeAT
mkdir $wd/bed
mkdir $wd/bimbam
bimbam2=$wd/bimbam/sweden_200p
mkdir $wd/bimbam2

##remove accessions we don't need (keep the first 200 genotypes)

plink --bfile $wd/sweden --make-bed --allow-no-sex --keep $wd/keep.txt --missing-code N,-9,0,NA --no-pheno -out $wd/sweden_200 --keep-allele-order

##filter out rare SNPs

plink --bfile $wd/sweden_200 --allow-no-sex --maf 0.1 --missing-code N,-9,0,NA --make-bed --out $wd/sweden_200_MAF10_filt1 --keep-allele-order

##remove poorly genotyped SNPs

plink --bfile $wd/sweden_200_MAF10_filt1 --allow-no-sex --geno 0.05 --missing-code N,-9,0,NA --make-bed --out $wd/sweden_200_MAF10_filt2 --keep-allele-order

##filter out rare SNPs again write bed (for bayesR) and recode (for gemma)

plink --bfile $wd/sweden_200_MAF10_filt2 --allow-no-sex -maf 0.10 --missing-code N,-9,0,NA --no-pheno --make-bed --set-missing-var-ids @_# --out $outbed --keep-allele-order

##make a freq file

plink --bfile $outbed --freq --out $outbed --keep-allele-order

plink --bfile $wd/sweden_200_MAF10_filt2 --allow-no-sex -maf 0.10 --missing-code N,-9,0,NA --no-pheno -recode A-transpose --set-missing-var-ids @_# --keep-allele-order --out $recodeAT 

plink --bfile $wd/sweden_200_MAF10_filt2 --allow-no-sex -maf 0.10 --missing-code N,-9,0,NA --no-pheno -recode bimbam --keep-allele-order --out $bimbam2

##to make it into a mean genotype file, it looks like I just need to remove the centimorgan column

##make an accession list to make sure we get the order right
head $bimbam2.recode.geno.txt -n 3 | tail -n 1 | sed -e 's/IND,//g' > $bimbam2.acclist.txt

cp $bimbam2* /home/benjamin/data/sweden/adaptation/GWA/snps

##to make it into a mean genotype file, it looks like I just need to remove the centimorgan column

tail -n +2 $recodeAT.traw | awk -F "\t" '{print $1"_"$4"\t"$4"\t"$1}'| sed 's/\t/,/g' > $bimbam.bimbam.map
cat $bimbam.bimbam.map | awk -F "," '{print $1}' > $wd/snps_list.temp
tail -n +2 $recodeAT.traw | cut -f5- | sed 's/\t/,/g' >  $wd/bimbam.temp

paste -d "," $wd/snps_list.temp $wd/bimbam.temp > $bimbam.bimbam.geno
rm $wd/snps_list.temp
rm $wd/bimbam.temp

head -n 1 $recodeAT.traw | cut -f7- | sed 's/[_][0-9]*/\n/g' > $bimbam.bimbam.acclist

##check the results
head $bimbam.bimbam.geno
head $bimbam.bimbam.map
head $bimbam.bimbam.acclist

## how many genotypes in the snps file:
head -n 1 $bimbam.bimbam.geno | awk -F "," '{print NF-3}' ##should say 200
## how many accessions in the list of accessions: should say the same + 1 because there is an empty line
wc -l $bimbam.bimbam.acclist ## 201








## over to the microbiota analysis folder.

#cp $bimbam.bimbam.map $adapt
#cp $bimbam.bimbam.geno $adapt
#cp $bimbam.bimbam.acclist $adapt


#cp $outbed* $adapt



# tail -n +2 $recodeAT.traw | awk -F "\t" '{print $1"_"$4"\t"$4"\t"$1}'| sed 's/\t/,/g' > $bimbam.bimbam.map
# cat $bimbam.bimbam.map | awk -F "," '{print $1}' > $wd/snps_list.temp
# tail -n +2 $recodeAT.traw | cut -f5- | sed 's/\t/,/g' >  $wd/bimbam.temp

# paste -d "," $wd/snps_list.temp $wd/bimbam.temp > $bimbam.bimbam.geno
# rm $wd/snps_list.temp
# rm $wd/bimbam.temp

# head -n 1 $recodeAT.traw | cut -f7- | sed 's/[_][0-9]*/\n/g' > $bimbam.bimbam.acclist

# ##check the results
# head $bimbam.bimbam.geno
# head $bimbam.bimbam.map
# head $bimbam.bimbam.acclist

# ## how many genotypes in the snps file:
# head -n 1 $bimbam.bimbam.geno | awk -F "," '{print NF-3}' ##should say 200
# ## how many accessions in the list of accessions: should say the same + 1 because there is an empty line
# wc -l $bimbam.bimbam.acclist ## 201

# ## over to the adaptation analysis folder.

# cp $bimbam.bimbam.map $adapt
# cp $bimbam.bimbam.geno $adapt
# cp $bimbam.bimbam.acclist $adapt

# cp $outbed* $adapt
# ##in the adaptation folder, make a backup fam file
# cp $adapt/sweden.fam $adapt/sweden.fam.backup

# ## make kinship matrices with all 200 genotypes
# for i in $(seq 1 200); do echo 1; done >  $wd/fakephen.txt

# cd $wd
# gemma -gk 1 -g $bimbam.bimbam.geno -p $wd/fakephen.txt -o K_200_MAF10
# gemma -gk 2 -g $bimbam.bimbam.geno -p $wd/fakephen.txt -o K_200_MAF10

# ##move the kinship matrix to the adaptation analysis folder
# cp $wd/output/K_200_MAF10* $adapt
