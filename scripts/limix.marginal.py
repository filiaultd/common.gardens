#!/usr/bin/env python
# coding: utf-8

# # Introduction
# Script to run marginal GWAS in limix with field experiment SNPs 
# DLF 11 June 20
# 

# ## Set up environment

# In[103]:


from limix.qtl import scan
import pandas as pd
import numpy as np
import h5py
from bisect import bisect
from limix import plot
import random
from sklearn.decomposition import PCA
import argparse
import os
from statsmodels.stats import multitest
import math
import limix


# ## Settings

# In[104]:


# phenotype
#pheno_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/fitness.BLUPS.with.Kgroup.txt'
#pheno_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/all.bioclim.v2.PC.txt'
pheno_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/slopeNS.blups.txt'
# genotype
geno_file = '/groups/nordborg/projects/field_experiments/001.common.reference.files/006.field.SNPs.Fernando/imputed/02_2.3M_200Swedes.biallelic.imputed.filtered.bed'

# kinship matrix
kin_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/K.matrix.200Swedes.labels.txt'
# minor allele frequency threshold, not set in this script.  This was set when generating the input bim file
#MAF_thrs = 0.1 
# output results
pheno_name = "mu2.sl"
output_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/limix/' + pheno_name +'.limix.results.csv'


# ## Read phenotypes

# Read in phenotype file and subset by phenotype of interest

# In[105]:


phenoAll = pd.read_csv(pheno_file, index_col = 1, sep=" ")
phenoAll.index = phenoAll.index.map(str) ## need to make index type match that of genotype index so can check if in same order

pheno = phenoAll[['id',pheno_name]]
# get indices of NA values to filter other data later
ind = np.where(pheno[pheno_name].isnull())

# remove NA values
pheno = pheno[np.isfinite(pheno)]
pheno = pheno.dropna()
pheno_ids = pheno[['id']]
pheno = pheno[[pheno_name]]
#print(pheno.shape)
#print(pheno.head(n=5))

#print(pheno.isnull().any().any()) #checks for NAs
Y=pheno.to_numpy()
print(Y.shape)


# ## Read genotypes

# In[106]:


import limix
(bim, fam, bed) = limix.io.plink.read(geno_file, verbose=False)


# In[107]:


print(bim.head())


# In[108]:


print(fam.head())


# In[109]:


print(bed.compute())
print(type(bed))


# In[110]:


G=bed.compute()
G=G.transpose()
#remove NA accessions noted with ind giving rownumber
G=np.delete(G, ind, axis=0)

print(G.shape)
print(G[1:10])
print(type(G))


# ## Read K matrix

# In[111]:


K = pd.read_csv(kin_file, index_col = 0, sep=" ")
K.index = K.index.map(str)
Kp = K
#print(K.head)
#print(type(K))
K = K.to_numpy()

#remove any NA accessions with index
#note, this only works because I already know that all these matrices are in the same order
#if writing for different data or troubleshooting, might need to address this
K=np.delete(K, ind, axis=0)
K=np.delete(K, ind, axis=1)
print(type(K))
print(K.shape)
print(K)


# In[112]:


# some random idea about calculating K matrices
# can ignore for now
#from numpy import dot
#K_all = dot(bed, bed.T)


# ## Make sure all datasets in same accession order

# In[113]:


# also remove NA row from fam, which gives information about genotypes
fam.drop(fam.index[ind], inplace=True)
# and from Kp, which gives K accession info
Kp.drop(Kp.index[ind], inplace=True)
Kp.drop(Kp.columns[ind], axis=1, inplace=True)


GP_check = fam.index == pheno_ids.index
PK_check = pheno_ids.index == Kp.index
if np.count_nonzero(GP_check == False) != 0 and np.count_nonzero(PK_check == False) != 0 : print("CAUTION: Not all data files are in the same order!")
if np.count_nonzero(GP_check == False) == 0 and np.count_nonzero(PK_check == False) == 0 : print("All input files in same order")


# ## Run one marginal GWAS

# In[114]:


r = scan(G=G, Y=Y, K = K, lik = 'normal', M = None, verbose = True)


# ## Output results

# In[115]:


chrom = bim[['chrom']]
pos = bim[['pos']]
#extract pvals
pvalues = r.stats.pv20.tolist()
#extract effect sizes
effsizes = r.effsizes['h2']['effsize'][r.effsizes['h2']['effect_type'] == 'candidate'].tolist()

#gwas_results = pd.DataFrame(list(zip(chrom, pos, pvalues, effsizes)), columns = ['chr', 'pos', 'pvalue', 'GVE'])
gwas_results = np.c_[chrom, pos, pvalues, effsizes]
gwas_results = pd.DataFrame(data=gwas_results, index=None, columns=["chrom", "pos", "pv", "GVE"])
#print(gwas_results.head)
print(type(gwas_results))
print(gwas_results[1:10])
gwas_results.to_csv(output_file, index = False)


# ## Manhattan plot

# In[116]:


gwas_results.dtypes
#df.astype('int64').dtypes
gwas_results["chrom"] = gwas_results['chrom'].astype('int')
gwas_results["pos"] = gwas_results['pos'].astype('int')
gwas_results["pv"] = gwas_results['pv'].astype('float')
gwas_results["GVE"] = gwas_results['GVE'].astype('float')
gwas_results.dtypes


# In[117]:


Bonferroni = multitest.multipletests(pvalues, alpha = 0.05, method = 'fdr_bh')[3]
plot.manhattan(gwas_results)
plt = plot.get_pyplot() 
_ = plt.axhline(-math.log(Bonferroni, 10), color='red')  
plt.savefig(f'{output_file}.manhattan.png')
plt.close()


# ## QQplot

# In[118]:


plot.qqplot(gwas_results.pv)
plt = plot.get_pyplot()
plt.savefig(f'{output_file}.QQplot.png')
plt.close()


# ## Calculate heritability

# In[119]:


from numpy.random import RandomState
from limix.her import estimate

print(estimate(y=Y, K=K, lik="normal", verbose=False))  


# Heritabilities for fitness phenotypes.  
# Did this manually, but could for-loop it when writing a proper script
# 

# ADA_2011 0.4928714721332532  
# RAM_2011 0.6061879022675384  
# ULL_2011 0.8039236546423144  
# RAT_2011 0.6902180797117754  
# ADA_2012 0.6399118475365775  
# RAM_2012 0.6529507897596929  
# ULL_2012 0.7554280346125117  
# RAT_2012 0.2978431659879658  

# Heritabilities for climate phenotypes  
# PC1 0.9999999999999952  
# PC2 0.9999999999999997   
# PC3 0.9999999999999997   
# PC4 0.9999999999999999  
# BIO12 FAIL!!  

# Heritabilities for slope phenotypes  
# ar1.sl 0.5792442458886554  
# ar2.sl 0.38996558249896895  
# au1.sl 0.4424553627725288  
# au2.sl 0.5373931581307294  
# mr1.sl 0.573535965255268  
# mr2.sl 0.36429838945393733  
# mu1.sl 0.585551103248642  
# mu2.sl 0.6234137171713826

# In[ ]:




