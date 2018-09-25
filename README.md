# UKB commands
## This repository is for running GWAS in the UKB

### Step 1: Dowloading files from the UKB

Ensure the key is downloaded in the file and saved as .ukbkey

```bash
for i in {1..22}; do ./ukbgene imp -c$i; done  #this downloads the BGEN files

for i in {1..22}; do ./ukbgene imp -m -c$i; done

```



### Step 2: Converting BGEN to plink files

```bash
for i in {1..21}; do ./plink2 --bgen ukb_imp_chr${i}_v3.bgen --sample ukb20904_imp_chr${i}_v3_s487334.sample --make-bed -out ukbchr${i} --maf 0.01 --geno 0.05 --threads 10 --hwe 0.000001 --mind 0.05; done

```

### Step 3: Creating the required phenotype and covariate files in R

Please see generic script below. Specific scripts are attached.

```R
# R program ukb22918.tab created 2018-09-17 by ukb2r.cpp Mar 14 2018 14:22:05
library(data.table)
library(plyr)


bd <- fread("~/UKB_v2/ukb22918.tab", header=TRUE, sep="\t")
bd$f.53.0.0 <- as.Date(bd$f.53.0.0)
bd$f.53.1.0 <- as.Date(bd$f.53.1.0)
bd$f.53.2.0 <- as.Date(bd$f.53.2.0)
bd$f.20400.0.0 <- as.Date(bd$f.20400.0.0) # date of cancer diagnosis
bd$f.40005.0.0 <- as.Date(bd$f.40005.0.0)
bd$f.40005.1.0 <- as.Date(bd$f.40005.1.0)
bd$f.40005.2.0 <- as.Date(bd$f.40005.2.0)
bd$f.40005.3.0 <- as.Date(bd$f.40005.3.0)
bd$f.40005.4.0 <- as.Date(bd$f.40005.4.0)
bd$f.40005.5.0 <- as.Date(bd$f.40005.5.0)
bd$f.40005.6.0 <- as.Date(bd$f.40005.6.0)
bd$f.40005.7.0 <- as.Date(bd$f.40005.7.0)
bd$f.40005.8.0 <- as.Date(bd$f.40005.8.0)
bd$f.40005.9.0 <- as.Date(bd$f.40005.9.0)
bd$f.40005.10.0 <- as.Date(bd$f.40005.10.0)
bd$f.40005.11.0 <- as.Date(bd$f.40005.11.0)
bd$f.40005.12.0 <- as.Date(bd$f.40005.12.0)
bd$f.40005.13.0 <- as.Date(bd$f.40005.13.0)
bd$f.40005.14.0 <- as.Date(bd$f.40005.14.0)
bd$f.40005.15.0 <- as.Date(bd$f.40005.15.0)
bd$f.40005.16.0 <- as.Date(bd$f.40005.16.0)
bd$f.40005.17.0 <- as.Date(bd$f.40005.17.0)
bd$f.40005.18.0 <- as.Date(bd$f.40005.18.0)
bd$f.40005.19.0 <- as.Date(bd$f.40005.19.0)
bd$f.40005.20.0 <- as.Date(bd$f.40005.20.0)
bd$f.40005.21.0 <- as.Date(bd$f.40005.21.0)
bd$f.40005.22.0 <- as.Date(bd$f.40005.22.0)
bd$f.40005.23.0 <- as.Date(bd$f.40005.23.0)
bd$f.40005.24.0 <- as.Date(bd$f.40005.24.0)
bd$f.40005.25.0 <- as.Date(bd$f.40005.25.0)
bd$f.40005.26.0 <- as.Date(bd$f.40005.26.0)
bd$f.40005.27.0 <- as.Date(bd$f.40005.27.0)
bd$f.40005.28.0 <- as.Date(bd$f.40005.28.0)
bd$f.40005.29.0 <- as.Date(bd$f.40005.29.0)
bd$f.40005.30.0 <- as.Date(bd$f.40005.30.0)
bd$f.40005.31.0 <- as.Date(bd$f.40005.31.0)
bd$f.42006.0.0 <- as.Date(bd$f.42006.0.0) # first trike outcome
bd$f.42008.0.0 <- as.Date(bd$f.42008.0.0) # ischemic stroke outcome
bd$f.42010.0.0 <- as.Date(bd$f.42010.0.0) # introcerebral heamorrhage
bd$f.42012.0.0 <- as.Date(bd$f.42012.0.0) # subarachnoid heamorrhage

#How to create a file for analysis
pheno  = bd[,c("f.eid", "f.25741.2.0", "f.31.0.0", "f.22009.0.1",
                    "f.22009.0.2", "f.22009.0.3", "f.22009.0.4", "f.22009.0.5",
                    "f.22009.0.6", "f.22009.0.7", "f.22009.0.8", "f.22009.0.9",
                    "f.22009.0.10", "f.22009.0.11", "f.22009.0.12", "f.22009.0.13",
                    "f.22009.0.14", "f.22009.0.15", "f.22009.0.16", "f.22009.0.17",
                    "f.22009.0.18", "f.22009.0.19", "f.22009.0.20", "f.22009.0.21",
                    "f.22009.0.22", "f.22009.0.23", "f.22009.0.24", "f.22009.0.25",
                    "f.22009.0.26", "f.22009.0.27", "f.22009.0.28", "f.22009.0.29",
                    "f.22009.0.30", "f.22009.0.31", "f.22009.0.32", "f.22009.0.33",
                    "f.22009.0.34", "f.22009.0.35", "f.22009.0.36", "f.22009.0.37",
                    "f.22009.0.38", "f.22009.0.39", "f.22009.0.40", "f.34.0.0", 
                    "f.22006.0.0","f.22000.0.0", "f.22007.0.0", "f.22027.0.0", 
                    "f.22001.0.0", "f.22021.0.0")]

pheno2 = pheno[!is.na(pheno$f.25741.2.0),] #remove items that are na in the pheno
pheno2 = pheno2[!is.na(pheno2$f.22006.0.0),] #remove non-europeans
pheno2$checksex = ifelse(pheno2_2$f.31.0.0 == pheno2$f.22001.0.0, "correct", "incorrect")
pheno2 = subset(pheno2, pheno2 == "correct") # remove sex mismatches
pheno2 = pheno2[is.na(pheno2$f.22027.0.0),] #remove excessive heterozygosity

#The dreaded removal of related individuals. Not needed if using BOLT-LMM.
rel = fread("ukb20904_rel_s488302.dat")
rel2 = subset(rel, Kinship >  0.0884) #equivalent to 3rd degree relatives
alpha = count(rel2, vars = "ID1")
setnames(alpha, 2, "ID1_freq")
rel2 = merge(rel2, alpha)
one = subset(rel2, ID1_freq > 1)
pheno2 = pheno2[!(pheno2$f.eid %in% one$ID1),]
two = subset(rel2, ID1_freq < 2)
pheno2 = pheno2[!(pheno2$f.eid %in% two$ID2),]



headmotion_pheno = pheno2[,c("f.eid", "f.eid", "f.25741.2.0")]
headmotioncovar = pheno2[,c("f.eid", "f.eid", "f.22009.0.1",
"f.22009.0.2", "f.22009.0.3", "f.22009.0.4", "f.22009.0.5",
"f.22009.0.6", "f.22009.0.7", "f.22009.0.8", "f.22009.0.9",
"f.22009.0.10", "f.22009.0.11", "f.22009.0.12", "f.22009.0.13",
"f.22009.0.14", "f.22009.0.15", "f.22009.0.16", "f.22009.0.17",
"f.22009.0.18", "f.22009.0.19", "f.22009.0.20", "f.22009.0.21",
"f.22009.0.22", "f.22009.0.23", "f.22009.0.24", "f.22009.0.25",
"f.22009.0.26", "f.22009.0.27", "f.22009.0.28", "f.22009.0.29",
"f.22009.0.30", "f.22009.0.31", "f.22009.0.32", "f.22009.0.33",
"f.22009.0.34", "f.22009.0.35", "f.22009.0.36", "f.22009.0.37",
"f.22009.0.38", "f.22009.0.39", "f.22009.0.40", "f.22001.0.0",
"f.34.0.0", "f.22000.0.0")]

write.table(headmotion_pheno, file = "headmotionpheno.txt", row.names = F, col.names = F, quote = F)

write.table(headmotioncovar, file = "headmotioncovar.txt", row.names = F, col.names = F, quote = F)

## f.25741.2.0 is the head motion variable. 
##22009 - PCs; 22006 - Ethnic grouping

```



### Step 4: Running GWAS in plink

```bash
{1..21}; do ./plink2 --bfile ukbchr$i --linear --pheno phenotypefile --covar covariatefile --out outfilename$i --threads 20; done

```



### Step 5: Combining the files and producing a reasonable output

Next, let's combine all the individual chromosome files and produce a single GWAS file, and a variation of that for LDSC.
```{R}

library(data.table)

setwd("/mnt/b2/home4/arc/vw260/UKB_v2/headmotion_rs/")
data1 = fread("headmotionukb1.PHENO1.glm.linear")
data1 = subset(data1, TEST == "ADD")

for (i in 2:22){
  a = fread(paste0("headmotionukb", i, ".PHENO1.glm.linear"), header = T)
  a2 = subset(a, TEST == "ADD")
  data1 = rbind(data1, a2)
}

save(data1, file = "headmotion_rs_fullGWAS.RData")

data1 = data1[,c( "ID", "ALT", "REF", "BETA", "SE", "P", "OBS_CT")]

setnames(data1, "ALT", "A1")
setnames(data1, "REF", "A2")
setnames(data1, "OBS_CT", "N")

write.table(data1, file = "headmotion_rs_forLDSC.txt", row.names = F, col.names = T, quote = F )

rm(list = ls())

```

### Step 5: Running genetic correlations and heritability using LDSC

```bash

```


### Step 6: Running GWAS in BOLT-LMM


```bash
./bolt --bed=ukbchr22.bed --bim=ukbchr22.bim --fam=ukbchr22.fam --phenoFile=headmotionphenobolt.txt --phenoCol=rs_headmotion --covarFile=headmotioncovarbolt.txt --covarCol=f.22000.0.0 --covarCol=f.22001.0.0 --qCovarCol=f.22009.0.{1..20} --qCovarCol=f.34.0.0 --covarMaxLevels=200 --lmm --LDscoresFile=LDSCORE.1000G_EUR.tab.gz --geneticMapFile=genetic_map_hg19_withX.txt.gz --lmmForceNonInf --numThreads=20 --statsFile=UKB_chr21headmotion.txt

```


### Step 7: Combining the files and creating manhattan and QQ plots



### Step 7: Running heritability estimates in BOLT-LMM 


