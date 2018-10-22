# UKB commands
## This repository is for running GWAS in the UKB

### Step 1: Dowloading files from the UKB

First, read this file for downloading UKB files: http://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=100314

Ensure the key is downloaded in the file and saved as .ukbkey

```bash
for i in {1..22}; do ./ukbgene imp -c$i; done  #this downloads the BGEN files

for i in {1..22}; do ./ukbgene imp -m -c$i; done

./ukbgene imp -cX
./ukbgene imp -cXY

./ukbgene imp -m -cX
./ukbgene imp -m -cXY
```



### Step 2: Converting BGEN to plink files

```bash
for i in {1..21}; do ./plink2 --bgen ukb_imp_chr${i}_v3.bgen --sample ukb20904_imp_chr${i}_v3_s487334.sample --make-bed -out ukbchr${i} --maf 0.01 --geno 0.05 --threads 10 --hwe 0.000001 --mind 0.05; done

```

### Step 3: Creating the required phenotype and covariate files in R

Please see specific scripts in the attached folder.



### Step 4: Running GWAS in plink

```bash
{1..21}; do ./plink2 --bfile ukbchr$i --linear --pheno phenotypefile --covar covariatefile --out outfilename$i --threads 20; done

```



### Step 5: Combining the files and producing a reasonable output

Next, let's combine all the individual chromosome files and produce a single GWAS file, and a variation of that for LDSC.
```{R}

setwd("/mnt/b2/home4/arc/vw260/UKB_v2/Plink_files/")
library(data.table)
data1 = fread("tfMRIheadmotionukb1.PHENO1.glm.linear")
data1 = subset(data1, TEST == "ADD")

for (i in 2:22){
  a = fread(paste0("tfMRIheadmotionukb", i, ".PHENO1.glm.linear"), header = T)
  a2 = subset(a, TEST == "ADD")
  data1 = rbind(data1, a2)
}

save(data1, file = "headmotion_tfMRI_fullGWAS.RData")

data1 = data1[,c( "ID", "ALT", "REF", "BETA", "SE", "P", "OBS_CT")]

setnames(data1, "ID", "SNP")
setnames(data1, "ALT", "A1")
setnames(data1, "REF", "A2")
setnames(data1, "OBS_CT", "N")

write.table(data1, file = "headmotion_tfMRI_forLDSC.txt", row.names = F, col.names = T, quote = F )

rm(list = ls())

```

### Step 5: Running genetic correlations and heritability using LDSC

First,  munge the data in R
```R
munge(c("headmotion_tfMRI_forLDSC.txt"), "w_hm3.snplist",trait.names="headmotion_tfMRI", c(9966), info.filter = 0.9, maf.filter = 0.01)
```

Then run the scripts in LDSC to calculate heritability and gen cor. On the clusters, use Head1

```bash
./ldsc.py \
--h2 ./sumstats/headmotion_tfMRI.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--out ./gencorresults/headmotion_rs \
--w-ld-chr eur_w_ld_chr/ 

./ldsc.py \
--ref-ld-chr eur_w_ld_chr/ \
--out ./gencorresults/headmotion_tfMRI_rg \
--rg ./sumstats/headmotion_tfMRI.sumstats.gz,./sumstats/headmotion_rs.sumstats.gz,./sumstats/scz3.sumstats,./sumstats/adhd2.sumstats,./sumstats/autsim_ipsych.sumstats,./sumstats/sniekerscognition.sumstats,./sumstats/chronotype.sumstats,./sumstats/anorexia3.sumstats,./sumstats/edu2.sumstats,./sumstats/BPD.sumstats,./sumstats/SWB.sumstats,./sumstats/ICV.sumstats,./sumstats/neuroticism.sumstats,./sumstats/alzheimers.sumstats,./anxietycc.sumstats \
--w-ld-chr eur_w_ld_chr/ 

```


### Step 6: Running GWAS in BOLT-LMM

First lets produce a modelSNPs file
```R
library(data.table)
library(dplyr)

data1 = fread("~/UKB_v2/ukb_mfi_chr1_v3.txt")

for (i in 2:22){
  data2  = fread(paste0("~/UKB_v2/ukb_mfi_chr", i, "_v3.txt"))
  data1 = rbind(data1, data2)
  data1 = subset(data1, V6 > 0.05) #keep only SNPs with MAF > 0.2
  data1 = subset(data1, V8 > 0.95) #keep only SNPs with info > 0.99
  rm(data2)
}


data2 = subset(data1, V8 > 0.99)
data2 = sample_n(data2, 1000000)

data3 = data2[,c("V2")]

write.table(data3, file = "~/UKB_v2/Plink_files/modelSNPs.txt", row.names = F, col.names = T, quote = F)
```

We run two BOLT-LMMs, the first only for the autosomes, and the second including the X-chromosome as specified by the authors of BOLT
```bash
##resting state
./bolt --bed=ukbchr{1:22}.bed --bim=ukbchr{1:22}.bim --fam=chr21bolt.fam --phenoFile=headmotionphenobolt.txt --phenoCol=rs_headmotion --covarFile=headmotioncovarbolt.txt --covarCol=f.22000.0.0 --covarCol=f.22001.0.0 --qCovarCol=f.22009.0.{1..20} --qCovarCol=f.21003.2.0 --covarCol=f.54.2.0 --covarMaxLevels=200 --lmm --LDscoresFile=LDSCORE.1000G_EUR.tab.gz --geneticMapFile=genetic_map_hg19_withX.txt.gz --lmmForceNonInf --numThreads=10 --statsFile=UKB_rsheadmotionbolt.txt --remove=headmotionremovebolt.txt --modelSnps=modelSNPs.txt


##tfMRI
./bolt --bed=ukbchr{1:22}.bed --bim=ukbchr{1:22}.bim --fam=chr21bolt.fam --phenoFile=taskheadmotionphenobolt.txt --phenoCol=t_headmotion --covarFile=taskheadmotioncovarbolt.txt --covarCol=f.22000.0.0 --covarCol=f.22001.0.0 --qCovarCol=f.22009.0.{1..20} --qCovarCol=f.21003.2.0 --covarCol=f.54.2.0 --covarMaxLevels=200 --lmm --LDscoresFile=LDSCORE.1000G_EUR.tab.gz --geneticMapFile=genetic_map_hg19_withX.txt.gz --lmmForceNonInf --numThreads=20 --statsFile=UKB_taskheadmotionbolt_v2.txt --remove=taskheadmotionremovebolt.txt --modelSnps=modelSNPs.txt



./bolt --bed=ukbchr{1:22}.bed --bim=ukbchr{1:22}.bim --fam=chr21bolt.fam --phenoFile=headmotionphenobolt.txt --phenoCol=rs_headmotion --covarFile=headmotioncovarbolt.txt --covarCol=f.22000.0.0 --covarCol=f.22001.0.0 --qCovarCol=f.22009.0.{1..20} --qCovarCol=f.21003.2.0 --covarCol=f.54.2.0 --covarMaxLevels=200 --lmm --LDscoresFile=LDSCORE.1000G_EUR.tab.gz --geneticMapFile=genetic_map_hg19_withX.txt.gz --lmmForceNonInf --numThreads=10 --statsFile=UKB_rsheadmotionbolt.txt --remove=headmotionremovebolt.txt --modelSnps=modelSNPs.txt





./bolt --bed=ukbchr{1:23}.bed --bim=ukbchr{1:23}.bim --fam=chr21bolt.fam --phenoFile=headmotionphenobolt.txt --phenoCol=rs_headmotion --covarFile=headmotioncovarbolt.txt --covarCol=f.22000.0.0 --covarCol=f.22001.0.0 --qCovarCol=f.22009.0.{1..20} --qCovarCol=f.34.0.0 --covarMaxLevels=200 --lmm --LDscoresFile=LDSCORE.1000G_EUR.tab.gz --geneticMapFile=genetic_map_hg19_withX.txt.gz --lmmForceNonInf --numThreads=10 --statsFile=UKB_rsheadmotionbolt_X.txt --remove=headmotionremovebolt.txt --modelSnps=modelSNPs.txt
```


### Step 8: Let's MTAG this up!

First, read the files in R, and create the necessary files for MTAG

```R
library(data.table)

taskheadmotion = fread("~/UKB_v2/Plink_files/UKB_taskheadmotionbolt.txt")
rheadmotion = fread("~/UKB_v2/Plink_files/UKB_rsheadmotionbolt.txt")

setnames(taskheadmotion, 1, "snpid")
setnames(taskheadmotion, 2, "chr")
setnames(taskheadmotion, 3, "bpos")
setnames(taskheadmotion, 5, "a1")
setnames(taskheadmotion, 6, "a2")
setnames(taskheadmotion, 7, "freq")
setnames(taskheadmotion, 11, "pval")
taskheadmotion$z = taskheadmotion$BETA/taskheadmotion$SE
taskheadmotion$n = "10000"
head(taskheadmotion)

taskheadmotionmtag = taskheadmotion[,c("snpid", "chr", "bpos", "a1", "a2", "freq", "z","pval", "n")]

write.table(taskheadmotionmtag, file = "/mnt/b2/home4/arc/vw260/UKB_v1/Plink_files/mtag/mtag/taskheadmotionmtag.txt",
            row.names =F, col.names = T, quote = F)


setnames(rheadmotion, 1, "snpid")
setnames(rheadmotion, 2, "chr")
setnames(rheadmotion, 3, "bpos")
setnames(rheadmotion, 5, "a1")
setnames(rheadmotion, 6, "a2")
setnames(rheadmotion, 7, "freq")
setnames(rheadmotion, 11, "pval")
rheadmotion$z = rheadmotion$BETA/rheadmotion$SE
rheadmotion$n = "11000"

rheadmotionmtag = rheadmotion[,c("snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n")]

write.table(rheadmotionmtag, file = "/mnt/b2/home4/arc/vw260/UKB_v1/Plink_files/mtag/mtag/rheadmotionmtag.txt",
            row.names =F, col.names = T, quote = F)

```

Next, let's MTAG this

```bash
 python ~/UKB_v1/Plink_files/mtag/mtag/mtag.py --sumstats ~/UKB_v1/Plink_files/mtag/mtag/rheadmotionmtag.txt,~/UKB_v1/Plink_files/mtag/mtag/taskheadmotionmtag.txt --out ~/UKB_v1/Plink_files/mtag/mtag/rntheadmotionboltmtag --n_min 0.0 --stream_stdout &
```

### Step 8: Read the files, and upload them onto FUMA to get some functional insights

```


```

