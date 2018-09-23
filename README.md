# UKB commands
## This repository is for running GWAS in the UKB

### Step 1: Dowloading files from the UKB



### Step 2: Converting BGEN to plink files


### Step 3: Running GWAS in plink


### Step 4: Combining the files and producing a reasonable output

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


### Step 6: Running GWAS in BOLT-LMM


### Step 7: Combining the files and creating manhattan and QQ plots



### Step 7: Running heritability estimates in BOLT-LMM 


