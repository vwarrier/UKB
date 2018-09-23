# UKB
UKB commands


This repository is for running GWAS in the UKB



``{R}

library(data.table)

setwd("/mnt/b2/home4/arc/vw260/UKB_v2/headmotion_rs/")
data1 = fread("headmotionukb1.PHENO1.glm.linear")
data1 = subset(data1, TEST == "ADD")

for (i in 2:22){
  a = fread(paste0("headmotionukb", i, ".PHENO1.glm.linear"), header = T)
  a2 = subset(a, TEST == "ADD")
  data1 = rbind(data1, a2)
}

```
