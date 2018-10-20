## This script now generates files for plink and BOLT for the headmotion statistics in the UK Biobank

# R program ukb22918.tab created 2018-09-17 by ukb2r.cpp Mar 14 2018 14:22:05
library(data.table)
library(plyr)


bd <- fread("~/UKB_v2/ukb24196.tab", header=TRUE, sep="\t")
#bd$f.53.0.0 <- as.Date(bd$f.53.0.0)
#bd$f.53.1.0 <- as.Date(bd$f.53.1.0)
#bd$f.53.2.0 <- as.Date(bd$f.53.2.0)
#bd$f.20400.0.0 <- as.Date(bd$f.20400.0.0) # date of cancer diagnosis
#bd$f.40005.0.0 <- as.Date(bd$f.40005.0.0)
#bd$f.40005.1.0 <- as.Date(bd$f.40005.1.0)
#bd$f.40005.2.0 <- as.Date(bd$f.40005.2.0)
#bd$f.40005.3.0 <- as.Date(bd$f.40005.3.0)
#bd$f.40005.4.0 <- as.Date(bd$f.40005.4.0)
#bd$f.40005.5.0 <- as.Date(bd$f.40005.5.0)
#bd$f.40005.6.0 <- as.Date(bd$f.40005.6.0)
#bd$f.40005.7.0 <- as.Date(bd$f.40005.7.0)
#bd$f.40005.8.0 <- as.Date(bd$f.40005.8.0)
#bd$f.40005.9.0 <- as.Date(bd$f.40005.9.0)
#bd$f.40005.10.0 <- as.Date(bd$f.40005.10.0)
#bd$f.40005.11.0 <- as.Date(bd$f.40005.11.0)
#bd$f.40005.12.0 <- as.Date(bd$f.40005.12.0)
#bd$f.40005.13.0 <- as.Date(bd$f.40005.13.0)
#bd$f.40005.14.0 <- as.Date(bd$f.40005.14.0)
#bd$f.40005.15.0 <- as.Date(bd$f.40005.15.0)
#bd$f.40005.16.0 <- as.Date(bd$f.40005.16.0)
#bd$f.40005.17.0 <- as.Date(bd$f.40005.17.0)
#bd$f.40005.18.0 <- as.Date(bd$f.40005.18.0)
#bd$f.40005.19.0 <- as.Date(bd$f.40005.19.0)
#bd$f.40005.20.0 <- as.Date(bd$f.40005.20.0)
#bd$f.40005.21.0 <- as.Date(bd$f.40005.21.0)
#bd$f.40005.22.0 <- as.Date(bd$f.40005.22.0)
#bd$f.40005.23.0 <- as.Date(bd$f.40005.23.0)
#bd$f.40005.24.0 <- as.Date(bd$f.40005.24.0)
#bd$f.40005.25.0 <- as.Date(bd$f.40005.25.0)
#bd$f.40005.26.0 <- as.Date(bd$f.40005.26.0)
#bd$f.40005.27.0 <- as.Date(bd$f.40005.27.0)
#bd$f.40005.28.0 <- as.Date(bd$f.40005.28.0)
#bd$f.40005.29.0 <- as.Date(bd$f.40005.29.0)
#bd$f.40005.30.0 <- as.Date(bd$f.40005.30.0)
#bd$f.40005.31.0 <- as.Date(bd$f.40005.31.0)
#bd$f.42006.0.0 <- as.Date(bd$f.42006.0.0) # first trike outcome
#bd$f.42008.0.0 <- as.Date(bd$f.42008.0.0) # ischemic stroke outcome
#bd$f.42010.0.0 <- as.Date(bd$f.42010.0.0) # introcerebral heamorrhage
#bd$f.42012.0.0 <- as.Date(bd$f.42012.0.0) # subarachnoid heamorrhage

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
pheno2$checksex = ifelse(pheno2$f.31.0.0 == pheno2$f.22001.0.0, "correct", "incorrect")
pheno2 = subset(pheno2, checksex == "correct") # remove sex mismatches
pheno2 = pheno2[is.na(pheno2$f.22027.0.0),] #remove excessive heterozygosity


####FILES FOR PLINK#####

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



####FILES FOR BOLT-LMM###

headmotion_pheno = pheno2[,c("f.eid", "f.eid", "f.25741.2.0")]
setnames(headmotion_pheno, 1, "FID")
setnames(headmotion_pheno, 2, "IID")
setnames(headmotion_pheno, 3, "rs_headmotion")

write.table(headmotion_pheno, file = "~/UKB_v2/headmotion_rs/headmotionphenobolt.txt", row.names = F, col.names = T, quote = F)


headmotion_covar = pheno2[,c("f.eid", "f.eid", "f.22009.0.1",
                              "f.22009.0.2", "f.22009.0.3", "f.22009.0.4", "f.22009.0.5",
                              "f.22009.0.6", "f.22009.0.7", "f.22009.0.8", "f.22009.0.9",
                              "f.22009.0.10", "f.22009.0.11", "f.22009.0.12", "f.22009.0.13",
                              "f.22009.0.14", "f.22009.0.15", "f.22009.0.16", "f.22009.0.17",
                              "f.22009.0.18", "f.22009.0.19", "f.22009.0.20", "f.22009.0.21",
                              "f.22009.0.22", "f.22009.0.23", "f.22009.0.24", "f.22009.0.25",
                              "f.22009.0.26", "f.22009.0.27", "f.22009.0.28", "f.22009.0.29",
                              "f.22009.0.30", "f.22009.0.31", "f.22009.0.32", "f.22009.0.33",
                              "f.22009.0.34", "f.22009.0.35", "f.22009.0.36", "f.22009.0.37",
                              "f.22009.0.38", "f.22009.0.39", "f.22009.0.40", "f.34.0.0", "f.22000.0.0",
                              "f.22001.0.0")]

setnames(headmotion_covar, 1, "FID")
setnames(headmotion_covar, 2, "IID")


write.table(headmotion_covar, file = "~/UKB_v2/headmotion_rs/headmotioncovarbolt.txt", row.names = F, col.names = T, quote = F)


remove_file = bd[!(bd$f.eid %in% headmotion_pheno$FID),]

remove_file = remove_file[,c("f.eid", "f.eid")]
setnames(remove_file, 1, "FID")
setnames(remove_file, 2, "IID")

fam = fread("~/UKB_v2/Plink_files/ukbchr21.fam")
remove3 = remove_file[remove_file$FID %in% fam$V1,]

write.table(remove3, file = "~/UKB_v2/Plink_files/headmotionremovebolt.txt", row.names = F, col.names = T, quote = F)

## f.25741.2.0 is the head motion variable. 
##22009 - PCs; 22006 - Ethnic grouping; 22001 - genetic sex; 22000 - batch; 34 year of birth
