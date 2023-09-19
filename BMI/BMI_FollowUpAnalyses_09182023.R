###------------------------###
### BMI Follow-up Analyses ###
###------------------------###

# Updated 09/18/2023

#### Inputs ####

# 1.  QC'd genotype matrix (genomatrix). You can use the same matrix used for CNV main effect and CNVxPRS interaction tests.
#     The necessary columns are:
#     a. sampleID
#     b. res_bmi_inv
#     c. sex
#     d. age
#     e. PC1, PC2, PC3... PC10
#     f. PRS_BMI_zscore (Yengo for EstBB, Geisinger, MVP. Locke for UKBB)
#     g. CNV genotype for all testable CNVs
# Place your genotype matrix into this directory for analysis

# 2. By chromosome PRS (Yengo for EstBB, Geisinger, MVP. Locke for UKBB)
#    After running provided script, put 22 output files into PRS_BMI_by_chromosome

# 3. tab-separated medication info formatted like so (1 for medication, 0 for no medication)
#   # sampleID medication
#     sample1  0
#     sample2  1
#     sample3  0
#     sample4  1

#### Outputs ####
#    Running this script on your data will write outputs to the output folder in this directory.
#    After running, share contents of this folder with msacks@ucsd.edu and jsebat@health.ucsd.edu

#### Document outline ####

# Collaborators should only need to edit section 0!

# 0. Import data and set up working environment (edit this section for your data)
# 1. Create additional columns for analysis 
# 2. Create additional dataframes for analysis
# 3. PRS quadratic
# 4. PRS odd/even
# 5. CNV_effect_size x PRS
# 6. PRS effect by CNV genotype
# 7. PRS effect by CNV group
# 8. CNV effect by PRS quantile
# 9. CNV group effect by PRS quantile
# 10. Sex x PRS interaction
# 11. Sex x CNV genotype interaction
# 12. Sex x CNV group interaction
# 13. PRS effect by sex
# 14. CNV genotype effect by sex
# 15. CNV group effect by sex
# 16. Medication x CNV group interaction
# 17. Medication x PRS interaction
# 18. Medication effect by CNV group
# 19. Medication effect by PRS group


###-----------------------------------------------###
### 0. Set up working directory and import data   ###
###  Sites should ONLY need to edit THIS section  ###
###-----------------------------------------------###

library(dplyr)
library(purrr)

####  Set working directory to this directory

# Fill in path to this directory
setwd('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByMolly/UKBB/BMI') 
dir.create('output')

#### Read in CNV main effects
main_effects <- read.csv('CNVmainEffects_4cohorts_allCNVs_BMI_20230531.csv')
testable <- main_effects[which((main_effects$cohort == 'Meta-analysis') & (main_effects$testable_CNV_n225 == 'yes')),]

#### Read in your genomatrix

## Read in your genomatrix (fill in path to your genomatrix)
my_genomatrix <- read.csv('UKBB_BMI_PRS-CS_20230824.csv')

## initialize blank genomatrix for this script
genomatrix <- data.frame(matrix(nrow=nrow(my_genomatrix), ncol=0))

## Rename columns (if necessary)
# covariates/ sample info
genomatrix$sampleID <- as.character(my_genomatrix$identifier) # whatever your sampleID column is named
genomatrix$sex <- my_genomatrix$sex_famfile_2022 # whatever your sex column is named
genomatrix$sex <- as.factor(genomatrix$sex)
genomatrix$age <- my_genomatrix$age # whatever your age (in years) column is named
genomatrix$res_bmi_inv <- my_genomatrix$res_bmi_inv # whatever your normalized BMI column is called
genomatrix$PC1 <- my_genomatrix$PC1 # PCs
genomatrix$PC2 <- my_genomatrix$PC2
genomatrix$PC3 <- my_genomatrix$PC3
genomatrix$PC4 <- my_genomatrix$PC4
genomatrix$PC5 <- my_genomatrix$PC5
genomatrix$PC6 <- my_genomatrix$PC6
genomatrix$PC7 <- my_genomatrix$PC7
genomatrix$PC8 <- my_genomatrix$PC8
genomatrix$PC9 <- my_genomatrix$PC9
genomatrix$PC10 <- my_genomatrix$PC10

# prs
genomatrix$PRS_BMI_zscore <- my_genomatrix$BMI_SCORE1_SUM_zscore # whatever your z-scored PRS_SCORE1_SUM column is called

# testable CNVs (rename based on what yours are named)
genomatrix$locus_13q12.12 <- my_genomatrix$locus_13q12.12_allele_specific
genomatrix$locus_15q11.2__BP1_BP2_CYFIP <- my_genomatrix$locus_15q11.2__BP1_BP2_CYFIP_allele_specific
genomatrix$locus_15q13.3_inclCHRNA7 <- my_genomatrix$locus_15q13.3_inclCHRNA7_allele_specific
genomatrix$locus_16p11.2_distal <- my_genomatrix$locus_16p11.2_distal_allele_specific
genomatrix$locus_16p11.2_proximal <- my_genomatrix$locus_16p11.2_proximal_allele_specific
genomatrix$locus_16p12.1_BP1_BP2 <- my_genomatrix$locus_16p12.1_BP1_BP2_allele_specific
genomatrix$locus_16p12.1_BP2_BP3 <- my_genomatrix$locus_16p12.1_BP2_BP3_allele_specific
genomatrix$locus_16p13.11_BP1_BP2 <- my_genomatrix$locus_16p13.11_BP1_BP2_allele_specific
genomatrix$locus_16p13.11_BP1_BP3 <- my_genomatrix$locus_16p13.11_BP1_BP3_allele_specific
genomatrix$locus_17p12 <- my_genomatrix$locus_17p12_allele_specific
genomatrix$locus_17q12 <- my_genomatrix$locus_17q12_allele_specific
genomatrix$locus_1q21.1 <- my_genomatrix$locus_1q21.1_allele_specific
genomatrix$locus_1q21.1_TAR_BP1_BP2 <- my_genomatrix$locus_1q21.1_TAR_BP1_BP2_allele_specific
genomatrix$locus_1q21.1_TAR_BP1_BP3 <- my_genomatrix$locus_1q21.1_TAR_BP1_BP3_allele_specific
genomatrix$locus_1q21.1_TAR_BP2_BP3 <- my_genomatrix$locus_1q21.1_TAR_BP2_BP3_allele_specific
genomatrix$locus_22q11.2_BP1_BP2_A_B <- my_genomatrix$locus_22q11.2_BP1_BP2_A_B_allele_specific
genomatrix$locus_22q11.2_BP1_BP4_A_D <- my_genomatrix$locus_22q11.2_BP1_BP4_A_D_allele_specific
genomatrix$locus_22q11.2_BP2_BP4_B_D <- my_genomatrix$locus_22q11.2_BP2_BP4_B_D_allele_specific
genomatrix$locus_22q11.2_BP3_BP4_C_D <- my_genomatrix$locus_22q11.2_BP3_BP4_C_D_allele_specific
genomatrix$locus_22q11.2_BP6_BP8_F_H <- my_genomatrix$locus_22q11.2_BP6_BP8_F_H_allele_specific
genomatrix$locus_22q13 <- my_genomatrix$locus_22q13_allele_specific
genomatrix$locus_2q13_NPHP1 <- my_genomatrix$locus_2q13_NPHP1_allele_specific

# Read in by chromosome PLINK output
# sampleIDs should match sampleIDs in genomatrix

prefix <- 'PRS_BMI_by_chromosome/UKBB_bmi_locke_pst_eff_a1_b0.5_phi1e-02.txt.chr'
suffix <- '.sscore'
names <- c('FID', 'IID', 'ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'SCORE1_AVG', 'SCORE1_SUM')

for(i in 1:22){
  # read in data
  PRS_i <- read.table(paste(prefix, as.character(i), suffix, sep=''), header = FALSE, col.names = names)
  rownames(PRS_i) <- as.character(PRS_i$IID)
  
  # filter to rows that are in genomatrix
  PRS_i <- PRS_i[genomatrix$sampleID, ]
  
  # add to genomatrix
  genomatrix[,paste('chr', as.character(i), 'PRS_BMI', sep='')] <- PRS_i$SCORE1_SUM
}

# Add medication information (if applicable)
genomatrix$medication <- my_genomatrix$medication # add medication column

# Additional QC (for UKBB, there were some PRS outliers, so we filter those out)
genomatrix <- genomatrix[which((genomatrix$PRS_BMI_zscore > -3) & (genomatrix$PRS_BMI_zscore < 3)),] # 

# The rest of this script should be able to run now. If there are errors, check:
# 1. Sex is coded as 1=Male, 2=Female
# 2. CNVs are coded 1,2,3 for del,no_cnv,dup
# 3. Sample IDs match between files
# 4. Columns consisting of ints/ floats are not represented as characters

###---------------------------------------------------###
### 1.  Create additional columns for analysis        ###
###---------------------------------------------------###

### Create odd/ even PRS columns

genomatrix$odd_chr_PRS_BMI_sum <- 
  genomatrix$chr1PRS_BMI +
  genomatrix$chr3PRS_BMI +
  genomatrix$chr5PRS_BMI +
  genomatrix$chr7PRS_BMI +
  genomatrix$chr9PRS_BMI +
  genomatrix$chr11PRS_BMI +
  genomatrix$chr13PRS_BMI +
  genomatrix$chr15PRS_BMI +
  genomatrix$chr17PRS_BMI +
  genomatrix$chr19PRS_BMI +
  genomatrix$chr21PRS_BMI 

genomatrix$odd_chr_PRS_BMI_zscore <- (genomatrix$odd_chr_PRS_BMI_sum - mean(genomatrix$odd_chr_PRS_BMI_sum))/sd(genomatrix$odd_chr_PRS_BMI_sum)

genomatrix$even_chr_PRS_BMI_sum <- 
  genomatrix$chr2PRS_BMI +
  genomatrix$chr4PRS_BMI +
  genomatrix$chr6PRS_BMI +
  genomatrix$chr8PRS_BMI +
  genomatrix$chr10PRS_BMI +
  genomatrix$chr12PRS_BMI +
  genomatrix$chr14PRS_BMI +
  genomatrix$chr16PRS_BMI +
  genomatrix$chr18PRS_BMI +
  genomatrix$chr20PRS_BMI +
  genomatrix$chr22PRS_BMI 

genomatrix$even_chr_PRS_BMI_zscore <- (genomatrix$even_chr_PRS_BMI_sum - mean(genomatrix$even_chr_PRS_BMI_sum))/sd(genomatrix$even_chr_PRS_BMI_sum)

### Create PRS squared column
genomatrix$PRS_BMI_zscore_squared <- genomatrix$PRS_BMI_zscore * genomatrix$PRS_BMI_zscore


### Create PRS quantile columns
breaks <- .25
genomatrix$odd_chr_PRS_BMI_quantile <- as.numeric(cut(genomatrix$odd_chr_PRS_BMI_zscore, 
                                                      breaks = quantile(genomatrix$odd_chr_PRS_BMI_zscore, 
                                                                        probs = seq(0, 1, by = breaks)), 
                                                      labels = FALSE))
genomatrix$even_chr_PRS_BMI_quantile <- as.numeric(cut(genomatrix$even_chr_PRS_BMI_zscore, 
                                                       breaks = quantile(genomatrix$even_chr_PRS_BMI_zscore, 
                                                                         probs = seq(0, 1, by = breaks)), 
                                                       labels = FALSE))
genomatrix$PRS_BMI_quantile <- as.numeric(cut(genomatrix$PRS_BMI_zscore, 
                                                      breaks = quantile(genomatrix$PRS_BMI_zscore, 
                                                                        probs = seq(0, 1, by = breaks)), 
                                                      labels = FALSE))

### Create CNV genotype column for all testable CNVs (allows for individuals with up to 4 CNVs)
testable_CNVs <- testable$locus
rownames(testable) <- testable$locus
keys <- testable$locus

genomatrix$CNV_genotype_1 <- NA
genomatrix$CNV_genotype_2 <- NA
genomatrix$CNV_genotype_3 <- NA
genomatrix$CNV_genotype_4 <- NA

for(i in seq_along(testable_CNVs)){
  CNV <- testable_CNVs[i] # get name of CNV
  loc <- paste("locus", substr(CNV, 1, nchar(CNV) - 4), sep='_') # get locus
  cn <- ifelse(testable[CNV, "genotype"] == 'deletion', 1, 3) # get copy number
  
  # get rows with CNV that don't already have another CNV
  rows1 <- which((genomatrix[[loc]] == cn) & (is.na(genomatrix$CNV_genotype_1)))  
  
  # get rows with CNV that already have 1 CNV
  rows2 <- which((genomatrix[[loc]] == cn) & (!is.na(genomatrix$CNV_genotype_1) & (is.na(genomatrix$CNV_genotype_2)))) 
  
  # get rows with CNV that already have 2 CNV
  rows3 <- which((genomatrix[[loc]] == cn) & (!is.na(genomatrix$CNV_genotype_2) & (is.na(genomatrix$CNV_genotype_3)))) 
  
  # get rows with CNV that already have 3 CNVs
  rows4 <- which((genomatrix[[loc]] == cn) & (!is.na(genomatrix$CNV_genotype_3) & (is.na(genomatrix$CNV_genotype_4)))) 
  
  
  # Fill in genotypes
  genomatrix[rows4, "CNV_genotype_4"] <- CNV
  genomatrix[rows3, "CNV_genotype_3"] <- CNV
  genomatrix[rows2, "CNV_genotype_2"] <- CNV
  genomatrix[rows1, "CNV_genotype_1"] <- CNV
  
}


### Create CNV group column
# individuals with CNVs in multiple groups are NA for this
genomatrix$group1 <- NA
genomatrix$group1 <- testable[genomatrix$CNV_genotype_1, 'CNV_effect_label']
genomatrix$group2 <- NA
genomatrix$group2 <- testable[genomatrix$CNV_genotype_2, 'CNV_effect_label']
genomatrix$group3 <- NA
genomatrix$group3 <- testable[genomatrix$CNV_genotype_3, 'CNV_effect_label']
genomatrix$group4 <- NA
genomatrix$group4 <- testable[genomatrix$CNV_genotype_4, 'CNV_effect_label']

genomatrix <- genomatrix %>%
   mutate(group = 
            ifelse((is.na(group1)) & (is.na(group2)) & (is.na(group3)) & (is.na(group4)), NA, 
                   ifelse(is.na(group2), group1,
                          ifelse(is.na(group3), ifelse(group1 == group2, group2, NA), 
                                 ifelse(is.na(group4), ifelse((group1 == group2) & (group2==group3), group3, NA), NA)))))

# switch suggestive positive to no effect
genomatrix[which(genomatrix$group == 'suggestive positive'), 'group'] <- 'no effect'

# set no effect as reference
genomatrix$group <- relevel(as.factor(genomatrix$group), ref='no effect')

### Create CNV effect size column for only significant positive/negative CNVs
# add effect size column that is 0 for CNVs with no effect or suggestive positive
testable$effect_size <- 0
testable[which(testable$CNV_effect_label %in% c('negative', 'positive')), 'effect_size'] <- 
  testable[which(testable$CNV_effect_label %in% c('negative', 'positive')), 'estimate_CNV_main']
testable['NA',] <- c(0,0,0,0,0,0,0,0,0,0)

# For individuals with multiple positive or negative CNVs, add effect sizes
genomatrix$CNV_effect_continuous <- 0
genomatrix$CNV_effect_continuous <- genomatrix$CNV_effect_continuous + testable[as.character(genomatrix$CNV_genotype_1),'effect_size']
genomatrix$CNV_effect_continuous <- genomatrix$CNV_effect_continuous + testable[as.character(genomatrix$CNV_genotype_2),'effect_size']
genomatrix$CNV_effect_continuous <- genomatrix$CNV_effect_continuous + testable[as.character(genomatrix$CNV_genotype_3),'effect_size']
genomatrix$CNV_effect_continuous <- genomatrix$CNV_effect_continuous + testable[as.character(genomatrix$CNV_genotype_4),'effect_size']


###---------------------------------------------------###
### 2.   Create additional dataframes for analysis    ###
###---------------------------------------------------###

# Sex
genomatrix_male <- genomatrix[which(genomatrix$sex == 1),]
genomatrix_female <- genomatrix[which(genomatrix$sex == 2),]

# CNV group
genomatrix_no_CNV <- genomatrix[which(is.na(genomatrix$group)),]
genomatrix_no_effect <- genomatrix[which(genomatrix$group == 'no effect'),]
genomatrix_neg <- genomatrix[which(genomatrix$group == 'negative'),]
genomatrix_positive <- genomatrix[which(genomatrix$group == 'positive'),]


###--------------------------------------------###
###  3. PRS quadratic                          ###
###--------------------------------------------###

dir.create('output/3_PRS_quadratic')

model_quadratic <- lm(res_bmi_inv ~ PRS_BMI_zscore + PRS_BMI_zscore_squared + age + sex + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

results_quadratic <- data.frame(summary(model_quadratic)$coefficients)
write.csv(results_quadratic, 'output/3_PRS_quadratic/PRS_quadratic.csv')


model_linear <- lm(res_bmi_inv ~ PRS_BMI_zscore + age + sex + PC1 + PC2 + PC3 + PC4 + 
                        PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      data = genomatrix, na.action=na.exclude)

results_linear <- data.frame(summary(model_linear)$coefficients)
write.csv(results_quadratic, 'output/3_PRS_quadratic/PRS_linear.csv')

###--------------------------------------------###
###  4. PRS odd/even                           ### 
###--------------------------------------------###

dir.create('output/4_PRS_odd_even')

model <- lm(res_bmi_inv ~ odd_chr_PRS_BMI_zscore*even_chr_PRS_BMI_zscore + age + sex + PC1 + PC2 + PC3 + PC4 + 
                     PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                   data = genomatrix, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/4_PRS_odd_even/PRS_oddeven.csv')

model <- lm(res_bmi_inv ~ odd_chr_PRS_BMI_zscore + age + sex + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/4_PRS_odd_even/PRS_odd.csv')


model <- lm(res_bmi_inv ~ even_chr_PRS_BMI_zscore + age + sex + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/4_PRS_odd_even/PRS_even.csv')


odd_by_even <- data.frame(matrix(nrow=0, ncol=7))
even_by_odd <- data.frame(matrix(nrow=0, ncol=7))

for(i in (1:4)){
  odd <- genomatrix[which(genomatrix$odd_chr_PRS_BMI_quantile == i),]
  model_odd <- lm(res_bmi_inv ~ 
                    even_chr_PRS_BMI_zscore +
                    age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                  data = odd, na.action=na.exclude)
  results <- summary(model_odd)$coefficients
  estimate <- as.numeric(results['even_chr_PRS_BMI_zscore', 1])
  se <- as.numeric(results['even_chr_PRS_BMI_zscore', 2])
  z <- as.numeric(results['even_chr_PRS_BMI_zscore', 3])
  p <- as.numeric(results['even_chr_PRS_BMI_zscore', 4])
  ci_lower <- estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  even_by_odd[nrow(even_by_odd) +1,] <- c(i, estimate, se, ci_lower, ci_upper, z, p)
  
  even <- genomatrix[which(genomatrix$even_chr_PRS_BMI_quantile == i),]
  model_even <- lm(res_bmi_inv ~ 
                     odd_chr_PRS_BMI_zscore +
                     age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                   data = even, na.action=na.exclude)
  results <- summary(model_even)$coefficients
  estimate <- as.numeric(results['odd_chr_PRS_BMI_zscore', 1])
  se <- as.numeric(results['odd_chr_PRS_BMI_zscore', 2])
  z <- as.numeric(results['odd_chr_PRS_BMI_zscore', 3])
  p <- as.numeric(results['odd_chr_PRS_BMI_zscore', 4])
  ci_lower <- estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  odd_by_even[nrow(odd_by_even) +1,] <- c(i, estimate, se, ci_lower, ci_upper, z, p)
  
}
colnames(odd_by_even) <- c('quantile', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')
colnames(even_by_odd) <- c('quantile', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')

write.csv(odd_by_even, 'output/4_PRS_odd_even/odd_by_even.csv')
write.csv(even_by_odd, 'output/4_PRS_odd_even/even_by_odd.csv')

###--------------------------------------------###
###  5. CNV_effect_size x PRS                  ### 
###--------------------------------------------###

dir.create('output/5_CNVxPRS_continuous')

model <- lm(res_bmi_inv ~ 
              PRS_BMI_zscore*CNV_effect_continuous +
              age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = odd, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/5_CNVxPRS_continuous/CNVxPRS_continuous.csv')


###--------------------------------------------###
###  6. PRS effect by CNV genotype             ### 
###--------------------------------------------###

dir.create('output/6_PRS_by_CNVgenotype')

out <- data.frame(matrix(nrow=0, ncol=8))

for(i in seq_along(keys)){
  genotype <- keys[i]
  df <- genomatrix[which(genomatrix$CNV_genotype_1 == genotype),]
  if(nrow(df) == 0){
    out[nrow(out) +1,] <- c(genotype,nrow(df), 0, 0, 0, 0, 0, 0)
    next
  }
  model <- lm(res_bmi_inv ~ 
                PRS_BMI_zscore +
                age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = df, na.action=na.exclude)
  results <- summary(model)$coefficients
  estimate <- as.numeric(results['PRS_BMI_zscore', 1])
  se <- as.numeric(results['PRS_BMI_zscore', 2])
  z <- as.numeric(results['PRS_BMI_zscore', 3])
  p <- as.numeric(results['PRS_BMI_zscore', 4])
  ci_lower <- estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  out[nrow(out) +1,] <- c(genotype, nrow(df), estimate, se, ci_lower, ci_upper, z, p)
}

colnames(out) <- c('locus', 'count', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')

write.csv(out, 'output/6_PRS_by_CNVgenotype/PRS_by_CNVgenotype.csv')


###--------------------------------------------###
###  7. PRS effect by CNV group                ### 
###--------------------------------------------###

dir.create('output/7_PRS_by_CNVgroup')

genotypes <- c(NA, 'negative', 'no effect', 'positive')
out <- data.frame(matrix(nrow=0, ncol=7))

for(i in seq_along(genotypes)){
  gen <- genotypes[i]
  df <- genomatrix[which(genomatrix$group== gen),]
  if(is.na(gen)){
    df <- genomatrix[which(is.na(genomatrix$group)),]
  }
  model <- lm(res_bmi_inv ~ PRS_BMI_zscore + age + sex + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = df, na.action=na.exclude)
  res <- summary(model)$coefficients
  estimate <- as.numeric(res['PRS_BMI_zscore', 1])
  se <- as.numeric(res['PRS_BMI_zscore', 2])
  ci_lower <-estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  z <- as.numeric(res['PRS_BMI_zscore', 3])
  p <- as.numeric(res['PRS_BMI_zscore', 4])
  out[nrow(out) + 1,] <- c(gen, estimate, se, ci_lower, ci_upper, z, p)
}


colnames(out) <- c('CNV_group', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')
write.csv(out, 'output/7_PRS_by_CNVgroup/PRS_by_CNVgroup.csv')


###--------------------------------------------###
###  8. CNV effect by PRS quantile             ### 
###--------------------------------------------###

dir.create('output/8_CNV_by_PRS')

out <- data.frame(matrix(nrow=0, ncol=9))

for(i in seq_along(keys)){
  genotype <- keys[i]
  geno <- paste('locus', substr(genotype, 1, nchar(genotype) - 4),  sep='_')
  cn <- ifelse(substr(genotype, nchar(genotype) - 3, nchar(genotype)) == '_DUP', 3, 1)
  for(j in 1:4){
    df <- genomatrix[which((genomatrix$PRS_BMI_quantile == j)),]
    if(nrow(df) == 0){
      out[nrow(out) +1,] <- c(genotype,j,nrow(df), 0, 0, 0, 0, 0, 0)
      next
    }
    carriers <- df[which(df$CNV_genotype_1 == genotype),]
    if(nrow(carriers) == 0){
      out[nrow(out) +1,] <- c(genotype,j,nrow(carriers), 0, 0, 0, 0, 0, 0)
      next
    }
    df$gen <- df[[geno]]
    df$gen <- relevel(as.factor(df$gen), ref="2")
    model <- lm(res_bmi_inv ~ gen + age + sex + PC1 + PC2 + PC3 + PC4 + 
                  PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data = df, na.action=na.exclude)
    res <- summary(model)$coefficients
    cov <- paste('gen', as.character(cn),sep='')
    estimate <- as.numeric(res[cov, 1])
    se <- as.numeric(res[cov, 2])
    ci_lower <-estimate - (1.96*se)
    ci_upper <- estimate + (1.96*se)
    z <- as.numeric(res[cov, 3])
    p <- as.numeric(res[cov, 4])
    out[nrow(out) + 1,] <- c(genotype, j, nrow(carriers), estimate, se, ci_lower, ci_upper, z, p)
  }
}


colnames(out) <- c('locus', 'PRS_quantile', 'count', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')
write.csv(out, 'output/8_CNV_by_PRS/CNV_by_PRS.csv')



###--------------------------------------------###
###  9. CNV group effect by PRS quantile       ### 
###--------------------------------------------###

dir.create('output/9_CNVgroup_by_PRS')

genotypes <- c(NA, 'negative', 'no effect', 'positive')
out <- data.frame(matrix(nrow=0, ncol=8))

genomatrix$group <- relevel(as.factor(genomatrix$group), ref='no effect')

for(j in 1:4){
  df <- genomatrix[which(genomatrix$PRS_BMI_quantile == j),]
  model <- lm(res_bmi_inv ~ group + age + sex + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = df, na.action=na.exclude)
  res <- summary(model)$coefficients
  estimate <- as.numeric(res['groupnegative', 1])
  se <- as.numeric(res['groupnegative', 2])
  ci_lower <-estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  z <- as.numeric(res['groupnegative', 3])
  p <- as.numeric(res['groupnegative', 4])
  out[nrow(out) + 1,] <- c('negative', j, estimate, se, ci_lower, ci_upper, z, p)
  
  estimate <- as.numeric(res['grouppositive', 1])
  se <- as.numeric(res['grouppositive', 2])
  ci_lower <-estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  z <- as.numeric(res['grouppositive', 3])
  p <- as.numeric(res['grouppositive', 4])
  out[nrow(out) + 1,] <- c('positive', j, estimate, se, ci_lower, ci_upper, z, p)
}

colnames(out) <- c('CNV_group', 'PRS_quantile', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')
write.csv(out, 'output/9_CNVgroup_by_PRS/CNVgroup_by_PRS.csv')


###--------------------------------------------###
###  10. PRS x Sex interaction                 ### 
###--------------------------------------------###

dir.create('output/10_PRSxSex')

model <- lm(res_bmi_inv ~ 
              PRS_BMI_zscore*sex +
              age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/10_PRSxSex/PRSxSex.csv')


###--------------------------------------------###
###  11. CNV genotype x Sex interaction        ### 
###--------------------------------------------###
dir.create('output/11_CNVxSex')

out <- data.frame(matrix(nrow=0, ncol=7))

for(i in seq_along(keys)){
  genotype <- keys[i]
  geno <- paste('locus', substr(genotype, 1, nchar(genotype) - 4),  sep='_')
  cn <- ifelse(substr(genotype, nchar(genotype) - 3, nchar(genotype)) == '_DUP', 3, 1)
  carriers <- genomatrix[which(genomatrix$CNV_genotype_1 == genotype),]
  if(length(table(carriers$sex)) == 0){
    out[nrow(out) + 1,] <- c(genotype, NA, NA, NA, NA, NA, NA)
    next
  }
  if(nrow(carriers) == 0){
    out[nrow(out) + 1,] <- c(genotype, NA, NA, NA, NA, NA, NA)
    next
  }
  df$gen <- df[[geno]]
  df$gen <- relevel(as.factor(df$gen), ref="2")
  model <- lm(res_bmi_inv ~ gen*sex + age  + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = df, na.action=na.exclude)
  res <- summary(model)$coefficients
  cov <- paste('gen', as.character(cn), ":sex2", sep='')
  estimate <- as.numeric(res[cov, 1])
  se <- as.numeric(res[cov, 2])
  ci_lower <-estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  z <- as.numeric(res[cov, 3])
  p <- as.numeric(res[cov, 4])
  out[nrow(out) + 1,] <- c(genotype, estimate, se, ci_lower, ci_upper, z, p)
  }


colnames(out) <- c('locus', 'CNVxSex_estimate', 'CNVxSex_se', 'CNVxSex_ci_lower', 'CNVxSex_ci_upper', 'CNVxSex_z', 'CNVxSex_p')
write.csv(out, 'output/11_CNVxSex/CNVxSex.csv')


###--------------------------------------------###
###  12. CNV group x Sex interaction           ### 
###--------------------------------------------###
dir.create('output/12_CNVgroupxSex')

out <- data.frame(matrix(nrow=0, ncol=7))

model <- lm(res_bmi_inv ~ group*sex + age + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix, na.action=na.exclude)
res <- summary(model)$coefficients

estimate <- as.numeric(res['groupnegative:sex2', 1])
se <- as.numeric(res['groupnegative:sex2', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['groupnegative:sex2', 3])
p <- as.numeric(res['groupnegative:sex2', 4])
out[nrow(out) + 1,] <- c('negative', estimate, se, ci_lower, ci_upper, z, p)

estimate <- as.numeric(res['grouppositive:sex2', 1])
se <- as.numeric(res['grouppositive:sex2', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['grouppositive:sex2', 3])
p <- as.numeric(res['grouppositive:sex2', 4])
out[nrow(out) + 1,] <- c('positive', estimate, se, ci_lower, ci_upper, z, p)



colnames(out) <- c('group', 'CNV_groupxSex_estimate', 'CNV_groupxSex_se', 'CNV_groupxSex_ci_lower', 'CNV_groupxSex_ci_upper', 'CNV_groupxSex_z', 'CNV_groupxSex_p')
write.csv(out, 'output/12_CNVgroupxSex/CNVgroupxSex.csv')


###--------------------------------------------###
###  13. PRS effect by sex                     ### 
###--------------------------------------------###

dir.create('output/13_PRS_by_Sex')

model <- lm(res_bmi_inv ~ 
              PRS_BMI_zscore +
              age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix_male, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/13_PRS_by_Sex/PRSeffect_male.csv')

model <- lm(res_bmi_inv ~ 
              PRS_BMI_zscore +
              age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix_female, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/13_PRS_by_Sex/PRSeffect_female.csv')


###--------------------------------------------###
###  14. CNV genotype effect by sex            ### 
###--------------------------------------------###

dir.create('output/14_CNV_by_Sex')

out <- data.frame(matrix(nrow=0, ncol=8))

for(i in seq_along(keys)){
  genotype <- keys[i]
  geno <- paste('locus', substr(genotype, 1, nchar(genotype) - 4),  sep='_')
  cn <- ifelse(substr(genotype, nchar(genotype) - 3, nchar(genotype)) == '_DUP', 3, 1)
  carriers1 <- genomatrix_male[which(genomatrix_male$CNV_genotype_1 == genotype),]
  carriers2 <- genomatrix_female[which(genomatrix_female$CNV_genotype_1 == genotype),]
  if(nrow(carriers1) == 0 | nrow(carriers2) == 0){
    next
  }
  genomatrix_male$gen <- genomatrix_male[[geno]]
  genomatrix_male$gen <- relevel(as.factor(genomatrix_male$gen), ref="2")
  
  model_male <- lm(res_bmi_inv ~ gen + age + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = genomatrix_male, na.action=na.exclude)
  res <- summary(model_male)$coefficients
  cov <- paste('gen', as.character(cn), sep='')
  estimate <- as.numeric(res[cov, 1])
  se <- as.numeric(res[cov, 2])
  ci_lower <-estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  z <- as.numeric(res[cov, 3])
  p <- as.numeric(res[cov, 4])
  out[nrow(out) + 1,] <- c(genotype, 'male', estimate, se, ci_lower, ci_upper, z, p)
  
  genomatrix_female$gen <- genomatrix_female[[geno]]
  genomatrix_female$gen <- relevel(as.factor(genomatrix_female$gen), ref="2")
  
  model_female <- lm(res_bmi_inv ~ gen + age + PC1 + PC2 + PC3 + PC4 + 
                     PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                   data = genomatrix_female, na.action=na.exclude)
  res <- summary(model_female)$coefficients
  cov <- paste('gen', as.character(cn), sep='')
  estimate <- as.numeric(res[cov, 1])
  se <- as.numeric(res[cov, 2])
  ci_lower <-estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  z <- as.numeric(res[cov, 3])
  p <- as.numeric(res[cov, 4])
  out[nrow(out) + 1,] <- c(genotype, 'female', estimate, se, ci_lower, ci_upper, z, p)
}

colnames(out) <- c('CNV', 'sex', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')

write.csv(out, 'output/14_CNV_by_Sex/CNV_by_Sex.csv')


###--------------------------------------------###
###  15. CNV group effect by sex               ### 
###--------------------------------------------###

dir.create('output/15_CNVgroup_by_Sex')


model_male <- lm(res_bmi_inv ~ group + age + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix_male, na.action=na.exclude)
res <- summary(model_male)$coefficients
write.csv(res, 'output/15_CNVgroup_by_Sex/CNVgroup_by_Sex_male.csv')

model_female <- lm(res_bmi_inv ~ group + age + PC1 + PC2 + PC3 + PC4 + 
                   PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix_female, na.action=na.exclude)
res <- summary(model_female)$coefficients
write.csv(res, 'output/15_CNVgroup_by_Sex/CNVgroup_by_Sex_female.csv')

###--------------------------------------------###
###  16. Medication x CNV group interaction    ### 
###--------------------------------------------###

if(length(table(genomatrix$medication)) < 2){
  stop("Analysis finished (not testing GeneXDrug)")
}

dir.create('output/16_CNVgroupxMedication')

out <- data.frame(matrix(nrow=0, ncol=7))

model <- lm(res_bmi_inv ~ group*medication + sex + age + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)
res <- summary(model)$coefficients

estimate <- as.numeric(res['groupnegative:medication', 1])
se <- as.numeric(res['groupnegative:medication', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['groupnegative:medication', 3])
p <- as.numeric(res['groupnegative:medication', 4])
out[nrow(out) + 1,] <- c('negative', estimate, se, ci_lower, ci_upper, z, p)

estimate <- as.numeric(res['grouppositive:medication', 1])
se <- as.numeric(res['grouppositive:medication', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['grouppositive:medication', 3])
p <- as.numeric(res['grouppositive:medication', 4])
out[nrow(out) + 1,] <- c('positive', estimate, se, ci_lower, ci_upper, z, p)



colnames(out) <- c('group', 'CNV_groupxMedication_estimate', 'CNV_groupxMedication_se', 'CNV_groupxMedication_ci_lower', 'CNV_groupxMedication_ci_upper', 'CNV_groupxMedication_z', 'CNV_groupxMedication_p')
write.csv(out, 'output/16_CNVgroupxMedication/CNVgroupxMedication.csv')


###--------------------------------------------###
###  17. Medication x PRS interaction          ### 
###--------------------------------------------###

dir.create('output/17_PRSxMedication')

model <- lm(res_bmi_inv ~ 
              PRS_BMI_zscore*medication +
              age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix, na.action=na.exclude)

results <- data.frame(summary(model)$coefficients)
write.csv(results, 'output/17_PRSxMedication/PRSxMedication.csv')


###--------------------------------------------###
###  18. Medication effect by CNV group        ### 
###--------------------------------------------###

dir.create('output/18_Medication_by_CNVgroup')

out <- data.frame(matrix(nrow=0, ncol=7))

model <- lm(res_bmi_inv ~ medication + age + sex + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix_no_CNV, na.action=na.exclude)
res <- summary(model)$coefficients
estimate <- as.numeric(res['medication', 1])
se <- as.numeric(res['medication', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['medication', 3])
p <- as.numeric(res['medication', 4])
out[nrow(out) + 1,] <- c('No CNV', estimate, se, ci_lower, ci_upper, z, p)

model <- lm(res_bmi_inv ~ medication + age + sex + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix_neg, na.action=na.exclude)
res <- summary(model)$coefficients
estimate <- as.numeric(res['medication', 1])
se <- as.numeric(res['medication', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['medication', 3])
p <- as.numeric(res['medication', 4])
out[nrow(out) + 1,] <- c('Negative', estimate, se, ci_lower, ci_upper, z, p)

model <- lm(res_bmi_inv ~ medication + age + sex + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix_no_effect, na.action=na.exclude)
res <- summary(model)$coefficients
estimate <- as.numeric(res['medication', 1])
se <- as.numeric(res['medication', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['medication', 3])
p <- as.numeric(res['medication', 4])
out[nrow(out) + 1,] <- c('No effect', estimate, se, ci_lower, ci_upper, z, p)

model <- lm(res_bmi_inv ~ medication + age + sex + PC1 + PC2 + PC3 + PC4 + 
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            data = genomatrix_positive, na.action=na.exclude)
res <- summary(model)$coefficients
estimate <- as.numeric(res['medication', 1])
se <- as.numeric(res['medication', 2])
ci_lower <-estimate - (1.96*se)
ci_upper <- estimate + (1.96*se)
z <- as.numeric(res['medication', 3])
p <- as.numeric(res['medication', 4])
out[nrow(out) + 1,] <- c('Positive', estimate, se, ci_lower, ci_upper, z, p)


colnames(out) <- c('CNV_group', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')
write.csv(out, 'output/18_Medication_by_CNVgroup/Medication_by_CNVgroup.csv')


###--------------------------------------------###
###  19. Medication effect by PRS group        ###    
###--------------------------------------------###


dir.create('output/19_Medication_by_PRS')

genotypes <- c(NA, 'negative', 'no effect', 'positive')
out <- data.frame(matrix(nrow=0, ncol=7))

for(j in 1:4){
  df <- genomatrix[which(genomatrix$PRS_BMI_quantile == j),]
  model <- lm(res_bmi_inv ~ medication + age + sex + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data = df, na.action=na.exclude)
  res <- summary(model)$coefficients
  estimate <- as.numeric(res['medication', 1])
  se <- as.numeric(res['medication', 2])
  ci_lower <-estimate - (1.96*se)
  ci_upper <- estimate + (1.96*se)
  z <- as.numeric(res['medication', 3])
  p <- as.numeric(res['medication', 4])
  out[nrow(out) + 1,] <- c(j, estimate, se, ci_lower, ci_upper, z, p)
}

colnames(out) <- c('PRS_quantile', 'estimate', 'se', 'ci_lower', 'ci_upper', 'z', 'p')
write.csv(out, 'output/19_Medication_by_PRS/Medication_by_PRS.csv')




model_meds <- lm(res_bmi_inv ~ PRS_BMI_zscore + medication + age + sex + PC1 + PC2 + PC3 + PC4 + 
                   PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix, na.action=na.exclude)

model_no_meds <- lm(res_bmi_inv ~  PRS_BMI_zscore + age + sex + PC1 + PC2 + PC3 + PC4 + 
                   PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                 data = genomatrix, na.action=na.exclude)


anova(model_no_meds, model_meds)
