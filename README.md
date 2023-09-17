# README.md

This repository contains code to perform follow-up analyses for our Biobank Meta-analysis of BMI and Height.


## Instructions

### 1. Clone this directory to a working directory
`cd your_working_directory`
`git clone 

### 2. Calculate by-chromosome polygenic risk scores for BMI and Height
The script `calculate_PRS_by_chromosome.sh` is provided for you to calculate by chromosome polygenic risk scores. 
Edit the script as needed and run on your data.

Required input:
- PRS-CS sumstats for BMI (Yengo et al) and Height
- PLINK fileset (.bed, .bim, .fam). This script will run faster if these are split by chromosome, but can also be run on the whole genome
### 3. Copy by chromosome PRSs for BMI to your_working_directory/Biobank_MetaAnalysis_FollowUp/BMI/PRS_BMI_by_chromosome

### 4. Copy by chromosome PRSs for Height to your_working_directory/Biobank_MetaAnalysis_FollowUp/Height/PRS_HEIGHT_by_chromosome

### 4. Upload genotype matrices to BMI and Height folder


### 5. Run follow-up analyses for Height
