This repository contains code to perform follow-up analyses for our Biobank Meta-analysis of BMI and Height.

# Instructions

### 1. Download this repo by clicking Code > Download .zip

### 2. Calculate by-chromosome polygenic risk scores for BMI and Height
The script `calculate_PRS_by_chromosome.sh` is provided for you to calculate by chromosome polygenic risk scores. 
Edit the script as needed and run on your data.

Required input:
- PRS-CS sumstats for BMI (Yengo et al) and Height
- PLINK fileset (.bed, .bim, .fam). This script will run faster if these are split by chromosome, but can also be run on the whole genome
### 3. Copy by chromosome PRSs for BMI to `BMI/PRS_BMI_by_chromosome`
There should be 22 files
### 4. Copy by chromosome PRSs for Height to `Height/PRS_HEIGHT_by_chromosome`
There should be 22 files
### 5. Upload genotype matrices to `BMI/` and `Height/`
The genotype matrices used to calculate CNVxPRS interaction effect for both Height and BMI should be sufficient here. These files should contain the columns:
- res_bmi_inv (or res_height_inv)
- sex, coded as 1 for male 2 for female
- `age` in years
- Ancestry PCs 1-10
- CNV genotypes for each of the testable CNVs (coded as 1 for del, 2 for no CNV, and 3 for dup)
- PRS-BMI zscore (or PRS-Height zscore)
- medication (if applicable, BMI only)
Column names aren't important here as long as all the data is present. We will edit those in step 6!

### 6. Edit section 0 of Height/Height_FollowUpAnalyses_09182023.R and /BMI/BMI_FollowUpAnalyses_09182023.R
This section reads in your data and compiles it into a genotype matrix so that the statistical tests can be run automatically. 
For example, in my genotype matrix, sex was in the column `my_genomatrix$sex_famfile_2022`. This script creates a new genotype matrix where the sex column is called `genomatrix$sex`. 

### 7. Run Height/Height_FollowUpAnalyses_09182023.R and /BMI/BMI_FollowUpAnalyses_09182023.R
After editing section 0 to properly read in your data, the rest of the script should run without any additional user input.
Output will be written to the `Height/output` and `BMI/output` directories.
For UKBB (~250000 subjects), these scripts each took less than 5 minutes to run on my laptop. 

### 8. Send `Height/output` and `BMI/output` to msacks@ucsd.edu and jsebat@health.ucsd.edu
For any troubleshooting, feel free to email msacks@ucsd.edu. Thank you for your cooperation!

