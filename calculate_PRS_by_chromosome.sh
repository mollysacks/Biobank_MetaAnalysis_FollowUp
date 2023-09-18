#!/bin/bash

#SBATCH --account=your_account
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --job-name=PRS_by_chr
#SBATCH --output=PRSbychr.%j.%N.out
#SBATCH --partition=your_partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH -t 24:00:00
#SBATCH --mem=100G
module load cpu/0.15.4
module load gcc/10.2.0
module load plink

mkdir /your/output/directory
mkdir /your/output/directory/sumstats_by_chr

# for each set of sumstats
for sumstats in /your/sumstats/directory/*; do
  filename=$(basename "$sumstats")
  # for each chromosome
  for i in {1..22}; do
    # Create sumstats for chromosome i
    sumstats_bychr="/your/output/directory/sumstats_by_chr/${filename}.chr${i}.txt"
    grep "^${i}" "${sumstats}" > "${sumstats_bychr}"

    # calculate PRS for chromosome i (if your bed files aren't split by chromosome that's ok)
    plink2 \
    --bfile /your/plink/fileset/your_dataset_chr${i} \
    --score ${sumstats_bychr} 2 4 6 cols=+scoresums \
    --out /your/output/directory/your_dataset_${filename}.chr${i} &
  done
done

wait