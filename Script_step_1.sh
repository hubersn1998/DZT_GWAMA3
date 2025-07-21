#! /bin/bash

# Step 1: making a normal and sparse GRM
# Input files:
# Phenotyped_GT_set.bim/bed/fam": Plink files containing the SNPs for the GRM for all participants with phenotype data, obtained through script step 0.

#Begin script step 1

# Define the genotype plink files
GTPHENO="./Path/To/Directory/Output/Phenotyped_GT_set"

# Make the GRM
OUTPUT0=Pheno_
./king -b $GTPHENO".bed" --unrelated --degree 2 --prefix ${OUTPUT0}
OUTPUT1=./Path/To/Directory/Output/Prune
./plink2 --bfile $GTPHENO --keep ${OUTPUT0}unrelated.txt --indep 50 5 2 --out ${OUTPUT1}
OUTPUT2=./Path/To/Directory/Output/Prune_GT_set
./plink2 --bfile $GTPHENO --extract ${OUTPUT1}.prune.in --make-bed --out ${OUTPUT2}
OUTPUT3=./Path/To/Directory/Output/Phenotyped_GRM
./gcta64 --bfile ${OUTPUT2} --make-grm --autosome --thread-num 12 --out ${OUTPUT3}
#^  Note threadnum can be adjusted for the number of cores in your system.
OUTPUT4=./Path/To/Directory/Output/Phenotyped_GRM_sparse
./gcta64 --grm ${OUTPUT3} --make-bK-sparse 0.05 --out ${OUTPUT4}

# End script step 1
# essential output files:
# Phenotyped_GRM_sparse: File containing the GRM, required for script step 2
