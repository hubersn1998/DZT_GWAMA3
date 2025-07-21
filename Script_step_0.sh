# Step 0: Select SNPs for the GRM

# Input files required for all 3 scripts:

#1 -Phenotype.dat: a tab or space text file containing the samples with three columns: Family_ID, Individual_ID and Phenotype.
#  The phenotype should be coded as 0= controls and 1 = case. The file should not contain a header. 

#2 - Qcovar.dat: a tab or space delimited text file containing the samples with the columns: Family_ID, Individual_ID, and then 10 columns 
#  with the 1000 genomes projected PCS and if needed other study specific quantitative covariates.

#3 - Covar.dat: a tab or space delimited text file containing the samples with the columns: Family_ID, Individual_ID, and then SEX and 
# any lab, batch, or other study specific (dummy-)categorical covariates. 

#4 HRC imputed chromosome files 22 for the autosomes and one for chromosome X. e.g. chr1.dose.vf

# Make sure the HM3_snplist.dat is also in this location.
# We assume that the SNPID column in the VCF files contains RSids if this is not the case we can also provide ch:bp locations.
 
# Begin script step 0

# NOTE all tools need to be made executable.
chmod +x ./plink2
chmod +x ./gcta
chmod +x ./bgenix
chmod +x ./tabix
chmod +x ./king
chmod +x ./qctool220

if [ ! -e Output ]; then
    mkdir Output
fi

#Define the phenotype file
PHENO="/Path/To/Directory/Phenotype_File.dat"

# Select the SNPs for the GRM for the autosomal Chromosomes

for chr in $(seq 1 22); do
INPUTF1=/Path/To/Directory/HRC_chr$chr.vcf.gz
OUTPUT1=${OUTPUT0}${chr}
# change --vcf to the appropriate option for your input genotype file type, e.g "--bfile" for plink format.
./plink2 --vcf $INPUTF1 --extract HM3_snplist.dat --make-bed --out $OUTPUT1
done

# Merge all the chromosomes together
ls ${OUTPUT0}*.bim | sed 's/.bim//' > mergelist.dat
OUTPUT2=./Output/HM3_SNPS
./plink2 --merge-list mergelist.dat --make-bed --out ${OUTPUT2}

# Select the fam and Id id from the phenotype file
# Assuming NA is mising value and all people are also having the covariates.
cat $PHENO | grep -v "NA" | awk '{print $1" "$2}' > WithPheno.dat		
# Select only the participants with phenotype data from the plink files containing the SNPs for the GRM
# Note make sure that the FAM & IDs of the VCF file match the IDs of the phenotype and covariate files.
# Also make sure that these IDs match in the generated plink files.
 
OUTPUT3=./Output/Phenotyped_GT_set
./plink2 --bfile ${OUTPUT2} --keep WithPheno.dat --maf 0.01 --hwe 0.0001 --make-bed --out ${OUTPUT3}

# End script step 0

# Output files:
# HM3_SNPS.bim/bed/fam: Plink file containing the SNPs for the GRM for all participants of your cohort
# Phenotyped_GT_set.bim/bed/fam": Plink files containing the SNPs for the GRM for all participants with phenotype data, required for script step 1
