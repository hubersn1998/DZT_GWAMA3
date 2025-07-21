# Step 2: The GWAS
# Definitions:
# change variable COHORTNAME to the one of your cohort or study
# change variable GWASTYPE to DZ if your phenotype was based on DZ twins instead of mothers of a DZ twin
yymmdd=$(date "+%y%m%d")
COHORTNAME=NAME_OF_YOUR_COHORT
GWASTYPE=MODZ
NTHREAD=YOUR_SETTINGS
# Input files:

# -Phenotype_File.dat: a text file containing the samples with three columns: Family_ID, Individual_ID and Phenotype. 
# The phenotype should be coded as 0= controls and 1 = case. The file should not contain a header. 
# -Covar_File.dat: A text file containing the columns with the Family_ID, Individual_ID, and all categorial covariates
# -QCovar_File.dat: A text file containing the columns with the Family_ID, Individual_ID, and all quantitative covariates
# -Phenotyped_GT_set.bim/bed/fam": Plink files containing the SNPs for the GRM for all participants with phenotype data, obtained through script step 0.
# -Phenotyped_GRM_sparse: File containing the GRM, output of step1
# Remember that sample names need to match in all files! 
# If needed to check and fix the sample names, remove the #'s starting at "ORGFam" and end after "sample"}'"

#ORGFam="{your data.fam}"	
#^ <- This needs to be a PLink family file which has sex defined in there e.g. from the original genotypes prior to imputation.
#^ <- in case no sex information is available you could also use a famele dummy, only the logs male / female counts won't be ok then.

#for chr in $(seq 1 22); do
#ORGsample="/Path/To/directory/"YOUR_FILE_NAME"$chr".sample"
#NEWsample="NEW_NAME"$chr".sample"
#head -2 $ORGsample | sed 's/ID_1 ID_2 missing/ID_1 ID_2 missing sex/' | sed s'/0 0 0/0 0 0 D/' > $NEWsample
#awk '{print $1"_"$2" "$5}' $ORGFam | sort -k1 > _tmp2.dat
#awk '{print $1" "$3" "NR}' $ORGsample | sed '1,2d' | sort -k1 > _tmp1.dat
#join -1 1 -2 1 _tmp1.dat _tmp2.dat | sort -n -k3 | awk '{print $1" "$2" "$4}' | sed 's/_/ /' >> $NEWsample
# check
#paste $ORGsample $NEWsample | tr '\t' ' ' | sed '1,2d' | awk '{print $1" "$4"_"$5}' | awk '{if ($1!=$2) print "Check files fam and sample"}'
#done

#chr=X
#ORGsample="/Path/To/directory/"YOUR_FILE_NAME"$chr".sample"
#NEWsample="NEW_NAME"$chr".sample"
#head -2 $ORGsample | sed 's/ID_1 ID_2 missing/ID_1 ID_2 missing sex/' | sed s'/0 0 0/0 0 0 D/' > $NEWsample
#awk '{print $1"_"$2" "$5}' $ORGFam | sort -k1 > _tmp2.dat
#awk '{print $1" "$3" "NR}' $ORGsample | sed '1,2d' | sort -k1 > _tmp1.dat
#join -1 1 -2 1 _tmp1.dat _tmp2.dat | sort -n -k3 | awk '{print $1" "$2" "$4}' | sed 's/_/ /' >> $NEWsample
# check
#paste $ORGsample $NEWsample | tr '\t' ' ' | sed '1,2d' | awk '{print $1" "$4"_"$5}' | awk '{if ($1!=$2) print "Check files fam and sample"}'



# Input genotypes
# -Genotype$chr.bgen/sample= Bgen files of the imputed snps for all autosomal chromosomes and the X chromosome
# and /or
# -Genotype$chr.bgen/sample= Bgen files of the imputed snps for all autosomal chromosomes and the X chromosome

# Convert VCF to BGEN if needed, if so remove the #'s starting at "for", ending after "./qctool220.
#for chr in $(seq 1 22); do
#INPUTF1="HRC_CHR"$chr".vcf.gz"
#OUTPUT2="HRC_CHR"$chr
#./tabix -p vcf $INPUTF1
#./qctool220 -g $INPUTF1 -vcf-genotype-field GP -bgen-permitted-input-rounding-error 0.01 -bgen-bits 8 -ofiletype bgen_v1.2 -og $OUTPUT2".bgen" -os $OUTPUT2".sample"
#./bgenix -g $OUTPUT2".bgen" -index
#done
#chr=X
#INPUTF1="HRC_CHR"$chr".vcf.gz"
#OUTPUT2="HRC_CHR"$chr
#./qctool220 -g $INPUTF1 -vcf-genotype-field GP -bgen-permitted-input-rounding-error 0.01 -bgen-bits 8 -ofiletype bgen_v1.2 -og $OUTPUT2".bgen" -os $OUTPUT2".sample"



# Begin script step 2


# GWAS
GTbase="/Path/To/Directory/Phenotyped_GT_set"

# Null-model using independent genotyped SNPs
PHENO=./Path/To/Directory/PhenoFast.dat
QCOV=./Path/To/Directory/QCovFast.dat
COV=./Path/To/Directory/CovFast.dat
OUT0=./Path/To/Directory/Twin_PCs_MD
GRM=./Path/To/Directory/Phenotyped_GRM_sparse
./gcta64 --bfile $GTbase --grm-sparse ${GRM} --est-vg HE --model-only --fastGWA-mlm-binary --pheno ${PHENO} --qcovar ${QCOV} --covar ${COV} --thread-num ${NTHREAD} --maf 0.0005 --out ${OUT0}

#GWAS on autosomal chromosomes
for chr in $(seq 1 22); do
BGEN="HRC_CHR"$chr".bgen"
SAMPLE="HRC_CHR"$chr".sample"
OUT1="TWIN_PCs_CHR"$chr
./gcta64 --bgen $BGEN --sample $SAMPLE --load-model ${OUT0}.fastGWA --thread-num ${NTHREAD} --out $OUT1
done

#GWAS on X chromosomes
chr=X
BGEN="HRC_CHR"$chr".bgen"
SAMPLE="HRC_CHR"$chr".sample"
OUT1="TWIN_PCs_CHR"$chr
./gcta64 --bgen $BGEN --sample $SAMPLE --load-model ${OUT0}.fastGWA --thread-num ${NTHREAD} --out $OUT1

# Combine all the output files 
c=1
head -1 ./Output/TWIN_PCs_CHR${c}.fastGWA >./Output/TWIN_AllCHR_intermediate.dat
for c in $(echo $(seq 1 22) X); do
    tail -n +2 ./Output/TWIN_PCs_CHR${c}.fastGWA
done | sort | uniq >>./Output/TWIN_AllCHR_intermediate.dat
mv ./Output/TWIN_AllCHR_intermediate.dat ./Output/${yymmdd}_TWIN_AllCHR_${COHORTNAME}_${GWASTYPE}.dat
gzip -f ./Output/${yymmdd}_TWIN_AllCHR_${COHORTNAME}_${GWASTYPE}.dat

# End script step 2

#Output:
#./Output/${yymmdd}_TWIN_AllCHR_${COHORTNAME}_${GWASTYPE}.dat.gz": GWAS results to be send to the consortium

# Not if you have both MODZT as well as a DZT set, then repeat all scripts twice using different
# Phenotype and covariate files.

# The following sections can then be commented out with "#" or removed from the scripts.

# STEP 0:
#for chr in $(seq 1 22); do
#INPUTF1=/Path/To/Directory/HRC_chr$chr.vcf.gz
#OUTPUT1=${OUTPUT0}${chr}
# change --bfile to the appropriate option for your input genotype file type, e.g "--vcf" for VCF format.
#./plink2 --bfile $INPUTF1 --extract HM3_snplist.dat --make-bed --out $OUTPUT1
#done

# Merge all the chromosomes together
#ls ${OUTPUT0}*.bim | sed 's/.bim//' > mergelist.dat
#OUTPUT2=./Output/HM3_SNPS
#./plink2 --merge-list mergelist.dat --make-bed --out ${OUTPUT2}

# STEP 2 once generated the BGEN files do not need to be regenerated.

