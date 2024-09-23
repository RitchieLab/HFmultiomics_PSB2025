#!/bin/bash

############################################
##SET UP ENVIRONMENT

module load singularity

cohorts=("invnorm_max_lv_indexed" "invnorm_min_lv_indexed" "invnorm_lvef" "Khurshid_2023_v19_seg_lvmi_adjusted.bolt.imputed.filtered")

##CREATE VARIABLES FOR LOOP
# LSB_JOBINDEX starts with 1 so you will need to subtract one to get the array
# indexes to work correctly
IDX=$((LSB_JOBINDEX-1))
FILE_IDX=$((IDX))

tissues=('ARIC_EA_hg38' 'ARIC_AA_hg38')

METAXCAN=filepath/software 
DATA=filepath/MVP_cardio 
RESULTS=filepath/PWAS/output
GWAS=filepath/MVP_cardio

export METAXCAN=filepath/software
export DATA=filepath/MVP_cardio
export RESULTS=filepath/PWAS/output
export GWAS=filepath/MVP_cardio

printf "Run association\n\n"

module load python

#run TWAS
for tissue in "${tissues[@]}"; do singularity exec -B /project /project/ritchie/projects/CodeWorks_Projects/metaxcan.sif /opt/conda/bin/python /app/MetaXcan/software/SPrediXcan.py \
--model_db_path ~/group/personal/rasika/Cardio_WAS/PWAS/ARIC_Models/ARIC-selected/$tissue.db \
--covariance ~/group/personal/rasika/Cardio_WAS/PWAS/ARIC_Models/ARIC-selected/$tissue.txt.gz \
--model_db_snp_key varID \
--gwas_file $DATA/harmonized_gwas/MRI_${cohorts[$FILE_IDX]}.harmonized.txt.gz \
--snp_column panel_variant_id \
--chromosome_column chromosome \
--position_column position \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--pvalue_column pvalue \
--keep_non_rsid \
--additional_output \
--throw \
--output_file $RESULTS/MRI_${cohorts[$FILE_IDX]}.$tissue.spredixcan.txt; done

