#!/bin/bash

######################################################################
##SET UP ENVIRONMENT

module load singularity

cohorts=("invnorm_max_lv_indexed" "invnorm_min_lv_indexed" "invnorm_lvef" "Khurshid_2023_v19_seg_lvmi_adjusted.bolt.imputed.filtered")
sample_n=(40000 40000 40000 43230) 

##CREATE VARIABLES FOR LOOP
# LSB_JOBINDEX starts with 1 so you will need to subtract one to get the array
# indexes to work correctly
IDX=$((LSB_JOBINDEX-1))
FILE_IDX=$((IDX))

tissues=('Artery_Aorta' 'Artery_Coronary' 'Artery_Tibial' 'Heart_Atrial_Appendage' 'Whole_Blood' 'Heart_Left_Ventricle' 'Adipose_Visceral_Omentum' 'Adipose_Subcutaneous' 'Liver' 'Kidney_Cortex')

METAXCAN=filepath/MetaXcan/software 
DATA=filepath/MVP_cardio 
RESULTS=filepath/MetaXcan_out
GWAS=filepath/MVP_cardio

export METAXCAN=filepath/software
export DATA=filepath/MVP_cardio
export RESULTS=filepath/MetaXcan_out
export GWAS=filepath/MVP_cardio

printf "Run association\n\n"

module load python

#run harmonization
python $DATA/summary-gwas-imputation/src/gwas_parsing.py \
-gwas_file $GWAS/HF/MRI_traits/${cohorts[$FILE_IDX]}.tsv.gz \
-snp_reference_metadata $DATA/summary-gwas-imputation/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map SNP variant_id \
-output_column_map ALLELE1 non_effect_allele \
-output_column_map ALLELE0 effect_allele \
-output_column_map BETA effect_size \
-output_column_map SE standard_error \
-output_column_map P_LINREG pvalue \
-output_column_map CHR chromosome \
--chromosome_format \
-output_column_map BP position \
-output_column_map A0FREQ frequency \
--insert_value sample_size ${sample_n[$FILE_IDX]} \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue effect_size standard_error sample_size \
-output $DATA/harmonized_gwas/MRI_${cohorts[$FILE_IDX]}.harmonized.txt.gz
 

#run TWAS
for tissue in "${tissues[@]}"; do singularity exec -B /project ./metaxcan.sif /opt/conda/bin/python /app/MetaXcan/software/SPrediXcan.py \
--model_db_path $DATA/eqtl/mashr/mashr_$tissue.db \
--covariance $DATA/eqtl/mashr/mashr_$tissue.txt.gz \
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

