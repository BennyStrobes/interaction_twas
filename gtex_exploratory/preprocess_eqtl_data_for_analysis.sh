#!/bin/bash
#SBATCH -t 0-7:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 



tissue_names_file="${1}"
gtex_v10_pc_genes_gtf="${2}"
gtex_tpm_expression="${3}"
gtex_subject_attributes_file="${4}"
gtex_sample_attributes_file="${5}"
processed_genotype_data_dir="${6}"
gtex_per_tissue_covariate_dir="${7}"
processed_eqtl_data_dir="${8}"


source ~/.bashrc
conda activate plink_env


python preprocess_eqtl_data_for_analysis.py $tissue_names_file $gtex_v10_pc_genes_gtf $gtex_tpm_expression $gtex_subject_attributes_file $gtex_sample_attributes_file $processed_genotype_data_dir $gtex_per_tissue_covariate_dir $processed_eqtl_data_dir