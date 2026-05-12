#!/bin/bash
#SBATCH -t 0-44:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 



chrom_num="${1}"
expr_pc_file="${2}"
expression_file="${3}"
genotype_indices_file="${4}"
processed_genotype_data_dir="${5}"
onek_genomes_plink_files="${6}"
predicted_expression_dir="${7}"
num_pcs="${8}"

source ~/.bashrc
conda activate plink_env

date

python generate_predicted_expression_across_contexts.py $chrom_num $expr_pc_file $expression_file $genotype_indices_file $processed_genotype_data_dir $onek_genomes_plink_files $predicted_expression_dir ${num_pcs}

date