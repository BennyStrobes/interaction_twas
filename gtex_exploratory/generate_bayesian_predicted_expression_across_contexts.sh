#!/bin/bash
#SBATCH -t 0-60:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=40GB 



expr_pc_file="${1}"
expression_file="${2}"
genotype_indices_file="${3}"
processed_genotype_data_dir="${4}"
onek_genomes_plink_files="${5}"
predicted_expression_dir="${6}"
num_pcs="${7}"

source ~/.bashrc
conda activate plink_env

date

python generate_bayesian_predicted_expression_across_contexts.py $expr_pc_file $expression_file $genotype_indices_file $processed_genotype_data_dir $onek_genomes_plink_files $predicted_expression_dir ${num_pcs}

date