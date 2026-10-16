#!/bin/bash
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=30GB 


expr_pc_file="${1}"
expression_file="${2}"
genotype_indices_file="${3}"
training_sample_indices_file="${4}"
processed_genotype_data_dir="${5}"
onek_genomes_plink_files="${6}"
predicted_expression_dir="${7}"

echo ${predicted_expression_dir}

source ~/.bashrc
conda activate plink_env

date


python generate_bayesian_predicted_expression.py $expr_pc_file $expression_file $genotype_indices_file $training_sample_indices_file $processed_genotype_data_dir $onek_genomes_plink_files $predicted_expression_dir

date