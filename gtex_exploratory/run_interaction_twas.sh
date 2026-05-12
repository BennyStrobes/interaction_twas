#!/bin/bash
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=25GB 



pred_expresssion_output_stem="${1}"
trait_sumstat_file="${2}"
onek_genomes_plink_files="${3}"
interaction_twas_output_stem="${4}"
inf_version="${5}"
tau="${6:-1.0}"

source ~/.bashrc
conda activate plink_env

date

echo ${interaction_twas_output_stem}
echo ${tau}

python run_interaction_twas.py $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem $inf_version $tau


date
