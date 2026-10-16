########################
# Input data
########################


# Processed expression directory
processed_eqtl_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/interaction_twas/processed_eqtl_data/"

# Directory containing genotype data
processed_eqtl_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Sumstats
sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/sumstats/"

# Plink files for hg38 one-thousands genomes snps
onek_genomes_plink_files="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/"

# GTEx per tissue covariate file
gtex_per_tissue_covariate_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/per_tissue_covariates/"

# Gotten from https://github.com/BennyStrobes/qtl_mediated_h2/tree/main/old/gtex_downsampling/tglr
tglr_code_dir="/lab-share/CHIP-Strober-e2/Public/ben/tglr/"

# LDSC code dir (Potetnially want up to date version)
ldsc_code_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/ldsc/"

# Hapmap3 snps used for ldsc regression snps
hapmap3_rsid_file="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/w_hm3.noMHC.snplist"

# BaselineLD Dir
baselineLD_anno_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"



#################
# Output data
#################

# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/genome_wide_expression_prediction/"

# Directory containing down-sampling summary
down_sampling_summary_dir=${output_root}"down_sampling_summary/"

# Results of Bayesian Inference analysis
bayesian_model_inference_dir=${output_root}"bayesian_model_inference/"



#################
# Code
#################

# Focus on Muscle-Skeletal tissue for primary analysis
# tissue-specific input data
tissue_name="muscle_skeletal"
expression_file=${processed_eqtl_data_dir}${tissue_name}"_inverse_normal_transformed_expression.txt"
expr_pc_file=${processed_eqtl_data_dir}${tissue_name}"_expression_pcs_top_50.txt"
genotype_indices_file=${processed_eqtl_data_dir}${tissue_name}"_genotype_mapping_to_expression_samples.txt"





####################
# Organize down-sampling and folds
n_folds="5"
downsampling_fold_sample_summary_file=$down_sampling_summary_dir${tissue_name}"_downsampling_and_fold_summary.txt"
if false; then
source ~/.bashrc
conda activate plink_env
python organize_downsampling_and_folds_for_single_tissue.py $genotype_to_expression_mapping_file $n_folds $downsampling_fold_sample_summary_file $down_sampling_summary_dir${tissue_name}"_"
fi



####################
# Predict expression using Bayesian approach
if false; then
tail -n +2 $downsampling_fold_sample_summary_file | while IFS=$'\t' read -r downsample_fraction test_fold_iter n_training_samples n_test_samples training_sample_indices_file test_sample_indices_file
do
	if [[ $test_fold_iter != "0" ]]; then
		continue
	fi

	output_stem=$bayesian_model_inference_dir${tissue_name}"_downsample_"${downsample_fraction}"_testfold_"${test_fold_iter}

	sbatch generate_bayesian_predicted_expression.sh $expr_pc_file $expression_file $genotype_indices_file $training_sample_indices_file $processed_eqtl_genotype_data_dir $onek_genomes_plink_files $output_stem
done

fi





























