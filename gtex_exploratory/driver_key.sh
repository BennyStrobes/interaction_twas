#################
# Input data
#################

# Directory containing genotype data
processed_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"

# Gtex cross tissue tpm expression
gtex_tpm_expression="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"

# Gtex subject attributes (has ancestry)
gtex_subject_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/genotype_dbgap_download/phs000424.v11.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"

# GTEx sample attributes (has mapping from sample ID to tissue name)
gtex_sample_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/gtex_sample_attributes/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"

# Tissue names
cross_tissue_tissue_names_file="/lab-share/CHIP-Strober-e2/Public/ben/interaction_twas/input_data/10_tissue_list.txt"
whole_blood_tissue_names_file="/lab-share/CHIP-Strober-e2/Public/ben/interaction_twas/input_data/whole_blood_tissue_list.txt"


sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/sumstats/"

# Plink files for hg38 one-thousands genomes snps
onek_genomes_plink_files="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/"

# GTEx per tissue covariate file
gtex_per_tissue_covariate_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/per_tissue_covariates/"


#################
# Output data
#################

# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/interaction_twas/"

# Processed expression directory
processed_eqtl_data_dir=${output_root}"processed_eqtl_data/"

# Processed expression directory
predicted_expression_dir=${output_root}"predicted_expression/"

# Processed expression directory
bayesian_predicted_expression_dir=${output_root}"bayesian_predicted_expression/"


# interaction twas results directory
interaction_twas_results_dir=${output_root}"interaction_twas/"


#################
# Code
#################
# Process expression data
if false; then
tissue_version="ten_tissue"
tissue_names_file=${cross_tissue_tissue_names_file}
specific_processed_eqtl_output_stem=$processed_eqtl_data_dir${tissue_version}
sbatch preprocess_eqtl_data_for_analysis.sh $tissue_names_file $gtex_v10_pc_genes_gtf $gtex_tpm_expression $gtex_subject_attributes_file $gtex_sample_attributes_file $processed_genotype_data_dir $gtex_per_tissue_covariate_dir $specific_processed_eqtl_output_stem

tissue_version="whole_blood"
tissue_names_file=${whole_blood_tissue_names_file}
specific_processed_eqtl_output_stem=$processed_eqtl_data_dir${tissue_version}
sbatch preprocess_eqtl_data_for_analysis.sh $tissue_names_file $gtex_v10_pc_genes_gtf $gtex_tpm_expression $gtex_subject_attributes_file $gtex_sample_attributes_file $processed_genotype_data_dir $gtex_per_tissue_covariate_dir $specific_processed_eqtl_output_stem
fi



# Run model to predict expression from interaction model
if false; then
for chrom_num in {1..22}
do
for tissue_version in "whole_blood" "ten_tissue"
do
num_pcs="5"
expr_pc_file=${processed_eqtl_data_dir}${tissue_version}"_expression_pcs_top_50.txt"
expression_file=${processed_eqtl_data_dir}${tissue_version}"_inverse_normal_transformed_expression.txt"
genotype_indices_file=${processed_eqtl_data_dir}${tissue_version}"_genotype_mapping_to_expression_samples.txt"
output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
sbatch generate_predicted_expression_across_contexts.sh $chrom_num $expr_pc_file $expression_file $genotype_indices_file $processed_genotype_data_dir $onek_genomes_plink_files $output_stem $num_pcs
done
done
fi


# Run model to predict expression from interaction model
if false; then
tissue_version="ten_tissue"

num_pcs="5"
expr_pc_file=${processed_eqtl_data_dir}${tissue_version}"_expression_pcs_top_50.txt"
expression_file=${processed_eqtl_data_dir}${tissue_version}"_inverse_normal_transformed_expression.txt"
genotype_indices_file=${processed_eqtl_data_dir}${tissue_version}"_genotype_mapping_to_expression_samples.txt"
output_stem=$bayesian_predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
sbatch generate_bayesian_predicted_expression_across_contexts.sh $expr_pc_file $expression_file $genotype_indices_file $processed_genotype_data_dir $onek_genomes_plink_files $output_stem $num_pcs

tissue_version="whole_blood"

num_pcs="5"
expr_pc_file=${processed_eqtl_data_dir}${tissue_version}"_expression_pcs_top_50.txt"
expression_file=${processed_eqtl_data_dir}${tissue_version}"_inverse_normal_transformed_expression.txt"
genotype_indices_file=${processed_eqtl_data_dir}${tissue_version}"_genotype_mapping_to_expression_samples.txt"
output_stem=$bayesian_predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
sbatch generate_bayesian_predicted_expression_across_contexts.sh $expr_pc_file $expression_file $genotype_indices_file $processed_genotype_data_dir $onek_genomes_plink_files $output_stem $num_pcs
fi






if false; then
# Interaction TWAS model
trait_name="UKB_460K.blood_MONOCYTE_COUNT"
tissue_version="ten_tissue"
num_pcs="5"
inf_version="joint_genes"
tau="2.0"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"
interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_"${inf_version}"_"${tau}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem $inf_version $tau

trait_name="UKB_460K.blood_MONOCYTE_COUNT"
tissue_version="ten_tissue"
num_pcs="5"
inf_version="independent_genes"
tau="2.0"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"
interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_"${inf_version}"_"${tau}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem $inf_version $tau




# Interaction TWAS model
trait_name="UKB_460K.blood_MONOCYTE_COUNT"
tissue_version="ten_tissue"
num_pcs="5"
inf_version="joint_genes"
tau="1.5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"
interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_"${inf_version}"_"${tau}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem $inf_version $tau

trait_name="UKB_460K.blood_MONOCYTE_COUNT"
tissue_version="ten_tissue"
num_pcs="5"
inf_version="independent_genes"
tau="1.5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"
interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_"${inf_version}"_"${tau}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem $inf_version $tau


# Interaction TWAS model
trait_name="UKB_460K.blood_MONOCYTE_COUNT"
tissue_version="ten_tissue"
num_pcs="5"
inf_version="joint_genes"
tau="3.0"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"
interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_"${inf_version}"_"${tau}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem $inf_version $tau

trait_name="UKB_460K.blood_MONOCYTE_COUNT"
tissue_version="ten_tissue"
num_pcs="5"
inf_version="independent_genes"
tau="3.0"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"
interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_"${inf_version}"_"${tau}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem $inf_version $tau
fi




if false; then

trait_name="UKB_460K.blood_MONOCYTE_COUNT"
tissue_version="whole_blood"
num_pcs="5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"

interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem


trait_name="UKB_460K.biochemistry_LDLdirect"
tissue_version="ten_tissue"
num_pcs="5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"

interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem



trait_name="UKB_460K.disease_AID_ALL"
tissue_version="ten_tissue"
num_pcs="5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"

interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem

trait_name="UKB_460K.body_WHRadjBMIz"
tissue_version="ten_tissue"
num_pcs="5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"

interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem


trait_name="UKB_460K.lung_FEV1FVCzSMOKE"
tissue_version="ten_tissue"
num_pcs="5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"

interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem


trait_name="UKB_460K.disease_HYPERTENSION_DIAGNOSED"
tissue_version="ten_tissue"
num_pcs="5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"

interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem


trait_name="UKB_460K.biochemistry_HbA1c"
tissue_version="ten_tissue"
num_pcs="5"

pred_expresssion_output_stem=$predicted_expression_dir${tissue_version}"_"${num_pcs}"_pcs"
trait_sumstat_file=${sumstats_dir}${trait_name}".sumstats"

interaction_twas_output_stem=${interaction_twas_results_dir}${trait_name}"_"${tissue_version}"_"${num_pcs}"_interaction_twas_results"

sbatch run_interaction_twas.sh $pred_expresssion_output_stem $trait_sumstat_file $onek_genomes_plink_files $interaction_twas_output_stem
fi
