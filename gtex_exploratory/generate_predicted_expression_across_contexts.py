import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink
import time
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_predict, KFold
from sklearn.metrics import r2_score
import scipy.stats as stats



def extract_dictionary_list_of_1_kg_variants(filer):
	dicti = {}
	f = open(filer)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		variant_name_1 = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5] + '_b38'
		variant_name_2 = 'chr' + data[0] + '_' + data[3] + '_' + data[5] + '_' + data[4] + '_b38'

		dicti[variant_name_1] = 1
		dicti[variant_name_2] = 1

	f.close()

	return dicti



def load_in_expression_pcs(expr_pc_file):
	tmp = np.loadtxt(expr_pc_file, dtype=str, delimiter='\t')
	sample_names = tmp[1:,0]
	expr_pcs = tmp[1:,1:].astype(float)
	return sample_names, expr_pcs

def extract_expression_sample_names(expression_file):
	f = open(expression_file)
	header = f.readline().rstrip().split('\t')
	f.close()
	return np.asarray(header[6:])

def check_sample_order_and_dimensions(expression_sample_names, pc_sample_names, expression_pcs, genotype_indices):
	if np.array_equal(expression_sample_names, pc_sample_names) == False:
		print('Expression sample names and expression PC sample names are not ordered the same')
		pdb.set_trace()
	if len(expression_sample_names) != expression_pcs.shape[0]:
		print('Number of expression samples does not match number of expression PC rows')
		pdb.set_trace()
	if len(expression_sample_names) != len(genotype_indices):
		print('Number of expression samples does not match number of genotype indices')
		pdb.set_trace()
	return

def filter_and_standardize_snp_genotype_matrix(geno_mat_donor, gene_bim, genotype_indices):
	unique_genotype_indices = np.unique(genotype_indices)
	geno_mat_observed_donors = geno_mat_donor[:, unique_genotype_indices]
	row_means = np.nanmean(geno_mat_observed_donors, axis=1)
	valid_snp_rows = np.isfinite(row_means)
	geno_mat_donor = geno_mat_donor[valid_snp_rows, :]
	geno_mat_observed_donors = geno_mat_observed_donors[valid_snp_rows, :]
	gene_bim = gene_bim[valid_snp_rows, :]
	row_means = row_means[valid_snp_rows]

	nan_rows, nan_cols = np.where(np.isnan(geno_mat_donor))
	geno_mat_donor[nan_rows, nan_cols] = row_means[nan_rows]
	nan_rows, nan_cols = np.where(np.isnan(geno_mat_observed_donors))
	geno_mat_observed_donors[nan_rows, nan_cols] = row_means[nan_rows]

	snp_stds = np.std(geno_mat_observed_donors, axis=1)
	valid_snp_rows = np.isfinite(snp_stds) & (snp_stds > 0)
	geno_mat_donor = geno_mat_donor[valid_snp_rows, :]
	gene_bim = gene_bim[valid_snp_rows, :]
	snp_stds = snp_stds[valid_snp_rows]
	row_means = row_means[valid_snp_rows]

	geno_mat_donor = (geno_mat_donor - row_means[:, None])/snp_stds[:, None]
	geno_mat = geno_mat_donor[:, genotype_indices]
	return geno_mat, gene_bim

def extract_boolean_vector_of_whether_gtex_snp_id_in_1kg(snp_ids, one_kg_variant_dicti):
	boolers = []
	for snp_id in snp_ids:
		if snp_id in one_kg_variant_dicti:
			boolers.append(True)
		else:
			boolers.append(False)
	return np.asarray(boolers)

def create_donor_aware_folds(expression_sample_names, n_folds):
	donor_ids = []
	for sample_name in expression_sample_names:
		donor_ids.append(sample_name.split(':')[0])
	donor_ids = np.asarray(donor_ids)
	unique_donor_ids = np.unique(donor_ids)
	unique_donor_ids = np.random.RandomState(1).permutation(unique_donor_ids)
	donor_folds = np.array_split(unique_donor_ids, n_folds)

	sample_folds = []
	for donor_fold in donor_folds:
		sample_folds.append(np.isin(donor_ids, donor_fold))
	return sample_folds

def construct_context_interaction_design_matrices(geno_mat, expression_pcs):
	genotype_main_effects = np.transpose(geno_mat)
	interaction_terms = []
	for pc_num in range(expression_pcs.shape[1]):
		interaction_terms.append(genotype_main_effects*expression_pcs[:, pc_num][:, None])
	interaction_terms = np.hstack(interaction_terms)
	return genotype_main_effects, interaction_terms

def fit_ridge_model_with_two_penalties(genotype_main_effects, interaction_terms, expr_vec, lambda_main, lambda_interaction, train_indices, test_indices):
	penalized_design = np.hstack((genotype_main_effects, interaction_terms))
	train_stds = np.std(penalized_design[train_indices, :], axis=0)
	valid_penalized_columns = np.isfinite(train_stds) & (train_stds > 0)
	penalized_design = penalized_design[:, valid_penalized_columns]

	num_main_effects = genotype_main_effects.shape[1]
	raw_penalty_weights = np.hstack((np.ones(num_main_effects)*lambda_main, np.ones(interaction_terms.shape[1])*lambda_interaction))
	penalty_weights = raw_penalty_weights[valid_penalized_columns]

	train_design = penalized_design[train_indices, :]
	test_design = penalized_design[test_indices, :]
	xtx = np.dot(np.transpose(train_design), train_design)
	xty = np.dot(np.transpose(train_design), expr_vec[train_indices])
	beta = np.linalg.solve(xtx + np.diag(penalty_weights), xty)
	genetic_pred = np.dot(test_design, beta)
	return genetic_pred

def fit_main_effect_ridge_on_full_data(genotype_main_effects, expr_vec, alpha):
	model = Ridge(alpha=alpha, fit_intercept=True)
	model.fit(genotype_main_effects, expr_vec)
	return model.intercept_, model.coef_

def calculate_correlation(observed, predicted):
	if np.std(observed) == 0 or np.std(predicted) == 0:
		return np.nan
	return np.corrcoef(observed, predicted)[0,1]

def select_main_effect_alpha_with_cv(genotype_main_effects, expr_vec, alphas, kf):
	best_alpha = None
	best_correlation = -np.inf
	for alpha in alphas:
		model = Ridge(alpha=alpha, fit_intercept=True)
		pred = cross_val_predict(model, genotype_main_effects, expr_vec, cv=kf)
		corry = calculate_correlation(expr_vec, pred)
		if corry > best_correlation:
			best_correlation = corry
			best_alpha = alpha
	return best_alpha, best_correlation

def select_two_penalty_alphas_with_cv(genotype_main_effects, interaction_terms, expr_vec, alphas, alpha_inters, kf):
	best_alpha_main = None
	best_alpha_interaction = None
	best_correlation = -np.inf
	for alpha_main in alphas:
		for alpha_interaction in alpha_inters:
			pred = np.zeros(len(expr_vec))
			for train_indices, test_indices in kf.split(expr_vec):
				pred[test_indices] = fit_ridge_model_with_two_penalties(genotype_main_effects, interaction_terms, expr_vec, alpha_main, alpha_interaction, train_indices, test_indices)
			corry = calculate_correlation(expr_vec, pred)
			if corry > best_correlation:
				best_correlation = corry
				best_alpha_main = alpha_main
				best_alpha_interaction = alpha_interaction
	return best_alpha_main, best_alpha_interaction, best_correlation

def estimate_holdout_prediction_performance(genotype_main_effects, interaction_terms, expr_vec, alphas, alpha_inters):
	outer_kf = KFold(n_splits=5, shuffle=True, random_state=2)
	outer_train_indices, outer_test_indices = next(outer_kf.split(expr_vec))
	inner_kf = KFold(n_splits=4, shuffle=True, random_state=1)

	outer_train_main_effects = genotype_main_effects[outer_train_indices, :]
	outer_train_interactions = interaction_terms[outer_train_indices, :]
	outer_train_expr = expr_vec[outer_train_indices]
	best_holdout_alpha, inner_main_correlation = select_main_effect_alpha_with_cv(outer_train_main_effects, outer_train_expr, alphas, inner_kf)
	best_holdout_alpha_main, best_holdout_alpha_interaction, inner_interaction_correlation = select_two_penalty_alphas_with_cv(outer_train_main_effects, outer_train_interactions, outer_train_expr, alphas, alpha_inters, inner_kf)

	main_model = Ridge(alpha=best_holdout_alpha, fit_intercept=True)
	main_model.fit(genotype_main_effects[outer_train_indices, :], expr_vec[outer_train_indices])
	main_holdout_pred = main_model.predict(genotype_main_effects[outer_test_indices, :])
	main_holdout_correlation = calculate_correlation(expr_vec[outer_test_indices], main_holdout_pred)

	interaction_holdout_pred = fit_ridge_model_with_two_penalties(genotype_main_effects, interaction_terms, expr_vec, best_holdout_alpha_main, best_holdout_alpha_interaction, outer_train_indices, outer_test_indices)
	interaction_holdout_correlation = calculate_correlation(expr_vec[outer_test_indices], interaction_holdout_pred)

	return main_holdout_correlation, interaction_holdout_correlation, best_holdout_alpha, best_holdout_alpha_main, best_holdout_alpha_interaction, inner_main_correlation, inner_interaction_correlation

def fit_ridge_model_with_two_penalties_on_full_data(genotype_main_effects, interaction_terms, expr_vec, lambda_main, lambda_interaction, num_pcs):
	penalized_design = np.hstack((genotype_main_effects, interaction_terms))
	design_stds = np.std(penalized_design, axis=0)
	valid_penalized_columns = np.isfinite(design_stds) & (design_stds > 0)
	valid_design = penalized_design[:, valid_penalized_columns]

	num_main_effects = genotype_main_effects.shape[1]
	raw_penalty_weights = np.hstack((np.ones(num_main_effects)*lambda_main, np.ones(interaction_terms.shape[1])*lambda_interaction))
	penalty_weights = raw_penalty_weights[valid_penalized_columns]

	xtx = np.dot(np.transpose(valid_design), valid_design)
	xty = np.dot(np.transpose(valid_design), expr_vec)
	valid_beta = np.linalg.solve(xtx + np.diag(penalty_weights), xty)
	beta = np.zeros(penalized_design.shape[1])
	beta[valid_penalized_columns] = valid_beta
	main_beta = beta[:num_main_effects]
	interaction_beta = beta[num_main_effects:].reshape((num_pcs, num_main_effects)).T
	return main_beta, interaction_beta

def create_gene_weight_file_stem(predicted_expression_dir, chrom_num, ensamble_id, gene_symbol):
	safe_ensamble_id = ensamble_id.replace(':', '_').replace('/', '_')
	safe_gene_symbol = gene_symbol.replace(':', '_').replace('/', '_')
	return predicted_expression_dir + 'gene_model_weights/chr' + str(chrom_num) + '_' + safe_ensamble_id + '_' + safe_gene_symbol

def extract_bim_alleles(bim):
	if 'a0' in bim.columns and 'a1' in bim.columns:
		return np.asarray(bim['a0']).astype(str), np.asarray(bim['a1']).astype(str)
	if 'allele1' in bim.columns and 'allele2' in bim.columns:
		return np.asarray(bim['allele1']).astype(str), np.asarray(bim['allele2']).astype(str)
	print('Could not find allele columns in BIM dataframe')
	pdb.set_trace()

def save_main_effect_model_weights(output_file, gene_bim, gene_chrom, main_beta):
	t = open(output_file, 'w')
	t.write('snp_id\tchrom\tposition\tbim_allele0\tbim_allele1\tmain_effect\n')
	for snp_num in range(gene_bim.shape[0]):
		t.write(gene_bim[snp_num, 0] + '\t' + gene_chrom + '\t' + gene_bim[snp_num, 1] + '\t' + gene_bim[snp_num, 2] + '\t' + gene_bim[snp_num, 3] + '\t' + str(main_beta[snp_num]) + '\n')
	t.close()
	return

def save_interaction_model_weights(output_file, gene_bim, gene_chrom, main_beta, interaction_beta):
	t = open(output_file, 'w')
	header = ['snp_id', 'chrom', 'position', 'bim_allele0', 'bim_allele1', 'main_effect']
	for pc_num in range(interaction_beta.shape[1]):
		header.append('interaction_effect_PC' + str(pc_num + 1))
	t.write('\t'.join(header) + '\n')
	for snp_num in range(gene_bim.shape[0]):
		data = [gene_bim[snp_num, 0], gene_chrom, gene_bim[snp_num, 1], gene_bim[snp_num, 2], gene_bim[snp_num, 3], str(main_beta[snp_num])]
		for pc_num in range(interaction_beta.shape[1]):
			data.append(str(interaction_beta[snp_num, pc_num]))
		t.write('\t'.join(data) + '\n')
	t.close()
	return


def run_standard_eqtl_calling_for_gene(expr_vec, expression_pcs, geno_mat):
	# geno_mat is SNPs x samples. Standard eQTL scans are marginal SNP tests with PCs as covariates.
	covariates = np.hstack((np.ones((len(expr_vec), 1)), expression_pcs))
	genotype_mat = np.transpose(geno_mat)
	num_samples = len(expr_vec)
	num_covariates = covariates.shape[1]
	degrees_of_freedom = num_samples - num_covariates - 1
	if degrees_of_freedom <= 0:
		print('Not enough samples to run standard eQTL calling')
		pdb.set_trace()

	covariate_xtx = np.dot(np.transpose(covariates), covariates)
	expr_covariate_beta = np.linalg.solve(covariate_xtx, np.dot(np.transpose(covariates), expr_vec))
	genotype_covariate_beta = np.linalg.solve(covariate_xtx, np.dot(np.transpose(covariates), genotype_mat))
	resid_expr = expr_vec - np.dot(covariates, expr_covariate_beta)
	resid_genotype = genotype_mat - np.dot(covariates, genotype_covariate_beta)

	xtx = np.sum(np.square(resid_genotype), axis=0)
	valid_snps = np.isfinite(xtx) & (xtx > 0)
	beta = np.zeros(genotype_mat.shape[1]) + np.nan
	beta_se = np.zeros(genotype_mat.shape[1]) + np.nan
	pvalue = np.zeros(genotype_mat.shape[1]) + np.nan

	beta[valid_snps] = np.dot(np.transpose(resid_genotype[:, valid_snps]), resid_expr)/xtx[valid_snps]
	resid = resid_expr[:, None] - (resid_genotype[:, valid_snps]*beta[valid_snps][None, :])
	residual_variance = np.sum(np.square(resid), axis=0)/degrees_of_freedom
	beta_se[valid_snps] = np.sqrt(residual_variance/xtx[valid_snps])
	t_stat = beta[valid_snps]/beta_se[valid_snps]
	pvalue[valid_snps] = 2.0*stats.t.sf(np.abs(t_stat), degrees_of_freedom)

	return beta, beta_se, pvalue, valid_snps

def run_interaction_qtl_calling_for_gene(expr_vec, expression_pcs, geno_mat):
	# Marginal interaction tests: expr ~ intercept + all PCs + SNP + SNP:PC_k.
	genotype_mat = np.transpose(geno_mat)
	num_samples = len(expr_vec)
	num_pcs = expression_pcs.shape[1]
	interaction_beta = np.zeros((genotype_mat.shape[1], num_pcs)) + np.nan
	interaction_beta_se = np.zeros((genotype_mat.shape[1], num_pcs)) + np.nan
	interaction_pvalue = np.zeros((genotype_mat.shape[1], num_pcs)) + np.nan
	valid_tests = np.zeros((genotype_mat.shape[1], num_pcs), dtype=bool)

	for pc_num in range(num_pcs):
		interaction_terms = genotype_mat*expression_pcs[:, pc_num][:, None]
		for snp_num in range(genotype_mat.shape[1]):
			covariates = np.hstack((np.ones((num_samples, 1)), expression_pcs, genotype_mat[:, snp_num][:, None]))
			degrees_of_freedom = num_samples - covariates.shape[1] - 1
			if degrees_of_freedom <= 0:
				continue
			covariate_xtx = np.dot(np.transpose(covariates), covariates)
			expr_covariate_beta = np.linalg.solve(covariate_xtx, np.dot(np.transpose(covariates), expr_vec))
			interaction_covariate_beta = np.linalg.solve(covariate_xtx, np.dot(np.transpose(covariates), interaction_terms[:, snp_num]))
			resid_expr = expr_vec - np.dot(covariates, expr_covariate_beta)
			resid_interaction = interaction_terms[:, snp_num] - np.dot(covariates, interaction_covariate_beta)
			xtx = np.sum(np.square(resid_interaction))
			if np.isfinite(xtx) == False or xtx <= 0:
				continue
			beta = np.dot(resid_interaction, resid_expr)/xtx
			resid = resid_expr - (resid_interaction*beta)
			residual_variance = np.sum(np.square(resid))/degrees_of_freedom
			beta_se = np.sqrt(residual_variance/xtx)
			t_stat = beta/beta_se
			interaction_beta[snp_num, pc_num] = beta
			interaction_beta_se[snp_num, pc_num] = beta_se
			interaction_pvalue[snp_num, pc_num] = 2.0*stats.t.sf(np.abs(t_stat), degrees_of_freedom)
			valid_tests[snp_num, pc_num] = True

	return interaction_beta, interaction_beta_se, interaction_pvalue, valid_tests

def residualize_and_standardize_expression(expr_vec, expression_pcs):
	covariates = np.hstack((np.ones((len(expr_vec), 1)), expression_pcs))
	covariate_beta = np.linalg.solve(np.dot(np.transpose(covariates), covariates), np.dot(np.transpose(covariates), expr_vec))
	residual_expr_vec = expr_vec - np.dot(covariates, covariate_beta)
	residual_std = np.std(residual_expr_vec)
	if residual_std == 0:
		print('Residualized expression has zero variance')
		pdb.set_trace()
	residual_expr_vec = (residual_expr_vec - np.mean(residual_expr_vec))/residual_std
	return residual_expr_vec

######################
# Command line args
######################
chrom_num = sys.argv[1]
expr_pc_file = sys.argv[2]
expression_file = sys.argv[3]
genotype_indices_file = sys.argv[4]
processed_genotype_data_dir = sys.argv[5]
onek_genomes_plink_files = sys.argv[6]
predicted_expression_dir = sys.argv[7]
num_pcs = int(sys.argv[8])
if predicted_expression_dir.endswith('/') == False:
	predicted_expression_dir = predicted_expression_dir + '/'

# Pick a distance threshold
distance_threshold = 50000.0
min_snps_per_gene=50
n_folds = 5

os.makedirs(predicted_expression_dir, exist_ok=True)
os.makedirs(predicted_expression_dir + 'gene_model_weights/', exist_ok=True)
model_summary_file = predicted_expression_dir + 'chr' + str(chrom_num) + '_gene_model_summary.txt'
model_summary_handle = open(model_summary_file, 'w')
model_summary_handle.write('ensamble_id\tgene_symbol\tchrom\ttss\tn_snps\tbest_eqtl_pvalue\tbest_interaction_qtl_pvalue\tmain_model_holdout_correlation\tinteraction_model_holdout_correlation\tmain_model_alpha\tinteraction_model_alpha_main\tinteraction_model_alpha_interaction\tmain_model_weight_file\tinteraction_model_weight_file\n')

# Extact dictionary list of 1KG variants
one_kg_variant_dicti = extract_dictionary_list_of_1_kg_variants(onek_genomes_plink_files + '1000G.EUR.hg38.' + chrom_num + '.bim')

# Load in expression pcs
pc_sample_names, expression_pcs = load_in_expression_pcs(expr_pc_file)
full_expression_pcs = np.copy(expression_pcs)
expression_pcs = expression_pcs[:, :num_pcs]


# Load in genotype indices
genotype_indices = np.loadtxt(genotype_indices_file).astype(int)

# Check sample ordering and dimensions across expression, PCs, and genotype mapping
expression_sample_names = extract_expression_sample_names(expression_file)
check_sample_order_and_dimensions(expression_sample_names, pc_sample_names, expression_pcs, genotype_indices)

# Load in genotype data
# string of chromosome name
chrom_string = 'chr' + str(chrom_num)
# Load in chromosome plink data
(bim, fam, G) = read_plink(processed_genotype_data_dir + 'gtex_v9_eqtl_chr' + str(chrom_num))
snp_ids = np.asarray(bim['snp'])
snp_pos = np.asarray(bim['pos'])
bim_ref_alleles, bim_alt_alleles = extract_bim_alleles(bim)
bim_mat = np.transpose(np.vstack((snp_ids.astype(str), snp_pos.astype(str), bim_ref_alleles.astype(str), bim_alt_alleles.astype(str))))

# Extract boolean vector of whether gtex snp_ids in 1kg
gtex_snp_id_in_1kg = extract_boolean_vector_of_whether_gtex_snp_id_in_1kg(snp_ids, one_kg_variant_dicti)

# Now loop through genes
f = open(expression_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Extract info for a given gene
	gene_descriptor_fields = np.asarray(data[:6])
	ensamble_id = gene_descriptor_fields[0]
	gene_symbol = gene_descriptor_fields[1]
	gene_chrom = gene_descriptor_fields[2]
	# Skip genes not on desired chrom
	if gene_chrom != chrom_string:
		continue
	if gene_descriptor_fields[5] == '+':
		tss = int(gene_descriptor_fields[3])
	else:
		tss = int(gene_descriptor_fields[4])
	expr_vec = np.asarray(data[6:]).astype(float)


	# Extract valid cis snps for the gene
	valid_cis_snps = (np.abs(snp_pos - tss) <= distance_threshold) & (gtex_snp_id_in_1kg)

	# Filter out genes with too few snps
	if np.sum(valid_cis_snps) < min_snps_per_gene:
		continue

	# Load in genotype data for gene
	geno_mat_donor = G[valid_cis_snps,:].compute()
	geno_mat_donor = 2.0 - geno_mat_donor
	# Extract snp info for gene
	gene_bim = bim_mat[valid_cis_snps,:]
	geno_mat, gene_bim = filter_and_standardize_snp_genotype_matrix(geno_mat_donor, gene_bim, genotype_indices)
	if geno_mat.shape[0] < min_snps_per_gene:
		continue

	# Run standard marginal eQTL calling for this gene with expression PCs as covariates
	beta, beta_se, pvalue, valid_eqtl_snps = run_standard_eqtl_calling_for_gene(expr_vec, expression_pcs, geno_mat)
	interaction_beta, interaction_beta_se, interaction_pvalue, valid_interaction_tests = run_interaction_qtl_calling_for_gene(expr_vec, expression_pcs, geno_mat)
	best_eqtl_pvalue = np.nanmin(pvalue)
	best_interaction_qtl_pvalue = np.nanmin(interaction_pvalue)

	# Residualize expression on expression PCs and standardize residual expression
	residual_expr_vec = residualize_and_standardize_expression(expr_vec, full_expression_pcs)
	genotype_main_effects, interaction_terms = construct_context_interaction_design_matrices(geno_mat, expression_pcs)


	alphas = [10, 100, 1000, 10000]
	alpha_inters = [1000, 10000, 100000, 1000000]
	kf = KFold(n_splits=4, shuffle=True, random_state=1)

	main_holdout_correlation, interaction_holdout_correlation, best_holdout_alpha, best_holdout_alpha_main, best_holdout_alpha_interaction, inner_main_correlation, inner_interaction_correlation = estimate_holdout_prediction_performance(genotype_main_effects, interaction_terms, residual_expr_vec, alphas, alpha_inters)

	best_alpha, best_r2 = select_main_effect_alpha_with_cv(genotype_main_effects, residual_expr_vec, alphas, kf)
	best_alpha_main3, best_alpha_interaction3, best_r23 = select_two_penalty_alphas_with_cv(genotype_main_effects, interaction_terms, residual_expr_vec, alphas, alpha_inters, kf)

	main_model_intercept, main_model_beta = fit_main_effect_ridge_on_full_data(genotype_main_effects, residual_expr_vec, best_alpha)
	interaction_model_main_beta, interaction_model_interaction_beta = fit_ridge_model_with_two_penalties_on_full_data(genotype_main_effects, interaction_terms, residual_expr_vec, best_alpha_main3, best_alpha_interaction3, expression_pcs.shape[1])

	gene_weight_file_stem = create_gene_weight_file_stem(predicted_expression_dir, chrom_num, ensamble_id, gene_symbol)
	main_model_weight_file = gene_weight_file_stem + '_main_effect_ridge_weights.txt'
	interaction_model_weight_file = gene_weight_file_stem + '_interaction_ridge_weights.txt'
	save_main_effect_model_weights(main_model_weight_file, gene_bim, gene_chrom, main_model_beta)
	save_interaction_model_weights(interaction_model_weight_file, gene_bim, gene_chrom, interaction_model_main_beta, interaction_model_interaction_beta)

	model_summary_handle.write(ensamble_id + '\t' + gene_symbol + '\t' + gene_chrom + '\t' + str(tss) + '\t' + str(gene_bim.shape[0]) + '\t' + str(best_eqtl_pvalue) + '\t' + str(best_interaction_qtl_pvalue) + '\t' + str(main_holdout_correlation) + '\t' + str(interaction_holdout_correlation) + '\t' + str(best_alpha) + '\t' + str(best_alpha_main3) + '\t' + str(best_alpha_interaction3) + '\t' + main_model_weight_file + '\t' + interaction_model_weight_file + '\n')
	model_summary_handle.flush()

	print(ensamble_id)
	print(main_holdout_correlation, interaction_holdout_correlation)


		


f.close()
model_summary_handle.close()
