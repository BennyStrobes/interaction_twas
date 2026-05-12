import numpy as np
import os
import sys
import pdb
import pickle
from pandas_plink import read_plink
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_predict, KFold
from sklearn.metrics import r2_score
import scipy.stats as stats
from numba import njit
import time


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

def extract_gene_metadata_from_expression_file(expression_file):
	gene_metadata = {}
	f = open(expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		gene_descriptor_fields = np.asarray(data[:6])
		ensamble_id = gene_descriptor_fields[0]
		if gene_descriptor_fields[5] == '+':
			tss = int(gene_descriptor_fields[3])
		else:
			tss = int(gene_descriptor_fields[4])

		gene_metadata[ensamble_id] = {}
		gene_metadata[ensamble_id]['gene_symbol'] = gene_descriptor_fields[1]
		gene_metadata[ensamble_id]['gene_chrom'] = gene_descriptor_fields[2]
		gene_metadata[ensamble_id]['tss'] = tss
	f.close()
	return gene_metadata

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
	return geno_mat_donor, geno_mat, gene_bim

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


@njit(fastmath=True)
def compute_gene_effect_xtxs(gene_donor_genotype, genotype_indices, expression_pcs):
	n_snps = gene_donor_genotype.shape[0]
	n_samples = len(genotype_indices)
	n_pcs = expression_pcs.shape[1]
	main_effect_xtx = np.zeros(n_snps)
	interaction_effect_xtx = np.zeros((n_snps, n_pcs))

	for snp_iter in range(n_snps):
		for sample_iter in range(n_samples):
			genotype_val = gene_donor_genotype[snp_iter, genotype_indices[sample_iter]]
			genotype_sq = genotype_val*genotype_val
			main_effect_xtx[snp_iter] = main_effect_xtx[snp_iter] + genotype_sq
			for pc_iter in range(n_pcs):
				interaction_effect_xtx[snp_iter, pc_iter] = interaction_effect_xtx[snp_iter, pc_iter] + (genotype_sq*expression_pcs[sample_iter, pc_iter]*expression_pcs[sample_iter, pc_iter])

	return main_effect_xtx, interaction_effect_xtx


def initialize_model_parameters(gene_training_data, genotype_indices, expression_pcs):
	n_genes = len(gene_training_data)
	ordered_genes = np.sort(np.asarray([*gene_training_data]))
	n_pcs = expression_pcs.shape[1]


	main_effects = []
	interaction_effects = []
	resid_expression = []
	main_effect_xtxs = []
	interaction_effect_xtxs = []

	for gene_id in ordered_genes:
		gene_donor_genotype = np.ascontiguousarray(gene_training_data[gene_id]['donor_genotype'])
		n_gene_snps = gene_donor_genotype.shape[0]
		gene_main_effect_xtx, gene_interaction_effect_xtx = compute_gene_effect_xtxs(gene_donor_genotype, genotype_indices, expression_pcs)
		resid_expression.append(np.copy(gene_training_data[gene_id]['resid_expr']))
		main_effects.append(np.zeros(n_gene_snps))
		interaction_effects.append(np.zeros((n_gene_snps, n_pcs)))
		main_effect_xtxs.append(gene_main_effect_xtx)
		interaction_effect_xtxs.append(gene_interaction_effect_xtx)

	# Hyperparameters
	resid_variances = np.ones(n_genes)
	main_effect_spike_prob = 0.5
	interaction_effects_spike_probs = np.ones(n_pcs)*0.5
	main_effect_slab_var = .1
	interaction_effects_slab_vars = np.ones(n_pcs)*.1
	resid_variance_prior_shape = 1e-6
	resid_variance_prior_scale = 1e-6
	spike_prob_prior_alpha = 0.1
	spike_prob_prior_beta = 0.1
	slab_variance_prior_shape = 1e-6
	slab_variance_prior_scale = 1e-6

	model_parameters = {}
	model_parameters['resid_variances'] = resid_variances
	model_parameters['resid_variance_prior_shape'] = resid_variance_prior_shape
	model_parameters['resid_variance_prior_scale'] = resid_variance_prior_scale
	model_parameters['spike_prob_prior_alpha'] = spike_prob_prior_alpha
	model_parameters['spike_prob_prior_beta'] = spike_prob_prior_beta
	model_parameters['slab_variance_prior_shape'] = slab_variance_prior_shape
	model_parameters['slab_variance_prior_scale'] = slab_variance_prior_scale
	model_parameters['main_effect_spike_prob'] = main_effect_spike_prob
	model_parameters['interaction_effects_spike_probs'] = interaction_effects_spike_probs
	model_parameters['main_effect_slab_var'] = main_effect_slab_var
	model_parameters['interaction_effects_slab_vars'] = interaction_effects_slab_vars
	model_parameters['resid_expression'] = resid_expression
	model_parameters['main_effects'] = main_effects
	model_parameters['interaction_effects'] = interaction_effects
	model_parameters['main_effect_xtxs'] = main_effect_xtxs
	model_parameters['interaction_effect_xtxs'] = interaction_effect_xtxs
	model_parameters['ordered_genes'] = ordered_genes
	model_parameters['n_pcs'] = n_pcs

	return model_parameters


@njit(fastmath=True)
def _sample_spike_slab_beta_from_sufficient_stats(xtx, xty, resid_var, spike_prob, slab_var):
	posterior_var = 1.0/((xtx/resid_var) + (1.0/slab_var))
	posterior_mean = posterior_var*xty/resid_var

	if spike_prob < 1e-12:
		spike_prob = 1e-12
	elif spike_prob > (1.0 - 1e-12):
		spike_prob = 1.0 - 1e-12

	log_bayes_factor = 0.5*(np.log(posterior_var/slab_var) + ((posterior_mean*posterior_mean)/posterior_var))
	log_slab_weight = np.log(1.0 - spike_prob) + log_bayes_factor
	log_spike_weight = np.log(spike_prob)
	if log_slab_weight > log_spike_weight:
		slab_posterior_prob = 1.0/(1.0 + np.exp(log_spike_weight - log_slab_weight))
	else:
		weight_ratio = np.exp(log_slab_weight - log_spike_weight)
		slab_posterior_prob = weight_ratio/(1.0 + weight_ratio)

	if np.random.random() < slab_posterior_prob:
		return np.random.normal(posterior_mean, np.sqrt(posterior_var))
	return 0.0


@njit(fastmath=True)
def update_causal_effect_sizes_for_a_single_gene(gene_resid_expression, gene_main_effects, gene_interaction_effects, gene_donor_genotype, genotype_indices, expression_pcs, gene_main_effect_xtx, gene_interaction_effect_xtx, gene_resid_var, main_effect_spike_prob, main_effect_slab_var, interaction_effects_spike_probs, interaction_effects_slab_vars):
	n_snps = gene_donor_genotype.shape[0]
	n_samples = len(genotype_indices)
	n_pcs = expression_pcs.shape[1]

	for snp_iter in range(n_snps):
		old_beta = gene_main_effects[snp_iter]
		xtx = gene_main_effect_xtx[snp_iter]

		if xtx > 0.0:
			xty = 0.0
			for sample_iter in range(n_samples):
				covariate_val = gene_donor_genotype[snp_iter, genotype_indices[sample_iter]]
				xty = xty + (covariate_val*gene_resid_expression[sample_iter])
			xty = xty + (old_beta*xtx)
			new_beta = _sample_spike_slab_beta_from_sufficient_stats(xtx, xty, gene_resid_var, main_effect_spike_prob, main_effect_slab_var)
		else:
			new_beta = 0.0
		gene_main_effects[snp_iter] = new_beta
		delta_beta = old_beta - new_beta
		if delta_beta != 0.0:
			for sample_iter in range(n_samples):
				covariate_val = gene_donor_genotype[snp_iter, genotype_indices[sample_iter]]
				gene_resid_expression[sample_iter] = gene_resid_expression[sample_iter] + (covariate_val*delta_beta)

		for pc_iter in range(n_pcs):
			old_beta = gene_interaction_effects[snp_iter, pc_iter]
			xtx = gene_interaction_effect_xtx[snp_iter, pc_iter]

			if xtx > 0.0:
				xty = 0.0
				for sample_iter in range(n_samples):
					covariate_val = gene_donor_genotype[snp_iter, genotype_indices[sample_iter]]*expression_pcs[sample_iter, pc_iter]
					xty = xty + (covariate_val*gene_resid_expression[sample_iter])
				xty = xty + (old_beta*xtx)
				new_beta = _sample_spike_slab_beta_from_sufficient_stats(xtx, xty, gene_resid_var, interaction_effects_spike_probs[pc_iter], interaction_effects_slab_vars[pc_iter])
			else:
				new_beta = 0.0
			gene_interaction_effects[snp_iter, pc_iter] = new_beta
			delta_beta = old_beta - new_beta
			if delta_beta != 0.0:
				for sample_iter in range(n_samples):
					covariate_val = gene_donor_genotype[snp_iter, genotype_indices[sample_iter]]*expression_pcs[sample_iter, pc_iter]
					gene_resid_expression[sample_iter] = gene_resid_expression[sample_iter] + (covariate_val*delta_beta)

	return gene_main_effects, gene_interaction_effects, gene_resid_expression


def update_residual_variance_for_a_single_gene(gene_resid_expression, resid_variance_prior_shape, resid_variance_prior_scale):
	posterior_shape = resid_variance_prior_shape + (len(gene_resid_expression)/2.0)
	posterior_scale = resid_variance_prior_scale + (np.sum(np.square(gene_resid_expression))/2.0)
	return 1.0/np.random.gamma(posterior_shape, scale=(1.0/posterior_scale))


def update_causal_effect_sizes_and_residual_variances(gene_training_data, model_parameters, genotype_indices, expression_pcs):
	# Loop through genes
	for gene_iter, gene_name in enumerate(model_parameters['ordered_genes']):

		# update causal effects for a single gene
		updated_gene_main_effects, updated_gene_interaction_effects, updated_gene_resid_expression = update_causal_effect_sizes_for_a_single_gene(model_parameters['resid_expression'][gene_iter], model_parameters['main_effects'][gene_iter], model_parameters['interaction_effects'][gene_iter], gene_training_data[gene_name]['donor_genotype'], genotype_indices, expression_pcs, model_parameters['main_effect_xtxs'][gene_iter], model_parameters['interaction_effect_xtxs'][gene_iter], model_parameters['resid_variances'][gene_iter], model_parameters['main_effect_spike_prob'], model_parameters['main_effect_slab_var'], model_parameters['interaction_effects_spike_probs'], model_parameters['interaction_effects_slab_vars'])
		model_parameters['resid_expression'][gene_iter] = updated_gene_resid_expression
		model_parameters['main_effects'][gene_iter] = updated_gene_main_effects
		model_parameters['interaction_effects'][gene_iter] = updated_gene_interaction_effects

		# Update resid variance for this gene
		model_parameters['resid_variances'][gene_iter] = update_residual_variance_for_a_single_gene(model_parameters['resid_expression'][gene_iter], model_parameters['resid_variance_prior_shape'], model_parameters['resid_variance_prior_scale'])


	return model_parameters



def update_priors(gene_training_data, model_parameters):
	spike_alpha = model_parameters['spike_prob_prior_alpha']
	spike_beta = model_parameters['spike_prob_prior_beta']
	slab_shape = model_parameters['slab_variance_prior_shape']
	slab_scale = model_parameters['slab_variance_prior_scale']
	n_pcs = model_parameters['n_pcs']

	main_n_spike = 0
	main_n_slab = 0
	main_slab_sum_sq = 0.0
	interaction_n_spike = np.zeros(n_pcs)
	interaction_n_slab = np.zeros(n_pcs)
	interaction_slab_sum_sq = np.zeros(n_pcs)

	for gene_iter, gene_name in enumerate(model_parameters['ordered_genes']):
		main_effects = model_parameters['main_effects'][gene_iter]
		main_slab_indices = main_effects != 0.0
		main_n_slab = main_n_slab + np.sum(main_slab_indices)
		main_n_spike = main_n_spike + np.sum(main_slab_indices == False)
		main_slab_sum_sq = main_slab_sum_sq + np.sum(np.square(main_effects[main_slab_indices]))

		interaction_effects = model_parameters['interaction_effects'][gene_iter]
		for pc_iter in range(n_pcs):
			pc_effects = interaction_effects[:, pc_iter]
			pc_slab_indices = pc_effects != 0.0
			interaction_n_slab[pc_iter] = interaction_n_slab[pc_iter] + np.sum(pc_slab_indices)
			interaction_n_spike[pc_iter] = interaction_n_spike[pc_iter] + np.sum(pc_slab_indices == False)
			interaction_slab_sum_sq[pc_iter] = interaction_slab_sum_sq[pc_iter] + np.sum(np.square(pc_effects[pc_slab_indices]))

	model_parameters['main_effect_spike_prob'] = np.random.beta(spike_alpha + main_n_spike, spike_beta + main_n_slab)
	model_parameters['main_effect_slab_var'] = 1.0/np.random.gamma(slab_shape + (main_n_slab/2.0), scale=(1.0/(slab_scale + (main_slab_sum_sq/2.0))))

	for pc_iter in range(n_pcs):
		model_parameters['interaction_effects_spike_probs'][pc_iter] = np.random.beta(spike_alpha + interaction_n_spike[pc_iter], spike_beta + interaction_n_slab[pc_iter])
		model_parameters['interaction_effects_slab_vars'][pc_iter] = 1.0/np.random.gamma(slab_shape + (interaction_n_slab[pc_iter]/2.0), scale=(1.0/(slab_scale + (interaction_slab_sum_sq[pc_iter]/2.0))))

	return model_parameters

def predict_expression_levels_using_bayesian_inference(gene_training_data, genotype_indices, expression_pcs, max_iters=1000, burn_in_iters=700):
	if burn_in_iters >= max_iters:
		print('Burn-in iterations must be less than max iterations')
		pdb.set_trace()

	model_parameters = initialize_model_parameters(gene_training_data, genotype_indices, expression_pcs)
	posterior_main_effect_sums = []
	posterior_interaction_effect_sums = []
	for gene_iter, gene_name in enumerate(model_parameters['ordered_genes']):
		posterior_main_effect_sums.append(np.zeros(model_parameters['main_effects'][gene_iter].shape))
		posterior_interaction_effect_sums.append(np.zeros(model_parameters['interaction_effects'][gene_iter].shape))
	n_posterior_samples = 0

	#####
	t1 = time.time()
	# Begin iterative algorithm
	for itera in range(max_iters):
		# Update causal effect sizes and residual variances
		model_parameters = update_causal_effect_sizes_and_residual_variances(gene_training_data, model_parameters, genotype_indices, expression_pcs)

		# Update priors
		model_parameters = update_priors(gene_training_data, model_parameters)
		t2 = time.time()
		print('##########################')
		print('Iteration: ' + str(itera),flush=True)
		print('time: ' + str(t2-t1),flush=True)
		print(model_parameters['main_effect_spike_prob'],flush=True)
		print(model_parameters['main_effect_slab_var'],flush=True)
		print(model_parameters['interaction_effects_spike_probs'],flush=True)
		print(model_parameters['interaction_effects_slab_vars'],flush=True)
		t1 = time.time()

		if itera >= burn_in_iters:
			for gene_iter, gene_name in enumerate(model_parameters['ordered_genes']):
				posterior_main_effect_sums[gene_iter] = posterior_main_effect_sums[gene_iter] + model_parameters['main_effects'][gene_iter]
				posterior_interaction_effect_sums[gene_iter] = posterior_interaction_effect_sums[gene_iter] + model_parameters['interaction_effects'][gene_iter]
			n_posterior_samples = n_posterior_samples + 1

	posterior_mean_effects = {}
	posterior_mean_effects['ordered_genes'] = model_parameters['ordered_genes']
	posterior_mean_effects['main_effects'] = []
	posterior_mean_effects['interaction_effects'] = []
	for gene_iter, gene_name in enumerate(model_parameters['ordered_genes']):
		posterior_mean_effects['main_effects'].append(posterior_main_effect_sums[gene_iter]/n_posterior_samples)
		posterior_mean_effects['interaction_effects'].append(posterior_interaction_effect_sums[gene_iter]/n_posterior_samples)
	posterior_mean_effects['n_posterior_samples'] = n_posterior_samples
	return posterior_mean_effects

######################
# Command line args
######################
expr_pc_file = sys.argv[1]
expression_file = sys.argv[2]
genotype_indices_file = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
onek_genomes_plink_files = sys.argv[5]
predicted_expression_dir = sys.argv[6]
num_pcs = int(sys.argv[7])
if predicted_expression_dir.endswith('/') == False:
	predicted_expression_dir = predicted_expression_dir + '/'

# Pick a distance threshold
distance_threshold = 50000.0
min_snps_per_gene=50
n_folds = 5

os.makedirs(predicted_expression_dir, exist_ok=True)
os.makedirs(predicted_expression_dir + 'gene_model_weights/', exist_ok=True)


#####################
# (1) Load in data
######################

# Load in expression pcs
pc_sample_names, expression_pcs = load_in_expression_pcs(expr_pc_file)
full_expression_pcs = np.copy(expression_pcs)
expression_pcs = np.ascontiguousarray(expression_pcs[:, :num_pcs])

# Load in genotype indices
genotype_indices = np.ascontiguousarray(np.loadtxt(genotype_indices_file).astype(np.int64))

gene_training_data = {}

for chrom_num_int in range(1, 23):
	chrom_num = str(chrom_num_int)

	# Extact dictionary list of 1KG variants
	one_kg_variant_dicti = extract_dictionary_list_of_1_kg_variants(onek_genomes_plink_files + '1000G.EUR.hg38.' + chrom_num + '.bim')

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
		std_geno_mat_donor, std_geno_mat, gene_bim = filter_and_standardize_snp_genotype_matrix(geno_mat_donor, gene_bim, genotype_indices)
		if std_geno_mat.shape[0] < min_snps_per_gene:
			continue

		# Run standard marginal eQTL calling for this gene with expression PCs as covariates
		beta, beta_se, pvalue, valid_eqtl_snps = run_standard_eqtl_calling_for_gene(expr_vec, expression_pcs, std_geno_mat)
		interaction_beta, interaction_beta_se, interaction_pvalue, valid_interaction_tests = run_interaction_qtl_calling_for_gene(expr_vec, expression_pcs, std_geno_mat)
		best_eqtl_pvalue = np.nanmin(pvalue)
		best_interaction_qtl_pvalue = np.nanmin(interaction_pvalue)

		# Residualize expression on expression PCs and standardize residual expression
		residual_expr_vec = residualize_and_standardize_expression(expr_vec, full_expression_pcs)
		#genotype_main_effects, interaction_terms = construct_context_interaction_design_matrices(std_geno_mat, expression_pcs)

		if np.array_equal(std_geno_mat_donor[:, genotype_indices], std_geno_mat) == False:
			print('assumption erororo')
			pdb.set_trace()


		if ensamble_id in gene_training_data:
			print('assumption erororo')
			pdb.set_trace()

		gene_training_data[ensamble_id] = {}
		gene_training_data[ensamble_id]['gene_symbol'] = gene_symbol
		gene_training_data[ensamble_id]['gene_chrom'] = gene_chrom
		gene_training_data[ensamble_id]['tss'] = tss
		gene_training_data[ensamble_id]['resid_expr'] = residual_expr_vec
		gene_training_data[ensamble_id]['donor_genotype'] = std_geno_mat_donor
		gene_training_data[ensamble_id]['bim'] = gene_bim
		gene_training_data[ensamble_id]['best_eqtl_pvalue'] = best_eqtl_pvalue
		gene_training_data[ensamble_id]['best_interaction_eqtl_pvalue'] = best_interaction_qtl_pvalue



	f.close()


'''
gene_training_data_file = predicted_expression_dir + 'gene_training_data.pkl'
with open(gene_training_data_file, 'wb') as output_handle:
	pickle.dump(gene_training_data, output_handle)

gene_training_data_file = predicted_expression_dir + 'gene_training_data.pkl'
with open(gene_training_data_file, 'rb') as input_handle:
	gene_training_data = pickle.load(input_handle)
'''






#####################
# (2) Run Inference
######################
posterior_mean_effects = predict_expression_levels_using_bayesian_inference(gene_training_data, genotype_indices, expression_pcs, max_iters=500, burn_in_iters=400)


#####################
# (3) Print to output
######################
for chrom_num_int in range(1, 23):
	chrom_num = str(chrom_num_int)
	chrom_string = 'chr' + chrom_num

	model_summary_file = predicted_expression_dir + 'chr' + str(chrom_num) + '_gene_model_summary.txt'
	model_summary_handle = open(model_summary_file, 'w')
	model_summary_handle.write('ensamble_id\tgene_symbol\tchrom\ttss\tn_snps\tbest_eqtl_pvalue\tbest_interaction_qtl_pvalue\tmain_model_holdout_correlation\tinteraction_model_holdout_correlation\tmain_model_alpha\tinteraction_model_alpha_main\tinteraction_model_alpha_interaction\tmain_model_weight_file\tinteraction_model_weight_file\n')

	gene_metadata = extract_gene_metadata_from_expression_file(expression_file)
	nan_string = 'nan'

	for gene_iter, ensamble_id in enumerate(posterior_mean_effects['ordered_genes']):
		gene_data = gene_training_data[ensamble_id]
		if 'gene_chrom' in gene_data:
			gene_chrom = gene_data['gene_chrom']
			gene_symbol = gene_data['gene_symbol']
			tss = gene_data['tss']
		else:
			gene_chrom = gene_metadata[ensamble_id]['gene_chrom']
			gene_symbol = gene_metadata[ensamble_id]['gene_symbol']
			tss = gene_metadata[ensamble_id]['tss']

		if gene_chrom != chrom_string:
			continue

		gene_bim = gene_data['bim']
		main_model_beta = posterior_mean_effects['main_effects'][gene_iter]
		interaction_model_main_beta = posterior_mean_effects['main_effects'][gene_iter]
		interaction_model_interaction_beta = posterior_mean_effects['interaction_effects'][gene_iter]

		gene_weight_file_stem = create_gene_weight_file_stem(predicted_expression_dir, chrom_num, ensamble_id, gene_symbol)
		main_model_weight_file = gene_weight_file_stem + '_main_effect_bayesian_weights.txt'
		interaction_model_weight_file = gene_weight_file_stem + '_interaction_bayesian_weights.txt'
		save_main_effect_model_weights(main_model_weight_file, gene_bim, gene_chrom, main_model_beta)
		save_interaction_model_weights(interaction_model_weight_file, gene_bim, gene_chrom, interaction_model_main_beta, interaction_model_interaction_beta)

		model_summary_handle.write(ensamble_id + '\t' + gene_symbol + '\t' + gene_chrom + '\t' + str(tss) + '\t' + str(gene_bim.shape[0]) + '\t' + str(gene_data['best_eqtl_pvalue']) + '\t' + str(gene_data['best_interaction_eqtl_pvalue']) + '\t' + nan_string + '\t' + nan_string + '\t' + nan_string + '\t' + nan_string + '\t' + nan_string + '\t' + main_model_weight_file + '\t' + interaction_model_weight_file + '\n')
		model_summary_handle.flush()


	model_summary_handle.close()




