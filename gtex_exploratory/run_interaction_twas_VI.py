import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink
import time
import pickle
from numba import njit



def create_mapping_from_gtex_id_to_rsids(bim):
	bim_arr = np.asarray(bim)
	mapping = {}
	mapping2 = {}
	for row_iter in range(bim_arr.shape[0]):
		rsid = bim_arr[row_iter, 1]

		gtex_id1 = 'chr' + bim_arr[row_iter,0] + '_' + str(bim_arr[row_iter, 3]) + '_' + bim_arr[row_iter, 4] + '_' + bim_arr[row_iter, 5] + '_b38'
		gtex_id2 = 'chr' + bim_arr[row_iter,0] + '_' + str(bim_arr[row_iter, 3]) + '_' + bim_arr[row_iter, 5] + '_' + bim_arr[row_iter,4] + '_b38'

		if gtex_id1 in mapping or gtex_id2 in mapping:
			continue
		mapping[gtex_id1] = rsid
		mapping[gtex_id2] = rsid

		if gtex_id1 in mapping2 or gtex_id2 in mapping2:
			continue

		mapping2[gtex_id1] = (row_iter, 1.0)
		mapping2[gtex_id2] = (row_iter, -1.0)

	return mapping, mapping2

def load_in_interaction_effect_weights(interaction_ridge_weights_file):
	g = open(interaction_ridge_weights_file)
	main_effects = []
	interaction_effects = []

	variant_ids = []

	head_counter = 0
	for line in g:
		line = line.rstrip()
		data = line.split('\t')
		if head_counter == 0:
			head_counter = head_counter + 1
			continue
		gtex_id = data[0]
		if gtex_id.split('_')[2] != data[4]:
			print('assumption erroorro')
			pdb.set_trace()
		if gtex_id.split('_')[3] != data[3]:
			print('assumption eroror')
			pdb.set_trace()
		main_effects.append(float(data[5]))
		interaction_effects.append(np.asarray(data[6:]).astype(float))
		variant_ids.append(gtex_id)

	g.close()

	main_effects = np.asarray(main_effects)
	interaction_effects = np.asarray(interaction_effects)
	variant_ids = np.asarray(variant_ids)
	return main_effects, interaction_effects, variant_ids

def load_in_genotype_indices_plus_sign_flip_info(gene_gtex_ids, gtex_id_to_geno_index):
	gene_genotype_indices = []
	gene_genotype_sign_flips = []

	for gtex_id in gene_gtex_ids:

		geno_index, sign_flipper = gtex_id_to_geno_index[gtex_id]
		gene_genotype_indices.append(geno_index)
		gene_genotype_sign_flips.append(sign_flipper)

	return np.asarray(gene_genotype_indices), np.asarray(gene_genotype_sign_flips)

def load_in_gwas_sumstats(trait_sumstat_file):
	mapping = {}
	head_count = 0
	f = open(trait_sumstat_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 5:
			print('assumptioneornronro')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			continue

		rsid = data[0]
		alleles = (data[1], data[2])
		N = float(data[3])
		zed = float(data[4])

		if rsid in mapping:
			print('asssumptioneornro')
			pdb.set_trace()

		mapping[rsid] = (zed, N, alleles)

	f.close()

	return mapping

def impute_missing_gwas_z_scores(gwas_zeds, LD, ridge=1e-3):
	gwas_zeds = np.asarray(gwas_zeds).astype(float)
	observed = np.isfinite(gwas_zeds)
	missing = observed == False
	imputed_gwas_zeds = np.copy(gwas_zeds)

	if np.sum(missing) == 0:
		return imputed_gwas_zeds
	if np.sum(observed) == 0:
		print('No observed GWAS z-scores')
		pdb.set_trace()

	LD_obs_obs = LD[observed, :][:, observed]
	LD_miss_obs = LD[missing, :][:, observed]
	regularized_LD_obs_obs = LD_obs_obs + ridge*np.eye(LD_obs_obs.shape[0])
	imputed_gwas_zeds[missing] = np.dot(LD_miss_obs, np.linalg.solve(regularized_LD_obs_obs, gwas_zeds[observed]))
	return imputed_gwas_zeds


def run_held_out_twas(valid_testing_chroms, pred_expresssion_output_stem, trait_sumstat_file, onek_genomes_plink_files, inferred_ee, held_out_twas_analysis_output_file):
	t = open(held_out_twas_analysis_output_file,'w')
	t.write('gene_id\tshared_twas_z\tinteraction_twas_z\n')
	for chrom_num in valid_testing_chroms:

		# Load in plink data for this chromosome
		genotype_stem = onek_genomes_plink_files + '1000G.EUR.hg38.' + str(chrom_num) 
		(bim, fam, G) = read_plink(genotype_stem)

		# Create mapping from rsid to gtex id
		gtex_id_to_rsid, gtex_id_to_geno_index = create_mapping_from_gtex_id_to_rsids(bim)

		# Now loop through genes
		gene_pred_summary_file = pred_expresssion_output_stem + '/' + 'chr' + str(chrom_num) + '_gene_model_summary.txt'

		# Load in gwas sumstats
		rsid_to_gwas_zed_and_N = load_in_gwas_sumstats(trait_sumstat_file)


		f = open(gene_pred_summary_file)
		head_count = 0
		misses = 0
		successes = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			ens_id = data[0]
			min_pval = float(data[5])
			min_interaction_qtl_pval = float(data[6])
			if min_pval > .001 and min_interaction_qtl_pval > .001:
				continue 


			# Load in stuff from gene prediction file
			interaction_ridge_weights_file = data[13]
			main_effects, interaction_effects, gene_gtex_ids = load_in_interaction_effect_weights(interaction_ridge_weights_file)

			main_weights_file = data[12]
			only_main_effects = np.loadtxt(main_weights_file,dtype=str,delimiter='\t')[1:,-1].astype(float)

			# LOAD IN LD
			# Extract ordered genotype indices, and sign flip info
			gene_genotype_indices, gene_genotype_sign_flips = load_in_genotype_indices_plus_sign_flip_info(gene_gtex_ids, gtex_id_to_geno_index)
			geno_mat_donor = G[gene_genotype_indices,:].compute()
			geno_mat_donor = 2.0 - geno_mat_donor
			flip_indices = gene_genotype_sign_flips == -1.0
			geno_mat_donor[flip_indices, :] = 2.0 - geno_mat_donor[flip_indices, :]
			if np.sum(np.sum(np.isnan(geno_mat_donor))) > 0:
				print('assumpriontoineron')
				pdb.set_trace()
			LD = np.corrcoef(geno_mat_donor)

			# Load in GWAS
			gwas_zeds= []
			gwas_NNs = []
			for gtex_id in gene_gtex_ids:
				rsid = gtex_id_to_rsid[gtex_id]

				if rsid in rsid_to_gwas_zed_and_N:
					zed, NN, alleles = rsid_to_gwas_zed_and_N[rsid]

					if alleles[0] == gtex_id.split('_')[3] and alleles[1] == gtex_id.split('_')[2]:
						zed = zed*1.0
					elif alleles[0] == gtex_id.split('_')[2] and alleles[1] == gtex_id.split('_')[3]:
						zed = zed*-1.0
					else:
						print('assumption erororo')
						pdb.set_trace()				

					gwas_zeds.append(zed)
					gwas_NNs.append(NN)
				else:
					gwas_zeds.append(np.nan)
					gwas_NNs.append(np.nan)


			gwas_zeds = np.asarray(gwas_zeds).astype(float)
			gwas_NNs = np.asarray(gwas_NNs).astype(float)
			obs_fraction = np.sum(np.isnan(gwas_zeds)==False)/len(gwas_zeds)
			if obs_fraction < .05:
				continue

			if len(np.unique(gwas_NNs[np.isnan(gwas_NNs)==False])) != 1:
				print('assumption erroror')
				pdb.set_trace()

			gwas_sample_size = np.unique(gwas_NNs[np.isnan(gwas_NNs)==False])[0]

			# Impute missing zeds
			imputed_gwas_zeds = impute_missing_gwas_z_scores(gwas_zeds, LD)

			# Pred gene deltas
			pred_gene_deltas = main_effects + np.dot(interaction_effects, inferred_ee)

			interaction_twas_z = pred_gene_deltas @ imputed_gwas_zeds / np.sqrt(pred_gene_deltas @ LD @ pred_gene_deltas)
			standard_twas_z = only_main_effects @ imputed_gwas_zeds / np.sqrt(only_main_effects @ LD @ only_main_effects)


			t.write(ens_id + '\t' + str(standard_twas_z) + '\t' + str(interaction_twas_z) + '\n')


		f.close()

	t.close()

	return 

def load_in_training_data(valid_chroms, pred_expresssion_output_stem, trait_sumstat_file, onek_genomes_plink_files):
	training_data_obj = {}

	for chrom_num in valid_chroms:

		# Load in plink data for this chromosome
		genotype_stem = onek_genomes_plink_files + '1000G.EUR.hg38.' + str(chrom_num) 
		(bim, fam, G) = read_plink(genotype_stem)

		# Create mapping from rsid to gtex id
		gtex_id_to_rsid, gtex_id_to_geno_index = create_mapping_from_gtex_id_to_rsids(bim)

		# Now loop through genes
		gene_pred_summary_file = pred_expresssion_output_stem + '/' + 'chr' + str(chrom_num) + '_gene_model_summary.txt'

		# Load in gwas sumstats
		rsid_to_gwas_zed_and_N = load_in_gwas_sumstats(trait_sumstat_file)


		f = open(gene_pred_summary_file)
		head_count = 0
		misses = 0
		successes = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			ens_id = data[0]
			min_pval = float(data[5])
			min_interaction_qtl_pval = float(data[6])
			if min_pval > .001 and min_interaction_qtl_pval > .001:
				continue 


			# Load in stuff from gene prediction file
			interaction_ridge_weights_file = data[13]
			main_effects, interaction_effects, gene_gtex_ids = load_in_interaction_effect_weights(interaction_ridge_weights_file)

			# LOAD IN LD
			# Extract ordered genotype indices, and sign flip info
			gene_genotype_indices, gene_genotype_sign_flips = load_in_genotype_indices_plus_sign_flip_info(gene_gtex_ids, gtex_id_to_geno_index)
			geno_mat_donor = G[gene_genotype_indices,:].compute()
			geno_mat_donor = 2.0 - geno_mat_donor
			flip_indices = gene_genotype_sign_flips == -1.0
			geno_mat_donor[flip_indices, :] = 2.0 - geno_mat_donor[flip_indices, :]
			if np.sum(np.sum(np.isnan(geno_mat_donor))) > 0:
				print('assumpriontoineron')
				pdb.set_trace()
			LD = np.corrcoef(geno_mat_donor)

			# Load in GWAS
			gwas_zeds= []
			gwas_NNs = []
			for gtex_id in gene_gtex_ids:
				rsid = gtex_id_to_rsid[gtex_id]

				if rsid in rsid_to_gwas_zed_and_N:
					zed, NN, alleles = rsid_to_gwas_zed_and_N[rsid]

					if alleles[0] == gtex_id.split('_')[3] and alleles[1] == gtex_id.split('_')[2]:
						zed = zed*1.0
					elif alleles[0] == gtex_id.split('_')[2] and alleles[1] == gtex_id.split('_')[3]:
						zed = zed*-1.0
					else:
						print('assumption erororo')
						pdb.set_trace()				

					gwas_zeds.append(zed)
					gwas_NNs.append(NN)
				else:
					gwas_zeds.append(np.nan)
					gwas_NNs.append(np.nan)


			gwas_zeds = np.asarray(gwas_zeds).astype(float)
			gwas_NNs = np.asarray(gwas_NNs).astype(float)
			obs_fraction = np.sum(np.isnan(gwas_zeds)==False)/len(gwas_zeds)
			if obs_fraction < .05:
				continue

			if len(np.unique(gwas_NNs[np.isnan(gwas_NNs)==False])) != 1:
				print('assumption erroror')
				pdb.set_trace()

			gwas_sample_size = np.unique(gwas_NNs[np.isnan(gwas_NNs)==False])[0]

			# Impute missing zeds
			imputed_gwas_zeds = impute_missing_gwas_z_scores(gwas_zeds, LD)

			if ens_id in training_data_obj:
				print('assumptionerororo')
				pdb.set_trace()

			training_data_obj[ens_id] = {}
			training_data_obj[ens_id]['LD'] = LD
			training_data_obj[ens_id]['main_effects'] = main_effects
			training_data_obj[ens_id]['interaction_effects'] = interaction_effects
			training_data_obj[ens_id]['GWAS_Z'] = imputed_gwas_zeds
			training_data_obj[ens_id]['GWAS_N'] = gwas_sample_size
			training_data_obj[ens_id]['variant_ids'] = gene_gtex_ids

		f.close()



	return training_data_obj




def initialize_data(training_data_obj,LL):
	n_genes = len(training_data_obj)
	ordered_genes = np.sort([*training_data_obj])

	n_components = training_data_obj[ordered_genes[0]]['interaction_effects'].shape[1]
	gene_alphas = np.zeros((n_genes, LL))
	EEs = np.zeros((n_components, LL))

	mixture_probs = np.ones(LL)*.1
	mixture_vars = np.ones(LL)*.00001

	gene_vars = []
	for gene_id in ordered_genes:
		gene_var = np.dot(np.dot(training_data_obj[gene_id]['main_effects'], training_data_obj[gene_id]['LD']), training_data_obj[gene_id]['main_effects'])
		if np.isfinite(gene_var) == False or gene_var <= 0:
			print('Gene variance must be positive')
			pdb.set_trace()
		gene_vars.append(gene_var)
	gene_vars = np.asarray(gene_vars)


	return ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars


def logistic_function(log_odds):
	if log_odds >= 0.0:
		return 1.0/(1.0 + np.exp(-log_odds))
	exp_log_odds = np.exp(log_odds)
	return exp_log_odds/(1.0 + exp_log_odds)

def compute_alpha_variational_params(gene_delta, residual_gwas_zed, LD, gwas_sample_size, mixture_prob, mixture_var, tau):
	if mixture_var <= 0:
		print('Mixture variance must be positive')
		pdb.set_trace()
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()
	mixture_prob = np.clip(mixture_prob, 1e-12, 1.0 - 1e-12)

	delta_R_delta = np.dot(gene_delta, np.dot(LD, gene_delta))
	if np.isfinite(delta_R_delta) == False or delta_R_delta <= 0:
		return 0.0, 0.0, 0.0

	likelihood_precision = (gwas_sample_size/tau)*delta_R_delta
	likelihood_mean_numer = (np.sqrt(gwas_sample_size)/tau)*np.dot(gene_delta, residual_gwas_zed)
	posterior_var = 1.0/(likelihood_precision + (1.0/mixture_var))
	posterior_mean = posterior_var*likelihood_mean_numer

	log_bayes_factor = 0.5*(np.log(posterior_var) - np.log(mixture_var)) + 0.5*(posterior_mean*posterior_mean/posterior_var)
	log_prior_odds = np.log(mixture_prob) - np.log(1.0 - mixture_prob)
	log_posterior_odds = log_prior_odds + log_bayes_factor
	posterior_inclusion_prob = logistic_function(log_posterior_odds)

	alpha_mean = posterior_inclusion_prob*posterior_mean
	alpha_second_moment = posterior_inclusion_prob*(posterior_var + (posterior_mean*posterior_mean))
	return alpha_mean, posterior_inclusion_prob, alpha_second_moment

def initialize_ashr_gene_prior(n_factors, initial_inclusion_prob=.1):
	ashr_grid_sds = np.asarray([0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
	ashr_grid_vars = ashr_grid_sds*ashr_grid_sds
	ashr_weights = np.zeros((n_factors, len(ashr_grid_vars)))
	ashr_weights[:, 0] = 1.0 - initial_inclusion_prob
	ashr_weights[:, 1:] = initial_inclusion_prob/(len(ashr_grid_vars) - 1.0)
	return ashr_grid_vars, ashr_weights

def compute_alpha_variational_params_ashr(gene_delta, residual_gwas_zed, LD, gwas_sample_size, ashr_component_weights, ashr_grid_vars, tau):
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()
	if len(ashr_component_weights) != len(ashr_grid_vars):
		print('Ashr grid weights and variances must have same length')
		pdb.set_trace()

	ashr_component_weights = np.asarray(ashr_component_weights).astype(float)
	ashr_grid_vars = np.asarray(ashr_grid_vars).astype(float)
	if np.sum(ashr_component_weights) <= 0:
		print('Ashr component weights must sum to a positive value')
		pdb.set_trace()
	ashr_component_weights = ashr_component_weights/np.sum(ashr_component_weights)
	ashr_component_weights = np.clip(ashr_component_weights, 1e-300, 1.0)
	ashr_component_weights = ashr_component_weights/np.sum(ashr_component_weights)

	component_log_weights = np.zeros(len(ashr_grid_vars))
	component_means = np.zeros(len(ashr_grid_vars))
	component_second_moments = np.zeros(len(ashr_grid_vars))

	delta_R_delta = np.dot(gene_delta, np.dot(LD, gene_delta))
	if np.isfinite(delta_R_delta) == False or delta_R_delta <= 0:
		component_probs = np.zeros(len(ashr_grid_vars))
		component_probs[0] = 1.0
		return 0.0, 0.0, 0.0, component_probs

	likelihood_precision = (gwas_sample_size/tau)*delta_R_delta
	likelihood_mean_numer = (np.sqrt(gwas_sample_size)/tau)*np.dot(gene_delta, residual_gwas_zed)

	for component_iter, component_var in enumerate(ashr_grid_vars):
		component_log_weights[component_iter] = np.log(ashr_component_weights[component_iter])
		if component_var <= 0.0:
			continue
		posterior_var = 1.0/(likelihood_precision + (1.0/component_var))
		posterior_mean = posterior_var*likelihood_mean_numer
		component_log_weights[component_iter] = component_log_weights[component_iter] + (0.5*(np.log(posterior_var) - np.log(component_var)) + 0.5*(posterior_mean*posterior_mean/posterior_var))
		component_means[component_iter] = posterior_mean
		component_second_moments[component_iter] = posterior_var + (posterior_mean*posterior_mean)

	max_log_weight = np.max(component_log_weights)
	component_probs = np.exp(component_log_weights - max_log_weight)
	component_probs = component_probs/np.sum(component_probs)

	alpha_mean = np.sum(component_probs*component_means)
	alpha_second_moment = np.sum(component_probs*component_second_moments)
	alpha_inclusion_prob = np.sum(component_probs[ashr_grid_vars > 0.0])
	return alpha_mean, alpha_inclusion_prob, alpha_second_moment, component_probs

def update_ashr_gene_prior_weights_vi(gene_alpha_component_probs, component_weight_prior=.1):
	n_factors = gene_alpha_component_probs.shape[1]
	n_components = gene_alpha_component_probs.shape[2]
	ashr_weights = np.zeros((n_factors, n_components))

	for factor_iter in range(n_factors):
		component_counts = np.sum(gene_alpha_component_probs[:, factor_iter, :], axis=0)
		ashr_weights[factor_iter, :] = component_counts + component_weight_prior
		ashr_weights[factor_iter, :] = ashr_weights[factor_iter, :]/np.sum(ashr_weights[factor_iter, :])

	return ashr_weights

def update_gene_alphas_vi(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, tau):
	n_genes = gene_alphas.shape[0]
	n_factors = gene_alphas.shape[1]


	for gene_iter, gene_name in enumerate(ordered_genes):
		LD = training_data_obj[gene_name]['LD']
		gwas_zed = training_data_obj[gene_name]['GWAS_Z']
		gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
		main_effects = training_data_obj[gene_name]['main_effects']
		interaction_effects = training_data_obj[gene_name]['interaction_effects']
		gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])
		gene_deltas = []
		for factor_iter in range(n_factors):
			gene_deltas.append(gene_scale*(main_effects + np.dot(interaction_effects, EEs[:, factor_iter])))
		gene_deltas = np.asarray(gene_deltas)

		for factor_iter in range(n_factors):
			residual_gwas_zed = np.copy(gwas_zed)
			for other_factor_iter in range(n_factors):
				if other_factor_iter == factor_iter:
					continue
				residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*gene_alphas[gene_iter, other_factor_iter]*np.dot(LD, gene_deltas[other_factor_iter, :]))

			alpha_mean, alpha_inclusion_prob, alpha_second_moment = compute_alpha_variational_params(gene_deltas[factor_iter, :], residual_gwas_zed, LD, gwas_sample_size, mixture_probs[factor_iter], mixture_vars[factor_iter], tau)
			gene_alphas[gene_iter, factor_iter] = alpha_mean
			gene_alpha_inclusion_probs[gene_iter, factor_iter] = alpha_inclusion_prob
			gene_alpha_second_moments[gene_iter, factor_iter] = alpha_second_moment

	return gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments

def update_EEs_vi(gene_alphas, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, e_var, gene_vars, tau):
	if e_var <= 0:
		print('E prior variance must be positive')
		pdb.set_trace()
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()

	n_components = EEs.shape[0]
	n_factors = EEs.shape[1]

	for factor_iter in range(n_factors):
		posterior_precision = np.eye(n_components)/e_var
		posterior_mean_numer = np.zeros(n_components)

		for gene_iter, gene_name in enumerate(ordered_genes):
			LD = training_data_obj[gene_name]['LD']
			gwas_zed = training_data_obj[gene_name]['GWAS_Z']
			gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
			main_effects = training_data_obj[gene_name]['main_effects']
			interaction_effects = training_data_obj[gene_name]['interaction_effects']
			alpha = gene_alphas[gene_iter, factor_iter]
			alpha_second_moment = gene_alpha_second_moments[gene_iter, factor_iter]
			gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])

			residual_gwas_zed = np.copy(gwas_zed)
			for other_factor_iter in range(n_factors):
				if other_factor_iter == factor_iter:
					continue
				other_delta = gene_scale*(main_effects + np.dot(interaction_effects, EEs[:, other_factor_iter]))
				residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*gene_alphas[gene_iter, other_factor_iter]*np.dot(LD, other_delta))
			residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*alpha*gene_scale*np.dot(LD, main_effects))

			posterior_precision = posterior_precision + ((gwas_sample_size/tau)*alpha_second_moment*gene_scale*gene_scale*np.dot(np.transpose(interaction_effects), np.dot(LD, interaction_effects)))
			posterior_mean_numer = posterior_mean_numer + ((np.sqrt(gwas_sample_size)/tau)*alpha*gene_scale*np.dot(np.transpose(interaction_effects), residual_gwas_zed))

		posterior_mean = np.linalg.solve(posterior_precision, posterior_mean_numer)
		EEs[:, factor_iter] = posterior_mean

	return EEs

def compute_gene_delta(main_effects, interaction_effects, EE, gene_var):
	gene_scale = 1.0/np.sqrt(gene_var)
	return gene_scale*(main_effects + np.dot(interaction_effects, EE))

def compute_gene_predicted_z(LD, gwas_sample_size, gene_alpha, gene_delta):
	return np.sqrt(gwas_sample_size)*gene_alpha*np.dot(LD, gene_delta)

def initialize_nm_variant_effects(training_data_obj, ordered_genes):
	nm_variant_effects = []
	nm_variant_inclusion_probs = []
	nm_variant_second_moments = []
	for gene_name in ordered_genes:
		n_variants = len(training_data_obj[gene_name]['GWAS_Z'])
		nm_variant_effects.append(np.zeros(n_variants))
		nm_variant_inclusion_probs.append(np.zeros(n_variants))
		nm_variant_second_moments.append(np.zeros(n_variants))
	return nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments

def compute_gene_mediated_predicted_z(gene_alpha, training_data_obj, gene_name, EEs, gene_var):
	LD = training_data_obj[gene_name]['LD']
	gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
	main_effects = training_data_obj[gene_name]['main_effects']
	interaction_effects = training_data_obj[gene_name]['interaction_effects']

	predicted_z = np.zeros(LD.shape[0])
	for factor_iter in range(gene_alpha.shape[0]):
		gene_delta = compute_gene_delta(main_effects, interaction_effects, EEs[:, factor_iter], gene_var)
		predicted_z = predicted_z + compute_gene_predicted_z(LD, gwas_sample_size, gene_alpha[factor_iter], gene_delta)
	return predicted_z

def compute_nm_variant_predicted_z(LD, gwas_sample_size, nm_variant_effects):
	return np.sqrt(gwas_sample_size)*np.dot(LD, nm_variant_effects)

def compute_direct_variant_variational_params(variant_index, residual_gwas_zed, LD, gwas_sample_size, mixture_prob, mixture_var, tau):
	if mixture_var <= 0:
		print('Non-mediated variant mixture variance must be positive')
		pdb.set_trace()
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()
	mixture_prob = np.clip(mixture_prob, 1e-12, 1.0 - 1e-12)

	variant_ld_var = LD[variant_index, variant_index]
	if np.isfinite(variant_ld_var) == False or variant_ld_var <= 0:
		return 0.0, 0.0, 0.0

	likelihood_precision = (gwas_sample_size/tau)*variant_ld_var
	likelihood_mean_numer = (np.sqrt(gwas_sample_size)/tau)*residual_gwas_zed[variant_index]
	posterior_var = 1.0/(likelihood_precision + (1.0/mixture_var))
	posterior_mean = posterior_var*likelihood_mean_numer

	log_bayes_factor = 0.5*(np.log(posterior_var) - np.log(mixture_var)) + 0.5*(posterior_mean*posterior_mean/posterior_var)
	log_prior_odds = np.log(mixture_prob) - np.log(1.0 - mixture_prob)
	posterior_inclusion_prob = logistic_function(log_prior_odds + log_bayes_factor)

	variant_effect_mean = posterior_inclusion_prob*posterior_mean
	variant_effect_second_moment = posterior_inclusion_prob*(posterior_var + (posterior_mean*posterior_mean))
	return variant_effect_mean, posterior_inclusion_prob, variant_effect_second_moment

def update_gene_alphas_vi_with_nm_variants(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, nm_variant_effects, mixture_probs, mixture_vars, gene_vars, tau):
	n_factors = gene_alphas.shape[1]

	for gene_iter, gene_name in enumerate(ordered_genes):
		LD = training_data_obj[gene_name]['LD']
		gwas_zed = training_data_obj[gene_name]['GWAS_Z']
		gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
		main_effects = training_data_obj[gene_name]['main_effects']
		interaction_effects = training_data_obj[gene_name]['interaction_effects']
		gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])
		nm_predicted_z = compute_nm_variant_predicted_z(LD, gwas_sample_size, nm_variant_effects[gene_iter])

		gene_deltas = []
		for factor_iter in range(n_factors):
			gene_deltas.append(gene_scale*(main_effects + np.dot(interaction_effects, EEs[:, factor_iter])))
		gene_deltas = np.asarray(gene_deltas)

		for factor_iter in range(n_factors):
			residual_gwas_zed = np.copy(gwas_zed) - nm_predicted_z
			for other_factor_iter in range(n_factors):
				if other_factor_iter == factor_iter:
					continue
				residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*gene_alphas[gene_iter, other_factor_iter]*np.dot(LD, gene_deltas[other_factor_iter, :]))

			alpha_mean, alpha_inclusion_prob, alpha_second_moment = compute_alpha_variational_params(gene_deltas[factor_iter, :], residual_gwas_zed, LD, gwas_sample_size, mixture_probs[factor_iter], mixture_vars[factor_iter], tau)
			gene_alphas[gene_iter, factor_iter] = alpha_mean
			gene_alpha_inclusion_probs[gene_iter, factor_iter] = alpha_inclusion_prob
			gene_alpha_second_moments[gene_iter, factor_iter] = alpha_second_moment

	return gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments

def update_gene_alphas_vi_with_nm_variants_ashr(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, gene_alpha_component_probs, training_data_obj, ordered_genes, EEs, nm_variant_effects, ashr_weights, ashr_grid_vars, gene_vars, tau):
	n_factors = gene_alphas.shape[1]

	for gene_iter, gene_name in enumerate(ordered_genes):
		LD = training_data_obj[gene_name]['LD']
		gwas_zed = training_data_obj[gene_name]['GWAS_Z']
		gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
		main_effects = training_data_obj[gene_name]['main_effects']
		interaction_effects = training_data_obj[gene_name]['interaction_effects']
		gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])
		nm_predicted_z = compute_nm_variant_predicted_z(LD, gwas_sample_size, nm_variant_effects[gene_iter])

		gene_deltas = []
		for factor_iter in range(n_factors):
			gene_deltas.append(gene_scale*(main_effects + np.dot(interaction_effects, EEs[:, factor_iter])))
		gene_deltas = np.asarray(gene_deltas)

		for factor_iter in range(n_factors):
			residual_gwas_zed = np.copy(gwas_zed) - nm_predicted_z
			for other_factor_iter in range(n_factors):
				if other_factor_iter == factor_iter:
					continue
				residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*gene_alphas[gene_iter, other_factor_iter]*np.dot(LD, gene_deltas[other_factor_iter, :]))

			alpha_mean, alpha_inclusion_prob, alpha_second_moment, component_probs = compute_alpha_variational_params_ashr(gene_deltas[factor_iter, :], residual_gwas_zed, LD, gwas_sample_size, ashr_weights[factor_iter, :], ashr_grid_vars, tau)
			gene_alphas[gene_iter, factor_iter] = alpha_mean
			gene_alpha_inclusion_probs[gene_iter, factor_iter] = alpha_inclusion_prob
			gene_alpha_second_moments[gene_iter, factor_iter] = alpha_second_moment
			gene_alpha_component_probs[gene_iter, factor_iter, :] = component_probs

	return gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, gene_alpha_component_probs

def update_EEs_vi_with_nm_variants(gene_alphas, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, nm_variant_effects, e_var, gene_vars, tau):
	if e_var <= 0:
		print('E prior variance must be positive')
		pdb.set_trace()
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()

	n_components = EEs.shape[0]
	n_factors = EEs.shape[1]

	for factor_iter in range(n_factors):
		posterior_precision = np.eye(n_components)/e_var
		posterior_mean_numer = np.zeros(n_components)

		for gene_iter, gene_name in enumerate(ordered_genes):
			LD = training_data_obj[gene_name]['LD']
			gwas_zed = training_data_obj[gene_name]['GWAS_Z']
			gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
			main_effects = training_data_obj[gene_name]['main_effects']
			interaction_effects = training_data_obj[gene_name]['interaction_effects']
			alpha = gene_alphas[gene_iter, factor_iter]
			alpha_second_moment = gene_alpha_second_moments[gene_iter, factor_iter]
			gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])
			nm_predicted_z = compute_nm_variant_predicted_z(LD, gwas_sample_size, nm_variant_effects[gene_iter])

			residual_gwas_zed = np.copy(gwas_zed) - nm_predicted_z
			for other_factor_iter in range(n_factors):
				if other_factor_iter == factor_iter:
					continue
				other_delta = gene_scale*(main_effects + np.dot(interaction_effects, EEs[:, other_factor_iter]))
				residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*gene_alphas[gene_iter, other_factor_iter]*np.dot(LD, other_delta))
			residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*alpha*gene_scale*np.dot(LD, main_effects))

			posterior_precision = posterior_precision + ((gwas_sample_size/tau)*alpha_second_moment*gene_scale*gene_scale*np.dot(np.transpose(interaction_effects), np.dot(LD, interaction_effects)))
			posterior_mean_numer = posterior_mean_numer + ((np.sqrt(gwas_sample_size)/tau)*alpha*gene_scale*np.dot(np.transpose(interaction_effects), residual_gwas_zed))

		posterior_mean = np.linalg.solve(posterior_precision, posterior_mean_numer)
		EEs[:, factor_iter] = posterior_mean

	return EEs

@njit(fastmath=True)
def _logistic_function_numba(log_odds):
	if log_odds >= 0.0:
		return 1.0/(1.0 + np.exp(-log_odds))
	exp_log_odds = np.exp(log_odds)
	return exp_log_odds/(1.0 + exp_log_odds)

@njit(fastmath=True)
def _update_nm_variant_effects_for_single_gene_numba(LD, gwas_zed, sqrt_gwas_sample_size, gene_mediated_predicted_z, gene_nm_variant_effects, gene_nm_variant_inclusion_probs, gene_nm_variant_second_moments, nm_variant_mixture_prob, nm_variant_mixture_var, tau):
	n_variants = len(gene_nm_variant_effects)
	nm_predicted_z = np.zeros(n_variants)
	for variant_iter in range(n_variants):
		variant_effect = gene_nm_variant_effects[variant_iter]
		if variant_effect != 0.0:
			for variant_iter2 in range(n_variants):
				nm_predicted_z[variant_iter2] = nm_predicted_z[variant_iter2] + (sqrt_gwas_sample_size*variant_effect*LD[variant_iter2, variant_iter])

	if nm_variant_mixture_prob < 1e-12:
		nm_variant_mixture_prob = 1e-12
	elif nm_variant_mixture_prob > (1.0 - 1e-12):
		nm_variant_mixture_prob = 1.0 - 1e-12

	if nm_variant_mixture_var < 1e-12:
		nm_variant_mixture_var = 1e-12
	if tau < 1e-12:
		tau = 1e-12

	log_prior_odds = np.log(nm_variant_mixture_prob) - np.log(1.0 - nm_variant_mixture_prob)

	for variant_iter in range(n_variants):
		variant_ld_var = LD[variant_iter, variant_iter]
		if np.isfinite(variant_ld_var) == False or variant_ld_var <= 0.0:
			new_variant_effect = 0.0
			variant_inclusion_prob = 0.0
			variant_second_moment = 0.0
		else:
			old_variant_effect = gene_nm_variant_effects[variant_iter]
			old_variant_predicted_z_at_variant = sqrt_gwas_sample_size*old_variant_effect*variant_ld_var
			residual_gwas_zed_at_variant = gwas_zed[variant_iter] - gene_mediated_predicted_z[variant_iter] - nm_predicted_z[variant_iter] + old_variant_predicted_z_at_variant

			likelihood_precision = ((sqrt_gwas_sample_size*sqrt_gwas_sample_size)/tau)*variant_ld_var
			likelihood_mean_numer = (sqrt_gwas_sample_size/tau)*residual_gwas_zed_at_variant
			posterior_var = 1.0/(likelihood_precision + (1.0/nm_variant_mixture_var))
			posterior_mean = posterior_var*likelihood_mean_numer

			log_bayes_factor = 0.5*(np.log(posterior_var) - np.log(nm_variant_mixture_var)) + 0.5*(posterior_mean*posterior_mean/posterior_var)
			variant_inclusion_prob = _logistic_function_numba(log_prior_odds + log_bayes_factor)
			new_variant_effect = variant_inclusion_prob*posterior_mean
			variant_second_moment = variant_inclusion_prob*(posterior_var + (posterior_mean*posterior_mean))

		old_variant_effect = gene_nm_variant_effects[variant_iter]
		if old_variant_effect != new_variant_effect:
			delta_effect = new_variant_effect - old_variant_effect
			for variant_iter2 in range(n_variants):
				nm_predicted_z[variant_iter2] = nm_predicted_z[variant_iter2] + (sqrt_gwas_sample_size*delta_effect*LD[variant_iter2, variant_iter])

		gene_nm_variant_effects[variant_iter] = new_variant_effect
		gene_nm_variant_inclusion_probs[variant_iter] = variant_inclusion_prob
		gene_nm_variant_second_moments[variant_iter] = variant_second_moment

	return gene_nm_variant_effects, gene_nm_variant_inclusion_probs, gene_nm_variant_second_moments

def update_nm_variant_effects_vi(gene_alphas, training_data_obj, ordered_genes, EEs, gene_vars, nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments, nm_variant_mixture_prob, nm_variant_mixture_var, tau):
	for gene_iter, gene_name in enumerate(ordered_genes):
		LD = training_data_obj[gene_name]['LD']
		gwas_zed = training_data_obj[gene_name]['GWAS_Z']
		gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
		gene_mediated_predicted_z = compute_gene_mediated_predicted_z(gene_alphas[gene_iter, :], training_data_obj, gene_name, EEs, gene_vars[gene_iter])
		updated_nm_variant_effects, updated_nm_variant_inclusion_probs, updated_nm_variant_second_moments = _update_nm_variant_effects_for_single_gene_numba(np.ascontiguousarray(LD), np.ascontiguousarray(gwas_zed), np.sqrt(gwas_sample_size), np.ascontiguousarray(gene_mediated_predicted_z), nm_variant_effects[gene_iter], nm_variant_inclusion_probs[gene_iter], nm_variant_second_moments[gene_iter], nm_variant_mixture_prob, nm_variant_mixture_var, tau)
		nm_variant_effects[gene_iter] = updated_nm_variant_effects
		nm_variant_inclusion_probs[gene_iter] = updated_nm_variant_inclusion_probs
		nm_variant_second_moments[gene_iter] = updated_nm_variant_second_moments

	return nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments

def update_nm_variant_mixture_prior_vi(nm_variant_inclusion_probs, nm_variant_second_moments, pi_prior_a=.1, pi_prior_b=.1, variance_prior_a=1e-6, variance_prior_b=1e-6):
	total_variants = 0
	expected_num_included = 0.0
	expected_effect_sum_sq = 0.0
	for gene_iter in range(len(nm_variant_inclusion_probs)):
		total_variants = total_variants + len(nm_variant_inclusion_probs[gene_iter])
		expected_num_included = expected_num_included + np.sum(nm_variant_inclusion_probs[gene_iter])
		expected_effect_sum_sq = expected_effect_sum_sq + np.sum(nm_variant_second_moments[gene_iter])

	if total_variants == 0:
		print('No non-mediated variants available')
		pdb.set_trace()

	nm_variant_mixture_prob = (pi_prior_a + expected_num_included)/(pi_prior_a + pi_prior_b + total_variants)
	posterior_a = variance_prior_a + (expected_num_included/2.0)
	posterior_b = variance_prior_b + (expected_effect_sum_sq/2.0)
	nm_variant_mixture_var = posterior_b/(posterior_a + 1.0)
	if nm_variant_mixture_var <= 1e-12:
		nm_variant_mixture_var = 1e-12
	return nm_variant_mixture_prob, nm_variant_mixture_var

def update_gene_alphas_joint_genes_vi(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, resid_zed_array, per_gene_global_indices, tau):
	n_factors = gene_alphas.shape[1]

	for gene_iter, gene_name in enumerate(ordered_genes):
		LD = training_data_obj[gene_name]['LD']
		gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
		main_effects = training_data_obj[gene_name]['main_effects']
		interaction_effects = training_data_obj[gene_name]['interaction_effects']
		gene_indices = per_gene_global_indices[gene_iter]
		gene_deltas = []
		for factor_iter in range(n_factors):
			gene_deltas.append(compute_gene_delta(main_effects, interaction_effects, EEs[:, factor_iter], gene_vars[gene_iter]))
		gene_deltas = np.asarray(gene_deltas)

		for factor_iter in range(n_factors):
			old_alpha = gene_alphas[gene_iter, factor_iter]
			old_predicted_z = compute_gene_predicted_z(LD, gwas_sample_size, old_alpha, gene_deltas[factor_iter, :])
			residual_gwas_zed = resid_zed_array[gene_indices] + old_predicted_z
			new_alpha, alpha_inclusion_prob, alpha_second_moment = compute_alpha_variational_params(gene_deltas[factor_iter, :], residual_gwas_zed, LD, gwas_sample_size, mixture_probs[factor_iter], mixture_vars[factor_iter], tau)
			new_predicted_z = compute_gene_predicted_z(LD, gwas_sample_size, new_alpha, gene_deltas[factor_iter, :])
			resid_zed_array[gene_indices] = residual_gwas_zed - new_predicted_z
			gene_alphas[gene_iter, factor_iter] = new_alpha
			gene_alpha_inclusion_probs[gene_iter, factor_iter] = alpha_inclusion_prob
			gene_alpha_second_moments[gene_iter, factor_iter] = alpha_second_moment

	return gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, resid_zed_array

def update_EEs_joint_genes_vi(gene_alphas, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, e_var, gene_vars, resid_zed_array, per_gene_global_indices, tau):
	if e_var <= 0:
		print('E prior variance must be positive')
		pdb.set_trace()
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()

	n_components = EEs.shape[0]
	n_factors = EEs.shape[1]

	for factor_iter in range(n_factors):
		old_EE = np.copy(EEs[:, factor_iter])
		posterior_precision = np.eye(n_components)/e_var
		posterior_mean_numer = np.zeros(n_components)

		for gene_iter, gene_name in enumerate(ordered_genes):
			LD = training_data_obj[gene_name]['LD']
			gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
			main_effects = training_data_obj[gene_name]['main_effects']
			interaction_effects = training_data_obj[gene_name]['interaction_effects']
			gene_indices = per_gene_global_indices[gene_iter]
			alpha = gene_alphas[gene_iter, factor_iter]
			alpha_second_moment = gene_alpha_second_moments[gene_iter, factor_iter]
			gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])

			old_delta = compute_gene_delta(main_effects, interaction_effects, old_EE, gene_vars[gene_iter])
			old_predicted_z = compute_gene_predicted_z(LD, gwas_sample_size, alpha, old_delta)
			residual_gwas_zed = resid_zed_array[gene_indices] + old_predicted_z
			residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*alpha*gene_scale*np.dot(LD, main_effects))

			posterior_precision = posterior_precision + ((gwas_sample_size/tau)*alpha_second_moment*gene_scale*gene_scale*np.dot(np.transpose(interaction_effects), np.dot(LD, interaction_effects)))
			posterior_mean_numer = posterior_mean_numer + ((np.sqrt(gwas_sample_size)/tau)*alpha*gene_scale*np.dot(np.transpose(interaction_effects), residual_gwas_zed))

		posterior_mean = np.linalg.solve(posterior_precision, posterior_mean_numer)
		new_EE = posterior_mean

		for gene_iter, gene_name in enumerate(ordered_genes):
			LD = training_data_obj[gene_name]['LD']
			gwas_sample_size = training_data_obj[gene_name]['GWAS_N']
			main_effects = training_data_obj[gene_name]['main_effects']
			interaction_effects = training_data_obj[gene_name]['interaction_effects']
			gene_indices = per_gene_global_indices[gene_iter]
			alpha = gene_alphas[gene_iter, factor_iter]

			old_delta = compute_gene_delta(main_effects, interaction_effects, old_EE, gene_vars[gene_iter])
			new_delta = compute_gene_delta(main_effects, interaction_effects, new_EE, gene_vars[gene_iter])
			old_predicted_z = compute_gene_predicted_z(LD, gwas_sample_size, alpha, old_delta)
			new_predicted_z = compute_gene_predicted_z(LD, gwas_sample_size, alpha, new_delta)
			resid_zed_array[gene_indices] = resid_zed_array[gene_indices] + old_predicted_z - new_predicted_z

		EEs[:, factor_iter] = new_EE

	return EEs, resid_zed_array

def update_mixture_priors_vi(gene_alpha_inclusion_probs, gene_alpha_second_moments, pi_prior_a=.1, pi_prior_b=.1, variance_prior_a=1e-6, variance_prior_b=1e-6):
	n_genes = gene_alpha_inclusion_probs.shape[0]
	n_factors = gene_alpha_inclusion_probs.shape[1]
	mixture_probs = np.zeros(n_factors)
	mixture_vars = np.zeros(n_factors)

	for factor_iter in range(n_factors):
		expected_num_included = np.sum(gene_alpha_inclusion_probs[:, factor_iter])
		expected_alpha_sum_sq = np.sum(gene_alpha_second_moments[:, factor_iter])
		mixture_probs[factor_iter] = (pi_prior_a + expected_num_included)/(pi_prior_a + pi_prior_b + n_genes)

		posterior_a = variance_prior_a + (expected_num_included/2.0)
		posterior_b = variance_prior_b + (expected_alpha_sum_sq/2.0)
		mixture_vars[factor_iter] = posterior_b/(posterior_a + 1.0)
		if mixture_vars[factor_iter] <= 1e-12:
			mixture_vars[factor_iter] = 1e-12

	return mixture_probs, mixture_vars

def save_ee_posterior_summary(ee_samples, ee_summary_output_file):
	ee_samples = np.asarray(ee_samples)
	if ee_samples.shape[0] == 0:
		print('No post-burn-in EE samples available')
		pdb.set_trace()
	n_samples = ee_samples.shape[0]
	n_components = ee_samples.shape[1]
	n_factors = ee_samples.shape[2]

	t = open(ee_summary_output_file, 'w')
	t.write('factor_num\te_dimension\tn_samples\tsampled_mean\tsampled_variance\tci_lower_95\tci_upper_95\n')
	for factor_iter in range(n_factors):
		for component_iter in range(n_components):
			samples = ee_samples[:, component_iter, factor_iter]
			sampled_mean = np.mean(samples)
			if n_samples > 1:
				sampled_variance = np.var(samples, ddof=1)
			else:
				sampled_variance = 0.0
			ci_lower = np.percentile(samples, 2.5)
			ci_upper = np.percentile(samples, 97.5)
			t.write('factor_' + str(factor_iter + 1) + '\tE' + str(component_iter + 1) + '\t' + str(n_samples) + '\t' + str(sampled_mean) + '\t' + str(sampled_variance) + '\t' + str(ci_lower) + '\t' + str(ci_upper) + '\n')
	t.close()
	return

def save_ee_point_estimate_summary(EEs, ee_summary_output_file):
	t = open(ee_summary_output_file, 'w')
	t.write('factor_num\te_dimension\tn_samples\tsampled_mean\tsampled_variance\tci_lower_95\tci_upper_95\n')
	for factor_iter in range(EEs.shape[1]):
		for component_iter in range(EEs.shape[0]):
			ee_estimate = EEs[component_iter, factor_iter]
			t.write('factor_' + str(factor_iter + 1) + '\tE' + str(component_iter + 1) + '\t1\t' + str(ee_estimate) + '\t0.0\t' + str(ee_estimate) + '\t' + str(ee_estimate) + '\n')
	t.close()
	return

def save_ashr_prior_summary(ashr_grid_vars, ashr_weights, ashr_prior_output_file):
	t = open(ashr_prior_output_file, 'w')
	t.write('factor_num\tcomponent_num\tprior_sd\tprior_variance\tcomponent_weight\n')
	for factor_iter in range(ashr_weights.shape[0]):
		for component_iter in range(ashr_weights.shape[1]):
			prior_var = ashr_grid_vars[component_iter]
			t.write('factor_' + str(factor_iter + 1) + '\tcomponent_' + str(component_iter) + '\t' + str(np.sqrt(prior_var)) + '\t' + str(prior_var) + '\t' + str(ashr_weights[factor_iter, component_iter]) + '\n')
	t.close()
	return

def interaction_twas_vi_inference(training_data_obj, LL=1, e_var=0.05, burn_in_iters=500, max_iters=1000, ee_summary_output_file=None):
	tau = 1.0
	# Initialize data
	ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars = initialize_data(training_data_obj,LL)
	gene_alpha_inclusion_probs = np.zeros(gene_alphas.shape)
	gene_alpha_second_moments = np.zeros(gene_alphas.shape)
	ee_samples = []


	# Begin iterative algorithm
	for itera in range(max_iters):

		t1 = time.time()
		# Ok first update gene alphas
		gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments = update_gene_alphas_vi(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, tau)

		# Next update EEs
		EEs = update_EEs_vi(gene_alphas, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, e_var, gene_vars, tau)

		# Update variational empirical Bayes parameters
		mixture_probs, mixture_vars = update_mixture_priors_vi(gene_alpha_inclusion_probs, gene_alpha_second_moments)
		t2 = time.time()

		print('###########')
		print('itera ' + str(itera))
		print(EEs[:,0])
		print(mixture_probs[0])
		print(mixture_vars[0])
		print(np.mean(gene_alpha_inclusion_probs[:,0]))

		if itera >= burn_in_iters:
			ee_samples.append(np.copy(EEs))

	if ee_summary_output_file is not None:
		save_ee_posterior_summary(ee_samples, ee_summary_output_file)
	return gene_alphas, EEs, mixture_probs, mixture_vars, gene_alpha_inclusion_probs, gene_alpha_second_moments


def interaction_twas_vi_inference_with_nm_variants(training_data_obj, LL=1, e_var=0.05, max_iters=100, ee_summary_output_file=None):
	tau = 1.0
	# Initialize data
	ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars = initialize_data(training_data_obj,LL)
	gene_alpha_inclusion_probs = np.zeros(gene_alphas.shape)
	gene_alpha_second_moments = np.zeros(gene_alphas.shape)
	nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments = initialize_nm_variant_effects(training_data_obj, ordered_genes)
	nm_variant_mixture_prob = .1
	nm_variant_mixture_var = .00001


	# Begin iterative algorithm
	for itera in range(max_iters):

		t1 = time.time()
		# Update genetically mediated gene-trait effects conditional on non-mediated variants
		gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments = update_gene_alphas_vi_with_nm_variants(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, nm_variant_effects, mixture_probs, mixture_vars, gene_vars, tau)

		# Update EEs conditional on non-mediated variants
		EEs = update_EEs_vi_with_nm_variants(gene_alphas, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, nm_variant_effects, e_var, gene_vars, tau)

		# Update non-mediated variant effects
		nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments = update_nm_variant_effects_vi(gene_alphas, training_data_obj, ordered_genes, EEs, gene_vars, nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments, nm_variant_mixture_prob, nm_variant_mixture_var, tau)

		# Update variational empirical Bayes parameters
		mixture_probs, mixture_vars = update_mixture_priors_vi(gene_alpha_inclusion_probs, gene_alpha_second_moments)
		nm_variant_mixture_prob, nm_variant_mixture_var = update_nm_variant_mixture_prior_vi(nm_variant_inclusion_probs, nm_variant_second_moments)
		t2 = time.time()

		print('###########')
		print('itera ' + str(itera), flush=True)
		print(EEs[:,0], flush=True)
		print(mixture_probs[0], flush=True)
		print(mixture_vars[0], flush=True)
		print(nm_variant_mixture_prob, flush=True)
		print(nm_variant_mixture_var, flush=True)

	if ee_summary_output_file is not None:
		save_ee_point_estimate_summary(EEs, ee_summary_output_file)
	return gene_alphas, EEs, mixture_probs, mixture_vars, gene_alpha_inclusion_probs, gene_alpha_second_moments, nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments, nm_variant_mixture_prob, nm_variant_mixture_var


def interaction_twas_vi_inference_with_nm_variants_ashr_grid(training_data_obj, LL=1, e_var=0.05, max_iters=100, ee_summary_output_file=None, ashr_prior_output_file=None):
	tau = 1.0
	# Initialize data
	ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars = initialize_data(training_data_obj,LL)
	gene_alpha_inclusion_probs = np.zeros(gene_alphas.shape)
	gene_alpha_second_moments = np.zeros(gene_alphas.shape)
	ashr_grid_vars, ashr_weights = initialize_ashr_gene_prior(LL)
	gene_alpha_component_probs = np.zeros((gene_alphas.shape[0], gene_alphas.shape[1], len(ashr_grid_vars)))
	gene_alpha_component_probs[:, :, 0] = 1.0
	nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments = initialize_nm_variant_effects(training_data_obj, ordered_genes)
	nm_variant_mixture_prob = .1
	nm_variant_mixture_var = .00001


	# Begin iterative algorithm
	for itera in range(max_iters):

		t1 = time.time()
		# Update genetically mediated gene-trait effects conditional on non-mediated variants
		gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, gene_alpha_component_probs = update_gene_alphas_vi_with_nm_variants_ashr(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, gene_alpha_component_probs, training_data_obj, ordered_genes, EEs, nm_variant_effects, ashr_weights, ashr_grid_vars, gene_vars, tau)

		# Update EEs conditional on non-mediated variants
		EEs = update_EEs_vi_with_nm_variants(gene_alphas, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, nm_variant_effects, e_var, gene_vars, tau)

		# Update non-mediated variant effects
		nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments = update_nm_variant_effects_vi(gene_alphas, training_data_obj, ordered_genes, EEs, gene_vars, nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments, nm_variant_mixture_prob, nm_variant_mixture_var, tau)

		# Update variational empirical Bayes parameters
		ashr_weights = update_ashr_gene_prior_weights_vi(gene_alpha_component_probs)
		nm_variant_mixture_prob, nm_variant_mixture_var = update_nm_variant_mixture_prior_vi(nm_variant_inclusion_probs, nm_variant_second_moments)
		t2 = time.time()

		print('###########')
		print('itera ' + str(itera), flush=True)
		print(EEs[:,0], flush=True)
		print(np.sum(ashr_weights[0, ashr_grid_vars > 0.0]), flush=True)
		print(ashr_weights[0, :], flush=True)
		print(nm_variant_mixture_prob, flush=True)
		print(nm_variant_mixture_var, flush=True)

	if ee_summary_output_file is not None:
		save_ee_point_estimate_summary(EEs, ee_summary_output_file)
	if ashr_prior_output_file is not None:
		save_ashr_prior_summary(ashr_grid_vars, ashr_weights, ashr_prior_output_file)
	return gene_alphas, EEs, ashr_weights, ashr_grid_vars, gene_alpha_inclusion_probs, gene_alpha_second_moments, gene_alpha_component_probs, nm_variant_effects, nm_variant_inclusion_probs, nm_variant_second_moments, nm_variant_mixture_prob, nm_variant_mixture_var





def initialize_data_joint_genes(training_data_obj,LL):
	n_genes = len(training_data_obj)
	ordered_genes = np.sort([*training_data_obj])

	n_components = training_data_obj[ordered_genes[0]]['interaction_effects'].shape[1]
	gene_alphas = np.zeros((n_genes, LL))
	EEs = np.zeros((n_components, LL))

	if np.sum(gene_alphas!=0.0) != 0:
		print('need to update residualizations for non-zero mean')
		pdb.set_trace()

	mixture_probs = np.ones(LL)*.5
	mixture_vars = np.ones(LL)*.001

	gene_vars = []
	for gene_id in ordered_genes:
		gene_var = np.dot(np.dot(training_data_obj[gene_id]['main_effects'], training_data_obj[gene_id]['LD']), training_data_obj[gene_id]['main_effects'])
		if np.isfinite(gene_var) == False or gene_var <= 0:
			print('Gene variance must be positive')
			pdb.set_trace()
		gene_vars.append(gene_var)
	gene_vars = np.asarray(gene_vars)


	# Create mapping from variant id to gwas z
	var_to_gwas_z = {}
	for gene_id in ordered_genes:
		gwas_z = training_data_obj[gene_id]['GWAS_Z']
		variant_ids = training_data_obj[gene_id]['variant_ids']

		for ii,variant_id in enumerate(variant_ids):
			if variant_id in var_to_gwas_z:
				var_to_gwas_z[variant_id].append(gwas_z[ii])
			else:
				var_to_gwas_z[variant_id] = [gwas_z[ii]]

	var_to_index = {}
	variant_array = []
	zed_array = []

	for indexer, variant_id in enumerate([*var_to_gwas_z]):
		var_to_index[variant_id] = indexer
		zed_array.append(np.mean(var_to_gwas_z[variant_id]))
		variant_array.append(variant_id)
	resid_zed_array = np.asarray(zed_array)
	variant_array = np.asarray(variant_array)

	# Create 
	per_gene_global_indices = []
	for gene_id in ordered_genes:
		tmp_arr = []
		for variant_id in training_data_obj[gene_id]['variant_ids']:
			tmp_arr.append(var_to_index[variant_id])
		tmp_arr = np.asarray(tmp_arr)
		per_gene_global_indices.append(tmp_arr)



	return ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars, resid_zed_array, per_gene_global_indices


def interaction_twas_joint_genes_vi_inference(training_data_obj, LL=1, e_var=0.05, burn_in_iters=500, max_iters=1000, ee_summary_output_file=None):
	tau = 1.0
	# Initialize data
	ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars, resid_zed_array, per_gene_global_indices = initialize_data_joint_genes(training_data_obj,LL)
	gene_alpha_inclusion_probs = np.zeros(gene_alphas.shape)
	gene_alpha_second_moments = np.zeros(gene_alphas.shape)
	ee_samples = []


	# Begin iterative algorithm
	for itera in range(max_iters):

		t1 = time.time()
		# Ok first update gene alphas
		gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, resid_zed_array = update_gene_alphas_joint_genes_vi(gene_alphas, gene_alpha_inclusion_probs, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, resid_zed_array, per_gene_global_indices, tau)

		# Next update EEs
		EEs, resid_zed_array = update_EEs_joint_genes_vi(gene_alphas, gene_alpha_second_moments, training_data_obj, ordered_genes, EEs, e_var, gene_vars, resid_zed_array, per_gene_global_indices, tau)

		# Update variational empirical Bayes parameters
		mixture_probs, mixture_vars = update_mixture_priors_vi(gene_alpha_inclusion_probs, gene_alpha_second_moments)
		t2 = time.time()

		print('###########')
		print('itera ' + str(itera), flush=True)
		print(EEs[:,0], flush=True)
		print(mixture_probs[0], flush=True)
		print(mixture_vars[0], flush=True)
		print(np.mean(gene_alpha_inclusion_probs[:,0]), flush=True)

		if itera >= burn_in_iters:
			ee_samples.append(np.copy(EEs))

	if ee_summary_output_file is not None:
		save_ee_posterior_summary(ee_samples, ee_summary_output_file)
	return gene_alphas, EEs, mixture_probs, mixture_vars, gene_alpha_inclusion_probs, gene_alpha_second_moments




####################
# Command line args
####################
pred_expresssion_output_stem = sys.argv[1]
trait_sumstat_file = sys.argv[2]
onek_genomes_plink_files = sys.argv[3]
interaction_twas_output_stem = sys.argv[4]
inference_version = sys.argv[5]
tau = 1.0

training_data_pickle_file = interaction_twas_output_stem + '_training_data_obj.pkl'
valid_training_chroms = np.arange(5,8)
valid_testing_chroms = np.arange(1,5)

training_data_obj = load_in_training_data(valid_training_chroms, pred_expresssion_output_stem, trait_sumstat_file, onek_genomes_plink_files)
#pickle.dump(training_data_obj, open(training_data_pickle_file, 'wb'))
#training_data_obj = pickle.load(open(training_data_pickle_file, 'rb'))


ee_summary_output_file = interaction_twas_output_stem + '_ee_posterior_summary.txt'
if inference_version == 'independent_genes':
	interaction_twas_vi_inference(training_data_obj, LL=1, e_var=1.0, ee_summary_output_file=ee_summary_output_file, burn_in_iters=600, max_iters=1000)
elif inference_version == 'independent_genes_nm_variants':
	interaction_twas_vi_inference_with_nm_variants(training_data_obj, LL=1, e_var=100.0, ee_summary_output_file=ee_summary_output_file, max_iters=400)
elif inference_version == 'independent_genes_nm_variants_ashr_grid' or inference_version == 'independent_genes_nm_variants_ashr':
	ashr_prior_output_file = interaction_twas_output_stem + '_gene_ashr_prior_summary.txt'
	interaction_twas_vi_inference_with_nm_variants_ashr_grid(training_data_obj, LL=1, e_var=100.0, ee_summary_output_file=ee_summary_output_file, ashr_prior_output_file=ashr_prior_output_file, max_iters=400)
elif inference_version == 'joint_genes':
	interaction_twas_joint_genes_vi_inference(training_data_obj, LL=1, e_var=.1,ee_summary_output_file=ee_summary_output_file, burn_in_iters=600, max_iters=1000)

# 
held_out_twas_analysis_output_file = interaction_twas_output_stem + '_held_out_twas.txt'
inferred_ee = np.loadtxt(ee_summary_output_file, dtype=str, delimiter='\t')[1:,3].astype(float)
run_held_out_twas(valid_testing_chroms, pred_expresssion_output_stem, trait_sumstat_file, onek_genomes_plink_files, inferred_ee, held_out_twas_analysis_output_file)


print(held_out_twas_analysis_output_file)
