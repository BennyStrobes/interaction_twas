import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink
import time
import pickle



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


	return ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars


def sample_alpha_from_spike_and_slab(gene_delta, residual_gwas_zed, LD, gwas_sample_size, mixture_prob, mixture_var, tau):
	if mixture_var <= 0:
		print('Mixture variance must be positive')
		pdb.set_trace()
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()
	mixture_prob = np.clip(mixture_prob, 1e-12, 1.0 - 1e-12)

	delta_R_delta = np.dot(gene_delta, np.dot(LD, gene_delta))
	if np.isfinite(delta_R_delta) == False or delta_R_delta <= 0:
		return 0.0

	likelihood_precision = (gwas_sample_size/tau)*delta_R_delta
	likelihood_mean_numer = (np.sqrt(gwas_sample_size)/tau)*np.dot(gene_delta, residual_gwas_zed)
	posterior_var = 1.0/(likelihood_precision + (1.0/mixture_var))
	posterior_mean = posterior_var*likelihood_mean_numer

	log_bayes_factor = 0.5*(np.log(posterior_var) - np.log(mixture_var)) + 0.5*(posterior_mean*posterior_mean/posterior_var)
	log_prior_odds = np.log(mixture_prob) - np.log(1.0 - mixture_prob)
	log_posterior_odds = log_prior_odds + log_bayes_factor
	posterior_inclusion_prob = 1.0/(1.0 + np.exp(-log_posterior_odds))

	if np.random.rand() > posterior_inclusion_prob:
		return 0.0
	return np.random.normal(posterior_mean, np.sqrt(posterior_var))

def update_gene_alphas(gene_alphas, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, tau):
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

			gene_alphas[gene_iter, factor_iter] = sample_alpha_from_spike_and_slab(gene_deltas[factor_iter, :], residual_gwas_zed, LD, gwas_sample_size, mixture_probs[factor_iter], mixture_vars[factor_iter], tau)

	return gene_alphas

def update_EEs(gene_alphas, training_data_obj, ordered_genes, EEs, e_var, gene_vars, tau):
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
			gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])

			residual_gwas_zed = np.copy(gwas_zed)
			for other_factor_iter in range(n_factors):
				if other_factor_iter == factor_iter:
					continue
				other_delta = gene_scale*(main_effects + np.dot(interaction_effects, EEs[:, other_factor_iter]))
				residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*gene_alphas[gene_iter, other_factor_iter]*np.dot(LD, other_delta))
			residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*alpha*gene_scale*np.dot(LD, main_effects))

			posterior_precision = posterior_precision + ((gwas_sample_size/tau)*alpha*alpha*gene_scale*gene_scale*np.dot(np.transpose(interaction_effects), np.dot(LD, interaction_effects)))
			posterior_mean_numer = posterior_mean_numer + ((np.sqrt(gwas_sample_size)/tau)*alpha*gene_scale*np.dot(np.transpose(interaction_effects), residual_gwas_zed))

		posterior_mean = np.linalg.solve(posterior_precision, posterior_mean_numer)
		posterior_cov = np.linalg.inv(posterior_precision)
		posterior_cov = (posterior_cov + np.transpose(posterior_cov))/2.0
		EEs[:, factor_iter] = np.random.multivariate_normal(posterior_mean, posterior_cov)

	return EEs

def compute_gene_delta(main_effects, interaction_effects, EE, gene_var):
	gene_scale = 1.0/np.sqrt(gene_var)
	return gene_scale*(main_effects + np.dot(interaction_effects, EE))

def compute_gene_predicted_z(LD, gwas_sample_size, gene_alpha, gene_delta):
	return np.sqrt(gwas_sample_size)*gene_alpha*np.dot(LD, gene_delta)

def update_gene_alphas_joint_genes(gene_alphas, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, resid_zed_array, per_gene_global_indices, tau):
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
			new_alpha = sample_alpha_from_spike_and_slab(gene_deltas[factor_iter, :], residual_gwas_zed, LD, gwas_sample_size, mixture_probs[factor_iter], mixture_vars[factor_iter], tau)
			new_predicted_z = compute_gene_predicted_z(LD, gwas_sample_size, new_alpha, gene_deltas[factor_iter, :])
			resid_zed_array[gene_indices] = residual_gwas_zed - new_predicted_z
			gene_alphas[gene_iter, factor_iter] = new_alpha

	return gene_alphas, resid_zed_array

def update_EEs_joint_genes(gene_alphas, training_data_obj, ordered_genes, EEs, e_var, gene_vars, resid_zed_array, per_gene_global_indices, tau):
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
			gene_scale = 1.0/np.sqrt(gene_vars[gene_iter])

			old_delta = compute_gene_delta(main_effects, interaction_effects, old_EE, gene_vars[gene_iter])
			old_predicted_z = compute_gene_predicted_z(LD, gwas_sample_size, alpha, old_delta)
			residual_gwas_zed = resid_zed_array[gene_indices] + old_predicted_z
			residual_gwas_zed = residual_gwas_zed - (np.sqrt(gwas_sample_size)*alpha*gene_scale*np.dot(LD, main_effects))

			posterior_precision = posterior_precision + ((gwas_sample_size/tau)*alpha*alpha*gene_scale*gene_scale*np.dot(np.transpose(interaction_effects), np.dot(LD, interaction_effects)))
			posterior_mean_numer = posterior_mean_numer + ((np.sqrt(gwas_sample_size)/tau)*alpha*gene_scale*np.dot(np.transpose(interaction_effects), residual_gwas_zed))

		posterior_mean = np.linalg.solve(posterior_precision, posterior_mean_numer)
		posterior_cov = np.linalg.inv(posterior_precision)
		posterior_cov = (posterior_cov + np.transpose(posterior_cov))/2.0
		new_EE = np.random.multivariate_normal(posterior_mean, posterior_cov)

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

def update_mixture_priors(gene_alphas, pi_prior_a=.1, pi_prior_b=.1, variance_prior_a=1e-6, variance_prior_b=1e-6):
	n_genes = gene_alphas.shape[0]
	n_factors = gene_alphas.shape[1]
	mixture_probs = np.zeros(n_factors)
	mixture_vars = np.zeros(n_factors)

	for factor_iter in range(n_factors):
		included_genes = gene_alphas[:, factor_iter] != 0.0
		num_included = np.sum(included_genes)
		mixture_probs[factor_iter] = np.random.beta(pi_prior_a + num_included, pi_prior_b + n_genes - num_included)

		posterior_a = variance_prior_a + (num_included/2.0)
		posterior_b = variance_prior_b + (np.sum(np.square(gene_alphas[included_genes, factor_iter]))/2.0)
		mixture_vars[factor_iter] = 1.0/np.random.gamma(posterior_a, scale=(1.0/posterior_b))

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

def interaction_twas_gibbs_inference(training_data_obj, LL=1, e_var=0.05, tau=1.0, burn_in_iters=500, max_iters=1000, ee_summary_output_file=None):
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()
	# Initialize data
	ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars = initialize_data(training_data_obj,LL)
	ee_samples = []


	# Begin iterative algorithm
	for itera in range(max_iters):

		t1 = time.time()
		# Ok first update gene alphas
		gene_alphas = update_gene_alphas(gene_alphas, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, tau)

		# Next update EEs
		EEs = update_EEs(gene_alphas, training_data_obj, ordered_genes, EEs, e_var, gene_vars, tau)

		# Update spike-and-slab prior parameters
		mixture_probs, mixture_vars = update_mixture_priors(gene_alphas)
		t2 = time.time()

		print('###########')
		print('itera ' + str(itera))
		print(EEs[:,0])
		print(mixture_probs[0])
		print(mixture_vars[0])

		if itera >= burn_in_iters:
			ee_samples.append(np.copy(EEs))

	if ee_summary_output_file is not None:
		save_ee_posterior_summary(ee_samples, ee_summary_output_file)
	return gene_alphas, EEs, mixture_probs, mixture_vars





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


def interaction_twas_joint_genes_gibbs_inference(training_data_obj, LL=1, e_var=0.05, tau=1.0, burn_in_iters=500, max_iters=1000, ee_summary_output_file=None):
	if tau <= 0:
		print('Tau must be positive')
		pdb.set_trace()
	# Initialize data
	ordered_genes, gene_alphas, EEs, mixture_probs, mixture_vars, gene_vars, resid_zed_array, per_gene_global_indices = initialize_data_joint_genes(training_data_obj,LL)
	ee_samples = []


	# Begin iterative algorithm
	for itera in range(max_iters):

		t1 = time.time()
		# Ok first update gene alphas
		gene_alphas, resid_zed_array = update_gene_alphas_joint_genes(gene_alphas, training_data_obj, ordered_genes, EEs, mixture_probs, mixture_vars, gene_vars, resid_zed_array, per_gene_global_indices, tau)

		# Next update EEs
		EEs, resid_zed_array = update_EEs_joint_genes(gene_alphas, training_data_obj, ordered_genes, EEs, e_var, gene_vars, resid_zed_array, per_gene_global_indices, tau)

		# Update spike-and-slab prior parameters
		mixture_probs, mixture_vars = update_mixture_priors(gene_alphas)
		t2 = time.time()

		print('###########')
		print('itera ' + str(itera), flush=True)
		print(EEs[:,0], flush=True)
		print(mixture_probs[0], flush=True)
		print(mixture_vars[0], flush=True)

		if itera >= burn_in_iters:
			ee_samples.append(np.copy(EEs))

	if ee_summary_output_file is not None:
		save_ee_posterior_summary(ee_samples, ee_summary_output_file)
	return gene_alphas, EEs, mixture_probs, mixture_vars




####################
# Command line args
####################
pred_expresssion_output_stem = sys.argv[1]
trait_sumstat_file = sys.argv[2]
onek_genomes_plink_files = sys.argv[3]
interaction_twas_output_stem = sys.argv[4]
inference_version = sys.argv[5]
if len(sys.argv) > 6:
	tau = float(sys.argv[6])
else:
	tau = 1.0
if tau <= 0:
	print('Tau must be positive')
	pdb.set_trace()

training_data_pickle_file = interaction_twas_output_stem + '_training_data_obj.pkl'
valid_training_chroms = np.arange(8,23)
print('need to change chroms')
valid_testing_chroms = np.arange(1,8)

training_data_obj = load_in_training_data(valid_training_chroms, pred_expresssion_output_stem, trait_sumstat_file, onek_genomes_plink_files)
#pickle.dump(training_data_obj, open(training_data_pickle_file, 'wb'))
#training_data_obj = pickle.load(open(training_data_pickle_file, 'rb'))


ee_summary_output_file = interaction_twas_output_stem + '_ee_posterior_summary.txt'
if inference_version == 'independent_genes':
	interaction_twas_gibbs_inference(training_data_obj, LL=1, tau=tau, ee_summary_output_file=ee_summary_output_file, burn_in_iters=600, max_iters=1000)
elif inference_version == 'joint_genes':
	interaction_twas_joint_genes_gibbs_inference(training_data_obj, LL=1, tau=tau, ee_summary_output_file=ee_summary_output_file, burn_in_iters=600, max_iters=1000)

# 
held_out_twas_analysis_output_file = interaction_twas_output_stem + '_held_out_twas.txt'
inferred_ee = np.loadtxt(ee_summary_output_file, dtype=str, delimiter='\t')[1:,3].astype(float)
run_held_out_twas(valid_testing_chroms, pred_expresssion_output_stem, trait_sumstat_file, onek_genomes_plink_files, inferred_ee, held_out_twas_analysis_output_file)


print(held_out_twas_analysis_output_file)
