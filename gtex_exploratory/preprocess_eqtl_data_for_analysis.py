import numpy as np
import os
import sys
import pdb
import gzip
from statistics import NormalDist
import rnaseqnorm
import pandas as pd

def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf):
	valid_chroms = {}
	for chrom_num in range(1,23):
		valid_chroms['chr' + str(chrom_num)] = 1


	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in valid_chroms:
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = (data[0], data[3], data[4], data[6])

	f.close()

	return dicti


def load_in_ancestry_labels(gtex_subject_attributes_file):
	f = gzip.open(gtex_subject_attributes_file, 'rt')
	mapping = {}
	arr = []
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		if line == '':
			continue
		if line.startswith('dbGaP_S'):
			continue
		data = line.split('\t')
		'''
		if len(data) != 190:
			print('assumption eroror')
			pdb.set_trace()
		'''
		subject_id = data[1]
		race = data[5]
		if subject_id in mapping:
			print('asssumption erororo')
			pdb.set_trace()
		mapping[subject_id] = race
		arr.append(race)
	f.close()
	arr = np.asarray(arr)
	return mapping

def load_in_tissue_data(tissue_names_file):
	f = open(tissue_names_file)
	tissue_names = []
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 2:
			print('assumptioneornronr')
			pdb.set_trace()
		tissue_names.append(data[0])
		mapping[data[1]] = data[0]
	f.close()

	return np.asarray(tissue_names), mapping

def create_mapping_from_gtex_sample_id_to_ind_tissue_id(gtex_sample_attributes_file, alt_gtex_tissue_to_gtex_tissue, valid_ind_tissue_ids):
	f = open(gtex_sample_attributes_file)
	mapping = {}
	used = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0].startswith('GTEX') == False:
			continue
		if len(data) < 18:
			continue
		gtex_sample_id = data[0]
		alt_tissue_id = data[6]
		typer = data[17]

		if typer != 'RNASEQ':
			continue
		if alt_tissue_id not in alt_gtex_tissue_to_gtex_tissue:
			continue
		gtex_tissue = alt_gtex_tissue_to_gtex_tissue[alt_tissue_id]
		gtex_id = gtex_sample_id.split('-')[0] + '-' + gtex_sample_id.split('-')[1]
		ind_tissue_id = gtex_id + ':' + gtex_tissue
		if ind_tissue_id in used:
			print('assumptioneornor')
			pdb.set_trace()
		if gtex_sample_id in mapping:
			print('assumption oeroror')
			pdb.set_trace()
		if ind_tissue_id not in valid_ind_tissue_ids:
			continue
		used[ind_tissue_id] = 0
		mapping[gtex_sample_id] = ind_tissue_id

	f.close()

	return mapping, used

def eqtl_sample_checking(ordered_gtex_tissues, gtex_per_tissue_covariate_dir, ind_tissue_ids):
	for gtex_tissue in ordered_gtex_tissues:
		cov_file = gtex_per_tissue_covariate_dir + gtex_tissue + '.v10.covariates.txt'
		data_tmp = np.loadtxt(cov_file,dtype=str,delimiter='\t')
		ind_ids = data_tmp[0,1:]

		for ind_id in ind_ids:
			ind_tissue_id = ind_id + ':' + gtex_tissue
			if ind_tissue_id not in ind_tissue_ids:
				print('assumption error')
				pdb.set_trace()
			ind_tissue_ids[ind_tissue_id] = 1

	for ind_tissue_id in [*ind_tissue_ids]:
		if ind_tissue_ids[ind_tissue_id] == 0:
			print('erroror')
			pdb.set_trace()
	return

def extract_ind_tissue_ids_used_for_eqtl_analysis(ordered_gtex_tissues, gtex_per_tissue_covariate_dir):
	valid_ind_tissue_ids = {}
	for gtex_tissue in ordered_gtex_tissues:
		cov_file = gtex_per_tissue_covariate_dir + gtex_tissue + '.v10.covariates.txt'
		data_tmp = np.loadtxt(cov_file,dtype=str,delimiter='\t')
		ind_ids = data_tmp[0,1:]

		for ind_id in ind_ids:
			ind_tissue_id = ind_id + ':' + gtex_tissue

			if ind_tissue_id in valid_ind_tissue_ids:
				print('assumptioneorroro')
				pdb.set_trace()
			valid_ind_tissue_ids[ind_tissue_id] = 1
	return valid_ind_tissue_ids


def generate_filtered_log_tpm_expression(filtered_log_tpm_expression_file,gtex_tpm_expression, gtex_sample_id_to_ind_tissue_id, gtex_id_to_ancestry_labels, pc_genes ):
	t = open(filtered_log_tpm_expression_file,'w')
	f = gzip.open(gtex_tpm_expression,'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if line.startswith('#') or line.startswith('59033'):
			continue

		if head_count == 0:
			head_count = head_count + 1
			all_sample_names = np.asarray(data[2:])
			valid_samples = []
			new_sample_names = []

			for ii, sample_name in enumerate(all_sample_names):
				if sample_name not in gtex_sample_id_to_ind_tissue_id:
					continue
				ind_tissue_id = gtex_sample_id_to_ind_tissue_id[sample_name]
				ind_id = ind_tissue_id.split(':')[0]
				ancestry_label = gtex_id_to_ancestry_labels[ind_id]
				if ancestry_label != '3':
					continue
				valid_samples.append(ii)
				new_sample_names.append(ind_tissue_id)
			valid_samples = np.asarray(valid_samples)
			new_sample_names = np.asarray(new_sample_names)
			t.write('Ensamble_id\tgene_symbol_id\tgene_chrom\tgene_start\tgene_end\tgene_strand' + '\t' + '\t'.join(new_sample_names) + '\n')
			continue

		# Extract lines from field
		full_ensamble_id = data[0]
		gene_id = data[1]
		short_ensample_id = full_ensamble_id.split('.')[0]
		if short_ensample_id not in pc_genes:
			continue
		gene_info = pc_genes[short_ensample_id]
		gene_chrom, gene_tss, gene_tes, gene_strand = pc_genes[short_ensample_id]

		expression_string = np.asarray(data[2:])[valid_samples]
		expression = expression_string.astype(float)
		if np.sum(expression > .1)/len(expression) < .1:
			continue
		
		log_expression = np.log2(expression + 1)

		t.write(full_ensamble_id + '\t' + gene_id + '\t' + gene_chrom + '\t' + gene_tss + '\t' + gene_tes + '\t' + gene_strand + '\t' + '\t'.join(log_expression.astype(str)) + '\n')


	f.close()
	t.close()

	return

def load_expression_file(expression_file):
	f = open(expression_file)
	header = f.readline().rstrip().split('\t')
	gene_infos = []
	expression = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		gene_infos.append(data[:6])
		expression.append(data[6:])
	f.close()

	return header, np.asarray(gene_infos), np.asarray(expression).astype(float)

def print_expression_file(output_file, header, gene_infos, expression):
	t = open(output_file,'w')
	t.write('\t'.join(header) + '\n')
	for gene_info, expr in zip(gene_infos, expression):
		t.write('\t'.join(gene_info) + '\t' + '\t'.join(expr.astype(str)) + '\n')
	t.close()
	return

def convert_observed_values_to_rank_values(observed_values, rank_values):
	num_values = len(observed_values)
	ordered_indices = np.argsort(observed_values, kind='mergesort')
	ordered_values = observed_values[ordered_indices]
	new_values = np.zeros(num_values)

	start_index = 0
	while start_index < num_values:
		end_index = start_index + 1
		while end_index < num_values and ordered_values[end_index] == ordered_values[start_index]:
			end_index = end_index + 1
		new_values[ordered_indices[start_index:end_index]] = np.mean(rank_values[start_index:end_index])
		start_index = end_index

	return new_values

def quantile_normalize_expression(input_expression_file, quantile_normalized_expression_file):
	header, gene_infos, expression = load_expression_file(input_expression_file)

	ordered_expression = np.sort(expression, axis=0)
	mean_ordered_expression = np.mean(ordered_expression, axis=1)
	quantile_normalized_expression = np.zeros(expression.shape)

	for sample_num in range(expression.shape[1]):
		quantile_normalized_expression[:, sample_num] = convert_observed_values_to_rank_values(expression[:, sample_num], mean_ordered_expression)


	print_expression_file(quantile_normalized_expression_file, header, gene_infos, quantile_normalized_expression)
	return

def inverse_normal_transform_expression(input_expression_file, inverse_normal_transformed_expression_file):
	header, gene_infos, expression = load_expression_file(input_expression_file)

	normal_dist = NormalDist()
	num_samples = expression.shape[1]
	normal_scores = []
	for sample_num in range(num_samples):
		normal_scores.append(normal_dist.inv_cdf((sample_num + .5)/num_samples))
	normal_scores = np.asarray(normal_scores)

	inverse_normal_transformed_expression = np.zeros(expression.shape)
	for gene_num in range(expression.shape[0]):
		inverse_normal_transformed_expression[gene_num, :] = convert_observed_values_to_rank_values(expression[gene_num, :], normal_scores)

	print_expression_file(inverse_normal_transformed_expression_file, header, gene_infos, inverse_normal_transformed_expression)
	return

def compute_expression_pcs(input_expression_file, expression_pc_output_file, num_pcs):
	header, gene_infos, expression = load_expression_file(input_expression_file)
	sample_names = np.asarray(header[6:])

	# PCA is run with samples as observations and genes as features.
	sample_gene_expression = np.transpose(expression)
	centered_expression = sample_gene_expression - np.mean(sample_gene_expression, axis=0)
	u, s, vh = np.linalg.svd(centered_expression, full_matrices=False)
	pc_loadings = u[:, :num_pcs]*s[:num_pcs]
	variance_explained = (s*s)/np.sum(s*s)

	t = open(expression_pc_output_file,'w')
	t.write('sample_id')
	for pc_num in range(num_pcs):
		t.write('\tPC' + str(pc_num + 1))
	t.write('\n')
	for sample_num, sample_name in enumerate(sample_names):
		t.write(sample_name + '\t' + '\t'.join(pc_loadings[sample_num, :].astype(str)) + '\n')
	t.close()

	variance_output_file = expression_pc_output_file + '.variance_explained.txt'
	t = open(variance_output_file,'w')
	t.write('pc_num\tvariance_explained\n')
	for pc_num in range(num_pcs):
		t.write('PC' + str(pc_num + 1) + '\t' + str(variance_explained[pc_num]) + '\n')
	t.close()
	return

def compute_umap_from_expression_pcs(expression_pc_file, umap_output_file):
	import umap

	f = open(expression_pc_file)
	head_count = 0
	sample_names = []
	pc_values = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_names.append(data[0])
		pc_values.append(data[1:])
	f.close()

	sample_names = np.asarray(sample_names)
	pc_values = np.asarray(pc_values).astype(float)
	umap_model = umap.UMAP(n_components=2, n_neighbors=15, min_dist=.1, metric='euclidean', random_state=1)
	embedding = umap_model.fit_transform(pc_values)

	t = open(umap_output_file,'w')
	t.write('sample_id\tUMAP1\tUMAP2\n')
	for sample_num, sample_name in enumerate(sample_names):
		t.write(sample_name + '\t' + str(embedding[sample_num, 0]) + '\t' + str(embedding[sample_num, 1]) + '\n')
	t.close()
	return

def generate_umap_tissue_identity_scatterplot(umap_output_file, umap_plot_file):
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

	f = open(umap_output_file)
	head_count = 0
	tissue_names = []
	umap1 = []
	umap2 = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_id = data[0]
		tissue_name = sample_id.split(':', 1)[1]
		tissue_names.append(tissue_name)
		umap1.append(float(data[1]))
		umap2.append(float(data[2]))
	f.close()

	tissue_names = np.asarray(tissue_names)
	umap1 = np.asarray(umap1)
	umap2 = np.asarray(umap2)
	unique_tissue_names = np.sort(np.unique(tissue_names))
	tissue_to_number = {}
	for tissue_num, tissue_name in enumerate(unique_tissue_names):
		tissue_to_number[tissue_name] = tissue_num
	tissue_numbers = np.asarray([tissue_to_number[tissue_name] for tissue_name in tissue_names])

	plt.figure(figsize=(12, 8))
	scatter = plt.scatter(umap1, umap2, c=tissue_numbers, cmap='gist_ncar', s=8, alpha=.75, rasterized=True)
	plt.xlabel('UMAP1')
	plt.ylabel('UMAP2')
	handles, labels = scatter.legend_elements(num=len(unique_tissue_names))
	plt.legend(handles, unique_tissue_names, markerscale=2, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=6, frameon=False)
	plt.tight_layout()
	plt.savefig(umap_plot_file, dpi=200)
	plt.close()
	return


def create_mapping_from_genotype_indices_to_expression_sample(ordered_expression_samples, processed_genotype_data_dir, genotype_indices_to_expression_file):
	fam_file = processed_genotype_data_dir + 'gtex_v9_eqtl_chr6.fam'
	f = open(fam_file)
	ind_id_to_geno_index = {}
	geno_index = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		gtex_ind_id = data[1].split('-')[0] + '-' + data[1].split('-')[1]
		if gtex_ind_id in ind_id_to_geno_index:
			print('assumptioneornro')
			pdb.set_trace()
		ind_id_to_geno_index[gtex_ind_id] = geno_index
		geno_index = geno_index+1
	f.close()

	t = open(genotype_indices_to_expression_file,'w')

	for expression_sample in ordered_expression_samples:
		gtex_ind_id = expression_sample.split(':')[0]
		geno_index = ind_id_to_geno_index[gtex_ind_id]
		t.write(str(geno_index) + '\n')

	t.close()
	return

def limit_expression_pcs_to_top_n(expression_pc_output_file, limited_expression_pc_output_file, top_n):
	f = open(expression_pc_output_file)
	header = f.readline().rstrip().split('\t')
	if len(header) < (top_n + 1):
		print('PC file must have sample_id plus at least ' + str(top_n) + ' PC columns')
		pdb.set_trace()
	header = header[:(top_n + 1)]

	sample_names = []
	top_pcs = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) < (top_n + 1):
			print('PC file row must have sample_id plus at least ' + str(top_n) + ' PC columns')
			pdb.set_trace()
		sample_names.append(data[0])
		top_pcs.append(data[1:(top_n + 1)])
	f.close()

	sample_names = np.asarray(sample_names)
	top_pcs = np.asarray(top_pcs).astype(float)
	for pc_num in range(top_n):
		pc_mean = np.mean(top_pcs[:, pc_num])
		pc_std = np.std(top_pcs[:, pc_num])
		if pc_std == 0:
			print('PC has zero variance')
			pdb.set_trace()
		top_pcs[:, pc_num] = (top_pcs[:, pc_num] - pc_mean)/pc_std

	t = open(limited_expression_pc_output_file,'w')
	t.write('\t'.join(header) + '\n')
	for sample_num, sample_name in enumerate(sample_names):
		t.write(sample_name + '\t' + '\t'.join(top_pcs[sample_num, :].astype(str)) + '\n')
	t.close()
	return







#######################
# Commmand line args
########################
tissue_names_file = sys.argv[1]
gtex_v10_pc_genes_gtf = sys.argv[2]
gtex_tpm_expression = sys.argv[3]
gtex_subject_attributes_file = sys.argv[4]
gtex_sample_attributes_file = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
gtex_per_tissue_covariate_dir = sys.argv[7]
processed_eqtl_data_dir = sys.argv[8]

# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)

# Create mapping from gtex id to ancestry label
# European is '3''
gtex_id_to_ancestry_labels = load_in_ancestry_labels(gtex_subject_attributes_file)

# Extract gtex tissue names
ordered_gtex_tissues, alt_gtex_tissue_to_gtex_tissue = load_in_tissue_data(tissue_names_file)

# Extract gtex individual: tissue pairs that were used for eQTL analysis
valid_ind_tissue_ids = extract_ind_tissue_ids_used_for_eqtl_analysis(ordered_gtex_tissues, gtex_per_tissue_covariate_dir)

# Create mapping from gtex_sample_id to gtex_id:tissue name
gtex_sample_id_to_ind_tissue_id, ind_tissue_ids = create_mapping_from_gtex_sample_id_to_ind_tissue_id(gtex_sample_attributes_file, alt_gtex_tissue_to_gtex_tissue, valid_ind_tissue_ids)

# Double check to make sure we got eqtl samples
eqtl_sample_checking(ordered_gtex_tissues, gtex_per_tissue_covariate_dir, ind_tissue_ids)

###########################
# Processed gene expression
###########################

# First filter to PC genes and lowly expression. transoform to log scale
filtered_log_tpm_expression_file = processed_eqtl_data_dir + '_filtered_log_tpm_expression.txt'
generate_filtered_log_tpm_expression(filtered_log_tpm_expression_file,gtex_tpm_expression, gtex_sample_id_to_ind_tissue_id, gtex_id_to_ancestry_labels, pc_genes )

# Quantile normalize samples
quantile_normalized_expression_file = processed_eqtl_data_dir + '_quantile_normalized_expression.txt'
quantile_normalize_expression(filtered_log_tpm_expression_file, quantile_normalized_expression_file)

# Inverse-normal transform genes
inverse_normal_transformed_expression_file = processed_eqtl_data_dir + '_inverse_normal_transformed_expression.txt'
inverse_normal_transform_expression(quantile_normalized_expression_file, inverse_normal_transformed_expression_file)

# Compute cross-tissue expression PCs
expression_pc_output_file = processed_eqtl_data_dir + '_expression_pcs.txt'
num_expression_pcs = 55
compute_expression_pcs(inverse_normal_transformed_expression_file, expression_pc_output_file, num_expression_pcs)

# Compute 2D UMAP from expression PCs
umap_output_file = processed_eqtl_data_dir + '_expression_pc_umap_2d.txt'
compute_umap_from_expression_pcs(expression_pc_output_file, umap_output_file)

# Plot 2D UMAP colored by tissue identity
umap_plot_file = processed_eqtl_data_dir + '_expression_pc_umap_2d_by_tissue.png'
generate_umap_tissue_identity_scatterplot(umap_output_file, umap_plot_file)

# Generate indices that will map genotype indices to expression sample
genotype_indices_to_expression_file = processed_eqtl_data_dir + '_genotype_mapping_to_expression_samples.txt'
ordered_expression_samples = np.loadtxt(expression_pc_output_file, dtype=str,delimiter='\t')[1:,0]
create_mapping_from_genotype_indices_to_expression_sample(ordered_expression_samples, processed_genotype_data_dir, genotype_indices_to_expression_file)

# Limit expression PCs to first five
nn = 5
top_five_expression_pc_output_file = processed_eqtl_data_dir + '_expression_pcs_top_' + str(nn) + '.txt'
limit_expression_pcs_to_top_n(expression_pc_output_file, top_five_expression_pc_output_file, nn)

nn = 50
top_five_expression_pc_output_file = processed_eqtl_data_dir + '_expression_pcs_top_' + str(nn) + '.txt'
limit_expression_pcs_to_top_n(expression_pc_output_file, top_five_expression_pc_output_file, nn)




