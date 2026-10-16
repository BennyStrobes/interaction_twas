import numpy as np
import sys


def print_samples_to_output(output_file, arr):
	t = open(output_file,'w')
	for ele in arr:
		t.write(str(ele) + '\n')
	t.close()
	return




#######################
# Command line args
#######################
genotype_to_expression_mapping_file = sys.argv[1]
n_folds = int(sys.argv[2])
output_file = sys.argv[3]
output_root = sys.argv[4]

# Initialize downsampling grid
downsampling_grid = np.asarray([1.0, .8, .6, .4, .2])

# total number of samples
tot_samples = len(np.loadtxt(genotype_to_expression_mapping_file))
sample_indices_vec = np.arange(tot_samples)

# Split sampled indices into n_folds
rng = np.random.default_rng(seed=123)  # optional seed
shuffled = rng.permutation(sample_indices_vec)
fold_sample_indices = np.array_split(shuffled, n_folds)

t = open(output_file,'w')
t.write('training_downsample_fraction\ttest_fold_iter\tn_training_samples\tn_test_samples\ttraining_sample_indices_file\ttest_sample_indices_file\n')


for test_fold_iter in range(n_folds):
	test_fold_sample_indices = fold_sample_indices[test_fold_iter]
	train_fold_sample_indices = np.concatenate([fold_sample_indices[ii] for ii in range(n_folds) if ii != test_fold_iter])


	for downsample_fraction in downsampling_grid:

		# Downsample training (no need to down-sample test)
		downsampled_train_fold_sample_indices = rng.choice(train_fold_sample_indices, size=int(np.ceil(downsample_fraction*len(train_fold_sample_indices))), replace=False)

		training_samples_file = output_root + 'downsample_' + str(downsample_fraction) + '_fold_' + str(test_fold_iter) + '_training_samples.txt'
		test_samples_file = output_root + 'downsample_' + str(downsample_fraction) + '_fold_' + str(test_fold_iter) + '_test_samples.txt'

		print_samples_to_output(training_samples_file, np.sort(downsampled_train_fold_sample_indices))
		print_samples_to_output(test_samples_file, np.sort(test_fold_sample_indices))


		t.write(str(downsample_fraction) + '\t' + str(test_fold_iter) + '\t' + str(len(downsampled_train_fold_sample_indices)) + '\t' + str(len(test_fold_sample_indices)) + '\t' + training_samples_file + '\t' + test_samples_file + '\n')




t.close()