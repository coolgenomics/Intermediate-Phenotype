import numpy as np
import sys
import pickle
import matopr
from simulation_individual import simulate_individual


# Calculate the RAF empirically from the provided genotype sample data set
def get_freqs(genos):
	return np.array([np.mean(genos[:, i]) / 2 for i in range(genos.shape[1])], dtype=np.float64)


# Get the mean and std of the trait associated with the weights empirically
def get_pred_norm_factor(gen2expr_wgtmat, freqs, alpha, size=5000):
	inds = simulate_individual(size, freqs)
	# using blockwise_dot to prevent memory overflow
	expression = matopr.blockwise_dot(gen2expr_wgtmat, inds.T, max_elements=int(2**20)).T
	prediction = expression.dot(alpha)
	return np.mean(prediction), np.std(prediction)


# Get the relevant parameters for simulation given the weights and the set up
# This partitions the alpha weight accoding to different sets of genotype.
def get_params(freqs, gen2expr_wgtmat, alpha, perm, gene_props):
	# expression weights for each disease (sub-phenotype)
	# the elements in each row of subalpha that does not correspond to the
	# disease of that row will be set to 0
	subalpha = np.zeros([len(gene_props), len(alpha)])
	# mean and std score for each disease, with the last element as the mean
	# and std score of the complex trait
	means = np.zeros(len(gene_props) + 1)
	stds = np.zeros(len(gene_props) + 1)
	start = 0
	for i in range(len(gene_props)):
		end = int(start + gene_props[i] * len(alpha))
		subalpha[i, perm[start:end]] = alpha[perm[start:end]]
		means[i], stds[i] = get_pred_norm_factor(gen2expr_wgtmat, freqs, subalpha[i])
		start = end
	# statistics for the complex trait
	mean, std = get_pred_norm_factor(gen2expr_wgtmat, freqs, alpha)
	means[-1] = mean
	stds[-1] = std
	return subalpha, means, stds


if __name__ == "__main__":
	weight_path, genos_path, result_path  = sys.argv[1:]
	gene_props = list(map(float, sys.argv[4:]))
	gen2expr_wgtmat, expr2row, snp2col, snps_ref, alpha = pickle.load(open(weight_path, "rb"))
	sample, snps, genos = pickle.load(open(genos_path, "rb"))
	freqs = get_freqs(genos)
	perm = np.random.permutation(len(alpha))
	subalpha, means, stds = get_params(freqs, gen2expr_wgtmat, alpha, perm, gene_props)
	pickle.dump((subalpha, means, stds, freqs, perm), open(result_path, "wb"))
