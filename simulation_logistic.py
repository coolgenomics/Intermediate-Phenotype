import numpy as np
import sys
import pickle
import matopr
from add_heritability import add_heritability
from simulation_individual import simulate_individual


# simulate cases and controls using liability threshold and logistic model
# according to the specified parameters
def get_cases_and_control(case_sizes, control_size, freqs, gen2expr_wgtmat, subalpha, means, stds, heritability, thresh=1.8, batch=3000, max_elements=int(2**20)):
	cases = []
	for i in range(len(case_sizes)):
		cases.append(np.empty([case_sizes[i], gen2expr_wgtmat.shape[0]]))
	control = np.empty([control_size, gen2expr_wgtmat.shape[0]])

	all_count = 0
	cur_count = 0
	cur_case_size = [0] * len(case_sizes)
	for case_size in case_sizes:
		all_count += case_size

	while cur_count < all_count:
		individuals = simulate_individual(batch, freqs)
		expression = matopr.blockwise_dot(gen2expr_wgtmat, individuals.T, max_elements=max_elements).T
		risk = expression.dot(subalpha.T)
		for i in range(risk.shape[1]):
			# add noise to data to adjust for heritibility
			risk[:, i] = add_heritability(stds[i], heritability, risk[:, i])
		for i in range(risk.shape[0]):
			for j in range(risk.shape[1]):
				z_score = (risk[i, j] - means[j]) * heritability / stds[j]
				prob = 1 / (1 + np.exp(thresh - z_score))
				np.random.seed()
				if np.random.rand() < prob and cur_case_size[j] < case_sizes[j]:
					cases[j][cur_case_size[j], :] = expression[i]
					cur_count += 1
					cur_case_size[j] += 1
		print("[CASE] generated", cur_count, "out of", all_count)

	cur_count = 0
	while cur_count < control_size:
		additional_size = min(batch, control_size - cur_count)
		individuals = simulate_individual(additional_size, freqs)
		expression = matopr.blockwise_dot(gen2expr_wgtmat, individuals.T, max_elements=max_elements).T
		control[cur_count:cur_count + additional_size, :] = expression
		cur_count += batch
		print("[CONTROL] generated", cur_count, "out of", control_size)
	return cases, control

if __name__ == "__main__":
	weight_path, param_path, result_path = sys.argv[1:4]
	heritability = float(sys.argv[4])
	control_size = int(sys.argv[5])
	case_sizes = list(map(int, sys.argv[6:]))
	gen2expr_wgtmat, expr2row, snp2col, snps_ref, alpha = pickle.load(open(weight_path, 'rb'))
	subalpha, means, stds, freqs, perm = pickle.load(open(param_path, 'rb'))
	cases, control = get_cases_and_control(case_sizes, control_size, freqs, gen2expr_wgtmat, subalpha, means, stds, heritability)
	pickle.dump((perm, means, stds, cases, control), open(result_path, "wb"))
