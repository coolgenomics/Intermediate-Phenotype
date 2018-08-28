import pickle
import sys
import numpy as np
from buhmbox import BuhmBox_TWAS


# calculate a sequence of heterogeneity values for expression with 
# the associated sample sizes
def heterogeneity_seq(expr, control, clist, resolution=2000):
	np.random.shuffle(expr)
	BB = BuhmBox_TWAS()
	x, res = [], []
	for i in range(5,len(expr), resolution):
		BB.bb(expr[range(i+1)], control)
		sbb = BB.get_values(clist)
		res.append(sbb)
		x.append(i+1)
	return (x, res)

# partition the expression list to the different subphenotypes
# according to the proportions
def partition_by_expression(perm, props):
	num_expr = case1.shape[1]
	start1 = 0
	end1 = int(props[0] * num_expr)
	start2 = end1
	end2 = start2+int(props[1] * num_expr)
	expr_1 = perm[start1:end1]
	expr_2 = perm[start2:end2]
	expr_neither = perm[end2:]
	return expr_1, expr_2, expr_neither

# calculate the heterogeneity sequence for the specified mixing
# proportion pi. 
def get_pi_seq(perm, cases, control, props, pi, use_expr_1=True):
	expr_1, expr_2, expr_neither = partition_by_expression(perm, props):
	case1, case2 = cases
	cases_comb = np.vstack((case1[0:len(case1)*pi], case2[0:len(case2)*(1-pi)]))
	if use_expr_1:
		return heterogeneity_seq(cases_comb, control, expr_1)
	else:
		return heterogeneity_seq(cases_comb, control, expr_2)


if __name__ == "__main__":
	pi = float(sys.argv[1])
	data_path, result_path = sys.argv[2:4]
	perm, means, stds, cases, control, gene_props = pickle.load(open(data_path, "rb"))
	seq = get_pi_seq(perm, cases, control, gene_props, pi)
	pickle.dump(seq, open(result_path, "wb"))


	