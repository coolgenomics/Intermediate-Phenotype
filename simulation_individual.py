import numpy as np


# generate the SNP vectors of individuals according to the specified RAF
def simulate_individual(num, freqs):
	num = int(num)
	inds = np.empty([num, len(freqs)])
	for i, p in enumerate(freqs):
		sprobs = [(1-p)*(1-p), 2*p*(1-p), p*p]
		inds[:, i] = np.random.choice(3,size=num,p=sprobs)
	return inds
