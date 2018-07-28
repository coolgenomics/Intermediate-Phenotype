import numpy as np

# Add noise to the data to conform to variance explained
def add_heritability(std, h, data):
	dev = np.sqrt(((1 - h) / h) * std**2)
	noise = np.random.normal(0, dev, len(data))
	return data + noise
