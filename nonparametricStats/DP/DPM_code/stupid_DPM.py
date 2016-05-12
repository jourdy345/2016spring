from __future__ import division
import numpy as np
import scipy as sc
from scipy.stats import norm, beta
from matplotlib import pyplot as plt

def simulate_stupidDPM(iter_num, M):
	# Generate mixture sample
	N = 1000
	mu = [0.0, 10.0, 3.0]

	components = np.random.choice(range(3), size = N, replace = True, p = [0.3, 0.5, 0.2])
	samples = [norm.rvs(size = 1, loc = mu[components[i]], scale = 1)[0] for i in range(N)]

	## Sample G from DP(M, G0)
	v = beta.rvs(a = 1.0, b = M, size = N)
	prob_vector = np.append(np.array(v[0]), v[1:] * np.cumprod(1.0 - v[:-1]))
	thetas = norm.rvs(size = N, loc = 1.0, scale = 1.0)

	### Initialize thetas
	thetas = np.random.choice(thetas, size = N, replace = True, p = prob_vector)

	### Start MCMC chain
	for i in xrange(iter_num):
		for j in xrange(N):
			theta_temp = np.append(thetas[:j], thetas[j+1:])
			p = np.append(norm.pdf(samples[j], loc = theta_temp, scale = 1.0), M * norm.pdf(samples[j], loc = 1.0, scale = np.sqrt(2.0)))
			p = p / sum(p)
			temp = np.random.choice(np.append(theta_temp, N), size = 1, replace = True, p = p)
			if (temp == N):
				thetas[j] = norm.rvs(size = 1, loc = 0.5 * (samples[j] + 1), scale = np.sqrt(0.5))
			else:
				thetas[j] = temp
		print(thetas)
	return {"thetas": thetas, "y": samples}

res = simulate_stupidDPM(200, 0.4)
print(np.unique(res["thetas"]))

# plt.scatter(range(N), res["thetas"])