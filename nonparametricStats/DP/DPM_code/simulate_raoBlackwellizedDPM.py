from __future__ import division
import numpy as np
import scipy as sc
from scipy.stats import norm, beta
from matplotlib import pyplot as plt

# def simulate_raoBlackwellizedDPM(iter_num, M):
# Generate mixture sample
N = 1000
mu = [0.0, 10.0, 3.0]
M = 0.4
components = np.random.choice(range(3), size = N, replace = True, p = [0.3, 0.5, 0.2])
samples = [norm.rvs(size = 1, loc = mu[components[i]], scale = 1)[0] for i in range(N)]

## Sample G from DP(M, G0)
v = beta.rvs(a = 1.0, b = M, size = N)
prob_vector = np.append(np.array(v[0]), v[1:] * np.cumprod(1.0 - v[:-1]))
thetas = norm.rvs(size = N, loc = 1.0, scale = 1.0)

### Initialize thetas
thetas = np.random.choice(thetas, size = N, replace = True, p = prob_vector)
theta_star = np.unique(thetas)
count_category = [sum(x == theta_star[i] for x in thetas) for i in xrange(len(theta_star))]
