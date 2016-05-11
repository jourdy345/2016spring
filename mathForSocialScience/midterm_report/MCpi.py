from scipy.stats import uniform
import math
from __future__ import division
def MCpi(n):
  return 4.0 * sum(uniform.rvs(size = n)**2.0 + uniform.rvs(size = n)**2.0 <= 1.0) / n