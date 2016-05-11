import math
from __future__ import division
def CBpi(tau):
  k1 = 545140134.0
  k2 = 13591409.0
  k3 = 640320.0
  k4 = 100100025.0
  k5 = 327843840.0
  k6 = 53360.0
  cbpi = 0.0
  n = 0.0
  eps = tau + 1
  while eps > tau:
    eps = ((-1.0)**n) * math.factorial(6.0*n) * ( k2 + n * k1) / ((math.factorial(n)**3.0) * math.factorial(3.0 * n) * ((8.0 * k4 * k5)**n))
    cbpi = cbpi + eps
    n += 1
  return k6 * math.sqrt(k3) / cbpi