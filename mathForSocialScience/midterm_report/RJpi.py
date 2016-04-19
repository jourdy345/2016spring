import math
from __future__ import division
def RJpi(tau):
	n = 0.0
	k1 = 26390.0
	k2 = 1103.0
	k3 = 396.0
	eps = tau + 1.0
	rjpi = 0
	epsList = []
	while(eps > tau):
		eps = math.factorial(4.0*n)/((math.factorial(n))**4.0) * (k1 * n + k2)/(k3**(4.0*n))
		epsList.append(eps)
		rjpi += eps
		n += 1
	rjpi = 9801.0/(math.sqrt(8)*rjpi)
	return rjpi, epsList



