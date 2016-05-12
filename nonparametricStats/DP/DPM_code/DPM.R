simulate_DPM = function(iter_num, M) {
  N = 1000

  components = sample(1:3, prob = c(0.3, 0.5, 0.2), size = N, replace = TRUE)
  mu = c(0, 10, 3)
  samples = rnorm(n = N, mean = mu[components], sd = 1)

  ### Sample G from DP(M, G0)
  v = rbeta(N, 1, M)
  prob_vector = c(v[1], v[-1] * cumprod(1-v[-N]))
  thetas = rnorm(N, 1, 1)
  ### Initialize thetas ###
  thetas = sample(thetas, N, replace = TRUE, prob = prob_vector)
  ### MCMC for DP mixture ###
  for (i in 1:iter_num) {
    for (n in 1:N) {
      p = c(dnorm(samples[n], thetas[-n], 1), M * dnorm(samples[n], 1, sqrt(2)))
      p = p / sum(p)
      temp = sample(c(thetas[-n], N), size = 1, replace = TRUE, prob = p)
      if (temp == N) {
        thetas[n] = rnorm(1, 0.5 * (samples[n]+1), sqrt(0.5))
      } else {
        thetas[n] = temp
      }
    }
  }
  list(thetas = thetas, y = samples)
}

