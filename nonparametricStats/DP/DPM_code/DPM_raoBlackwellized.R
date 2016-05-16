library(data.table)
generate_sample = function(N, M, mu0, sigma0) {
  ### Generate simulation data
  components = sample(1:3, prob = c(0.3, 0.5, 0.2), size = N, replace = TRUE)
  mu = c(0, 10, 3)
  samples = rnorm(n = N, mean = mu[components], sd = 1)
  original = data.table(y = samples, category = components)

  ### Sample G from DP(M, G0)
  v = rbeta(N, 1, M)
  prob_vector = c(v[1], v[-1] * cumprod(1-v[-N]))
  thetas = rnorm(N, mu0, sigma0)
  
  ### Initialize thetas
  thetas = sample(thetas, N, replace = TRUE, prob = prob_vector)
  theta_star = unique(thetas)
  
  category_variable = sample(1:length(theta_star), size = N, replace = TRUE)
  category_unique = unique(category_variable)
  DT = unname(as.matrix(cbind(samples, category_variable)))
  cat(DT, "\n")

  category_unique = sort(unique(DT[,2] ))
  print(category_unique)  
  for (k in 1:length(category_unique)) {
    for (loop in 1:N) {
      print(DT[loop,2])
      if (DT[loop,2] == category_unique[k]) {
        DT[loop,2] = k
      }
    }
    cat('\n')
  }
  DT
}


DPM_raoBlackwellized = function(iter_num, M, N, mu0, sigma0, sigma) {

  ### Generate simulation data
  components = sample(1:4, prob = c(0.3, 0.4, 0.2, 0.1), size = N, replace = TRUE)
  mu = c(0, 10, 3, 5)
  samples = rnorm(n = N, mean = mu[components], sd = 1)
  original = data.table(y = samples, category = components)

  ### Sample G from DP(M, G0)
  v = rbeta(N, 1, M+1)
  prob_vector = c(v[1], v[-1] * cumprod(1-v[-N]))
  thetas = rnorm(N, mu0, sigma0)
  
  ### Initialize thetas
  thetas = sample(thetas, N, replace = TRUE, prob = prob_vector)
  theta_star = unique(thetas)
  
  category_variable = sample(1:length(theta_star), size = N, replace = TRUE)
  category_unique = sort(unique(category_variable))
  DT = data.table(y = samples, category = category_variable)
  
  ### 
  for (k in 1:length(category_unique)) {
  	for (loop in 1:N) {
  		if (DT[loop,]$category == category_unique[k]) {
  			DT[loop,]$category = k
  		}
  	}	
  }

  K = length(category_unique)
  ### Start MCMC chain
  for (j in 1:iter_num) {
  	for (i in 1:N) {
			p = c()
  		temp_y = DT[i,]
  		# temp_DT = DT[category != temp_y$category,]
  		temp_DT = DT[-i,]
  		temp_category = unique(temp_DT$category)

  		for (temp_cat in temp_category) {
	  		temp_DT2 = temp_DT[category == temp_cat,]
	  		n_c = nrow(temp_DT2)
	  		mu_theta = (sigma0 * sum(temp_DT2$y) + mu0 * sigma) / (n_c * sigma0)
	  		sigma_theta = (sigma * sigma0) / (n_c * sigma0 + sigma)
  			p = c(p, n_c / (N-1+M) * dnorm(temp_y$y, mu_theta, sqrt(sigma_theta + sigma)))
  		}
  		p = c(p, M / (N-1+M) * dnorm(temp_y$y, mu0, sqrt(sigma0 + sigma)))
  		p = p / sum(p)
  		cat('p:', p, '\n')
  		cat('length: ', length(p), '\n')
  		cat('K+1: ', K+1, '\n')
  		x_vec = temp_category
  		DT[i,]$category = sample(c(x_vec,K+1), size = 1, prob = p)

  		category_unique = sort(unique(DT$category))
  		for (k in 1:length(category_unique)) {
		  	for (loop in 1:N) {
		  		if (DT[loop,]$category == category_unique[k]) {
		  			DT[loop,]$category = k
		  		}
		  	}	
		  }
		  K = length(category_unique)
  	}
  }
  list(true = original, DPM_estimate = DT)
}
