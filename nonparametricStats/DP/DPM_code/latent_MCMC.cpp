// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
#include <random>
using namespace Rcpp; using namespace arma;

unsigned int my_sample(arma::uvec x, arma::vec p)
{
	std::random_device generator;
	std::uniform_real_distribution<double> distribution;
	double u = distribution(generator);
	arma::vec temp_p = arma::cumsum(p);
	int n = x.size();
	int j;
	for (int i = 0; i < n; ++i) {
		if (u < temp_p(i)) {
			j = i;
			break;
		}
	}
	return x(j);
}

double dnorm(double x, double m, double s)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}

arma::vec remove_ith_vec(arma::vec x, int i) {
  int N = x.size();
  if (i == 1) {
    return x.tail(N-1);
  } else if (i == N) {
    return x.head(N-1);
  } else {
    arma::vec temp = x.head(i-1);
    arma::vec temp2 = x.tail(N-i);
    int i1 = temp.size();
    int i2 = temp2.size();
    int i3 = i1 + i2;
    arma::vec temp3(i3);
    for (int j = 0; j < i1; j++) {
      temp3(j) = temp(j);
    }

    for (int k = 0; k < i2; k++) {
      temp3(k + i1) = temp2(k);
    }
    return temp3;
  }
}

arma::uvec remove_ith_uvec(arma::uvec x, int i) {
  // arma::uvec tempp(x.size());
  // for (int i = 0; i < x.size(); i++) {
  //   tempp(i) = static_cast<unsigned int>(x(i));
  // }
  int N = x.size();
  if (i == 1) {
    return x.tail(N-1);
  } else if (i == N) {
    return x.head(N-1);
  } else {
    arma::uvec temp = x.head(i-1);
    arma::uvec temp2 = x.tail(N-i);
    int i1 = temp.size();
    int i2 = temp2.size();
    int i3 = i1 + i2;
    arma::uvec temp3(i3);
    for (int j = 0; j < i1; j++) {
      temp3(j) = temp(j);
    }

    for (int k = 0; k < i2; k++) {
      temp3(k + i1) = temp2(k);
    }
    return temp3;
  }
}


// [[Rcpp::export]]
Rcpp::List latent_MCMC(arma::mat DT, const int iter_num, const double M, const double mu0, const double sigma0, const double sigma)
{
  std::random_device generator;
	int N = DT.n_rows;
	arma::vec cat2 = DT.col(1);
	arma::vec y = DT.col(0);
	arma::vec theta = DT.col(2);
  arma::uvec cat(DT.n_rows);
	for (int v = 0; v < cat.size(); v++) {
    cat(v) = static_cast<unsigned int>(cat2(v));
  }
  arma::uvec category_unique = arma::unique(cat);
  unsigned int K = category_unique.size();
	for (int i = 0; i < iter_num; ++i) {
		for (int j = 0; j < N; ++j) {
			// Generate probability vector p
      arma::vec tempy = remove_ith_vec(y, i+1);
      arma::uvec tempcat = remove_ith_uvec(cat, i+1);
      arma::vec temptheta = remove_ith_vec(theta, i+1);
      arma::uvec unique_category = arma::unique(tempcat);
      int p_size = unique_category.size();
      arma::colvec p(p_size + 1); p.zeros();

      // Update categories
      for (int l = 0; l < unique_category.size(); ++l) {
        arma::uvec lth_cluster = arma::find(tempcat == unique_category(l));
        arma::vec lth_y = tempy(lth_cluster);
        arma::uvec lth_cat = tempcat(lth_cluster);
        arma::vec lth_theta = temptheta(lth_cluster);
    		double  lth_theta_val = lth_theta(0); 

        double n_c = static_cast<double>(lth_y.size());
        p(l) = n_c * dnorm(y(i), lth_theta_val, std::sqrt(sigma));
      }
      p(p_size) = M * dnorm(y(i), mu0, sigma + sigma0);
      p /= arma::sum(p);
      std::cout << "p: " << p << std::endl;
      std::cout << "# of unique values: " << p_size << std::endl;
      std::cout << "# of cluster: " << p.size() - 1 << std::endl;
      arma::uvec x_vec(p_size + 1);
      for (int iter; iter < p_size; ++iter) {
        x_vec(iter) = unique_category(iter);
      }
      x_vec(p_size) = K+1;
      unsigned int new_val = my_sample(x_vec, p);
      cat(i) = new_val;

      // Update theta
      if (new_val == K + 1) {
        std::normal_distribution<double> rnorm((y(i) * sigma0 + mu0 * sigma) / (sigma0 + sigma), std::sqrt(sigma*sigma0 / (sigma + sigma0)));
        theta(i) = rnorm(generator);
      } else {
        arma::uvec new_cluster = arma::find(cat == new_val);
        arma::vec new_theta = theta(new_cluster);
        theta(i) = new_theta(0);
      }

      arma::uvec category_unique = arma::unique(cat);
      for (int k = 0; k < category_unique.size(); ++k) {
        for (int loop = 0; loop < N; ++loop) {
          if (cat(loop) == category_unique(k)) {
            cat(loop) = k+1;
          }
        }
      }
      K = static_cast<unsigned int>(category_unique.size());
		}
	}
  return Rcpp::List::create(Rcpp::Named("y") = y,
                            Rcpp::Named("categories") = cat,
                            Rcpp::Named("theta") = theta);
}