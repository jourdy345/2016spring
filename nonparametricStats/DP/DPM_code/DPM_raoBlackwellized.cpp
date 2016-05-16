// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
using namespace Rcpp; using namespace arma;

unsigned int my_sample(arma::uvec x, arma::vec p)
{
  std::random_device generator;
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  double u = distribution(generator);
  arma::vec temp_p = arma::cumsum(p);
  int n = p.size();
  int j;
  for (int i = 0; i < n; ++i) {
    if (u < temp_p(i)) {
      j = i;
      break;
    }
  }
  return x(j);
}


double my_dnorm(double x, double m, double s)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;
    return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}

arma::vec remove_ith_vec(arma::vec x, int i)
{
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
    for (int j = 0; j < i1; ++j) {
      temp3(j) = temp(j);
    }

    for (int k = 0; k < i2; ++k) {
      temp3(k + i1) = temp2(k);
    }
    return temp3;
  }
}

arma::uvec remove_ith_uvec(arma::uvec x, int i)
{
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
    for (int j = 0; j < i1; ++j) {
      temp3(j) = temp(j);
    }

    for (int k = 0; k < i2; ++k) {
      temp3(k + i1) = temp2(k);
    }
    return temp3;
  }
}

// [[Rcpp::export]]
arma::uvec DPM_raoBlackwellized(arma::mat DT, const int iter_num, const double M, const double mu0, const double sigma0, const double sigma)
{
  int N = DT.n_rows;
  arma::vec cat2 = DT.col(1);
  arma::uvec cat(DT.n_rows); 
  for (int v = 0; v < cat.size(); ++v) {
    cat(v) = static_cast<int>(cat2(v));
  }
  arma::vec y = DT.col(0);
  arma::uvec category_unique = arma::unique(cat);
  unsigned int K = category_unique.size();

  for (int j = 0; j < iter_num; ++j) {
    for (int i = 0; i < N; ++i) {
      // Generate probability vector p
      arma::vec tempy = remove_ith_vec(y, i+1);
      arma::uvec tempcat = remove_ith_uvec(cat, i+1);
      arma::uvec unique_category = arma::unique(tempcat);
      int p_size = unique_category.size();
      arma::vec p(p_size + 1);

      for (int l = 0; l < unique_category.size(); ++l) {
        arma::uvec lth_cluster = arma::find(tempcat == unique_category(l));
        arma::vec lth_y = tempy(lth_cluster);
        arma::uvec lth_cat = tempcat(lth_cluster);
    
        double n_c = static_cast<double>(lth_y.size());
        double mu_theta = (sigma0 * arma::accu(lth_y) + mu0 * sigma) / (n_c * sigma0);
        double sigma_theta = (sigma * sigma0) / (n_c * sigma0 + sigma);
        p(l) = n_c / (static_cast<double>(N) - 1.0 + M) * my_dnorm(y(i), mu_theta, std::sqrt(sigma_theta + sigma));
      }
      p(p_size) = M / (static_cast<double>(N) - 1.0 + M) * my_dnorm(y(i), mu0, std::sqrt(sigma0 + sigma));
      p /= arma::sum(p);

      // std::cout << "p: " << p << std::endl;
      std::cout << "# of unique values: " << p_size << std::endl;
      std::cout << "# of cluster: " << p.size() - 1 << std::endl;
      arma::uvec x_vec(p_size + 1);
      for (int iter; iter < p_size; ++iter) {
        x_vec(iter) = unique_category(iter);
      }
      x_vec(p_size) = K+1;
      cat(i) = my_sample(x_vec, p);

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
  return cat;
}


// arma::mat gibbs(const int npost, const int nthin)
// {
//   std::random_device generator;
//   arma::mat ps(npost, 2); ps.zeros();
//   double x = 0;
//   double y = 0;
//   for (int i = 0; i < npost; ++i) {
//     for (int j = 0; j < nthin; ++j) {
//       std::gamma_distribution<double> rgamma(3.0, 1.0 / (y * y + 4.0));
//       std::normal_distribution<double> rnorm(1.0 / (x + 1.0), 1.0 / (std::sqrt(2.0 * x + 2.0)));
//       x = rgamma(generator);
//       y = rnorm(generator);
//     }
//     ps(i,0) = x;
//     ps(i,1) = y;
//   }
//   return ps;
// }