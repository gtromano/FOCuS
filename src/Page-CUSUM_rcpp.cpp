#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <algorithm>

/* ------------------------------------------------------------

                OFFLINE VERSIONS OF PageCUSUM
 
-------------------------------------------------------------- */

std::vector<double> pageCUSUM_step (std::vector<double> Q, const std::vector<double>& grid, double new_point) {
  
  std::transform(Q.begin(), Q.end(), grid.begin(), Q.begin(), [&new_point] (auto q, auto mu) {
    auto out = std::max(0.0, q + mu * (new_point - mu * 0.5));
    return out;
    });
  return Q;
}

double get_max (std::vector<double> v) {
  auto result = std::max_element(v.begin(), v.end());
  auto out = std::distance(v.begin(), result);
  return v[out];
}

// [[Rcpp::export]]
List PageCUSUM_offline(NumericVector Y, const double thres, const double& mu0, std::vector<double>& grid) {

  long t = 0;
  long cp = -1;
  
  std::vector<double> Q(grid.size(), 0.0);
  std::list<double> max_at_time_t;

  for (auto& y:Y) {
      t += 1;
      Q = pageCUSUM_step(Q, grid, y - mu0);
      
      max_at_time_t.push_back(get_max(Q));
      
      if (max_at_time_t.back() >= thres) {
      cp = t;
      break;
    }
  }
  

  // std::cout << info.Q1.size() << std::endl;
  // std::cout << info.global_max << std::endl;
  
  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("maxs") = max_at_time_t);
}



/* ------------------------------------------------------------

                OFFLINE VERSIONS OF CUSUM

-------------------------------------------------------------- */


double CUSUM_step (double Q, double new_point, double mu0) {
  return Q + (new_point - mu0);
}


// [[Rcpp::export]]
List CUSUM_offline(NumericVector Y, const double thres, const double& mu0) {

  long t = 0;
  long cp = -1;

  auto Q = 0.0;
  std::list<double> max_at_time_t;

  for (auto& y:Y) {
      t += 1;
      Q = CUSUM_step(std::move(Q), y, mu0);

      max_at_time_t.push_back(std::abs(Q));

      if (max_at_time_t.back() >= thres) {
      cp = t;
      break;
    }
  }


  // std::cout << info.Q1.size() << std::endl;
  // std::cout << info.global_max << std::endl;

  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("maxs") = max_at_time_t);
}
