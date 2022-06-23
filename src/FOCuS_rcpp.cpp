#include <Rcpp.h>
using namespace Rcpp;
#include "FOCuS.h"

// This is the script for exporting the C++ functions to R.

// converting c++ list into R list
std::list<List> convert_output_to_R(const std::list<Quadratic>& c_obj) {
  std::list<List> output;
  
  //std::cout << "last Q1"<< std::endl;
  for (auto& q:c_obj) {
    // coversion of the data from c++ into R
    // print(q);
    
    std::list<List> ints(q.ints.size());
    std::transform(q.ints.begin(), q.ints.end(), ints.begin(),
                   [](const auto &i){
                     return List::create(Rcpp::Named("l") = i.l,
                                         Rcpp::Named("u") = i.u);
                   });
    
    auto l = List::create(Rcpp::Named("a") = q.a,
                          Rcpp::Named("b") = q.b,
                          Rcpp::Named("c") = q.c,
                          Rcpp::Named("ints") = ints,
                          Rcpp::Named("max") = q.max);
    l.attr("class") = "Quadratic";
    output.push_back(l);
  }
  return output;
}


/* ------------------------------------------------------------
 
                  Online version - Rcpp wrapper
 
 -------------------------------------------------------------- */

// [[Rcpp::export(.FoCUS)]]
List FOCuS (Rcpp::Function dataGen, const double thres, const double& mu0, std::list<double>& grid, const double& K) {
  if (!std::isnan(grid.front())) {
    grid.push_back(INFINITY);
    grid.push_front(-INFINITY);
  }
  
  long t = 0;
  long cp = -1;

  Quadratic Q0, q1;
  Info info = {Q0, {q1}, 0};  
  
  try {
    // if we don't know the pre-change mean then replace the f with a lambda  
    if (std::isnan(mu0)) {
      while(true) {
        t++;
        double y = Rcpp::as<double>(dataGen());
        info = FOCuS_step(std::move(info), y, grid, K);
        if (info.global_max >= thres) {
          cp = t;
          break;
        }
      }
    } else {
      while(true) {
        t++;
        double y = Rcpp::as<double>(dataGen());
        info = FOCuS_step_sim(std::move(info), y - mu0, grid, K);
        if (info.global_max >= thres) {
          cp = t;
          break;
        }
      }
    }
    
  }
  catch (std::bad_alloc &e) {
    Rcpp::stop("insufficient memory");
  }
  catch (...) {
    auto last_Q1 = convert_output_to_R(info.Q1);
    return List::create(Rcpp::Named("t") = cp,
                        Rcpp::Named("Q1") = last_Q1,
                        Rcpp::Named("warning_message") = "The procedure was interrupted or terminated unexpectedly. The output was successfully returned, however there is a possibility it can be possibly corrupted.");;
  }
  

  auto last_Q1 = convert_output_to_R(info.Q1);
  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("Q1") = last_Q1);
    
} 


/* ------------------------------------------------------------

                OFFLINE VERSIONS OF FOCuS
 
-------------------------------------------------------------- */

// [[Rcpp::export(.FoCUS_offline)]]
List FOCuS_offline(NumericVector Y, const double thres, const double& mu0, std::vector<double>& training_data, std::list<double>& grid, const double& K) {

  // checks if we have a grid, if so adds infinity on both ends to avoid deletions
  if (!std::isnan(grid.front())) {
    grid.push_back(INFINITY);
    grid.push_front(-INFINITY);
  }
    
  
  long t {0};
  long cp {-1};
  
  Quadratic Q0, q1;
  Info info = {Q0, {q1}, 0};
  std::list<double> max_at_time_t;
  
  try {
    // if we have previous training data for FOCuS pre-change-unknown, then updates the Q0 accordingly
    if (!std::isnan(training_data.front())) {
      for (auto& y_train:training_data) {
        info = FOCuS_training_step(std::move(info), y_train, grid, K);
      }
      info.Q1 = {Q0};
    }
    
    
    // pre-change mean not known 
    if (std::isnan(mu0)) {
      for (auto& y:Y) {
        t += 1;
        info = FOCuS_step(std::move(info), y, grid, K);
        //print(info.Q1.front());
        max_at_time_t.push_back(info.global_max);
        if (info.global_max >= thres) {
          cp = t;
          break;
        }
      }
    } else { // pre change mean known
      for (auto& y:Y) {
        t += 1;
        info = FOCuS_step_sim(std::move(info), y - mu0, grid, K);
        max_at_time_t.push_back(info.global_max);
        if (info.global_max >= thres) {
          cp = t;
          break;
        }
      }
    }
    
  }
  catch (std::bad_alloc &e) {
    Rcpp::stop("insufficient memory");
  }
  catch (...) {
    auto last_Q1 = convert_output_to_R(info.Q1);
    return List::create(Rcpp::Named("t") = cp,
                        Rcpp::Named("Q1") = last_Q1,
                        Rcpp::Named("maxs") = max_at_time_t,
                        Rcpp::Named("warning_message") = "The procedure was interrupted or terminated unexpectedly. The output was successfully returned, however there is a possibility it can be possibly corrupted.");;
  }
  


  auto last_Q1 = convert_output_to_R(info.Q1);
  
  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("Q1") = last_Q1,
                      Rcpp::Named("maxs") = max_at_time_t);
}

// [[Rcpp::export(.FOCuS_Melk)]]
List FOCuS_melk(NumericVector Y, const double thres, const double& mu0, std::list<double>& grid, const double& K) {
  
  if (!std::isnan(grid.front())) {
    grid.push_back(INFINITY);
    grid.push_front(-INFINITY);
  }
  
  
  auto rgrid = grid;
  auto lgrid = grid;
  rgrid.remove_if([](auto& x){
    return x < 0;
  });
  lgrid.remove_if([](auto& x){
    return x > 0;
  });
  
  
  long t = 0;
  long cp = -1;
  
  Quadratic q1;
  mInfo info = {{q1}, {q1}, 0};
  std::list<double> max_at_time_t;
  
  for (auto& y:Y) {
    t += 1;
    //std::cout << "\n time " << t << std::endl;
    info = FOCuS_step_melk(std::move(info), y - mu0, lgrid, rgrid, K);
    max_at_time_t.push_back(info.global_max);
    if (info.global_max >= thres) {
      cp = t;
      break;
    }
    //std::cout << "__________________________________" << std::endl;
  }
  
  
  // std::cout << info.Q1.size() << std::endl;
  // std::cout << info.global_max << std::endl;
  
  auto last_Qright = convert_output_to_R(info.Qright);
  auto last_Qleft = convert_output_to_R(info.Qleft);
  
  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("Qleft") = last_Qleft,
                      Rcpp::Named("Qright") = last_Qright,
                      Rcpp::Named("maxs") = max_at_time_t);
}


// // [[Rcpp::export]]
// List FOCuS_offline_sim(NumericVector Y, double thres, std::list<double>& grid, const double& K = INFINITY) {
// 
//   if (!std::isnan(grid.front())) {
//     grid.push_back(INFINITY);
//     grid.push_front(-INFINITY);
//   }
// 
//   long t = 0;
//   long cp = -1;
// 
//   Quadratic Q0, q1;
//   Info info = {Q0, {q1}, 0};
// 
//   for (auto& y:Y) {
//     t += 1;
//     info = FOCuS_step_sim(std::move(info), y, grid, K);
// 
//     if (info.global_max >= thres) {
//       cp = t;
//       break;
//     }
//   }
// 
//   auto last_Q1 = convert_output_to_R(info.Q1);
// 
//   return List::create(Rcpp::Named("t") = cp,
//                       Rcpp::Named("Q1") = last_Q1);
// }



/* ------------------------------------------------------------
 
 Offline multivariate FOCuS - Rcpp wrapper
 
 -------------------------------------------------------------- */




// [[Rcpp::export(.FoCUS_mult_offline)]]
List FOCuS_mult_offline(NumericMatrix Y, const double thres, const double& mu0, std::vector<double>& training_data, std::list<double>& grid, const double& K) {
  
  // checks if we have a grid, if so adds infinity on both ends to avoid deletions
  if (!std::isnan(grid.front())) {
    grid.push_back(INFINITY);
    grid.push_front(-INFINITY);
  }
  
  
  //long t {0};
  long cp {-1};
  
  
  auto nr = Y.nrow();
  auto nc = Y.ncol();

  std::list<double> max_at_time_t;
  
  //std::vector<std::list<double>> maxs_at_time_t(nr, std::list<double>());
  NumericMatrix maxs_at_time_t(nr,nc);
  
  Quadratic Q0, q1;
  
  Info temp = {Q0, {q1}, 0};
  std::vector<Info> m_info(nr, temp); // Initializing the storage for the independent FOCuS traces
  
  
  
  try {
    
    
    for (auto t = 0; t<nc; t++) { // t indexes time (the columns of the matrix) 
      
      std::multiset<double> f_stats;
      
      for (auto j = 0; j<nr; j++) {      // j indexes the different sequences (the rows of the matrix)
        
        
        
        if (std::isnan(mu0)) {
          
          
          std::cout << "y: " << Y(j ,t) << " j: " << j << " t: " << t << std::endl;
          
          m_info[j] = FOCuS_step(m_info[j], Y(j, t), grid, K);
          //print(info.Q1.front());
          
          f_stats.insert(m_info[j].global_max);
          
          
          maxs_at_time_t(j, t) = m_info[j].global_max;

        }
        
        // just some printing for debugging
        if (j == (nr - 1)) {
          std::cout << "time: " << t << " - The elements within the set are: ";
          for (auto pr_ = f_stats.rbegin(); pr_ != f_stats.rend(); pr_++)
            std::cout << *pr_ << " ";
          std::cout << std::endl << "_________________________________" << std::endl;
          
        }
        

      }
    }
    
    
    
    
  }
  catch (std::bad_alloc &e) {
    Rcpp::stop("insufficient memory");
  }
  catch (...) {
    // auto last_Q1 = convert_output_to_R(info.Q1);
    // return List::create(Rcpp::Named("t") = cp,
    //                     Rcpp::Named("Q1") = last_Q1,
    //                     Rcpp::Named("maxs") = max_at_time_t,
    //                     Rcpp::Named("warning_message") = "The procedure was interrupted or terminated unexpectedly. The output was successfully returned, however there is a possibility it can be possibly corrupted.");;
  }
  
  
  
  //auto last_Q1 = convert_output_to_R(info.Q1);
  
  return List::create(Rcpp::Named("t") = cp,
                      //Rcpp::Named("Q1") = last_Q1,
                      Rcpp::Named("maxs") = maxs_at_time_t
                      );
}
