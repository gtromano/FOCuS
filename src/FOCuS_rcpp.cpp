#include <Rcpp.h>
using namespace Rcpp;
#include "FOCuS.h"

// This is the script for exporting the C++ functions to R.

// converting c++ list into R list
std::list<List> convert_output_to_R(const auto& c_obj) {
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
                          Rcpp::Named("ints") = ints);
    l.attr("class") = "Quadratic";
    output.push_back(l);
  }
  return output;
}


// [[Rcpp::export]]
List FOCuS_offline(NumericVector Y, double thres) {
  long t = 0;
  long cp = -1;
  
  Quadratic Q0, q1;
  Info info = {Q0, {q1}, 0};
  
  for (auto& y:Y) {
    t += 1;
    info = FOCuS_step(std::move(info), y);
    
    if (info.global_max >= thres) {
      cp = t;
      break;
    }
  }
  // std::cout << info.Q1.size() << std::endl;
  // std::cout << info.global_max << std::endl;

  auto last_Q1 = convert_output_to_R(info.Q1);
  
  return List::create(Rcpp::Named("cp") = cp,
                      Rcpp::Named("Q1") = last_Q1);
}




// [[Rcpp::export]]
List FOCuS_offline_V1(NumericVector Y, double thres) {
  long t = 0;
  long cp = -1;
  
  Quadratic Q0, q1;
  Info info = {Q0, {q1}, 0};
  
  for (auto& y:Y) {
    t += 1;
    FOCuS_step_V1(info, y);
    
    if (info.global_max >= thres) {
      cp = t;
      break;
    }
  }
  // std::cout << info.Q1.size() << std::endl;
  // std::cout << info.global_max << std::endl;
  
  auto last_Q1 = convert_output_to_R(info.Q1);
  
  return List::create(Rcpp::Named("cp") = cp,
                      Rcpp::Named("Q1") = last_Q1);
}