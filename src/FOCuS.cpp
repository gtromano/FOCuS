#include "FOCuS.h"


// update a quadratic with a new observation. To be used with
// std::move
auto update_quad(Quadratic q, const double& new_point, const double& offset = 0.0) {
  q.a -= 0.5;
  q.b += new_point;
  q.c += offset;
  return q;
}


// takes the information from the past and updates it
// to be used with std::move to avoid copy
Info FOCuS_step(Info info, const double& new_point) {
  
  // update the quad for the null
  // add std::move after testing
  info.Q0 = update_quad(std::move(info.Q0), new_point);
    
  // find the max of Q0
  double Q0_max = -INFINITY;
  for(const auto& i:info.Q0.ints) {
    auto m = std::get<0>(get_minimum(info.Q0, i));
    if (m > Q0_max)
      Q0_max = m;
  }

  // update the quad for the alternative and lowering by the max of Q0
  std::for_each(info.Q1.begin(), info.Q1.end(), [&new_point, &Q0_max](auto &q){
    q = update_quad(std::move(q), new_point, -Q0_max);
  });


  // lowering Q0 by the max of Q0
  info.Q0.c -= Q0_max;
  
  // get the new line
  Quadratic line; // remember that this is initialized at 0 0 0, for (-inf, inf)
  
  // trimming with the new line // add std::move
  info.Q1 = get_max_of_cost(std::move(info.Q1), std::move(line));
  

  // getting the maximums for each piecewise quadratic
  double global_max = -INFINITY;
  std::for_each(info.Q1.begin(), info.Q1.end(), [&global_max](const auto& q){
    double max = -INFINITY;
    // iterating in the intervals of the piecewise quadratics
    for(const auto& i:q.ints) {
      auto m = std::get<0>(get_minimum(q, i));
      if (m > max)
        max = m;
    }
    if (max > global_max)
      global_max = max;
  });
  
  info.global_max = std::move(global_max);
  
  // and we're done!
  return info;
}

// here we probably need the pruning function, to implement the pruning




///////////////////////////////////////////////////////////////////

// takes the information from the past and updates it
// to be used with std::move to avoid copy
void FOCuS_step_V1(Info& info, const double& new_point) {
  
  
  // update the quad for the null
  // add std::move after testing
  info.Q0 = update_quad(std::move(info.Q0), new_point);
  
  
  // find the max of Q0
  auto Q0_max = std::get<0>(get_minimum(info.Q0, info.Q0.ints.front()));

  
  // update the quad for the alternative
  std::for_each(info.Q1.begin(), info.Q1.end(), [&new_point, &Q0_max](auto &q){
    q = update_quad(std::move(q), new_point, -Q0_max);
  });
  
  
  
  // lowering Q0 by the max of Q0
  info.Q0.c -= Q0_max;
  
  // get the new line
  Quadratic line; // remember that this is initialized at 0 0 0, for (-inf, inf)
  
  // trimming with the new line // add std::move
  get_max_of_cost_V1(info.Q1, line);
  
  
  // getting the maximums for each piecewise quadratic
  // std::vector<double> maxs(info.Q1.size());
  // std::transform(info.Q1.begin(), info.Q1.end(), maxs.begin(), [](const auto& q){
  //   double max = -INFINITY;
  //   // iterating in the intervals of the piecewise quadratics
  //   for(const auto& i:q.ints) {
  //     auto m = std::get<0>(get_minimum(q, i));
  //     if (m > max)
  //       max = m;
  //   }
  //   return max;
  // });
  double global_max = -INFINITY;
  std::for_each(info.Q1.begin(), info.Q1.end(), [&global_max](const auto& q){
    double max = -INFINITY;
    // iterating in the intervals of the piecewise quadratics
    for(const auto& i:q.ints) {
      auto m = std::get<0>(get_minimum(q, i));
      if (m > max)
        max = m;
    }
    if (max > global_max)
      global_max = max;
  });
  
  
  info.global_max = std::move(global_max);
  
  // and we're done!
}
