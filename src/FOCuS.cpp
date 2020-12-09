#include "FOCuS.h"


// update a quadratic with a new observation. To be used with
// std::move
void update_quad(Quadratic& q, const double& new_point, const double& offset = 0.0) {
  q.a -= 0.5;
  q.b += new_point;
  q.c += offset;
  
}


// takes the information from the past and updates it
// to be used with std::move to avoid copy
Info FOCuS_step(Info info, const double& new_point, const std::list<double>& grid) {
  
  // update the quad for the null
  // add std::move after testing
  update_quad(info.Q0, new_point);
    
  // find the max of Q0
  info.Q0.max = std::get<0>(get_minimum(info.Q0, info.Q0.ints.front()));

  // update the quad for the alternative and lowering by the max of Q0
  for (auto& q:info.Q1)
    update_quad(q, new_point, -info.Q0.max);
  
  // lowering Q0 by the max of Q0
  info.Q0.c -= info.Q0.max;
  
  // get the new line
  Quadratic line; // remember that this is initialized at 0 0 0, for (-inf, inf)
  
  // trimming with the new line // add std::move
  info.Q1 = get_max_of_cost(std::move(info.Q1), std::move(line));
  
  // grid approximation
  if (!std::isnan(grid.front()))
    approximation_grid(info.Q1, grid);

  // getting the maximums for each piecewise quadratic
  double global_max = -INFINITY;
  std::for_each(info.Q1.begin(), info.Q1.end(), [&global_max](auto& q){
    double m = -INFINITY;
    for(const auto& i:q.ints) {
      m = std::get<0>(get_minimum(q, i));
      if (m > q.max)
        q.max = m;
    }
    if (q.max > global_max)
      global_max = q.max;
  });
  
  info.global_max = std::move(global_max);
  
  // and we're done!
  return info;
}


// This is needed for the simulations, in the sense that it assumes that the 
// data are centered on zero under the null
// takes the information from the past and updates it
// to be used with std::move to avoid copy
Info FOCuS_step_sim(Info info, const double& new_point, const std::list<double>& grid) {
  
  for (auto& q:info.Q1)
    update_quad(q, new_point);

  // get the new line
  Quadratic line; // remember that this is initialized at 0 0 0, for (-inf, inf)
  
  // trimming with the new line // add std::move
  info.Q1 = get_max_of_cost(std::move(info.Q1), std::move(line));

  // grid approximation
  if (!std::isnan(grid.front()))
    approximation_grid(info.Q1, grid);
  
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

