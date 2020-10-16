#include "FOCuS.h"


// update a quadratic with a new observation. To be used with
// std::move
auto update_quad(Quadratic q, const double& new_point) {
  q.a -= 0.5;
  q.b += new_point;
  return q;
}


// takes the information from the past and updates it
// to be used with std::move to avoid copy
Info FOCuS_step(Info info, const double& new_point) {
  
  // std::cout << "Q1 at beginning "<< std::endl;
  // for (auto& q:info.Q1)
  //   print(q);
  // std::cout << std::endl;
  
  // std::cout << << std::endl;
  
  // std::cout << "q0 at beginning"<< std::endl;
  // print(info.Q0);
  // std::cout << std::endl;
  
  
  // update the quad for the null
  // add std::move after testing
  info.Q0 = update_quad(std::move(info.Q0), new_point);
  
  // std::cout << "q0 after update"<< std::endl;
  // print(info.Q0);
  // std::cout << std::endl;
  
  
  // std::cout << "q0 after update "<< std::endl;
  // print(info.Q0);
  
  // update the quad for the alternative
  std::for_each(info.Q1.begin(), info.Q1.end(), [&new_point](auto &q){
    q = update_quad(std::move(q), new_point);
  });

  // std::cout << "Q1 after update"<< std::endl;
  // for (auto& q:info.Q1)
  //   print(q);
  // std::cout << std::endl;
  
    
  // find the max of Q0
  double Q0_max = -INFINITY;
  for(const auto& i:info.Q0.ints) {
    auto m = std::get<0>(get_minimum(info.Q0, i));
    if (m > Q0_max)
      Q0_max = m;
  }
  
  // std::cout << "Q0_max is: "<< Q0_max << std::endl;
  
  
  // lowering Q1 by the max Q0
  std::for_each(info.Q1.begin(), info.Q1.end(), [&Q0_max](auto &q){
    q.c -= Q0_max;
  });
  
  // std::cout << "Q1 after lowering"<< std::endl;
  // for (auto& q:info.Q1)
  //   print(q);
  // std::cout << std::endl;
  
  
  // lowering Q0 by the max of Q0
  info.Q0.c -= Q0_max;
  
  
  // std::cout << "q0 after lowering"<< std::endl;
  // print(info.Q0);
  // std::cout << std::endl;
  
  
  // get the new line
  Quadratic line; // remember that this is initialized at 0 0 0, for (-inf, inf)
  
  // trimming with the new line // add std::move
  info.Q1 = get_max_of_cost(std::move(info.Q1), std::move(line));
  
  // std::cout << "Q1 after max search"<< std::endl;
  // for (auto& q:info.Q1)
  //   print(q);
  // std::cout << std::endl;
  
  
  // getting the maximums for each piecewise quadratic
  std::vector<double> maxs(info.Q1.size());
  std::transform(info.Q1.begin(), info.Q1.end(), maxs.begin(), [](const auto& q){
    double max = -INFINITY;
    // iterating in the intervals of the piecewise quadratics
    for(const auto& i:q.ints) {
      auto m = std::get<0>(get_minimum(q, i));
      if (m > max)
        max = m;
    }
    return max;
  });
  
  info.global_max = *std::max_element(maxs.begin(), maxs.end());
  
  // and we're done!
  return info;
}

// here we probably need the pruning function, to implement the pruning