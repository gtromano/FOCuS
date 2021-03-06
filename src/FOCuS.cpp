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
Info FOCuS_step(Info info, const double& new_point, const std::list<double>& grid, const double& K = INFINITY) {
  
  // update the quad for the null
  // add std::move after testing
  //if (std::isinf(K)) {
  if (false) {
    update_quad(info.Q0, new_point);
  } else {
    
    auto m0 = - (info.Q0.b + new_point) / (2.0 * (info.Q0.a - 0.5)); // argmax of updated quad
    if (-.5 * (m0 * m0 - 2 * m0 * new_point + new_point * new_point) > - 2 * K) // if max of new quad is greater than -2 * K
      update_quad(info.Q0, new_point);
    else
      update_quad(info.Q0, - 2 * K);
  }
    
  // find the max of Q0
  info.Q0.max = std::get<0>(get_minimum(info.Q0, info.Q0.ints.front()));
  
  //if (std::isinf(K)) {
  if (false) {
    // update the quad for the alternative and lowering by the max of Q0
    for (auto& q:info.Q1)
      update_quad(q, new_point, -info.Q0.max);
  } else {
    update_cost_biweight(info.Q1, new_point, K, std::get<1>(get_minimum(info.Q0, info.Q0.ints.front())));
  }
  
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

Info FOCuS_step_sim(Info info, const double& new_point, const std::list<double>& grid, const double& K = INFINITY) {
  
  //std::cout << "AAAA";
  
  if (std::isinf(K)) {
    for (auto& q:info.Q1)
      update_quad(q, new_point);
  } else {
    update_cost_biweight(info.Q1, new_point, K, 0.0);
  }
  
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



// finds the minimum of a quadratic on an interval
double get_only_the_minimum(const Quadratic& q, const Interval& i){
  
  if (q.a == 0) {
    return q.c;
  }
  
  auto at = - q.b / (2.0 * q.a);
  if (at <= i.l) {
    at = i.l;
  } else if (at >= i.u) {
    at = i.u;
  }
  
  return  q.a * (at * at) + q.b * at + q.c;
}


/*
mInfo FOCuS_step_melk(mInfo info, const double& new_point, const std::list<double>& lgrid, const std::list<double>& rgrid, const double& K = INFINITY) {
  
  // info should have Qright and Qleft
  // this is the update
  if (std::isinf(K)) {
    for (auto& q:info.Qleft)
      update_quad(q, new_point);
    for (auto& q:info.Qright)
      update_quad(q, new_point);
  } else {
    update_cost_biweight(info.Qleft, new_point, K, 0.0);
    update_cost_biweight(info.Qright, new_point, K, 0.0);
  }
  
    if (new_point > 0) {
      // get the new line
      Quadratic rline; // remember that this is initialized at 0 0 0, for (-inf, inf)
      rline.ints.front().l = 0;
      
      // trimming with the new line // add std::move
      get_max_of_cost_melk_right(info.Qright, std::move(rline));
      
      // grid approximation
      if (!std::isnan(rgrid.front()))
        approximation_grid(info.Qright, rgrid);
      
      info.global_max = -INFINITY;
      // getting the maximums for right
      std::for_each(info.Qright.begin(), info.Qright.end(), [&info](auto& q){
        double m = -INFINITY;
        for(const auto& i:q.ints) {
          m = std::get<0>(get_minimum(q, i));
          if (m > q.max)
            q.max = m;
        }
        if (q.max > info.global_max)
          info.global_max = q.max;
      });
      
    } else {
      // get the new line
      Quadratic lline; // remember that this is initialized at 0 0 0, for (-inf, inf)
      lline.ints.front().u = 0;
      
      // trimming with the new line // add std::move
      get_max_of_cost_melk_left(info.Qleft, std::move(lline));
      
      // grid approximation
      if (!std::isnan(lgrid.front()))
        approximation_grid(info.Qleft, lgrid);
      
      // getting the maximums for left
      info.global_max = -INFINITY;
      std::for_each(info.Qleft.begin(), info.Qleft.end(), [&info](auto& q){
        double m = -INFINITY;
        for(const auto& i:q.ints) {
          m = std::get<0>(get_minimum(q, i));
          if (m > q.max)
            q.max = m;
        }
        if (q.max > info.global_max)
          info.global_max = q.max;
      });
    }
    
  // and we're done!
  return info;
}
*/


mInfo FOCuS_step_melk(mInfo info, const double& new_point, const std::list<double>& lgrid, const std::list<double>& rgrid, const double& K = INFINITY) {
  
  // info should have Qright and Qleft
  // this is the update
  if (std::isinf(K)) {
    for (auto& q:info.Qleft)
      update_quad(q, new_point);
    for (auto& q:info.Qright)
      update_quad(q, new_point);
  } else {
    update_cost_biweight(info.Qleft, new_point, K, 0.0);
    update_cost_biweight(info.Qright, new_point, K, 0.0);
  }
  
  if (new_point > 0) {
    // get the new line
    Quadratic rline; // remember that this is initialized at 0 0 0, for (-inf, inf)
    rline.ints.front().l = 0;
    
    // trimming with the new line // add std::move
    get_max_of_cost_melk_right(info.Qright, std::move(rline));
    
    // grid approximation
    if (!std::isnan(rgrid.front()))
      approximation_grid(info.Qright, rgrid);
    
    info.global_max = -INFINITY;
    // getting the maximums for right
    std::for_each(info.Qright.begin(), info.Qright.end(), [&info](auto& q){
      double m = -INFINITY;
      for(const auto& i:q.ints) {
        m = get_only_the_minimum(q, i);
        if (m > q.max)
          q.max = m;
      }
      if (q.max > info.global_max)
        info.global_max = q.max;
    });
    
  } else {
    // get the new line
    Quadratic lline; // remember that this is initialized at 0 0 0, for (-inf, inf)
    lline.ints.front().u = 0;
    
    // trimming with the new line // add std::move
    get_max_of_cost_melk_left(info.Qleft, std::move(lline));
    
    // grid approximation
    if (!std::isnan(lgrid.front()))
      approximation_grid(info.Qleft, lgrid);
    
    // getting the maximums for left
    info.global_max = -INFINITY;
    std::for_each(info.Qleft.begin(), info.Qleft.end(), [&info](auto& q){
      double m = -INFINITY;
      for(const auto& i:q.ints) {
        m = get_only_the_minimum(q, i);
        if (m > q.max)
          q.max = m;
      }
      if (q.max > info.global_max)
        info.global_max = q.max;
    });
  }
  
  // and we're done!
  return info;
}

