#ifndef ___QUAD_H___
#define ___QUAD_H___

#include <iostream>
#include <list>
#include <tuple>
#include <numeric>      // std::iota
#include <cmath>
#include <algorithm>

// interval structure
// set of real values l leq x leq u 
typedef struct {
  double l;           // lower extreme
  double u;           // upper extreme
} Interval;

// quadratic structure
typedef struct {
  private:
    Interval I = {-INFINITY, INFINITY};
  public:
    double a = 0; // a coefficient
    double b = 0; // b coefficient
    double c = 0; // c coefficient
    std::list<Interval> ints = {I}; // intervals list
    double max = 0; // maximum of the quadratic

} Quadratic;

Interval I(const double&, const double&);
bool inRange(const double&, const Interval&);
std::tuple<double, double> get_intersections (const Quadratic&, const Quadratic&);
void get_min_of_two_quadratics (Quadratic& q1, Quadratic& q2);
void print(const Quadratic&);
std::tuple<double, double> get_minimum(const Quadratic&, const Interval&);
std::list<Quadratic> get_max_of_cost(std::list<Quadratic>, Quadratic);
void approximation_grid (std::list<Quadratic>& Q, const std::list<double>& grid);

#endif