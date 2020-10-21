#ifndef ___FOCuS_H___
#define ___FOCuS_H___

#include <vector>
#include <algorithm>
#include "quadratic.h"

typedef struct {
  Quadratic Q0;
  std::list<Quadratic> Q1;
  double global_max;
} Info;


Info FOCuS_step(Info, const double&);
void FOCuS_step_V1(Info&, const double&);

#endif