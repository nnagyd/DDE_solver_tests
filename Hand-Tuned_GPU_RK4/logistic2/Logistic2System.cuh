#ifndef LORENZ_SYSTEM_
#define LORENZ_SYSTEM_

#include <cmath>

//------------------------------Global constants----------------------------------------
const double PI = 3.14159265358979323846264338327950;

//------------------------------ ODE Function --------------------------------
__forceinline__ __device__ void f(double * xd, const double t, const double x, const double xDelay1, const double xDelay2, const double p)
{
	*xd =  x*xDelay2*(p-xDelay1);
}

//------------------------------ Initial functions --------------------------------
double x0(double t, int dir)
{
	if (t < -1.5 || (t == -1.5 && dir == -1))                           return std::cos(4 * PI * t);
  if (t < -std::sqrt(2.0) || (t == -std::sqrt(2.0) && dir == -1))     return t * t;
  if (t < -1.101 || (t == -1.101 && dir == -1))                       return std::exp(t);
  if (t < -0.5 || (t == -0.5 && dir == -1))                           return 0;
  return t + 0.5;
}
double xd0(double t, int dir)
{
	if (t < -1.5 || (t == -1.5 && dir == -1))                           return -4 * PI * std::sin(4 * PI * t);
  if (t < -std::sqrt(2.0) || (t == -std::sqrt(2.0) && dir == -1))     return 2 * t;
  if (t < -1.101 || (t == -1.101 && dir == -1))                       return std::exp(t);
  if (t < -0.5 || (t == -0.5 && dir == -1))                           return 0;
  return 1.0;
}



#endif
