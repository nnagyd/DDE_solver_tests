#ifndef LORENZ_SYSTEM_
#define LORENZ_SYSTEM_

#include <cmath>

//------------------------------Global constants----------------------------------------
const double PI = 3.14159265358979323846264338327950;

//------------------------------ ODE Function --------------------------------
__forceinline__ __device__ void f(double * xd, const double t, const double x, const double xDelay, const double p)
{
	*xd =  x*(p-xDelay);
}

//------------------------------ Initial functions --------------------------------
double x0(double t, int dir)
{
	return 1.5 - cos(t);
}
double xd0(double t, int dir)
{
	return sin(t);
}



#endif
