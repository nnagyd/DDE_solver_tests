#ifndef LORENZ_H
#define LORENZ_H


__forceinline__ __device__ void PerThread_OdeFunction(int tid, int NT,
			double*f, double*x, double t, \
			RegisterStruct r, SharedParametersStruct s)
{
	f[0] = 10.0 * (r.xdelay[0]-x[0]);
	f[1] = r.p[0] * x[0] - x[1] - x[0]*x[2];
	f[2] = x[0]*x[1] - (8.0/3.0) * x[2];
}

// ACCESSORIES
__forceinline__ __device__ void PerThread_ActionAfterSuccessfulTimeStep(\
	int tid, int NT, \
	RegisterStruct &r, SharedParametersStruct s)
{

}

__forceinline__ __device__ void PerThread_Initialization(\
	int tid, int NT, \
	RegisterStruct &r, SharedParametersStruct s)
{

}

__forceinline__ __device__ void PerThread_Finalization(\
			int tid, int NT, \
			RegisterStruct &r, SharedParametersStruct s)
{
	r.acc[0] = r.x[0];
}




#endif
