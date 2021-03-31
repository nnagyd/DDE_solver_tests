#ifndef LOGISTIC2_H
#define LOGISTIC2_H

// SYSTEM
__forceinline__ __device__ void PerThread_OdeFunction(int tid, int NT, \
			double*    f, double*    x, double     t, \
			RegisterStruct r, SharedParametersStruct s)
{
	f[0] = x[0]*r.xdelay[1]*(r.p[0]-r.xdelay[0]);
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

}

#endif
