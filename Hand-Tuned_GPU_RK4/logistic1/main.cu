/*
Compiler options:
-O3 --std=c++14 --ptxas-options=-v --gpu-architecture=sm_61 -lineinfo -maxrregcount=128 -w --resource-usage
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "DDEInit.h"
#include "GPUTimer.h"
#include "Logistic1System.cuh"
#include "GPUDDESolver.cuh"

int main(int argc, char const *argv[])
{
		// get device information
		int dev = 0;
		cudaDeviceProp deviceProp;
		CHECK(cudaGetDeviceProperties(&deviceProp, dev));
		printf("Using Device %d: %s\n", dev, deviceProp.name);
		CHECK(cudaSetDevice(dev));

		//constants
    const unsigned int nrOfInitialPoints = 50;
    const unsigned int nrOfSteps = 10000;
		const unsigned int nrOfPoints = nrOfSteps + 2 * nrOfInitialPoints;
		const unsigned int nrOfParameters = 1<<18;
		const unsigned int batchSize = 1<<13;
		const unsigned int nrOfBatches = (nrOfParameters + batchSize - 1)/batchSize;

		//memory sizes
		size_t tValsInitLen = nrOfInitialPoints;
		size_t xValsInitLen = nrOfInitialPoints * batchSize;
		size_t tValsLen = nrOfPoints;
		size_t xValsLen = nrOfPoints * batchSize;

		//fill integration settings struct
		integrationSettings intSettings;
		intSettings.nrOfInitialPoints = nrOfInitialPoints;
		intSettings.nrOfParameters = batchSize;
		intSettings.nrOfPoints = nrOfPoints;
		intSettings.nrOfSteps = nrOfSteps;
		intSettings.t0 = 1.0;
		intSettings.dt = 10. / (nrOfSteps);

		//kernel configuration
		const unsigned int blocksize = 64;
		const unsigned int gridsize = (batchSize + blocksize - 1) / blocksize;
		dim3 block(blocksize);
		dim3 grid(gridsize);

		//parameter stuff, initial CPU and GPU memory
		double * parameterListHost = linspace(0, 4, nrOfParameters);
		double * parameterListDevice;
		cudaMalloc((void**)&parameterListDevice,batchSize * sizeof(double));

		//discretize initial functions
		double * tInit = linspaceDisc(-1, 0.0, nrOfInitialPoints);
		double * x0Init = discretize(x0, tInit, nrOfInitialPoints);
		double * xd0Init = discretize(xd0, tInit, nrOfInitialPoints);

		//copy initial conditions to new bigger arrays
		double * xlist = new double[xValsInitLen];
		double * xdlist = new double[xValsInitLen];
		for (size_t i = 0; i < nrOfInitialPoints; i++)
		{
			for (size_t j = 0; j < batchSize; j++)
			{
				unsigned int idx = i*batchSize + j;
				xlist[idx] = x0Init[i];
				xdlist[idx] = xd0Init[i];
			}
		}

		//allocate GPU memory
		cudaMalloc((void**)&intSettings.tVals, tValsLen*sizeof(double));
		cudaMalloc((void**)&intSettings.xVals, xValsLen*sizeof(double));
		cudaMalloc((void**)&intSettings.xdVals, xValsLen*sizeof(double));

		//end values
		double * endvalsDev, *endvalsHost;
		cudaMalloc((void**)&endvalsDev, batchSize*sizeof(double));
		endvalsHost = new double[batchSize];

		//copy the initial values to gpu memory
		cudaMemcpy(intSettings.tVals, tInit,tValsInitLen*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.xVals, xlist ,xValsInitLen*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.xdVals, xdlist,xValsInitLen*sizeof(double),cudaMemcpyHostToDevice);

		//information about the run
		printf("Memory size: %zd MB \n",(xValsLen*4+tValsLen)*sizeof(double)/1024/1024);
		printf("Launching kernel with <<<%d,%d>>> in %d batches\n",gridsize,blocksize,nrOfBatches);

		//save to file
		std::ofstream ofs("C:/Users/nnagy/Documents/Egyetem/HDS/DDE/results/gpu_logistic1/envals.txt");
		int id = nrOfInitialPoints + nrOfSteps;

		//execution in batches
		double tStart = seconds();

		//execute in batches
		for (size_t k = 0; k < nrOfBatches; k++)
		{
			CHECK(cudaMemcpy(parameterListDevice,parameterListHost + k*batchSize,batchSize * sizeof(double),cudaMemcpyHostToDevice));

			//launch kernel
			solver<<<grid,block>>>(intSettings, parameterListDevice, endvalsDev);
			CHECK(cudaDeviceSynchronize());

			//copy back endvals global memory
			CHECK(cudaMemcpy(endvalsHost,endvalsDev,batchSize*sizeof(double),cudaMemcpyDeviceToHost));

			for (size_t i = 0; i < batchSize; i++)
			{
					double x = endvalsHost[i];
					ofs << parameterListHost[k*batchSize + i] <<"\t" << x << "\n";
			}
		}
		double tEnd = seconds();
		printf("Execution finished for p = %d parameters in t = %lf s \n", nrOfParameters, (tEnd - tStart) );
		ofs.flush();
		ofs.close();

		//free gpu memomory
		cudaFree(parameterListDevice);

		//delete cpu memory
		delete parameterListHost;
		delete tInit, x0Init, xd0Init, xlist, xdlist;
}
