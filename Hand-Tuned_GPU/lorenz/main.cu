/*
Compiler options:
nvcc -O3 --std=c++14 --ptxas-options=-v --gpu-architecture=sm_61 -lineinfo -maxrregcount=128 -w --resource-usage -o lorenz.exe main.cu
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "DDEInit.h"
#include "GPUTimer.h"
#include "LorenzSystem.cuh"
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
		const unsigned int nrOfParameters = 1<<17;
		const unsigned int batchSize = 8192;
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
		intSettings.t0 = 13.0/28.0;
		intSettings.dt = 10.0 / nrOfSteps;

		//kernel configuration
		const unsigned int blocksize = 128;
		const unsigned int gridsize = (batchSize + blocksize - 1) / blocksize;
		dim3 block(blocksize);
		dim3 grid(gridsize);

		//parameter stuff, initial CPU and GPU memory
		double * parameterListHost = linspace(0, 50, nrOfParameters);
		double * parameterListDevice;
		cudaMalloc((void**)&parameterListDevice,batchSize * sizeof(double));

		//discretize initial functions
		double * tInit = linspaceDisc(-13.0/28.0, 0.0, nrOfInitialPoints);
		double * y0Init = discretize(y0, tInit, nrOfInitialPoints);
		double * yd0Init = discretize(yd0, tInit, nrOfInitialPoints);

		//copy initial conditions to new bigger arrays
		double * y0 = new double[xValsInitLen];
		double * yd0 = new double[xValsInitLen];
		for (size_t i = 0; i < nrOfInitialPoints; i++)
		{
			for (size_t j = 0; j < batchSize; j++)
			{
				unsigned int idx = i*batchSize + j;
				y0[idx] = y0Init[i];
				yd0[idx] = yd0Init[i];
			}
		}

		//allocate GPU memory
		cudaMalloc((void**)&intSettings.tVals, tValsLen*sizeof(double));
		cudaMalloc((void**)&intSettings.yVals, xValsLen*sizeof(double));
		cudaMalloc((void**)&intSettings.ydVals, xValsLen*sizeof(double));

		//end values
		double * endvalsDev, *endvalsHost;
		cudaMalloc((void**)&endvalsDev, batchSize*sizeof(double));
		endvalsHost = new double[batchSize];

		//copy the initial values to gpu memory
		cudaMemcpy(intSettings.tVals, tInit,tValsInitLen*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.yVals, y0,xValsInitLen*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.ydVals, yd0,xValsInitLen*sizeof(double),cudaMemcpyHostToDevice);

		//information about the run
		printf("Memory size: %zd MB \n",(xValsLen*4+tValsLen)*sizeof(double)/1024/1024);
		printf("Launching kernel with <<<%d,%d>>> in %d batches\n",gridsize,blocksize,nrOfBatches);

		//save to file
		std::ofstream ofs("C:/Users/nnagy/Documents/Egyetem/HDS/DDE/results/gpu_lorenz/endvals.txt");
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
		delete tInit, y0Init, yd0Init;
}
