/*
Compiler options:
nvcc -O3 --std=c++14 --ptxas-options=-v --gpu-architecture=sm_61 -lineinfo -maxrregcount=128 -w --resource-usage -o logistic2.exe main.cu
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "DDEInit.h"
#include "GPUTimer.h"
#include "Logistic2System.cuh"
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
		const unsigned int discInit = 4;
    const unsigned int nrOfInitialPoints = 50 + 2*discInit;
    const unsigned int nrOfSteps = 10015;
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
		intSettings.t0[0] = 1.0;
		intSettings.t0[1] = 2.0;
		intSettings.dt = 10. / (10000);

		//kernel configuration
		const unsigned int blocksize = 64;
		const unsigned int gridsize = (batchSize + blocksize - 1) / blocksize;
		dim3 block(blocksize);
		dim3 grid(gridsize);

		//parameter stuff, initial CPU and GPU memory
		double * parameterListHost = linspace(0, 2, nrOfParameters);
		double * parameterListDevice;
		cudaMalloc((void**)&parameterListDevice,batchSize * sizeof(double));

		//create mesh
		double discMesh[discInit+1] = { -1.5, -std::sqrt(2), -1.101, -0.5, 0 };
		double doublePoints[4] = { 0.5, -std::sqrt(2) + 2, -1.101 + 2, 1.5 };
		unsigned int meshLen = 6*(discInit+1);
		unsigned int meshId = 0;
		double * mesh = new double[meshLen];
		for (size_t i = 0; i < discInit+1; i++)
		{
			for (size_t j = 0; j < 6; j++)
			{
				size_t linId = i * 6 + j;
				mesh[linId] = discMesh[i] + 1.0 +  double(j);
			}
		}
		mesh = filter(mesh, &meshLen, 0);
		sort(mesh, meshLen);
		int * meshType = new int[meshLen];

		for (size_t i = 0; i < meshLen; i++)
		{
			//assuming a simple point
			meshType[i] = 1;
			for (size_t j = 0; j < 4; j++)
			{
				if (std::abs(mesh[i] - doublePoints[j]) < 1e-10)
				{
					meshType[i] = 2;
				}
			}
		}


		//fill initial function
		double disc[discInit] = { -1.5, -std::sqrt(2), -1.101, -0.5 };
		double* tInit = linspaceDisc(-2.0, 0.0, nrOfInitialPoints - 2*discInit, disc, discInit);
		double* x0Init = discretize(x0, tInit, nrOfInitialPoints);
		double* xd0Init = discretize(xd0, tInit, nrOfInitialPoints);

		//print mesh
		/*for (size_t i = 0; i < meshLen; i++)
		{
			printf("id=%f\tt=%lf\ttype=%d\n",i,mesh[i],meshType[i]);
		}*/


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
		cudaMalloc((void**)&intSettings.mesh, meshLen*sizeof(double));
		cudaMalloc((void**)&intSettings.meshType, meshLen*sizeof(int));
		cudaMalloc((void**)&intSettings.tVals, tValsLen*sizeof(double));
		cudaMalloc((void**)&intSettings.xVals, xValsLen*sizeof(double));
		cudaMalloc((void**)&intSettings.xdVals, xValsLen*sizeof(double));

		//end values
		double * endvalsDev, *endvalsHost;
		cudaMalloc((void**)&endvalsDev, batchSize*sizeof(double));
		endvalsHost = new double[batchSize];

		//copy the initial values to gpu memory
		cudaMemcpy(intSettings.mesh, mesh, meshLen*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.meshType, meshType, meshLen*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.tVals, tInit,tValsInitLen*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.xVals, xlist ,xValsInitLen*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(intSettings.xdVals, xdlist,xValsInitLen*sizeof(double),cudaMemcpyHostToDevice);

		//information about the run
		printf("Memory size: %zd MB \n",(xValsLen*4+tValsLen)*sizeof(double)/1024/1024);
		printf("Launching kernel with <<<%d,%d>>> in %d batches\n",gridsize,blocksize,nrOfBatches);

		//save to file
		std::ofstream ofs("C:/Users/nnagy/Documents/Egyetem/HDS/DDE/results/gpu_logistic2/endvals.txt");
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
