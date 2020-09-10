#include <iostream>
#include <fstream>
#define MAX_VECTOR_SIZE 256
#include "vectorclass.h"
#include "DDEInit.h"

#define vecSize 4
#define PI 3.1415926535897932385

double y0(double t, int)
{
	return -8 + std::sin(2 * PI * t);
}

double yd0(double t, int)
{
	return std::cos(2 * PI * t) * 2 * PI;
}

template<unsigned int unroll>
inline void F(Vec4d* xd, Vec4d* x, Vec4d* xDelay, Vec4d* p)
{
	for (unsigned int i = 0, j = 0; j < unroll; i += 3, j++)
	{
		xd[i] = 10. * (xDelay[j] - x[i]);
		xd[i + 1] = p[j] * x[i] - x[i + 1] - x[i] * x[i + 2];
		xd[i + 2] = x[i] * x[i + 1] - (8. / 3.) * x[i + 2];
	}
}

template<unsigned int unroll, unsigned int memSize>
struct vars
{
	vars() { return; };
	Vec4d kSum[3*unroll], kAct[3 * unroll], x[3 * unroll], xDelay[unroll], p[unroll], xTmp[3 * unroll];
	Vec4d xb[unroll], xn[unroll], xdb[unroll], xdn[unroll];
	double tb, deltat, pdeltat;
	double t, dt, dtAct, t0;
	unsigned int lastIndex;
};

template<unsigned int unroll, unsigned int memSize>
struct history
{
	double** xVals, ** xdVals, * tVals;
	history()
	{
		xVals = new double* [memSize];
		xdVals = new double* [memSize];
		tVals = (double*)_aligned_malloc(memSize * sizeof(double), 32);

		for (size_t i = 0; i < memSize; i++)
		{
			xVals[i] = (double*)_aligned_malloc(unroll * vecSize * sizeof(double), 32);
			xdVals[i] = (double*)_aligned_malloc(unroll * vecSize * sizeof(double), 32);
		}
	}
};

template<unsigned int unroll, unsigned int memSize>
inline void denseOutput(Vec4d* xDelay, double tDelay, vars<unroll, memSize>* v, history<unroll, memSize> h)
{
	unsigned int id = (*v).lastIndex;
	if (tDelay >= h.tVals[id]) //new step id needed
	{
		while (tDelay >= h.tVals[id]) //find next good value
		{
			id++;
		}
		(*v).lastIndex = id--;

		//load next step from memory
		(*v).tb = h.tVals[id];
		(*v).deltat = h.tVals[id + 1] - (*v).tb;
		(*v).pdeltat = 1.0 / (*v).deltat;

		//load from memory
		for (size_t i = 0; i < unroll; i++)
		{
			(*v).xb[i].load_a(h.xVals[id] + i * vecSize);
			(*v).xn[i].load_a(h.xVals[id + 1] + i * vecSize);
			(*v).xdb[i].load_a(h.xdVals[id] + i * vecSize);
			(*v).xdn[i].load_a(h.xdVals[id + 1] + i * vecSize);
		}
	}

	//calculate dense output
	for (size_t i = 0; i < unroll; i++)
	{
		double theta = (tDelay - (*v).tb) * (*v).pdeltat;
		double thetaM1 = theta - 1;
		xDelay[i] = -thetaM1 * (*v).xb[i] + theta * ((*v).xn[i] + thetaM1 * ((1 - 2 * theta) * ((*v).xn[i] - (*v).xb[i]) + (*v).deltat * (thetaM1 * (*v).xdb[i] + theta * (*v).xdn[i])));
	}
}

int main()
{
	//settings
	const unsigned int nrOfProblems = 262144;
	const unsigned int unroll = 2; //amount of manual unrolling of outer for loop
	const unsigned int outerForLoopStep = vecSize * unroll;
	const unsigned int nrOfInitialPoints = 100;
	const unsigned int nrOfSteps = 10000;
	const unsigned int memSize = nrOfInitialPoints + nrOfSteps;

	//create structs
	struct vars<unroll, memSize> v;
	struct history<unroll, memSize> h;
	v.dt = 10.0 / double(nrOfSteps - 2);
	v.t0 = 13.0 / 28.0;
	v.lastIndex = 0;

	//fill initial function
	double* tValsInit = linspace(-0.5, 0, nrOfInitialPoints);
	double* xValsInit = discretize(y0, tValsInit, nrOfInitialPoints);
	double* xdValsInit = discretize(yd0, tValsInit, nrOfInitialPoints);
	for (size_t i = 0; i < nrOfInitialPoints; i++)
	{
		h.tVals[i] = tValsInit[i];
		for (size_t j = 0; j < unroll * vecSize; j++)
		{
			h.xVals[i][j] = xValsInit[i];
			h.xdVals[i][j] = xdValsInit[i];
		}
	}

	//create parameter list
	double* p_Parameters = linspaceAlligned(0.0, 50.0, nrOfProblems);

	//create mesh
	const unsigned int meshLen = 4;
	unsigned int meshId = 0;
	double mesh[meshLen] = { 1.0, 2.0, 3.0, 10.0 };

	//save end values
	std::ofstream ofsEnd("lorenz_endvals.txt");

	double tStart = seconds();
	for (size_t i = 0; i < nrOfProblems; i += outerForLoopStep) //parameter sweep loop
	{
		//step 1: Initialize dynamic variables t,x and p
		unsigned int memoryId = nrOfInitialPoints - 1;
		v.t = 0.0;
		v.lastIndex = 0;
		meshId = 0;
		for (size_t j = 0, offset = i; j < unroll; j++, offset += 4) //initial setup before integration
		{
			v.p[j].load_a(p_Parameters + offset); //loading parameters from alligned memory

		}
		for (size_t j = 0; j < unroll*3; j++)
		{
			v.x[j] = -8.0; //initial conditions are fix
		}

		//step 2: first double point
		memoryId++;
		h.tVals[memoryId] = v.t;
		for (size_t j = 0; j < unroll; j++) //initial setup before integration
		{
			v.x[j*3 + 1].store_a(h.xVals[memoryId] + vecSize * j);
		}

		//step 3: integration
		while (memoryId + 1 < memSize) //start of integration loop
		{
			//step 3.1: assuming a simple step
			v.dtAct = v.dt;

			//step 3.2: detecting mesh point
			if (meshId < meshLen && mesh[meshId] < v.t + v.dtAct)
			{
				v.dtAct = mesh[meshId] - v.t;
				meshId++;
			}

			//step 3.3: RK4 step
			//K1
			denseOutput<unroll, memSize>(v.xDelay, v.t - v.t0, &v, h);
			F<unroll>(v.kSum, v.x, v.xDelay, v.p);
			for (size_t j = 0; j < unroll; j++)
			{
				v.kSum[j * 3 + 1].store_a(h.xdVals[memoryId] + j * vecSize);
			}

			for (size_t j = 0; j < unroll * 3; j++)
			{
				v.xTmp[j] = v.x[j] + 0.5 * v.dtAct * v.kSum[j];
			}

			//K2
			double tTmp = v.t + 0.5 * v.dtAct;
			denseOutput<unroll,memSize>(v.xDelay, tTmp - v.t0, &v, h);
			F<unroll>(v.kAct, v.xTmp, v.xDelay, v.p);
			for (size_t j = 0; j < unroll * 3; j++)
			{
				v.kSum[j] += 2.0 * v.kAct[j];
				v.xTmp[j] = v.x[j] + 0.5 * v.dtAct * v.kAct[j];
			}

			//K3
			F<unroll>(v.kAct, v.xTmp, v.xDelay, v.p);
			for (size_t j = 0; j < unroll * 3; j++)
			{
				v.kSum[j] += 2.0 * v.kAct[j];
				v.xTmp[j] = v.x[j] + v.dtAct * v.kAct[j];
			}

			//K4
			tTmp = v.t + v.dtAct;
			denseOutput<unroll, memSize>(v.xDelay, tTmp - v.t0, &v, h);
			F<unroll>(v.kAct, v.xTmp, v.xDelay, v.p);
				for (size_t j = 0; j < unroll * 3; j++)
			{
				v.kSum[j] += v.kAct[j];
			}

			//STEP 3.4: calculate and save new values
			v.t += v.dtAct;
			memoryId++;
			h.tVals[memoryId] = v.t;

			for (size_t j = 0; j < unroll*3; j++)
			{
				v.x[j] += (1.0 / 6.0) * v.dtAct * v.kSum[j];
			}
			for (size_t j = 0; j < unroll; j++)
			{
				v.x[j * 3 + 1].store_a(h.xVals[memoryId] + vecSize * j);
			}
		}//end of integration loop

		//save values for reference
		/*std::ofstream ofs("C:/Users/nnagy/Documents/Egyetem/HDS/DDE/results/handtuned_lorenz/handtuned_lorenz_1600.txt");
		ofs << std::setprecision(18) << v.p[0][0] <<"\n";
		for (size_t j = 0; j < memSize; j++)
		{
			ofs << h.tVals[j] << "\t" << h.xVals[j][0] << "\t" << h.xdVals[j][0] << "\n";
		}
		ofs.flush();
		ofs.close();*/

		for (size_t j = 0; j < unroll; j++)
		{
			for (size_t k = 0; k < vecSize; k++)
			{
				ofsEnd << v.t << "\t" << v.p[j][k] << "\t" << v.x[j * 3][k] << "\n";
			}
		}
	}//end of parameter sweep loop

	ofsEnd.flush();
	ofsEnd.close();

	double runTime = seconds() - tStart;

	std::cout << "Parameters = " << nrOfProblems << "  Rollout = " << unroll << "  Runtime = " << runTime << std::endl;

}
