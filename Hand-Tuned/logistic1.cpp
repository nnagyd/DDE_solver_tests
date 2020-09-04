#include <iostream>
#include <fstream>
#define MAX_VECTOR_SIZE 256
#include "vectorclass.h"
#include "DDEInit.h"

#define vecSize 4

double x0(double t,int)
{
	return 1.5 - cos(t);
}

double xd0(double t,int)
{
	return sin(t);
}

inline void F(Vec4d* xd, Vec4d x, Vec4d xDelay, Vec4d p)
{
	*xd =  x * (p - xDelay);
}

template<unsigned int unroll,unsigned int memSize>
struct vars
{
	vars() { return; };
	Vec4d kSum[unroll], kAct[unroll], x[unroll], xDelay[unroll], p[unroll], xTmp[unroll];  //24 * unroll [double]
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
		tVals = new double [memSize];

		for (size_t i = 0; i < memSize; i++)
		{
			xVals[i] = (double*)_aligned_malloc(unroll * vecSize * sizeof(double), 32);
			xdVals[i] = (double*)_aligned_malloc(unroll * vecSize * sizeof(double), 32);
		}
	}
};

template<unsigned int unroll, unsigned int memSize>
inline void denseOutput(Vec4d * xDelay, double tDelay, vars<unroll, memSize> *v, history<unroll, memSize> h)
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
		xDelay[i] = -thetaM1 * (*v).xb[i] + theta * ((*v).xn[i] + thetaM1 * ((1 - 2 * theta) * ((*v).xn[i] - (*v).xb[i]) + (*v).deltat * ( thetaM1 *  (*v).xdb[i] + theta * (*v).xdn[i])));
	}


}

int main()
{
	//settings
	const unsigned int nrOfProblems = 18144;
	const unsigned int unroll = 8; //amount of manual unrolling of outer for loop
	const unsigned int outerForLoopStep = vecSize * unroll;
	const unsigned int nrOfInitialPoints = 100;
	const unsigned int nrOfSteps = 10000;
	const unsigned int memSize = nrOfInitialPoints + nrOfSteps;

	//create structs
	struct vars<unroll,memSize> v;
	struct history<unroll,memSize> h;
	v.dt = 10.0 / double(nrOfSteps-2);
	v.t0 = 1;
	v.lastIndex = 0;

	//fill initial function
	double* tValsInit = linspace(-1.0, 0, nrOfInitialPoints);
	double* xValsInit = discretize(x0, tValsInit, nrOfInitialPoints);
	double* xdValsInit = discretize(xd0, tValsInit, nrOfInitialPoints);
	for (size_t i = 0; i < nrOfInitialPoints; i++)
	{
		h.tVals[i] = tValsInit[i];
		for (size_t j = 0; j < unroll*vecSize; j++)
		{
			h.xVals[i][j] = xValsInit[i];
			h.xdVals[i][j] = xdValsInit[i];
		}
	}

	//create parameter list
	double* p_Parameters = linspaceAlligned(0.0, 4.0, nrOfProblems);

	//create mesh
	const unsigned int meshLen = 4;
	unsigned int meshId = 0;
	double mesh[meshLen] = { 1.0, 2.0, 3.0, 10.0 };

	//save end values
	std::ofstream ofsEnd("endvals.txt");

	double tStart = seconds();
	for (size_t i = 0; i < nrOfProblems; i+= outerForLoopStep) //parameter sweep loop
	{
		//step 1: Initialize dynamic variables t,x and p
		unsigned int memoryId = nrOfInitialPoints - 1;
		v.t = 0.0;
		v.lastIndex = 0;
		meshId = 0;
		for (size_t j = 0, offset = i; j < unroll; j++, offset += 4) //initial setup before integration
		{
			v.p[j].load_a(p_Parameters + offset); //loading parameters from alligned memory
			v.x[j].load_a(h.xVals[memoryId] + j*vecSize);
		}


		//step 2: first double point
		memoryId++;
		h.tVals[memoryId] = v.t;
		for (size_t j = 0; j < unroll; j++) //initial setup before integration
		{
			v.x[j].store_a(h.xVals[memoryId] + vecSize * j);
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
			denseOutput(v.xDelay, v.t - v.t0, &v, h);
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kSum[j]), v.x[j], v.xDelay[j], v.p[j]);
				v.kSum[j].store_a(h.xdVals[memoryId] + j * vecSize);
				v.xTmp[j] = v.x[j] + 0.5 * v.dtAct * v.kSum[j];
			}

			//K2
			double tTmp = v.t + 0.5*v.dtAct;
			denseOutput(v.xDelay, tTmp - v.t0, &v, h);
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kAct[j]), v.xTmp[j], v.xDelay[j], v.p[j]);
				v.kSum[j] += 2*v.kAct[j];
				v.xTmp[j] = v.x[j] + 0.5 * v.dtAct * v.kAct[j];
			}

			//K3
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kAct[j]), v.xTmp[j], v.xDelay[j], v.p[j]);
				v.kSum[j] += 2*v.kAct[j];
				v.xTmp[j] = v.x[j] + v.dtAct * v.kAct[j];
			}

			//K4
			tTmp = v.t + v.dtAct;
			denseOutput(v.xDelay, tTmp - v.t0, &v, h);
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kAct[j]), v.xTmp[j], v.xDelay[j], v.p[j]);
				v.kSum[j] += v.kAct[j];
			}

			//STEP 3.4: calculate and save new values
			v.t += v.dtAct;
			memoryId++;
			h.tVals[memoryId] = v.t;
			for (size_t j = 0; j < unroll; j++)
			{
				v.x[j] += (1.0 / 6.0) * v.dtAct * v.kSum[j];
				v.x[j].store_a(h.xVals[memoryId] + vecSize * j);
			}
		}//end of integration loop

		//save values for reference
		/*std::ofstream ofs("handtuned_logistic1_100.txt");
		ofs << std::setprecision(18) << v.p[0][1];
		for (size_t j = 0; j < memSize; j++)
		{
			ofs << h.tVals[j] << "\t" << h.xVals[j][1] << "\t" << h.xdVals[j][1] << "\n";
		}
		ofs.flush();
		ofs.close();*/

		for (size_t j = 0; j < unroll; j++)
		{
			for (size_t k = 0; k < vecSize; k++)
			{
				ofsEnd << v.t << "\t" << v.p[j][k] << "\t" << v.x[j][k] << "\n";
			}
		}

		if (i % 100 == 1) //i % 100 == 1)
		{
			for (size_t j = 0; j < unroll; j++)
			{
				for (size_t k = 0; k < vecSize; k++)
				{
					std::cout << "t = " << v.t << "\tp=" << v.p[j][k] << "\tx=" << v.x[j][k] << "\n";
				}
			}
		}
	}//end of parameter sweep loop

	ofsEnd.flush();
	ofsEnd.close();

	double runTime = seconds() - tStart;

	std::cout << "Parameters = " << nrOfProblems << "  Rollout = " << unroll << "  Runtime = " << runTime << std::endl;

}
