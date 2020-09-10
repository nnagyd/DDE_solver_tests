#include <iostream>
#include <fstream>
#include <iomanip>
#define MAX_VECTOR_SIZE 256
#include "vectorclass.h"
#include "DDEInit.h"

#define vecSize 4
#define PI 3.14159265358979323846

double x0(double t, int dir)
{
	if (t < -1.5 || (t == -1.5 && dir == -1)) return std::cos(4 * PI * t);
	if (t < -std::sqrt(2) || (t == -std::sqrt(2) && dir == -1)) return t*t;
	if (t < -1.101 || (t == -1.101 && dir == -1)) return std::exp(t);
	if (t < -0.5 || (t == -0.5 && dir == -1)) return 0;
	return t + 0.5;
}

double xd0(double t, int dir)
{
	if (t < -1.5 || (t == -1.5 && dir == -1)) return -std::sin(4 * PI * t) * 4 * PI;
	if (t < -std::sqrt(2) || (t == -std::sqrt(2) && dir == -1)) return 2.0 * t;
	if (t < -1.101 || (t == -1.101 && dir == -1)) return std::exp(t);
	if (t < -0.5 || (t == -0.5 && dir == -1)) return 0;
	return 1.0;
}

inline void F(Vec4d* xd, Vec4d x, Vec4d xDelay0, Vec4d xDelay1, Vec4d p)
{
	xd[0] = x * xDelay1 * (p - xDelay0);
}

template<unsigned int unroll, unsigned int memSize>
struct vars
{
	vars() { return; };
	Vec4d kSum[unroll], kAct[unroll], x[unroll], xDelay0[unroll], xDelay1[unroll], p[unroll], xTmp[unroll];  //24 * unroll [double]
	Vec4d xb[2][unroll], xn[2][unroll], xdb[2][unroll], xdn[2][unroll];
	double tb[2], deltat[2], pdeltat[2];
	double t , dt, dtAct, t0[2];
	unsigned int lastIndex[2];
};

template<unsigned int unroll, unsigned int memSize>
struct history
{
	double** xVals, ** xdVals, * tVals;
	history()
	{
		xVals = new double* [memSize];
		xdVals = new double* [memSize];
		tVals = (double*)_aligned_malloc(memSize* sizeof(double), 32);

		for (size_t i = 0; i < memSize; i++)
		{
			xVals[i] = (double*)_aligned_malloc(unroll * vecSize * sizeof(double), 32);
			xdVals[i] = (double*)_aligned_malloc(unroll * vecSize * sizeof(double), 32);
		}
	}
};

template<unsigned int unroll, unsigned int memSize>
inline void denseOutput(Vec4d  * xDelay,unsigned int delayId, double tDelay, vars<unroll, memSize>* v, history<unroll, memSize> h)
{
	unsigned int id = (*v).lastIndex[delayId];
	if (tDelay >= h.tVals[id]) //new step id needed
	{
		while (tDelay >= h.tVals[id]) //find next good value
		{
			id++;
		}
		(*v).lastIndex[delayId] = id--;

		//load next step from memory
		(*v).tb[delayId] = h.tVals[id];
		(*v).deltat[delayId] = h.tVals[id + 1] - (*v).tb[delayId];
		(*v).pdeltat[delayId] = 1.0 / (*v).deltat[delayId];

		//load from memory
		for (size_t i = 0; i < unroll; i++)
		{
			(*v).xb[delayId][i].load_a(h.xVals[id] + i * vecSize);
			(*v).xn[delayId][i].load_a(h.xVals[id + 1] + i * vecSize);
			(*v).xdb[delayId][i].load_a(h.xdVals[id] + i * vecSize);
			(*v).xdn[delayId][i].load_a(h.xdVals[id + 1] + i * vecSize);
		}
	}

	//calculate dense output
	for (size_t i = 0; i < unroll; i++)
	{
		double theta = (tDelay - (*v).tb[delayId]) * (*v).pdeltat[delayId];
		double thetaM1 = theta - 1;
		xDelay[i] = -thetaM1 * (*v).xb[delayId][i] + theta * ((*v).xn[delayId][i] + thetaM1 * ((1 - 2 * theta) * ((*v).xn[delayId][i] - (*v).xb[delayId][i]) + (*v).deltat[delayId] * (thetaM1 * (*v).xdb[delayId][i] + theta * (*v).xdn[delayId][i])));
	}
}

int main()
{
	//settings
	const unsigned int nrOfProblems = 18144; //181440;
	const unsigned int unroll = 8; //amount of manual unrolling of outer for loop
	const unsigned int outerForLoopStep = vecSize * unroll;
	const unsigned int discInit = 4;
	const unsigned int nrOfInitialPoints = 100;
	const unsigned int nrOfInitialPointsWithDisc = nrOfInitialPoints + 2 * discInit;
	const unsigned int nrOfSteps = 10000;
	const unsigned int memSize = nrOfInitialPointsWithDisc + nrOfSteps + 14;

	//create structs
	struct vars<unroll, memSize> v;
	struct history<unroll, memSize> h;
	v.dt = 10.0 / double(nrOfSteps - 2);
	v.t0[0] = 1.0;
	v.t0[1] = 2.0;
	v.lastIndex[0] = 0;
	v.lastIndex[1] = 0;

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
			if (std::abs(mesh[i] - doublePoints[j]) < 1e-12)
			{
				meshType[i] = 2;
			}
		}
	}


	//fill initial function
	double disc[discInit] = { -1.5, -std::sqrt(2), -1.101, -0.5 };
	double* tValsInit = linspaceDisc(-2.0, 0.0, nrOfInitialPoints, disc, discInit);
	double* xValsInit = discretize(x0, tValsInit, nrOfInitialPointsWithDisc);
	double* xdValsInit = discretize(xd0, tValsInit, nrOfInitialPointsWithDisc);

	for (size_t i = 0; i < nrOfInitialPointsWithDisc; i++)
	{
		h.tVals[i] = tValsInit[i];
		for (size_t j = 0; j < unroll * vecSize; j++)
		{
			h.xVals[i][j] = xValsInit[i];
			h.xdVals[i][j] = xdValsInit[i];
		}
	}

	//create parameter list
	double* p_Parameters = linspaceAlligned(0.0, 2.0, nrOfProblems);

	//save end values
	std::ofstream ofsEnd("endvals.txt");

	double tStart = seconds();
	for (size_t i = 0; i < nrOfProblems; i += outerForLoopStep) //parameter sweep loop
	{
		//step 1: Initialize dynamic variables t,x and p
		unsigned int memoryId = nrOfInitialPointsWithDisc - 1;
		v.t = 0.0;
		v.lastIndex[0] = 0;
		v.lastIndex[1] = 0;
		meshId = 0;
		for (size_t j = 0, offset = i; j < unroll; j++, offset += 4) //initial setup before integration
		{
			v.p[j].load_a(p_Parameters + offset); //loading parameters from alligned memory
			v.x[j].load_a(h.xVals[memoryId] + j * vecSize);
		}

		//step 2: first double point
		memoryId++;
		h.tVals[memoryId] = v.t;
		for (size_t j = 0; j < unroll; j++) //initial setup before integration
		{
			v.x[j].store_a(h.xVals[memoryId] + vecSize * j);
		}

		//step 3: integration
		int stepType = 0;
		while (memoryId + 1 < memSize) //start of integration loop
		{
			//step 3.1: if previous step was a double mesh point
			if (stepType == 2)
			{
				v.dtAct = 2 * (1e-12);
			}
			else
			{
				v.dtAct = v.dt;
			}

			//step 3.2: assuming a simple step
			stepType = 0;


			//step 3.3: detecting mesh point
			if (meshId < meshLen && mesh[meshId] < v.t + v.dtAct)
			{
				stepType = meshType[meshId];

				if (stepType == 1)
				{
					v.dtAct = mesh[meshId] - v.t;
					meshId++;
				}

				if (stepType == 2)
				{
					v.dtAct = mesh[meshId] - v.t - 1e-12;
					meshId++;
				}


			}

			//step 3.4: RK4 step
			//K1
			denseOutput(v.xDelay0,0, v.t - v.t0[0], &v, h);
			denseOutput(v.xDelay1,1, v.t - v.t0[1], &v, h);
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kSum[j]), v.x[j], v.xDelay0[j], v.xDelay1[j], v.p[j]);
				v.kSum[j].store_a(h.xdVals[memoryId] + j * vecSize);
				v.xTmp[j] = v.x[j] + 0.5 * v.dtAct * v.kSum[j];
			}

			//K2
			double tTmp = v.t + 0.5 * v.dtAct;
			denseOutput(v.xDelay0, 0, tTmp - v.t0[0], &v, h);
			denseOutput(v.xDelay1, 1, tTmp - v.t0[1], &v, h);
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kAct[j]), v.xTmp[j], v.xDelay0[j], v.xDelay1[j], v.p[j]);
				v.kSum[j] += 2 * v.kAct[j];
				v.xTmp[j] = v.x[j] + 0.5 * v.dtAct * v.kAct[j];
			}

			//K3
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kAct[j]), v.xTmp[j], v.xDelay0[j], v.xDelay1[j], v.p[j]);
				v.kSum[j] += 2 * v.kAct[j];
				v.xTmp[j] = v.x[j] + v.dtAct * v.kAct[j];
			}

			//K4
			tTmp = v.t + v.dtAct;
			denseOutput(v.xDelay0, 0, tTmp - v.t0[0], &v, h);
			denseOutput(v.xDelay1, 1, tTmp - v.t0[1], &v, h);
			for (size_t j = 0; j < unroll; j++)
			{
				F(&(v.kAct[j]), v.xTmp[j], v.xDelay0[j], v.xDelay1[j], v.p[j]);
				v.kSum[j] += v.kAct[j];
			}

			//STEP 3.5: calculate and save new values
			v.t += v.dtAct;
			memoryId++;
			h.tVals[memoryId] = v.t;
			for (size_t j = 0; j < unroll; j++)
			{
				v.x[j] += (1.0 / 6.0) * v.dtAct * v.kSum[j];
				v.x[j].store_a(h.xVals[memoryId] + vecSize * j);
			}
		}//end of integration loop

		/*//save values for reference
		std::ofstream ofs("handtuned_logistic2_68.txt");
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

	}//end of parameter sweep loop

	ofsEnd.flush();
	ofsEnd.close();

	double runTime = seconds() - tStart;

	std::cout << "Parameters = " << nrOfProblems << "  Rollout = " << unroll << "  Runtime = " << runTime << std::endl;

}
