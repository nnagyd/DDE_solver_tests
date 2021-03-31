#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

// Solver Configuration
#define __MPGOS_PERTHREAD_SOLVER_DDE4 //RK4 solver
#define __MPGOS_PERTHREAD_PRECISION double
#define __MPGOS_PERTHREAD_NT    		1024 		// NumberOfThreads
#define __MPGOS_PERTHREAD_SD    		1     // SystemDimension
#define __MPGOS_PERTHREAD_DOD    		1     // DenseDimension
#define __MPGOS_PERTHREAD_NDELAY    2     // NumberOfDelays
#define __MPGOS_PERTHREAD_NCP    		1			// ControlParameters
#define __MPGOS_PERTHREAD_NDO   		2005   // NumberOfPointsOfDenseOutput

#include "SingleSystem_PerThread_DataStructures.cuh"
#include "logistic2_SystemDefinition.cuh"
#include "SingleSystem_PerThread_Interface.cuh"

using namespace std;

#define PI 3.14159265358979323846

void Linspace(vector<double>&, double, double, int);
void Discretize(vector<double>&, double f(double), vector<double>);

double f0(double t)
{
    if (t < -1.5)                return std::cos(4 * PI * t);
    if (t < -std::sqrt(2.0))     return t * t;
    if (t < -1.101)              return std::exp(t);
    if (t < -0.5)                return 0;
    return t + 0.5;
}

double fd0(double t)
{
    if (t < -1.5)                return -4 * PI * std::sin(4 * PI * t);
    if (t < -std::sqrt(2.0))     return 2 * t;
    if (t < -1.101)              return std::exp(t);
    if (t < -0.5)                return 0;
    return 1.0;
}
int main()
{
	//run configuration
  int NumberOfProblems = __MPGOS_PERTHREAD_NT;
  int BlockSize = 32;
  int CUDAdevice = 0;
  PrintPropertiesOfSpecificDevice(CUDAdevice);

	//parameters and initial conditions
	double tmin = 		0.0;
	double tmax = 		10.0;
	double dt = 			0.001;
  double pmin = 	0;
  double pmax = 	2;
  vector<double> p(NumberOfProblems,0);
	Linspace(p,pmin,pmax, NumberOfProblems);

	//initial points
	int pointsInit = 200;
	double tInitMin = -2.0;
	double tInitMax = 0;
	vector<double> t0list(pointsInit,0);
	vector<double> x0list(pointsInit,0);
	vector<double> xd0list(pointsInit,0);
	Linspace(t0list,tInitMin,tInitMax, pointsInit);
	Discretize(x0list,f0,t0list);
	Discretize(xd0list,fd0,t0list);


	//initialize solver
	ProblemSolver Solver(CUDAdevice);

	//Set Solver options
	Solver.SolverOption(ThreadsPerBlock, BlockSize);
	Solver.SolverOption(InitialTimeStep, dt);
	Solver.SolverOption(ActiveNumberOfThreads, NumberOfProblems);

	Solver.SolverOption(DenseOutputTimeStep, -1);
	Solver.SolverOption(DenseOutputVariableIndex, 0, 0);

	Solver.SolverOption(Delay, 0, 0, 1.0);
	Solver.SolverOption(Delay, 1, 0, 2.0);


	//fill solver object
	for (int i = 0; i < NumberOfProblems; i++)
	{
		Solver.SetHost(i,TimeDomain,0,tmin);
		Solver.SetHost(i,TimeDomain,1,tmax);
		Solver.SetHost(i,ActualTime,tmin);

		Solver.SetHost(i,ActualState,0,f0(0.0));
		Solver.SetHost(i,ControlParameters,0,p[i]);

		//fill initial dense output
		Solver.SetHost(i,DenseIndex,pointsInit);
		for (int j = 0; j < pointsInit; j++)
		{
			Solver.SetHost(i,DenseTime,j,t0list[j]);
			Solver.SetHost(i,DenseState,0,j,x0list[j]);
			Solver.SetHost(i,DenseDerivative,0,j,xd0list[j]);

		}
	}

	//synchronize
	Solver.SynchroniseFromHostToDevice(All);
	Solver.InsertSynchronisationPoint();
	Solver.SynchroniseSolver();

	//solve
	clock_t SimulationStart = clock();
	Solver.Solve();
	Solver.InsertSynchronisationPoint();
	Solver.SynchroniseSolver();
	clock_t SimulationTime = clock()-SimulationStart;
	std::cout << "Simulation Time: "<< 1000*SimulationTime/CLOCKS_PER_SEC << " ms"<<std::endl;

	//write back to CPU
	Solver.SynchroniseFromDeviceToHost(All);
	Solver.InsertSynchronisationPoint();
	Solver.SynchroniseSolver();

	//write to file
	ofstream DataFile("logistic2.txt");
	DataFile.precision(8);
	DataFile.flags(ios::scientific);
	for (int i = 0; i < NumberOfProblems; i++)
	{
		DataFile.width(13); DataFile << p[i] << ',';
		DataFile.width(13); DataFile << Solver.GetHost<double>(i, ActualState, 0) << '\n';
	}
	DataFile.flush();
	DataFile.close();

	cout << "Test finished!" << endl;
}

void Linspace(vector<double>& x, double B, double E, int N)
{
    double Increment;

	x[0]   = B;

	if ( N>1 )
	{
		x[N-1] = E;
		Increment = (E-B)/(N-1);

		for (int i=1; i<N-1; i++)
		{
			x[i] = B + i*Increment;
		}
	}
}

void Discretize(vector<double>& y, double f(double), vector<double> x)
{
  double Increment;

	for (int i = 0; i < x.size(); i++)
	{
		y[i] = f(x[i]);
	}
}
