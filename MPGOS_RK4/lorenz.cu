#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

// Solver Configuration
#define __MPGOS_PERTHREAD_SOLVER_DDE4 //DDE4 solver
#define __MPGOS_PERTHREAD_PRECISION double
#define __MPGOS_PERTHREAD_NT    		18144 // NumberOfThreads
#define __MPGOS_PERTHREAD_SD    		3     // SystemDimension
#define __MPGOS_PERTHREAD_DOD    		1     // DenseDimension
#define __MPGOS_PERTHREAD_NDELAY    1     // NumberOfDelays
#define __MPGOS_PERTHREAD_NCP    		1			// ControlParameters
#define __MPGOS_PERTHREAD_NSP       0     // Shared parameters
#define __MPGOS_PERTHREAD_NA        1     // Accessories
#define __MPGOS_PERTHREAD_NDO   		500   // NumberOfPointsOfDenseOutput

#include "SingleSystem_PerThread_DataStructures.cuh"
#include "lorenz_SystemDefinition.cuh"
#include "SingleSystem_PerThread_Interface.cuh"

#define PI 3.14159265358979323846

using namespace std;

void Linspace(vector<double>&, double, double);

int main()
{
  //run configuration
  int NumberOfProblems = __MPGOS_PERTHREAD_NT;
	int InitialPoints = 50;
  int BlockSize = 32;
  int CUDAdevice = 0;
  PrintPropertiesOfSpecificDevice(CUDAdevice);

  //time domain
	double tmin = 0.0;
	double tmax = 10.0;
	double dt = 0.001;

	//initial conditions
	double x0 = -8;
  double y0 = -8;
  double z0 = -8;

	//parameters
	double tau = 13.0/28.0;
	vector <double> rho(NumberOfProblems); //rho
	Linspace(rho,0,50.0);


	//initialize solver
	ProblemSolver Solver(CUDAdevice);
	Solver.SolverOption(ThreadsPerBlock, BlockSize);
	Solver.SolverOption(InitialTimeStep, dt);
	Solver.SolverOption(ActiveNumberOfThreads, NumberOfProblems);

	Solver.SolverOption(DenseOutputTimeStep, -1); //save every point in the dense output
	Solver.SolverOption(DenseOutputVariableIndex, 0, 1); //0. dense output -> 1. system variable

	Solver.SolverOption(Delay, 0, 0, tau); //0. delay -> 0. dense output (-> 1. system variable)



	for (int i = 0; i < NumberOfProblems; i++)
	{
		Solver.SetHost(i,TimeDomain,0,tmin);
		Solver.SetHost(i,TimeDomain,1,tmax);
		Solver.SetHost(i,ActualTime,tmin);

		Solver.SetHost(i,ActualState,0,x0);
		Solver.SetHost(i,ActualState,1,y0);
		Solver.SetHost(i,ActualState,2,z0);
		Solver.SetHost(i,ControlParameters,0,rho[i]);

		//fill initial dense output
		Solver.SetHost(i,DenseIndex,InitialPoints);
		double tDelayAct = -tau;
		double dtDelay = tau/(InitialPoints-1);
		for (int j = 0; j < InitialPoints; j++) {
			Solver.SetHost(i,DenseTime,j,tDelayAct);
			Solver.SetHost(i,DenseState,0,j,-8 + sin(2*PI*tDelayAct));
			Solver.SetHost(i,DenseDerivative,0,j,2*PI*cos(2*PI*tDelayAct));
			tDelayAct += dtDelay;
		}

	}

	//synchronize
	Solver.SynchroniseFromHostToDevice(All);
	Solver.InsertSynchronisationPoint();
	Solver.SynchroniseSolver();

	//solve
	clock_t SimulationStart = clock();

	Solver.Solve();
	Solver.SynchroniseFromDeviceToHost(All);
	Solver.InsertSynchronisationPoint();
	Solver.SynchroniseSolver();

	clock_t SimulationTime = clock()-SimulationStart;
	std::cout << "Simulation Time: "<< 1000*SimulationTime/CLOCKS_PER_SEC << " ms"<<std::endl;

	//save data
	ofstream DataFile;
	DataFile.open ("lorenzDelay.txt");
	DataFile.precision(14);
	DataFile.flags(ios::scientific);

	for (size_t tid = 0; tid < __MPGOS_PERTHREAD_NT; tid++)
	{
		DataFile.width(20); DataFile << Solver.GetHost<double>(tid,ControlParameters,0) << ",";
		DataFile.width(20); DataFile << Solver.GetHost<double>(tid,Accessories,0) << "\n";
	}


	DataFile.flush();
	DataFile.close();

  return 0;
}

void Linspace(vector<double>& x, double B, double E)
{
  double Increment;
	int N = x.size();
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
