// Lorenz.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "source/ParallelDDE_DP5.h"

#define PI 3.1415926535897932385

double y0(double t, int)
{
    return -8 + sin(2 * PI * t);
}

void lorenz(Vec4d* xd, Vec4d t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    xd[0] = 10.0 * (delay[0] - x[0]);
    xd[1] = p[0] * x[0] - x[1] - x[2] * x[0];
    xd[2] = x[0] * x[1] - (8.0 / 3.0) * x[2];
}

int main()
{
    //global constants
    const unsigned nrOfSteps = 2000;
    const unsigned nrOfParameters = 262144;

    double* parameters = linspace(0, 50, nrOfParameters);

    //global integrator settings
    ParallelDDE_DP5<3, 1, 1, 1> solver;

    solver.setRange(0.0, 10.0);
    solver.setNrOfSteps(nrOfSteps);
    solver.setGlobalTol(1e-8);

    solver.addDelay(0, 13.0 / 28.0, 0);
    solver.setX0(-8.0);

    solver.addInitialFunction(y0, 1, 0);
    solver.calculateIntegrationMesh();

    std::ofstream ofs("endvals.txt");
    ofs << std::setprecision(16);

    double start = seconds();
    for (size_t i = 0; i < nrOfParameters; i += vecSize)
    {
        solver.setParameters(parameters + i, 0);

        solver.integrate(lorenz);

        double* endvals = solver.getEndValues(0);
        int* steps = solver.getStepCount();
        for (size_t j = 0; j < vecSize; j++)
        {
            ofs << parameters[i + j] << "\t" << endvals[j] << "\t" << steps[j] << std::endl;
        }
        delete steps, endvals;
    }
    double runtime = seconds() - start;

    std::cout << "Parameter number = " << nrOfParameters << " runtime = " << runtime << " s \n";

    ofs.flush();
    ofs.close();

    delete parameters;
}
