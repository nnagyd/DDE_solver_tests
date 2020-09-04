// Lorenz.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "source/ParallelDDE_RK4.h"

#define PI 3.1415926535897932385

double y0(double t, int)
{
    return - 8 + sin(2 * PI * t);
}

double yd0(double t, int)
{
    return  cos(2 * PI * t) * 2 * PI;
}

void lorenz(Vec4d* xd, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    xd[0] = 10.0 * (delay[0] - x[0]);
    xd[1] = p[0] * x[0] - x[1] - x[2] * x[0];
    xd[2] = x[0] * x[1] - (8.0 / 3.0) * x[2];
}

int main()
{
    //global constants
    const unsigned nrOfSteps = 10000;
    const unsigned nrOfInitialPoints = 100;
    const unsigned nrOfParameters = 18144;
    const unsigned unroll = 8;

    double* parameters = linspace(0, 50, nrOfParameters);

    //global integrator settings
    ParallelDDE_RK4<3,1,1,1,0,unroll> solver;

    solver.setRange(0.0, 10.0);
    solver.setNrOfSteps(nrOfSteps);
    solver.setNrOfInitialPoints(nrOfInitialPoints);

    solver.addDelay(0, 13.0 / 28.0, 0);
    solver.setX0(-8.0);

    solver.calculateInitialTvals(-0.5);
    solver.calculateInitialXvals(y0, yd0,1,0);
    solver.calculateIntegrationMesh();

    std::ofstream ofs("endvals.txt");
    ofs << std::setprecision(16);

    double start = seconds();
    for (size_t i = 0; i < nrOfParameters; i+= vecSize*unroll)
    {
        solver.setParameters(parameters + i, 0);

        solver.integrate(lorenz);

        double * endvals = solver.getEndValues(0);
        for (size_t j = 0; j < vecSize * unroll; j++)
        {
            ofs << parameters[i + j] << "\t" << endvals[j] << std::endl;
        }

        delete endvals;
    }
    double runtime = seconds() - start;

    std::cout << "Parameter number = " << nrOfParameters << " rollout = " <<unroll <<  " runtime = " << runtime << " s \n";

    ofs.flush();
    ofs.close();

    delete parameters;
}
