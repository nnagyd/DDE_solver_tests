#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "source/ParallelDDE_RK4.h"

void f(Vec4d* xd, double t, Vec4d* x, Vec4d* xdelay, Vec4d* p)
{
    xd[0] = -9.81;
    xd[1] = x[0] - p[0] * xdelay[0];
}

void eventLocation(Vec4d* lst, double t, Vec4d* x, Vec4d* xdelay, Vec4d* p)
{
    lst[0] = x[1];
}

void eventIntervention(int eventId, int eventDir, Vec4db mask, double t, Vec4d* x, Vec4d* xdelay, Vec4d* p)
{
    if (eventId == 0 && eventDir == -1) //event in negative direction
    {
        x[0] = select(mask, -0.95 * x[0], x[0]);
        p[1] = select(mask, p[1] + Vec4d(1.0), p[1]);
    }
}

double v0(double t, int)
{
    return 0;
}

double vd0(double t, int)
{
    return 0;
}

int main()
{
    const int steps = 10000;
    const int unroll = 2;
    const int N = 512;
    const int p_len = N;
    double* p = linspace(0.1, 2, p_len);
    const int tau_len = N;
    double* tau = linspace(0.1, 2, p_len);

    //basic setup
    ParallelDDE_RK4<2, 1, 1, 2, 1, unroll> solver;
    solver.setRange(0.0, 60.0);
    solver.setNrOfSteps(steps);
    solver.setNrOfExtraPoints(5000);

    //delay
    solver.addDelay(0, 0.2, 0);

    //set initial conditions
    solver.setX0(0.0, 0); //v0
    solver.setX0(20.0, 1); //x0

    //initial points
    solver.setNrOfInitialPoints(2);
    solver.calculateInitialTvals(-3);
    solver.calculateInitialXvals(v0, vd0, 0, 0);

    //events
    solver.setEventPrecision(1e-6);
    solver.setEventFlushLookback(50);
    solver.addEventLocationFunction(&eventLocation);
    solver.addEventInterventionFunction(&eventIntervention);

    std::ofstream ofs("bounces.txt");
    ofs << std::setprecision(16);

    double iStart = seconds();
    for (size_t k = 0; k < tau_len; k++)
    {
        //add delays
        solver.addDelay(0, tau[k], 0);

        //calculate mesh
        solver.calculateIntegrationMesh();

        for (size_t i = 0; i < p_len; i += unroll * vecSize)
        {
            //add parameters
            solver.setParameters(p + i, 0);
            solver.setParameters(0.0, 1);

            //integrate
            solver.integrate(f);

            //save end values
            double* endvals = solver.getEndValues(1);
            double* pars1 = solver.getParameters(1);

            for (size_t j = 0; j < unroll * vecSize; j++)
            {
                ofs << p[i + j] << "\t" << tau[k] << "\t" << pars1[j] << "\n";
            }

            delete endvals, pars1;
        }
    }

    double iTime = seconds() - iStart;

    std::cout << "Unroll= " << unroll << " Parameters= " << N*N << " Runtime= " << iTime << "s" << std::endl;

    ofs.flush();
    ofs.close();

}
