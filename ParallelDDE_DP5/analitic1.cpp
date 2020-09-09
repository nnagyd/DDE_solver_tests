#include <iostream>
#include "source/ParallelDDE_DP5.h"

void analitic1(Vec4d* dx, Vec4d t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    dx[0] = -p[0] * delay[0];
}


double x0(double t, int)
{
    return 1;
}

int main()
{
    //initialize solver
    const unsigned Vars = 1;
    const unsigned Delays = 1;
    const unsigned DenseVars = 1;
    const unsigned ParametersPerEquation = 1;
    ParallelDDE_DP5<Vars, Delays, DenseVars, ParametersPerEquation> solver;

    //basic setup
    solver.setRange(0.0, 10);
    solver.setNrOfSteps(1000);
    solver.setGlobalTol(1e-10);

    //delay
    solver.addDelay(0, 1, 0);
    solver.setStepSizeBounds(1e-10, 0.33);
    solver.setStepSize(1e-2);

    //initial function
    solver.addInitialFunction(x0, 0, 0);

    //mesh (only discontinous point is t = 0)
    solver.calculateIntegrationMesh();
    solver.printMesh();
    solver.printLookupTables();

    //initial conditions
    solver.setX0(1.0, 0);

    //parameter
    solver.setParameters(1.0, 0);

    //integrate
    solver.integrate(analitic1);

    //save
    solver.save("tol_10.txt");
}
