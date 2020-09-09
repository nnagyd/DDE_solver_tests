#include <iostream>
#include "source/ParallelDDE_DP5.h";


void analitic2(Vec4d* xd, Vec4d t, Vec4d* x, Vec4d* xDelay, Vec4d* p)
{
    xd[0] = -p[0] * xDelay[0];
}

double x0(double t, int)
{
    if (t < -2) return 0;
    if (t < -1) return 1;
    return 0;
}

int main()
{
    //solver object
    ParallelDDE_DP5<1, 1, 1, 1> solver;

    //steps
    double tmax = 25.0;
    int steps = 1000;
    solver.setRange(0.0, tmax);
    solver.setNrOfSteps(steps + 5);
    solver.setStepSize(tmax / steps);

    //adaptive settings
    solver.setGlobalTol(1e-10);
    solver.setStepSizeBounds(1e-6, 1.0);

    //delays
    solver.addDelay(0, 3.0, 0);

    //mesh
    double* initC0 = new double[2]{ -1,-2 };
    solver.setInitialDisc(NULL, 0, initC0, 2);
    solver.calculateIntegrationMesh();
    solver.printMesh();

    //initial conditions
    solver.setX0(0);
    solver.addInitialFunction(x0,0, 0);

    //set parameters
    solver.setParameters(1.0, 0);

    solver.integrate(analitic2);

    solver.save("C:/Users/nnagy/Documents/Egyetem/HDS/DDE/results/sajat_dp5_analitic2/tol_10.txt");
}