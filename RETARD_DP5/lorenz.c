#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "retard.h"

#define PI 3.1415926535897932385

#define nrs 3 //dimension of the system

/* initial function
i: component index
x: independent variable
*/

double * parameters;
int nrOfParameters, parId;

FILE* fp;

double phi (unsigned i, double x)
{
  return -8 + sin(2*PI*x);
}

//definition of DDE
void fcn (unsigned n, double x, double* y, double* f)
{
  double p = 50; //parameters[parId];
  double ydelay = ylag(1, x-(13.0/28.0), phi); //(component, point)
  f[0] = 10*(ydelay - y[0]);
  f[1] = p*y[0]-y[1]-y[2]*y[0];
  f[2] = y[0]*y[1]-(8.0/3.0)*y[2];
}

//save output with this method
void solout (long nr, double xold, double x, double* y, unsigned n, int* irtrn)
{
  double ydelay = ylag(1, x-(13.0/28.0), phi); //(component, point)
  fprintf (fp,"%17.20lf\t%17.20lf\n", x, y[0]);
}

int main(int argc, char const *argv[])
{
  //set tolerance
  double tol = 1e-8;
  fp = fopen("tol_8.txt", "w");

  /*
  ------------------------Integrator settings-------------------------------
  */
  //fcn
  double x = 0.0; //initial x value (independent var)
  double y[nrs];
  y[0] = -8.0; y[1] = -8.0; y[2] = -8.0; //initial y value (dependent var)
  double xend = 10.0; //initial x value (independent var)
  double rtoler = tol; //relative tolerance
  double atoler = tol; //absolute tolerance
  int itoler = 0; //0 -> tol is a single value
  //solout
  int iout_i = 1; //1 -> call solout
  FILE* fileout_i = NULL; // message stream
  double uround_i = 0;
  double safe = 0;
  double fac1 = 0;
  double fac2 = 0;
  double beta = 0;
  double hmax_i = 0.1; //max step size
  double h = 1e-2; //initial step size
  long nmax = 50000; //maximum number of steps
  int meth = 1; //method
  long nstiff = -1; //don't detect stiffness
  unsigned int maxbst = 50000; //maximal number of back steps
  unsigned int nrdens = 1; //nr of vars with dense output
  unsigned int icont[1]; icont[0] = 1;
  unsigned int licont = 1;
  unsigned int ngrid = 4;
  double * grid = (double*)malloc((ngrid+1)*sizeof(double));
  for (size_t i = 0; i < ngrid+1; i++)
  {
    grid[i] = (double)(i+1)*(13.0 / 28.0);
    printf("%lf\n",grid[i]);
  }

  retard(nrs,fcn,x,y,xend,&rtoler,&atoler,itoler,
    solout,iout_i,fileout_i,uround_i,safe,
    fac1,fac2,beta,hmax_i,h,nmax,meth,nstiff,
    maxbst,nrdens,icont,licont,ngrid,grid);

    printf ("x=xend  y=%12.10f\r\n", y[0]);
    printf ("rtol=%12.10f   fcn=%li   step=%li   accpt=%li   rejct=%li\r\n",
  	  tol, nfcnRead(), nstepRead(), naccptRead(), nrejctRead());

  return 0;
}
