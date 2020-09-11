#include <math.h>
#include <stdio.h>
#include "retard.h"

#define nrs 1 //dimension of the system

/* initial function
i: component index
x: independent variable
*/

FILE* fp;

double phi (unsigned i, double x)
{
  if(x <= -2) return 0;
  if(x <= -1) return 1;
  return 0;
}

//definition of DDE
void fcn (unsigned n, double x, double* y, double* f)
{
  double xdelay = ylag(0, x-3.0, phi); //(component, point)
  f[0] = -xdelay;
}

//save output with this method
void solout (long nr, double xold, double x, double* y, unsigned n, int* irtrn)
{
  fprintf (fp,"%17.20lf\t%17.20lf\n", x, y[0]);
}

int main(int argc, char const *argv[])
{
  //set tolerance
  double tol=1e-10;
  fp = fopen("tol_10.txt", "w");

  /*
  ------------------------Integrator settings-------------------------------
  */
  //fcn
  double x = 0.0; //initial x value (independent var)
  double y[nrs];
  y[0] = 0.0; //initial y value (dependent var)
  double xend = 25.0; //initial x value (independent var)
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
  double hmax_i = 1.0; //max step size
  double h = 1e-2; //initial step size
  long nmax = 10000; //maximum number of steps
  int meth = 1; //method
  long nstiff = -1; //don't detect stiffness
  unsigned int maxbst = 1000; //maximal number of back steps
  unsigned int nrdens = 1; //nr of vars with dense output
  unsigned int icont[1]; icont[0] = 0;
  unsigned int licont = 1;
  unsigned int ngrid = 15;
  double grid[] = { 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0 };

  retard(nrs,fcn,x,y,xend,&rtoler,&atoler,itoler,
    solout,iout_i,fileout_i,uround_i,safe,
    fac1,fac2,beta,hmax_i,h,nmax,meth,nstiff,
    maxbst,nrdens,icont,licont,ngrid,grid);

    printf ("x=xend  y=%12.10f\r\n", y[0]);
    printf ("rtol=%12.10f   fcn=%li   step=%li   accpt=%li   rejct=%li\r\n",
  	  tol, nfcnRead(), nstepRead(), naccptRead(), nrejctRead());

  return 0;
}
