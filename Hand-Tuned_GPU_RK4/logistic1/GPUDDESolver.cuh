#ifndef GPU_DDE_SOLVER_
#define GPU_DDE_SOLVER_


//------------------------------ GPU Functions and Stuff --------------------------------
struct integrationSettings
{
	//delays
	double t0, dt;

	//counters
	unsigned int nrOfSteps, nrOfInitialPoints, nrOfPoints, nrOfParameters;

	//memory
	double * __restrict__ xVals, * __restrict__ xdVals, * __restrict__ tVals;
};

struct threadVariables
{
	unsigned int threadId;
	unsigned int lastIndex; //index prediction

	//integration variables:    4*3 + 6 = 18 double
	double x;
	double xTmp;
	double kAct;
	double kSum;
	double p, xDelay;
	double t, tTmp, dt;

	//memory loads
	double tb,xb,xn,xdb,xdn,pdeltat;
};

__forceinline__ __device__ void loadValues(double t, threadVariables  * vars, integrationSettings intSettings)
{
	unsigned int i = (*vars).lastIndex;
	unsigned int step;
	if(t < intSettings.tVals[i]) //on correct step
	{
		step = i-1;
	}
	else //next step is needed
	{
		step = i;
		unsigned int linIdx = (step + 1) * intSettings.nrOfParameters + (*vars).threadId;
		(*vars).lastIndex++;
		(*vars).tb = intSettings.tVals[step];
		(*vars).xb = (*vars).xn;
		(*vars).xn = intSettings.xVals[linIdx];
		(*vars).xdb = (*vars).xdn;
		(*vars).xdn = intSettings.xdVals[linIdx];
	}

}

__forceinline__ __device__ double denseOutput(double t, threadVariables  * __restrict__ vars, integrationSettings intSettings)
{
	double theta = (t - (*vars).tb) * (*vars).pdeltat;
	double thetaM1 = theta - 1;
	double res = -1. * thetaM1 * (*vars).xb + theta * ((*vars).xn + thetaM1 * (1. - 2.*theta)*((*vars).xn-(*vars).xb)+(*vars).dt*(thetaM1*(*vars).xdb + theta*(*vars).xdn));
	return res;
}

//----------------------------------- integration  -------------------------------------
__global__ void solver(integrationSettings intSettings, const double * parameters, double * endvals)
{
	//Starting thread
	threadVariables vars;

	//calculate thread id
	vars.threadId = blockDim.x*blockIdx.x + threadIdx.x;

	//read parameter from global memory
	if(vars.threadId < intSettings.nrOfParameters) vars.p = parameters[vars.threadId];
	else printf("INDEX OUT OF MEMORY");

	//initialize thread variables
	vars.dt = intSettings.dt;
	vars.pdeltat = 1.0 / vars.dt;
	vars.t = 0;
	vars.lastIndex = 0;

	vars.tb = intSettings.tVals[0];
	vars.xb = intSettings.xVals[0];
	vars.xn = intSettings.xVals[1];
	vars.xdb= intSettings.xdVals[0];
	vars.xdn= intSettings.xdVals[1];

	//read initial values
	unsigned int idx,idxLin;
	idx = intSettings.nrOfInitialPoints-1;
	idxLin = (intSettings.nrOfInitialPoints-1) * intSettings.nrOfParameters + vars.threadId;
	vars.x = intSettings.xVals[idxLin];
	idxLin += intSettings.nrOfParameters;
	idx++;

	//set initial values to save derivative from positive direction
	intSettings.tVals[idx] = vars.t;
	intSettings.xVals[idxLin] = vars.x;

	//integrate
	for(size_t stepNumber = 0; stepNumber < intSettings.nrOfSteps; stepNumber++)
	{
		//------------------------ LOAD DENSE denseOutput VALUES -----------------------------
		vars.tTmp = vars.t + 0.5*vars.dt;
		loadValues(vars.tTmp-intSettings.t0,&vars,intSettings);

		//----------------------------- START OF RK4 STEP ------------------------------------
		//k1
		vars.xDelay = denseOutput(vars.t-intSettings.t0, &vars, intSettings);
		f(&vars.kAct, vars.t, vars.x, vars.xDelay, vars.p);

		//saving the derivative
		intSettings.xdVals[idxLin] = vars.kAct;

		//k2
		vars.xDelay = denseOutput(vars.tTmp-intSettings.t0, &vars, intSettings);
		vars.kSum = vars.kAct;
		vars.xTmp = vars.x + 0.5*vars.dt*vars.kAct;
		f(&vars.kAct, vars.tTmp, vars.xTmp, vars.xDelay, vars.p);

		//k3
		vars.kSum += 2*vars.kAct;
		vars.xTmp = vars.x + 0.5*vars.dt*vars.kAct;
		f(&vars.kAct, vars.tTmp, vars.xTmp, vars.xDelay, vars.p);

		//k4
		vars.tTmp = vars.t + vars.dt;
		vars.xDelay = denseOutput(vars.tTmp-intSettings.t0, &vars, intSettings);
		vars.kSum += 2*vars.kAct;
		vars.xTmp = vars.x + vars.dt*vars.kAct;
		f(&vars.kAct, vars.tTmp, vars.xTmp, vars.xDelay, vars.p);

		//result of step
		vars.kSum += vars.kAct;
		vars.x += (1. / 6.) * vars.dt * vars.kSum;
		vars.t += vars.dt;
		//-----------------------------  END  OF RK4 STEP ------------------------------------

		//----------------------------- SAVE T AND X TO GLOBAL MEMORY ------------------------------------
		idxLin += intSettings.nrOfParameters;
		idx++;
		intSettings.tVals[idx] = vars.t;
		intSettings.xVals[idxLin] = vars.x;

		//printf("t=%6.3lf  x=%6.3lf delay=%6.3lf\n",vars.t,vars.x,vars.xDelay);
	}
	endvals[vars.threadId] = vars.x;
}

#endif
