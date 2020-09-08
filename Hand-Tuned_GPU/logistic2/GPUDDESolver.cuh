#ifndef GPU_DDE_SOLVER_
#define GPU_DDE_SOLVER_


//------------------------------ GPU Functions and Stuff --------------------------------
struct integrationSettings
{
	//delays
	double t0[2], dt;

	//counters
	unsigned int nrOfSteps, nrOfInitialPoints, nrOfPoints, nrOfParameters;

	//memory
	double * __restrict__ xVals, * __restrict__ xdVals, * __restrict__ tVals;

	//mesh
	double * mesh;
	int * meshType;
	unsigned int meshLen;
};

struct threadVariables //25 double
{
	unsigned int threadId;
	unsigned int meshId;
	int stepType;
	unsigned int lastIndex[2]; //index prediction    5 double

	//integration variables:    10 double
	double x;
	double xTmp;
	double kAct;
	double kSum;
	double p, xDelay[2];
	double t, tTmp, dt, dtBase;

	//memory loads
	double tb[2],xb[2],xn[2],xdb[2],xdn[2],deltat[2]; //12 double
};


__forceinline__ __device__ unsigned int findIndex(double t, unsigned int id, threadVariables * __restrict__ vars, integrationSettings intSettings)
{
	for (unsigned int i = (*vars).lastIndex[id]; i < intSettings.nrOfPoints; i++)
	{
		if (t < intSettings.tVals[i])
		{
			if (i >= 1) (*vars).lastIndex[id] = i;
			else (*vars).lastIndex[id] = 0;
			return i;
		}
	}
	return (unsigned int)0;
}

__forceinline__ __device__ void loadValues(double t, unsigned int id, threadVariables  * __restrict__ vars, integrationSettings intSettings)
{
	unsigned int step = findIndex(t,id, vars, intSettings) - 1;
	(*vars).tb[id] = intSettings.tVals[step];
	(*vars).xb[id] = intSettings.xVals[step * intSettings.nrOfParameters + (*vars).threadId];
	(*vars).xn[id] = intSettings.xVals[(step + 1) * intSettings.nrOfParameters + (*vars).threadId];
	(*vars).xdb[id] = intSettings.xdVals[step * intSettings.nrOfParameters + (*vars).threadId];
	(*vars).xdn[id] = intSettings.xdVals[(step + 1) * intSettings.nrOfParameters + (*vars).threadId];
	(*vars).deltat[id] = intSettings.tVals[step+1] - (*vars).tb[id];
}


__forceinline__ __device__ double denseOutput(double t, unsigned int id, threadVariables  * __restrict__ vars, integrationSettings intSettings)
{
	double theta = (t - (*vars).tb[id]) / (*vars).deltat[id];
	double res = (1 - theta)*(*vars).xb[id] + theta * (*vars).xn[id] + theta * (theta - 1)*((1 - 2 * theta)*((*vars).xn[id] - (*vars).xb[id]) + (theta - 1)*(*vars).deltat[id]*(*vars).xdb[id] + theta * (*vars).deltat[id]*(*vars).xdn[id]);
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
	vars.dtBase = intSettings.dt;
	vars.t = 0;
	vars.lastIndex[0] = 0;
	vars.lastIndex[1] = 0;

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
	vars.stepType = 0;
	vars.meshId = 0;
	for(size_t stepNumber = 0; stepNumber < intSettings.nrOfSteps; stepNumber++)
	{
		//if previous step was a double mesh point
		if (vars.stepType == 2)
		{
			vars.dt = 2 * (1e-10);
		}
		else
		{
			vars.dt = vars.dtBase;
		}

		//assuming a simple step
		vars.stepType = 0;

		//detecting mesh point
		if (vars.meshId < intSettings.meshLen && intSettings.mesh[vars.meshId] < vars.t + vars.dt)
		{
			vars.stepType = intSettings.meshType[vars.meshId];
			if (vars.stepType == 1)
			{
				vars.dt = intSettings.mesh[vars.meshId] - vars.t;
				vars.meshId++;
			}

			if (vars.stepType == 2)
			{
				vars.dt = intSettings.mesh[vars.meshId] - vars.t - 1e-10;
				vars.meshId++;
			}
		}

		//------------------------ LOAD DENSE denseOutput VALUES -----------------------------
		loadValues(vars.t-intSettings.t0[0],0,&vars,intSettings);
		loadValues(vars.t-intSettings.t0[1],1,&vars,intSettings);

		//----------------------------- START OF RK4 STEP ------------------------------------
		//k1
		vars.xDelay[0] = denseOutput(vars.t-intSettings.t0[0],0, &vars, intSettings);
		vars.xDelay[1] = denseOutput(vars.t-intSettings.t0[1],1, &vars, intSettings);
		f(&vars.kAct, vars.t, vars.x, vars.xDelay[0], vars.xDelay[1], vars.p);

		//saving the derivative
		intSettings.xdVals[idxLin] = vars.kAct;

		//k2
		vars.tTmp = vars.t + 0.5*vars.dt;
		loadValues(vars.tTmp-intSettings.t0[0],0,&vars,intSettings);
		loadValues(vars.tTmp-intSettings.t0[1],1,&vars,intSettings);
		vars.xDelay[0] = denseOutput(vars.tTmp-intSettings.t0[0],0, &vars, intSettings);
		vars.xDelay[1] = denseOutput(vars.tTmp-intSettings.t0[1],1, &vars, intSettings);
		vars.kSum = vars.kAct;
		vars.xTmp = vars.x + 0.5*vars.dt*vars.kAct;
		f(&vars.kAct, vars.tTmp, vars.xTmp, vars.xDelay[0], vars.xDelay[1], vars.p);

		//k3
		vars.kSum += 2*vars.kAct;
		vars.xTmp = vars.x + 0.5*vars.dt*vars.kAct;
		f(&vars.kAct, vars.tTmp, vars.xTmp, vars.xDelay[0], vars.xDelay[1], vars.p);

		//k4
		vars.tTmp = vars.t + vars.dt;
		loadValues(vars.tTmp-intSettings.t0[0],0,&vars,intSettings);
		loadValues(vars.tTmp-intSettings.t0[1],1,&vars,intSettings);
		vars.xDelay[0] = denseOutput(vars.tTmp-intSettings.t0[0],0, &vars, intSettings);
		vars.xDelay[1] = denseOutput(vars.tTmp-intSettings.t0[1],1, &vars, intSettings);
		vars.kSum += 2*vars.kAct;
		vars.xTmp = vars.x + vars.dt*vars.kAct;
		f(&vars.kAct, vars.tTmp, vars.xTmp, vars.xDelay[0], vars.xDelay[1], vars.p);

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

		/*if(vars.threadId == 0)
			printf("t=%6.3lf  x=%6.3lf  dt=%6.3lf\n",vars.t,vars.x,vars.dt);*/
	}
	endvals[vars.threadId] = vars.x;
}

#endif
