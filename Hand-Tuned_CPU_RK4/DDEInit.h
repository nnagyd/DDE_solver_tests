#ifndef DDE_INIT_H
#define DDE_INIT_H

#include <chrono>


double * discretize(double f(double t, int dir), double * tVals, unsigned int nrOfPoints){
	double * lst = new double[nrOfPoints];

	for (size_t i = 0; i < nrOfPoints; i++)
	{
		int dir = -1;
		if (i >= 1 && tVals[i] == tVals[i - 1])
		{
			dir = 1;
		}
		lst[i] = f(tVals[i], dir);
	}
	return lst;
}

double * linspace(double t0, double t1, unsigned int nr){
	double * lst = new double[nr];
	double dt = (t1 - t0) / double(nr - 1);
	double t = t0;
	for (size_t i = 0; i < nr; i++)
	{
		lst[i] = t;
		t += dt;
	}
	return lst;
}

double* linspaceAlligned(double t0, double t1, unsigned int nr)
{
	double* lst = (double*)_aligned_malloc(nr * sizeof(double),32);
	double dt = (t1 - t0) / double(nr - 1);
	double t = t0;
	for (size_t i = 0; i < nr; i++)
	{
		lst[i] = t;
		t += dt;
	}
	return lst;
}

double * linspaceDisc(double t0, double t1, unsigned int nr, double * tDisc = NULL, unsigned int nrOfDisc = 0, double eps = 0.0){
	double * lst = new double[nr + 2 * nrOfDisc];
	int * discMask = new int[nrOfDisc];
	//set all element to 0, set to 1 if the i. discontinouity is included
	for (size_t i = 0; i < nrOfDisc; i++)
	{
		discMask[i] = 0;
	}

	double dt = (t1 - t0) / (nr - 1);
	double t = t0;
	for (size_t i = 0; i < nr + 2 * nrOfDisc; i++)
	{
		bool set = true;
		for (size_t j = 0; j < nrOfDisc; j++)
		{

			if (!discMask[j] && fabs(tDisc[j] - t) < eps) //discontinuity happens at a point
			{
				lst[i] = tDisc[j] - eps;
				lst[i + 1] = tDisc[j];
				lst[i + 2] = tDisc[j] + eps;
				set = false;
				discMask[j] = 1;
				i += 2;
			}
			else if (!discMask[j] && tDisc[j] < t) //discontinuity far from point
			{
				lst[i] = tDisc[j] - eps;
				lst[i + 1] = tDisc[j] + eps;
				discMask[j] = 1;
				i += 2;
			}
		}
		if (set) lst[i] = t;
		t += dt;
	}
	return lst;
}

double* filter(double* original, unsigned int* nr, double min, double* toRemove = NULL, unsigned int nrToRemove = 0, double tol = 1e-10)
{
	unsigned int count = *nr;
	double* unique = new double[count];
	unsigned int uniqueNr = 0;
	for (size_t i = 0; i < count; i++)
	{
		bool set = true;
		for (size_t j = 0; j < uniqueNr; j++)
		{
			if (abs(original[i] - unique[j]) < tol) //already in new list
			{
				set = false;
			}
		}

		for (size_t j = 0; j < nrToRemove; j++)
		{
			if (abs(original[i] - toRemove[j]) < tol) //already in concurrent list
			{
				set = false;
			}
		}

		if (set && original[i] > min) //original[i] not in new list yet
		{
			unique[uniqueNr] = original[i];
			uniqueNr++;
		}
	}

	double* filtered = new double[uniqueNr];
	for (size_t i = 0; i < uniqueNr; i++)
	{
		filtered[i] = unique[i];
	}
	*nr = uniqueNr;
	return filtered;
}

void sort(double* lst, unsigned int len)
{
	for (size_t i = 1; i < len; i++)
	{
		for (size_t j = 1; j < len; j++)
		{
			if (lst[j] < lst[j - 1]) //swap
			{
				double tmp = lst[j];
				lst[j] = lst[j - 1];
				lst[j - 1] = tmp;
			}
		}
	}
}

//------------------------------ Timer stuff -------------------------------------------
uint64_t micros()
{
	uint64_t us = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::
		now().time_since_epoch()).count();
	return us;
}

double seconds()
{
	return double(micros()) / 1000. / 1000.;
}

#endif
