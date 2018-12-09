#ifndef FUNCS_CPP
#define FUNCS_CPP

#include <cmath>
double polyfit(double *y, int n)
{
	double sx, sy, sxx, sxy;
	sx = sy = sxx = sxy = 0.0;
	for (int i = 0; i < n; ++i)
	{
		sx += i;
		sy += y[i];
		sxx += (i * i);
		sxy += (i * y[i]);
	}
	sx /= n;
	sy /= n;
	sxx /= n;
	sxy /= n;
	return (sxy - sx * sy) / (sxx - sx * sx);

}
double corr_length(double *corr, int corr_size, double stop_point = exp(-8.0))
{
	double beg = corr[0];
	int n = 0;
	while (n < corr_size && corr[n] > beg * stop_point)
		++n;
	double *corr_local = new double[n];
	for (int i = 0; i < n; ++i)
		corr_local[i] = log(corr[i]);
	return -1.0 / polyfit(corr_local, n);
}

#endif