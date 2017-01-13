#if !defined(TOOLS_H)
#define TOOLS_H

#include <stdlib.h>
#include <math.h>

typedef unsigned int boolean;
#define false 0
#define true (!false)

void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V);
void ktau_matrix(double *data, int *d, int *N, double *out);
double debyen(const double x, const int n);

#endif
