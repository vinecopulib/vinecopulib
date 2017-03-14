#if !defined(TOOLS_H)
#define TOOLS_H

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V);
void ktau_matrix(double *data, int *d, int *N, double *out);
double debyen(const double x, const int n);

#endif
