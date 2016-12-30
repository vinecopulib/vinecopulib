/*
** memoryhandling.c - C code of the package CDRVine
**
** with contributions from Carlos Almeida, Aleksey Min,
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
**
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication.
**
*/

#if !defined(TOOLS_H)
#define TOOLS_H

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

typedef unsigned int boolean;
#define false 0
#define true (!false)

void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V);
void ktau_matrix(double *data, int *d, int *N, double *out);
double StableGammaDivision(double x1, double x2);

#endif
