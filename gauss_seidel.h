#ifndef __GAUSS_SEIDEL_H__
#define __GAUSS_SEIDEL_H__

#include "edo.h"

int gaussSeidel(Tridiag *sl, real_t *Y, int maxIt, real_t tol, real_t *normal2, int *numIt);

#endif // __GAUSS_SEIDEL_H__