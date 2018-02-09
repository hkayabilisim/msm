#ifndef ewaldHernquist_h
#define ewaldHernquist_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "msmLibrary.h"

void ewaldHernQuistRunner(char *dataFile, double alpha, double cutoffDirect, double cutoffFourier);
double ewaldHernquist(double *r, double *q, int n, double L,double alpha, double cutoffDirect,double cutoffFourier);
#endif /* ewaldHernquist_h */
