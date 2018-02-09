#ifndef ewaldClassical_h
#define ewaldClassical_h

#include <stdio.h>
void ewaldClassicalRunner(char *dataFile, double alpha, double cutoffDirect, double cutoffFourier);
double ewaldClassical(double *r, double *q, int n, double Lx,double Ly,double Lz,double alpha, double cutoffDirect,double cutoffFourier);

#endif /* ewaldClassical_h */
