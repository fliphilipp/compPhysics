/*
func.h
*/

#ifndef _func_h
#define _func_h

extern gsl_rng* initialize_rng();
extern double dist(double [2][3]);
extern double energy(double [2][3], double);
extern double wfunc(double [2][3], double);
extern double cosAngle(double [2][3]);
extern double laplogwfunc(double, double);
extern double blockav(double *, int);
extern void dist_nuc(double [2][3], double []);
extern void autocorr(double *, int, double *);

#endif
