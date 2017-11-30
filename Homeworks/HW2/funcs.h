/*
func.h
*/

#ifndef _func_h
#define _func_h

extern gsl_rng* initialize_rng();
extern double dist(double [2][3]);
extern double energy(double [][3], double);
extern double wfunc(double [][3], double);
extern void dist_nuc(double [][3], double []);
extern void autocorr(double *, int, double *);
extern void blockav(double *, int);
extern double grad(double, double);
extern double rescale(double, double [], double [], double, int, int, double);

#endif
