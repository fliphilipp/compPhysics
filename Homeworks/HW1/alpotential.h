/*
alpotential.h
 
Created by Anders Lindman on 2013-03-15.
*/

#ifndef _alpotential_h
#define _alpotential_h

extern void get_forces_AL(double[][3] , double[][3], double, int);
extern double get_potEn_AL(double[][3], double, int);
extern double get_kinEn_AL(double[][3], double, int);
extern double get_virial_AL(double[][3], double, int);
extern void get_acc_AL(double[][3], double[][3], double, double, int);
extern void get_pairwise_distances(double[], double[][3], double, int);
extern void add_box(double[][3], double[][3], double, double, double, int, int *);
extern void get_distances_largebox(double[][3], double, int, int);

#endif
