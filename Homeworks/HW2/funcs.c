// functions to use in main

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#define PI 3.141592653589

// random number generator initialization
gsl_rng* initialize_rng() {
	const gsl_rng_type *rng_type;
	gsl_rng *my_rng ;
	gsl_rng_env_setup();
	rng_type = gsl_rng_default;
	my_rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(my_rng, time(NULL));
	return my_rng;
}

// helper function that calculates distance between two particles in 3D
double dist(double pos[2][3]){
	double sqDist = 0.0;
	int i;
	for(i = 0; i < 3; i++){
		sqDist += pow(pos[0][i] - pos[1][i], 2);
	}
	return sqrt(sqDist);
}

// get the value of the trial wave function specified in Eq. 6 in the assignment
double wfunc(double pos[][3], double alpha){
	double r1 = sqrt(pow(pos[0][0], 2) + pow(pos[0][1], 2) + pow(pos[0][2], 2));
	double r2 = sqrt(pow(pos[1][0], 2) + pow(pos[1][1], 2) + pow(pos[1][2], 2));
	double d = dist(pos);
	double psi = exp(-2.0 * r1) * exp(-2.0 * r2) * exp(d/(2.0 * (1.0 + alpha * d)));
	return psi;
}


// get distances to the nucleus for the two electrons
void dist_nuc(double pos[][3], double distNuc[]){
	distNuc[0] = sqrt(pow(pos[0][0], 2) + pow(pos[0][1], 2) + pow(pos[0][2], 2));
	distNuc[1] = sqrt(pow(pos[1][0], 2) + pow(pos[1][1], 2) + pow(pos[1][2], 2));
}


// get the local energy
double energy(double pos[][3], double alpha){

}  		 


// get auto-correlation for statistical inefficiency
void autocorr(double *A, int len, double *result){

}


// get block average for statistical inefficiency
void blockav(double *A, int len){

}

// get gradient of ln(wfunc)
double grad(double dist, double alpha){

}

// rescale the alpha
double rescale(double alpha, double energy_l[], double grad_ln_wave[], double distance, int iteration, int start_average, double beta){

}

