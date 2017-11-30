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
double wfunc(double pos[2][3], double alpha) {
	double r1 = sqrt(pow(pos[0][0], 2) + pow(pos[0][1], 2) + pow(pos[0][2], 2));
	double r2 = sqrt(pow(pos[1][0], 2) + pow(pos[1][1], 2) + pow(pos[1][2], 2));
	double d = dist(pos);
	double psi = exp(-2.0 * r1) * exp(-2.0 * r2) * exp(d/(2.0 * (1.0 + alpha * d)));
	return psi;
}

// get distances to the nucleus for the two electrons
void dist_nuc(double pos[2][3], double distNuc[]) {
	distNuc[0] = sqrt(pow(pos[0][0], 2) + pow(pos[0][1], 2) + pow(pos[0][2], 2));
	distNuc[1] = sqrt(pow(pos[1][0], 2) + pow(pos[1][1], 2) + pow(pos[1][2], 2));
}

// get the cosine of the angle between the two electrons
double cosAngle(double pos[2][3]) {
	double r1norm = sqrt(pow(pos[0][0], 2) + pow(pos[0][1], 2) + pow(pos[0][2], 2));
	double r2norm = sqrt(pow(pos[1][0], 2) + pow(pos[1][1], 2) + pow(pos[1][2], 2));
	double dotProd = pos[0][0] * pos[1][0] + pos[0][1] * pos[1][1] + pos[0][2] * pos[1][2];
	double cosTheta = dotProd / (r1norm * r2norm);
	return cosTheta;
}


// get the local energy based on Eq. 7 in the assignment
double energy(double pos[2][3], double alpha) {
	int i;
	double r1norm = sqrt(pow(pos[0][0], 2) + pow(pos[0][1], 2) + pow(pos[0][2], 2));
	double r2norm = sqrt(pow(pos[1][0], 2) + pow(pos[1][1], 2) + pow(pos[1][2], 2));
	double r12 = dist(pos);
	double ar12p1 = 1.0 + alpha * r12;
	double dotProd = 0.0;
	for (i=0; i<3; i++) {
		dotProd += (pos[0][i] / r1norm - pos[1][i] / r2norm) * (pos[0][i] - pos[1][i]);
	}
	double locEn = -4.0 + dotProd / (r12 * pow(ar12p1, 2)) - 1.0 / (r12 * pow(ar12p1, 3)) - 1.0 / (4.0 * pow(ar12p1, 4)) + 1.0 / r12;
	return locEn;
}		 


// get auto-correlation for statistical inefficiency
void autocorr(double *array, int array_length, double *result) {
	int steps = 300;
	int decay_time;
	double mean_fi = 0.0;
	double square_mean_fi = 0.0;
	double correlation_func[steps];
	double fik[steps];

	// Array of zeros
	for (int k = 0; k < steps; ++k){
		fik[k] = 0.0;
	}

	// Calculate <f_i> and <f_i>^2
	for (int i = 0; i < array_length - steps; ++i){
		mean_fi += array[i] / (array_length - steps);
		square_mean_fi += pow(array[i],2) / (array_length - steps);
	}

	// Calculate <f_(i+k)*f_i> 
	for (int i = 0; i < array_length - steps; ++i){
		for (int k = 0; k < steps; ++k){
			fik[k] += (array[i + k] * array[i]) / (array_length - steps);
		} 
	}

	// Calculate correlation function
	for (int k = 0; k < steps; ++k){
		correlation_func[k] = ( fik[k] - pow(mean_fi,2) ) / ( square_mean_fi - pow(mean_fi,2)); 
	}

	// Calculate the time when the function will have decayed by approximately 10%
	decay_time = 0;
	while(correlation_func[decay_time] >= exp(-2)){
		decay_time++;
	}

	result[0] = mean_fi;
	result[1] = sqrt(decay_time * (square_mean_fi - mean_fi*mean_fi) / array_length);
	result[2] = decay_time;
}

// get block average for statistical inefficiency
void blockav(double *array, int array_length) {
	// Declaration and initiation of variables
	int i, j;
	int block_size = 500;
	int nbr_blocks = array_length/block_size;
	double block_means[nbr_blocks];
	double square_block_means[nbr_blocks];
	double var_f[nbr_blocks], total_error;
	double mean = 0.0, square_mean = 0.0, var_F, s;
	double mean_var_f = 0.0, square_mean_var_f = 0.0;

	// Determine variance for the whole array
	for(i = 0; i < array_length; i++){
		mean += array[i] / array_length;
		square_mean += pow(array[i],2) / array_length;
	}

	var_F = square_mean - mean*mean;
	
	// Determine average in each block
	for(i = 0; i < nbr_blocks; i++){
		
		block_means[i] = 0.0;
		square_block_means[i] = 0.0;
		var_f[i] = 0.0;
		
		for(j = 0; j < block_size; j++){
			block_means[i] += array[( i * block_size + j )] / block_size;
			square_block_means[i] += pow(array[( i * block_size + j )],2) / block_size;
		}
		
		var_f[i] = square_block_means[i] - pow(block_means[i],2);

		mean_var_f += var_f[i] / (block_size * nbr_blocks);
		square_mean_var_f += pow(var_f[i],2) / (block_size * nbr_blocks);
	}

	s = block_size * var_F / mean_var_f;
	printf("Statistical inefficiency s = %f\n", s);

	total_error = sqrt((square_mean_var_f - pow(mean_var_f,2)) / (block_size * s));
	printf("Result: %.4f | Error margin: %.4e \n", mean, total_error);
	printf("*******************************************\n");
}

// get gradient of ln(wfunc)
double grad(double dist, double alpha) {
	return 0.0;
}

// rescale the alpha
double rescale(double alpha, double energy_l[], double grad_ln_wave[], double distance, int iteration, int start_average, double beta) {
	return 0.0;
}

