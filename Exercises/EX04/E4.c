// EX04

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include "fft_func.h"
#define PI 3.141592653589
// Boltzmann constant from Wikipedia: 1.3806488 x 10^-23 Joule / Kelvin
#define kB (1.3806488 * pow(10,-23)) // units are right because we use both micrometers and microseconds
#define dt 0.1 // Timestep
#define steps 20000

// UNITS USED:
// temperature: Kelvin
// time: microseconds 
// distance: micrometers

// THIS IMPLIES
// kB: Joule / Kelvin
// frequency f0: Megahertz

gsl_rng* initialize_rng();

void task2(){

	// Define and initialize variables
	int nbr_trajectories = 10000; // number of trajectories to average
	int nbr_traj_to_plot = 5; // number of trajectories to plot
	double t[2];
	double x[nbr_trajectories], v[nbr_trajectories], a[nbr_trajectories];
	// random number generator
	gsl_rng *my_rng = initialize_rng();
	double rand1, rand2, G1, G2, v_temp;
	double c0, eta;
	double T = 297.0; // Kelvin
	double f_0 = 0.003; // Vibrational frequency in 1 / microseconds (i.e. MHz)
	double omega = 2 * PI * f_0; // Converting vibrational frequency to angular, for Brownian particle

	double m;
	double diameter = 2.79; // Units: [micrometers]
	double density = 2.65; // Units: [grams / cm^3]
	// Calculate volume of sphere. Units: [micrometers^3]
	double volume = (4.0 / 3.0) * PI * pow(diameter / 2.0, 3.0);
	printf("Volume: %e micrometers^3\n", volume);
	// Calculate mass of particle. Units: [grams]
	m = density * pow(10.0, -12.0) * volume; // Multiply by 10^-12 for cm^3 -> micrometers^3
	m = m * pow(10.0, -3.0); // Multiply by 10^-3 for g -> kg
	printf("Mass:   %e kg / micrometers^3\n", m);

	// Calculate v_th value
	double v_th = sqrt(kB * T / m);
	printf("v_th:   %.20f meters per second\n", v_th);

	// The two relaxtation times considered in this problem
	t[0] = 48.5;
	t[1] = 147.3;

	// when to save distributions
	int sxa[5] = {(int)(32.4/dt), (int)(231.5/dt), (int)(434.1/dt), (int)(645.3/dt), steps-1};
	int sva[5] = {(int)(0.1/dt), (int)(94.7/dt), (int)(293.6/dt), (int)(499.6/dt), steps-1};

	int sxb[9] = {(int)(38.9/dt), (int)(208.4/dt), (int)(377.9/dt), (int)(547.7/dt), (int)(716.3/dt),
	              (int)(886.2/dt), (int)(1057.7/dt), (int)(1223.5/dt), steps-1};

	int svb[9] = {(int)(0.1/dt), (int)(114.0/dt), (int)(283.1/dt), (int)(452.6/dt), (int)(620.6/dt),
	              (int)(791.6/dt), (int)(959.5/dt), (int)(1135.6/dt), steps-1};

	// for averages and standard deviations
	double x_avg, v_avg;
	double x_square, v_square;
	double x_std, v_std;

	// allocate memory for power spectrum
	double *meanposA = malloc((steps) * sizeof (double));
	double *meanvelA = malloc((steps) * sizeof (double));
	double *meanposB = malloc((steps) * sizeof (double));
	double *meanvelB = malloc((steps) * sizeof (double));
	double *pspecXA = malloc((steps) * sizeof (double));
	double *pspecVA = malloc((steps) * sizeof (double));
	double *pspecXB = malloc((steps) * sizeof (double));
	double *pspecVB = malloc((steps) * sizeof (double));
	double *freq = malloc((steps) * sizeof (double));

	// Open files for writing
	FILE *timefile;
	FILE *trajAx;
	FILE *trajAv;
	FILE *trajBx;
	FILE *trajBv;
	FILE *trajA_Average;
	FILE *trajB_Average;
	FILE *xdistA;
	FILE *vdistA;
	FILE *xdistB;
	FILE *vdistB;
	FILE *pspec;
	timefile = fopen("time.data","w");
	trajAx = fopen("trajAx.data","w");
	trajAv = fopen("trajAv.data","w");
	trajBx = fopen("trajBx.data","w");
	trajBv = fopen("trajBv.data","w");
	trajA_Average = fopen("trajA_Average.data","w");
	trajB_Average = fopen("trajB_Average.data","w");
	xdistA = fopen("xdistA.data","w");
	vdistA = fopen("vdistA.data","w");
	xdistB = fopen("xdistB.data","w");
	vdistB = fopen("vdistB.data","w");

	// Loop for the two relaxation times
	for (int i = 0; i < 2; ++i){
		// Calculate parameters
		eta = 1.0/t[i];
		c0 = exp(- eta * dt);
		printf("c0: %.20f\n", c0);

		// loop through time
		for (int step = 0; step < steps; ++step){

			if (i == 0) {fprintf(timefile, "%e\t", step * dt);}

			// reset averages across trajectories at each time step
			x_avg = 0.0;
			x_square = 0.0;
			v_avg = 0.0;
			v_square = 0.0;

			// loop through trajectories
			for (int trajectory = 0; trajectory < nbr_trajectories; ++trajectory){

				// initial conditions
				if (step == 0) {
					x[trajectory] = 0.1;
					v[trajectory] = 2.0 / 1000; // Convert milliseconds to microseconds
					a[trajectory] = - pow(omega,2.0) * x[trajectory];
				}

				// get random gaussian
				G1 = gsl_ran_gaussian (my_rng , 1);
				G2 = gsl_ran_gaussian (my_rng , 1);

				// Temporary velocity t + 1
				v_temp = (1.0 / 2.0) * a[trajectory] * dt + sqrt(c0) * v[trajectory] + v_th * sqrt(1 - c0) * G1;

				// Position x + 1
				x[trajectory] = x[trajectory] + v_temp * dt;

				// New particle acceleration
				a[trajectory] = - pow(omega,2.0) * x[trajectory];

				// Final velocity t + 1
				v[trajectory] = (1.0 / 2.0) * sqrt(c0) * a[trajectory] * dt + sqrt(c0) * v_temp + v_th * sqrt(1 - c0) * G2;

				// add to means
				x_avg += x[trajectory] / nbr_trajectories;
				v_avg += v[trajectory] / nbr_trajectories;
				x_square += pow(x[trajectory], 2.0) / nbr_trajectories;
				v_square += pow(v[trajectory], 2.0) / nbr_trajectories;

				// save positions and velocities to file (only for the first few trajectories)
				if (i == 0 && trajectory < nbr_traj_to_plot) {
					fprintf(trajAx, "%e\t", x[trajectory]);
					fprintf(trajAv, "%e\t", v[trajectory]);
					if (trajectory == nbr_traj_to_plot - 1) {
						fprintf(trajAx, "\n");
						fprintf(trajAv, "\n");
					}
				}
				if (i == 1 && trajectory < nbr_traj_to_plot) {
					fprintf(trajBx, "%e\t", x[trajectory]);
					fprintf(trajBv, "%e\t", v[trajectory]);
					if (trajectory == nbr_traj_to_plot - 1) {
						fprintf(trajBx, "\n");
						fprintf(trajBv, "\n");
					}
				}

				// save distributions at each max/min and at end
				if (i == 0 && (step == sxa[0] || step == sxa[1] || step == sxa[2] || step == sxa[3] || step == sxa[4])) {
					fprintf(xdistA, "%e\t", x[trajectory]);
				}

				if (i == 0 && (step == sva[0] || step == sva[1] || step == sva[2] || step == sva[3] || step == sva[4])) {
					fprintf(vdistA, "%e\t", v[trajectory]);
				}

				if (i == 1 && (step == sxb[0] || step == sxb[1] || step == sxb[2] || step == sxb[3] || step == sxb[4]
											  || step == sxb[5] || step == sxb[6] || step == sxb[7] || step == sxb[8])) {
					fprintf(xdistB, "%e\t", x[trajectory]);
				}

				if (i == 1 && (step == svb[0] || step == svb[1] || step == svb[2] || step == svb[3] || step == svb[4]
											  || step == svb[5] || step == svb[6] || step == svb[7] || step == svb[8])) {
					fprintf(vdistB, "%e\t", v[trajectory]);
				}

			}// end loop through trajectories

			// new line for distribution 
			if (i == 0 && (step == sxa[0] || step == sxa[1] || step == sxa[2] || step == sxa[3] || step == sxa[4])) {
				fprintf(xdistA, "\n");
				printf("Position distribution A saved at %.3f microseconds.\n", step * dt);
			}

			if (i == 0 && (step == sva[0] || step == sva[1] || step == sva[2] || step == sva[3] || step == sva[4])) {
				fprintf(vdistA, "\n");
				printf("Velocity distribution A saved at %.3f microseconds.\n", step * dt);
			}

			if (i == 1 && (step == sxb[0] || step == sxb[1] || step == sxb[2] || step == sxb[3] || step == sxb[4]
											|| step == sxb[5] || step == sxb[6] || step == sxb[7] || step == sxb[8])) {
				fprintf(xdistB, "\n");
				printf("Position distribution B saved at %.3f microseconds.\n", step * dt);
			}

			if (i == 1 && (step == svb[0] || step == svb[1] || step == svb[2] || step == svb[3] || step == svb[4]
											|| step == svb[5] || step == svb[6] || step == svb[7] || step == svb[8])) {
				fprintf(vdistB, "\n");
				printf("Velocity distribution B saved at %.3f microseconds.\n", step * dt);
			}

			// get the standard deviation after each time step
			x_std = sqrt(x_square - pow(x_avg, 2.0));
			v_std = sqrt(v_square - pow(v_avg, 2.0));

			// save average and standard deviation to file
			if (i == 0) {
				fprintf(trajA_Average, "%e \t %e \t %e \t %e\n", x_avg, x_std, v_avg, v_std);
				meanposA[step] = x_avg;
				meanvelA[step] = v_avg;
			}
			if (i == 1) {
				fprintf(trajB_Average, "%e \t %e \t %e \t %e\n", x_avg, x_std, v_avg, v_std);
				meanposB[step] = x_avg;
				meanvelB[step] = v_avg;
			}

		}// end loop through time

	}// end loop for the two relaxation times


	// make FFT (powerspectrum)
	powerspectrum(meanposA, pspecXA, steps);
	powerspectrum(meanvelA, pspecVA, steps);
	powerspectrum(meanposB, pspecXB, steps);
	powerspectrum(meanvelB, pspecVB, steps);
	powerspectrum_shift(pspecXA, steps);
	powerspectrum_shift(pspecVA, steps);
	powerspectrum_shift(pspecXB, steps);
	powerspectrum_shift(pspecVB, steps);
	fft_freq_shift(freq, dt, steps);

	pspec = fopen("pspec.data","w");
	for (int i = 0; i < steps; i++) {
		fprintf(pspec, "%e \t %e \t %e \t %e \t %e \n", freq[i], pspecXA[i], pspecVA[i], pspecXB[i], pspecVB[i]);
	}
	fclose(pspec);

	// Close data files
	fclose(timefile);
	fclose(trajAx);
	fclose(trajAv);
	fclose(trajBx);
	fclose(trajBv);
	fclose(trajA_Average);
	fclose(trajB_Average);
	fclose(xdistA);
	fclose(vdistA);
	fclose(xdistB);
	fclose(vdistB);
}


int main(){

	task2();
	return 0;
}

// helper function to initialize the rng
gsl_rng* initialize_rng() {
  const gsl_rng_type *rng_type;
  gsl_rng *my_rng ;
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  my_rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(my_rng, time(NULL));
  return my_rng;
}


