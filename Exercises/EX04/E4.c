// EX04

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_randist.h>
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

	// for averages and standard deviations
	double x_avg, v_avg;
	double x_square, v_square;
	double x_std, v_std;

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

				// save distributions after each cycle of the harmonic potential well
				if (step % (int)(1.0 / (f_0 * dt)) == 0) {
					if (i == 0) {
						fprintf(xdistA, "%e\t", x[trajectory]);
						fprintf(vdistA, "%e\t", v[trajectory]);
					}
					if (i == 1) {
						fprintf(xdistB, "%e\t", x[trajectory]);
						fprintf(vdistB, "%e\t", v[trajectory]);
					}
					if(trajectory == 0) {
						printf("Distribution saved at %.3f microseconds.\n", step * dt);
					}
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

			}// end loop through trajectories

			// new line for distribution 
			if (step % (1 + (int)(1.0 / (f_0 * dt))) == 0) {
				if (i == 0) {
					fprintf(xdistA, "\n");
					fprintf(vdistA, "\n");
				}
				if (i == 1) {
					fprintf(xdistB, "\n");
					fprintf(vdistB, "\n");
				}
			}

			// get the standard deviation after each time step
			x_std = sqrt(x_square - pow(x_avg, 2.0));
			v_std = sqrt(v_square - pow(v_avg, 2.0));

			// save average and standard deviation to file
			if (i == 0) {
				fprintf(trajA_Average, "%e \t %e \t %e \t %e\n", x_avg, x_std, v_avg, v_std);
			}
			if (i == 1) {
				fprintf(trajB_Average, "%e \t %e \t %e \t %e\n", x_avg, x_std, v_avg, v_std);
			}

		}// end loop through time

	}// end loop for the two relaxation times

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

	/*
	/////////////////////////////////////////////////////////////
	// Run for multiple trajectories
	// Basically a copy/paste of the above code, but not really..
	/////////////////////////////////////////////////////////////
	double x_avg[steps], v_avg[steps];
	double x_square[steps], v_square[steps];
	double x_variance, v_variance;
	double x_mul[steps], v_mul[steps]; // Use one dimensional arrays for this task for x and v

	x_mul[0] = 0.1;
	v_mul[0] = 2.0 / 1000;	// Convert milliseconds to microseconds
	x_avg[0] = x_mul[0];
	v_avg[0] = v_mul[0];

	// Open files for writing
	FILE *trajA_Average;
	FILE *trajB_Average;
	trajA_Average = fopen("trajA_Average.data","w");
	trajB_Average = fopen("trajB_Average.data","w");

	// Loop for the two relaxation times
	for (int i = 0; i < 2; ++i){
		// Calculate parameters
		eta = 1.0/t[i];
		c0 = exp(- eta * dt);

		for (int trajectory = 0; trajectory < nbr_trajectories_to_avg; ++trajectory){
				
				// Particle acceleration 
			// https://en.wikipedia.org/wiki/Particle_acceleration
			a = - pow(omega,2.0) * x_mul[0];
			
			for (int step = 1; step < steps; ++step){

				// Get random samples
			rand1 = gsl_rng_uniform(my_rng);

			// Use Bow-Muller transform to get values from normal distribution
			// Recommended on: https://stackoverflow.com/questions/2325472/generate-random-numbers-following-a-normal-distribution-in-c-c
			G1 = sqrt(-2.0 * log(rand1)) * cos(2.0 * PI * rand1);

			// Temporary velocity t + 1
			v_temp = (1.0 / 2.0) * a * dt + sqrt(c0) * v_mul[step - 1] \
														+ v_th * sqrt(1 - c0) * G1;

			// Position x + 1
		x_mul[step] = x_mul[step - 1] + v_temp * dt;

		// New particle acceleration
		a = - pow(omega,2.0) * x_mul[step];

		// Get new random sample
		rand2 = gsl_rng_uniform(my_rng);
		G2 = sqrt(-2.0 * log(rand2)) * cos(2.0 * PI * rand2);

		// Final velocity t + 1
		v_mul[step] = (1.0 / 2.0) * sqrt(c0) * a * dt + sqrt(c0) * v_temp \
													+ v_th * sqrt(1 - c0) * G2;

		// Get averages and variance
		x_avg += x_mul[step] / nbr_trajectories_to_avg;
		v_avg += v_mul[step] / nbr_trajectories_to_avg;
		x_square[step] += pow(x_avg, 2.0) / nbr_trajectories_to_avg;
		v_square[step] += pow(v_avg, 2.0) / nbr_trajectories_to_avg;
			}
		}
			
		// Export values from averages
		for (int step = 0; step < steps; ++step){
			// Compute variance for x and v
			x_variance = sqrt(x_avg - x_square[step]);
			v_variance = sqrt(v_avg - v_square[step]);
			if( i == 0 ){
				fprintf(trajA_Average, "%.20f \t %.20f \t %.20f \t %.20f\n", x_avg, x_variance, v_avg, v_variance);
			} else if ( i == 1 ){
				fprintf(trajB_Average, "%.20f \t %.20f \t %.20f \t %.20f\n", x_avg, x_variance, v_avg, v_variance);
			}
		}
	}
	// Save & close
	fclose(trajA_Average);
	fclose(trajB_Average);
	*/

}

/*
void task3(int time_i){

	// Define and initialize variables
	int nbr_trajectories = pow(10,3); // trajectories to bin
	double t[2];
	double x, v;
	double X[nbr_trajectories][2], V[nbr_trajectories][2];
	// random number generator
	gsl_rng *my_rng = initialize_rng();
	double rand1, rand2, G1, G2, v_temp;
	double c0, eta, a;
	double T = 297.0; // Kelvin
	double f_0 = 3.0; // Vibrational frequency for the Brownian particle
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

	x = 0.1;
	v = 2.0 / 1000;	// Convert milliseconds to microseconds

	// Open files for writing
	FILE *task3;
	task3 = fopen("task3.data","w");

	// Loop for the two relaxation times
	for (int i = 0; i < 2; ++i){
		// Calculate parameters
		eta = 1.0/t[i];
		c0 = exp(- eta * dt);

		for (int trajectory = 0; trajectory < nbr_trajectories; ++trajectory){
				
				// Particle acceleration 
			// https://en.wikipedia.org/wiki/Particle_acceleration
			a = - pow(omega,2.0) * x;
			
			for (int step = 1; step < steps; ++step){

				// Get random samples
			rand1 = gsl_rng_uniform(my_rng);

			// Use Bow-Muller transform to get values from normal distribution
			// Recommended on: https://stackoverflow.com/questions/2325472/generate-random-numbers-following-a-normal-distribution-in-c-c
			G1 = sqrt(-2.0 * log(rand1)) * cos(2.0 * PI * rand1);

			// Temporary velocity t + 1
			v_temp = (1.0 / 2.0) * a * dt + sqrt(c0) * v \
														+ v_th * sqrt(1 - c0) * G1;

			// Position x + 1
		x = x + v_temp * dt;

		// New particle acceleration
		a = - pow(omega,2.0) * x;

		// Get new random sample
		rand2 = gsl_rng_uniform(my_rng);
		G2 = sqrt(-2.0 * log(rand2)) * cos(2.0 * PI * rand2);

		// Final velocity t + 1
		v = (1.0 / 2.0) * sqrt(c0) * a * dt + sqrt(c0) * v_temp \
													+ v_th * sqrt(1 - c0) * G2;


		if(step == (time_i * dt)){
			X[trajectory][i] = x;
			V[trajectory][i] = v;
		 	break;
		}
			}
		}		
	}

	// Save & close
	for (int trajectory = 0; trajectory < nbr_trajectories; ++trajectory){
		fprintf(task3, "%.20f \t %.20f \t %.20f \t %.20f\n", X[trajectory][0], V[trajectory][0], X[trajectory][1], V[trajectory][1]);
	}
	fclose(task3);

}
*/

// int main(int argc, char const *argv[]){

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


