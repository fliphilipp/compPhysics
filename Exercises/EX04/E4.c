// EX03

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#define PI 3.141592653589
// Boltzmann constant from Wikipedia
// Units: [Joule / Kelvin]
#define kB (1.3806488 * pow(10,-23))


gsl_rng* initialize_rng();

void task2(){

	// Define and initialize variables
	int steps = pow(10,5);
	int nbr_trajectories = 5;
	double dt = 0.01; // Timestep
	double t[2];
	double x[steps][nbr_trajectories], v[steps][nbr_trajectories];
	int f_0 = 3; // Vibrational frequency for the Brownian particle
	// random number generator
	gsl_rng *my_rng = initialize_rng();
	double rand1, rand2, G1, G2, v_temp;
	double c0, eta, a;
	double T = 297.0; // Kelvin
	double freq = 3.0; // Vibrational frequency for Brownian particle
	

	double m;
	double diameter = 2.79; // Units: [micrometers]
	double density = 2.65; // Units: [grams / cm^3]
	// Calculate volume of sphere. Units: [micrometers^3]
	double volume = (4.0 / 3.0) * PI * pow(diameter / 2.0, 3.0);
	// Calculate mass of particle. Units: [grams]
	m = density * pow(10.0, -4.0) * volume; // Multiply by 10^-4 for um -> cm
	m = m * pow(10.0, -4.0); // To convert from grams to kilos
	printf("Particle mass: %e\n", m);

	// Calculate v_th value
	double v_th = sqrt(kB * T / m);
	printf("v_th: %e\n", v_th);

	// Set initial conditions for the trajectories
	for (int trajectory = 0; trajectory < nbr_trajectories; ++trajectory){
		x[0][trajectory] = 0.1;
		v[0][trajectory] = 2.0;	
	}

	// The two relaxtation times considered in this problem
	t[0] = 48.5;
	t[1] = 147.3; 

	// Open files for writing
	FILE *trajAx;
	FILE *trajAv;
	FILE *trajBx;
	FILE *trajBv;
	trajAx = fopen("trajAx.data","w");
	trajAv = fopen("trajAv.data","w");
	trajBx = fopen("trajBx.data","w");
	trajBv = fopen("trajBv.data","w");

	// Save initial values
	fprintf(trajAx, "%e \t %e \t %e \t %e \t %e\n", x[0][0], x[0][1], x[0][2], x[0][3], x[0][4]);
	fprintf(trajAv, "%e \t %e \t %e \t %e \t %e\n", v[0][0], v[0][1], v[0][2], v[0][3], v[0][4]);
	fprintf(trajBx, "%e \t %e \t %e \t %e \t %e\n", x[0][0], x[0][1], x[0][2], x[0][3], x[0][4]);
	fprintf(trajBv, "%e \t %e \t %e \t %e \t %e\n", v[0][0], v[0][1], v[0][2], v[0][3], v[0][4]);

	// Loop for the two relaxation times
	for (int i = 0; i < 2; ++i){
		// Calculate parameters
		eta = 1.0/t[i];
		c0 = exp(- eta * dt);
		printf("c0: %e\n", c0);

		for (int trajectory = 0; trajectory < nbr_trajectories; ++trajectory){

			for (int step = 1; step < steps; ++step){

				// Get random samples
		    rand1 = gsl_rng_uniform(my_rng);
		    rand2 = gsl_rng_uniform(my_rng);

		    // Use Bow-Muller transform to get values from normal distribution
		    // Recommended on: https://stackoverflow.com/questions/2325472/generate-random-numbers-following-a-normal-distribution-in-c-c
		    G1 = sqrt(-2.0 * log(rand1)) * cos(2.0 * PI * rand1);
		    G2 = sqrt(-2.0 * log(rand2)) * cos(2.0 * PI * rand2);

		    // Particle acceleration 
		    // https://en.wikipedia.org/wiki/Particle_acceleration
		    a = - pow(freq,2.0) * x[step - 1][trajectory];

		    // Temporary velocity t + 1
		    v_temp = (1.0 / 2.0) * a * dt + sqrt(c0) * v[step - 1][trajectory] \
		    											+ v_th * sqrt(1 - c0) * G1;

		    // Position x + 1
        x[step][trajectory] = x[step - 1][trajectory] + v_temp * dt;

        // New particle acceleration
        a = - pow(freq,2.0) * x[step][trajectory];

        // Final velocity t + 1
        v[step][trajectory] = (1.0 / 2.0) * sqrt(c0) * a * dt + sqrt(c0) * v_temp \
        											+ v_th * sqrt(1 - c0) * G2;

			}
			// printf("a: %e v_temp: %e x: %e\n", a, v_temp, x[steps][trajectory]);

		}

		// Save A
		switch(i) {
			case 0 :
				for (int step = 0; step < steps; ++step){
					fprintf(trajAx, "%e \t %e \t %e \t %e \t %e\n", x[step][0], x[step][1], x[step][2], x[step][3], x[step][4]);
					fprintf(trajAv, "%e \t %e \t %e \t %e \t %e\n", v[step][0], v[step][1], v[step][2], v[step][3], v[step][4]);
				}
				break;
			case 1 :
				for (int step = 0; step < steps; ++step){
					fprintf(trajBx, "%e \t %e \t %e \t %e \t %e\n", x[step][0], x[step][1], x[step][2], x[step][3], x[step][4]);
					fprintf(trajBv, "%e \t %e \t %e \t %e \t %e\n", v[step][0], v[step][1], v[step][2], v[step][3], v[step][4]);
				}
				break;
			default :
				printf("This shouldn't happen!\n");
		}

	}

	// Close
	fclose(trajAx);
	fclose(trajAv);
	fclose(trajBx);
	fclose(trajBv);
}

void task3(){

}

int main(int argc, char const *argv[]){

	// Run task 2
	task2();

	// Run task 3
	task3();

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
