
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "funcs.h"
#define PI 3.141592653589

// Main program 
int main(){

	int task = 1;  // which task to run

	// parameters
	double alpha = 0.1;
	double delta = 1.7;
	int tSteps = pow(10,6);

	// initialize variables
	int i, t;
	double rand; 
	int randi; 
	double pos[2][3];  // positions of the two electrons
	double nextPos[2][3];  // position potential next state
	int whichParticle;  // which particle to change
	int whichDim;  // which coordinate to change
	double dist_e;  // distance of the electrons
	double psi;  // wave function value
	double distNuc[2];  // distances of both electrons from nucleus
	double prob, newProb, probRatio;  // probabilities for metropolis weight functions
	int nAccepted = 0;

	FILE *file_dist_nuc;

	// initialize random number generator
	gsl_rng *my_rng = initialize_rng(); 

	if (task == 1) {
		alpha = 0.1;
		tSteps = pow(10,6);
	}

	file_dist_nuc = fopen("file-dist-nuc.dat","w");

	// set some initial positions of the electrons
	pos[0][0] = 1.0;
	pos[0][1] = 0.0;
	pos[0][2] = 0.0;
	pos[1][0] = -1.0;
	pos[1][1] = 0.0;
	pos[1][2] = 0.0;

	// initialize next positions
	for (i=0; i<3; i++) {
		nextPos[0][i] = pos[0][i];
		nextPos[1][i] = pos[1][i];
	}

	// calculate initial trial wave function (assignment, Eq. 6)
	psi = wfunc(pos, alpha);
	printf("psi = %e \n", psi);

	// get non-normalized probability weight function for metropolis
	prob = pow(psi, 2);

	// get initial distances from origin, write to file
	dist_nuc(pos, distNuc);
	fprintf(file_dist_nuc,"%e\n%e\n", distNuc[0], distNuc[1]);

	// metropolis loop
	for (t=1; t<tSteps; t++) {

		// get potential new position by changing a random coordinate 
		randi = (int) (6.0 * gsl_rng_uniform(my_rng));
		whichParticle = randi / 3;
		whichDim = randi % 3;
		nextPos[whichParticle][whichDim] = pos[whichParticle][whichDim] + delta * (gsl_rng_uniform(my_rng) - 0.5);

		// get new wave function
		psi = wfunc(nextPos, alpha);

		// get new (non-normalized) probability weight function for metropolis
		newProb = pow(psi, 2);

		// get the ratio between the new and old probabilities
		probRatio = newProb / prob;

		// accept new state if probRatio > rand
		rand = gsl_rng_uniform(my_rng);
		if (probRatio > rand) {
			nAccepted ++;
			pos[whichParticle][whichDim] = nextPos[whichParticle][whichDim];
			prob = newProb;
		}

		// get distances from origin, write to file
		dist_nuc(pos, distNuc);
		fprintf(file_dist_nuc,"%e\n%e\n", distNuc[0], distNuc[1]);

	}// end metropolis loop

	printf("Percentage accepted: %.2f\n", 1.0 * nAccepted / tSteps);

	fclose(file_dist_nuc); 


}// end main

