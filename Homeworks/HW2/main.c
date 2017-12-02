
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "funcs.h"
#define PI 3.141592653589

// Main program 
int main(int argc, char **argv){

	// Get which task to run from the system call arguments
	char *p;
	int task = strtol(argv[1], &p, 10);
	printf("Running Task %i\n", task);

	// parameters
	double alpha;
	double thisAlpha;
	double delta = 1.7;
	int tSteps = pow(10,6);
	int endEquilibration = pow(10,3);
	double alphaMin = 0.05;
	double alphaMax = 0.025;
	double alphastep = 0.01;
	int nSimulations = 20;
	int measureEvery = 1;
	int alphaStart = 0;
	double beta;

	// Boot selected task
	switch (task) {
		case 1 :
			delta = 1.7;
			tSteps = pow(10,6);
			endEquilibration = pow(10,3);
			alphaMin = 0.1;
			alphaMax = 0.1;
			alphastep = 0.1;
			nSimulations = 1;
			measureEvery = 1;
			break;
		case 2 :
			delta = 1.7;
			tSteps = pow(10,6);
			endEquilibration = pow(10,3);
			alphaMin = 0.1;
			alphaMax = 0.1;
			alphastep = 0.1;
			nSimulations = 1;
			measureEvery = 1;
			break;
		case 3 :
			delta = 1.7;
			tSteps = pow(10,6);
			endEquilibration = pow(10,4);
			alphaMin = 0.06;
			alphaMax = 0.25;
			alphastep = 0.01;
			nSimulations = 100;
			measureEvery = 40;
			break;
		case 4 :
			delta = 1.7;
			tSteps = 5*pow(10,6);
			endEquilibration = pow(10,4);
			alphaStart = pow(10,6);
			alphaMin = 0.145;
			alphaMax = 0.145;
			alphastep = 0.1;
			nSimulations = 100;
			measureEvery = 40;
			beta = 0.95; // run for different values, gives different output files
			break;
		case 5 :
			delta = 1.7;
			tSteps = pow(10,8);
			endEquilibration = pow(10,4);
			alphaMin = 0.1486; // this is the optimized value from task 4
			alphaMax = 0.1486;
			alphastep = 0.1;
			nSimulations = 1;
			measureEvery = 40;
			break;
		default :
			printf("Please select task from 1-5 by appending it as an argument\n");
			printf("when you run the binary. Usage: ./go [1-5]\n");
			exit(0);
	}

	printf("Good choice!\n");


	// initialize variables
	int i, t;  // iterator variables
	int iSim;
	int nMeas;
	int alphaMeas;
	int optimizedAlphaMeas;
	double rand;  // random double
	int randi;  // random integer
	double pos[2][3];  // positions of the two electrons
	double nextPos[2][3];  // position potential next state
	int whichParticle;  // which particle to change
	int whichDim;  // which coordinate to change
	double psi;  // wave function value
	double distNuc[2];  // distances of both electrons from nucleus
	double prob, newProb, probRatio;  // probabilities for metropolis weight functions
	int nAccepted;  // count number of moves that are accepted
	double result[3]; // for statistical inefficiency calculation
	double avgLocalE, localEstdev, statIneff;
	double block_stat_ineff;
	double locE, meanLocalE, mean2LocalE, stdevLocalE;
	double simMeanLocalE, simStdevLocalE;
	double dist_e;
	double lapLogWFunc;
	double sumE_L, sumLapLogWFunc, theProductSum;
	double currentMeanE_L, currentMeanLapLogWFunc, theProductMean;
	double gamma;
	char filename[80];
	double alphaSum;
	double sumOptimizedAlpha;

	FILE *file_dist_nuc;
	FILE *file_cosAngle;
	FILE *file_locEn;
	FILE *file_alpha_Eloc;
	FILE *file_alphaval;
	FILE *file_optimizedalphas;

	// initialize random number generator
	gsl_rng *my_rng = initialize_rng(); 


	int nMeasurements = (tSteps - endEquilibration) / measureEvery;

	// allocate memory for large arrays
	double *localE = malloc((tSteps - endEquilibration) * sizeof(double));

	file_alpha_Eloc = fopen("file-alpha-eloc.dat","w");
	file_optimizedalphas = fopen("file-optimized-alphas.dat","w");

	// loop through multiple values of alpha
	for (alpha = alphaMin; alpha <= alphaMax; alpha += alphastep) {

		simMeanLocalE = 0.0;
		simStdevLocalE = 0.0;
		sumOptimizedAlpha = 0.0;
		optimizedAlphaMeas = 0;

		// make multiple independent runs
		for (iSim = 0; iSim < nSimulations; iSim++) {

			// we use thisAlpha when scaling it, so that it doesn't affect the loop for task 3
			thisAlpha = alpha; 

			// just here for seeing if the wfunc and energy functions work: THEY DO!
			/*
			pos[0][0] = - 1.0;
			pos[0][1] = 1.0;
			pos[0][2] = 2.0;
			pos[1][0] = - 1.0;
			pos[1][1] = -2.0;
			pos[1][2] = 0.0;
			thisAlpha = 0.1;
			printf("\n\nTEST: wave function = %.10f\n", wfunc(pos,alpha));
			printf("TEST: local energy = %.10f\n\n", energy(pos,alpha));
			*/

			if (task < 3) {
				file_dist_nuc = fopen("file-dist-nuc.dat","w");
				file_cosAngle = fopen("file-cosAngle.dat","w");
				file_locEn = fopen("file-locEn.dat","w");
			}

			if (task == 4) {
				sprintf(filename, "file-alphaval%d.dat", (int) (beta*100));
				puts(filename);
				file_alphaval = fopen(filename,"w");
				fprintf(file_alphaval, "%d\t%e\n", 0, thisAlpha);
			}

			// set some initial positions of the electrons
			pos[0][0] = 2.0 * gsl_rng_uniform(my_rng) - 1.0;
			pos[0][1] = 2.0 * gsl_rng_uniform(my_rng) - 1.0;
			pos[0][2] = 2.0 * gsl_rng_uniform(my_rng) - 1.0;
			pos[1][0] = - 2.0 * gsl_rng_uniform(my_rng) - 1.0;
			pos[1][1] = - 2.0 * gsl_rng_uniform(my_rng) - 1.0;
			pos[1][2] = - 2.0 * gsl_rng_uniform(my_rng) - 1.0;

			// initialize possible next positions
			for (i=0; i<3; i++) {
				nextPos[0][i] = pos[0][i];
				nextPos[1][i] = pos[1][i];
			}

			// calculate initial trial wave function (assignment, Eq. 6)
			psi = wfunc(pos, thisAlpha);

			// get non-normalized probability weight function for metropolis
			prob = pow(psi, 2);

			// set means to zero before loop
			meanLocalE = 0.0;
			mean2LocalE = 0.0;

			nMeas = 1;
			alphaMeas = 0;
			alphaSum = 0.0;

			sumE_L = 0.0;
			sumLapLogWFunc = 0.0;
			theProductSum = 0.0;

			nAccepted = 0;


			// metropolis loop
			for (t=1; t<tSteps; t++) {

				// get potential new position by changing a random coordinate 
				randi = (int) (6.0 * gsl_rng_uniform(my_rng));
				whichParticle = randi / 3;
				whichDim = randi % 3;
				nextPos[whichParticle][whichDim] = pos[whichParticle][whichDim] + delta * (gsl_rng_uniform(my_rng) - 0.5);

				// get new wave function
				psi = wfunc(nextPos, thisAlpha);

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

				// start sampling data after equilibrating the system
				if (t >= endEquilibration && t % measureEvery == 0) {

					if (task < 3) {
						// get distances from origin, write to file
						dist_nuc(pos, distNuc);
						fprintf(file_dist_nuc, "%e\n%e\n", distNuc[0], distNuc[1]);

						// get the cosine of the angle between the electrons and write to file
						fprintf(file_cosAngle, "%e\n", cosAngle(pos));
					}

					// get the local energy
					locE = energy(pos, thisAlpha);
					meanLocalE += locE / nMeasurements;
					mean2LocalE += pow(locE, 2) / nMeasurements;

					if (task == 2) {
						localE[nMeas] = locE;
					}

					// for rescaling alpha
					if (task == 4) {

						// get the laplacian of the logarithmized wave function
						dist_e = dist(pos);
						lapLogWFunc = laplogwfunc(dist_e, thisAlpha);

						// get current averages
						sumE_L += locE;
						currentMeanE_L = sumE_L / nMeas;
						sumLapLogWFunc += lapLogWFunc;
						currentMeanLapLogWFunc = sumLapLogWFunc / nMeas;
						theProductSum += locE * lapLogWFunc;
						theProductMean = theProductSum / nMeas;

						// calculate scaling factor gamma
						gamma = pow(nMeas, -beta);

						// rescale alpha
						thisAlpha = thisAlpha - gamma * 2.0 * (theProductMean - currentMeanE_L * currentMeanLapLogWFunc);
						fprintf(file_alphaval, "%d\t%e\n", nMeas*measureEvery, thisAlpha);

						// average alpha results
						if (t > alphaStart) {
							alphaSum += thisAlpha;
							alphaMeas++;
						}
					}

					nMeas++;

				}

				// save local energy to file 
				if (task < 3) {
					fprintf(file_locEn, "%e\n", energy(pos, thisAlpha));
				}

			}// end metropolis loop

			stdevLocalE = sqrt((mean2LocalE - pow(meanLocalE, 2)) / nMeasurements);

			if (task != 3) {
				printf("\nRun: %d\n", iSim);
				printf("Percentage accepted: %.2f%%\n", 100.0 * nAccepted / tSteps);
				printf("Average local energy: %.5f +/- %.5f\n", meanLocalE, stdevLocalE);
			}

			// statistical inefficiency from autocorrelation and block averaging
			if (task == 2) {
				autocorr(localE, tSteps - endEquilibration, result);
				avgLocalE = result[0];
				localEstdev = result[1];
				statIneff = result[2];
				printf("Statistical inefficiency from correlation: %d\n", (int) result[2]);

				block_stat_ineff = blockav(localE, tSteps - endEquilibration);
				printf("Statistical inefficiency from block averaging: %.2f\n", block_stat_ineff);
			}


			simMeanLocalE += meanLocalE / nSimulations;
			simStdevLocalE += stdevLocalE / nSimulations / sqrt(nSimulations);

			// print results for each alpha, and for each simulation
			// if (task == 3) {
			// 	printf("alpha = %.3f,   sim:%2d,   E_loc = %.3f +/- %.5f\n", thisAlpha, iSim, meanLocalE, stdevLocalE);
			// }

			// close data files
			if (task < 3) {
				fclose(file_dist_nuc); 
				fclose(file_cosAngle); 
				fclose(file_locEn); 
			}

			if (task == 4){
				fclose(file_alphaval);
				printf("Value of beta used for scaling decay: %.3f\n", beta);
				printf("Optimized value of alpha: %.10f\n", alphaSum / alphaMeas);
				fprintf(file_optimizedalphas,"%e\n",alphaSum / alphaMeas);

				// count only count if something sensible comes out (sometimes it goes nuts...)
				if (alphaSum / alphaMeas > 0.13 && alphaSum / alphaMeas < 0.16) {
					sumOptimizedAlpha += alphaSum / alphaMeas;
					optimizedAlphaMeas ++;
				}
			}
		}// end multiple independent runs

		if (task == 3) {
			printf("alpha = %.3f,   average   E_loc = %.3f +/- %.5f\n\n", thisAlpha, simMeanLocalE, simStdevLocalE);
		}

		if (task == 4) {
			printf("Mean optimized alpha: %.10f\n", sumOptimizedAlpha / optimizedAlphaMeas);
		}

		// write
		fprintf(file_alpha_Eloc, "%e\t%e\t%e\t%e\n", thisAlpha, simMeanLocalE, simMeanLocalE - simStdevLocalE, simMeanLocalE + simStdevLocalE);

	}// end loop over alpha values

	// free memory
	free(localE); localE = NULL;

	printf("\n");
	return 0; // Indicator that everything went well.
}// end main

