#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include "initfcc.h"
#include "alpotential.h"

/* 
FUNCTIONS THAT WE CAN USE:
void init_fcc(double positions[][3], int N, double lattice_param)
double splineEval(double x, const double *table,int m)
double splineEvalDiff(double x, const double *table, int m)
void get_forces_AL(double forces[][3], double positions[][3], double cell_length, int nbr_atoms)
double get_potEn_AL(double positions[][3], double cell_length, int nbr_atoms)
double get_virial_AL(double positions[][3], double cell_length, int nbr_atoms)
*/

gsl_rng* initialize_rng();

int main()
{	
	// parameters
	int ndim = 4;  // grid size in each dimension  
	int tsteps = pow(10,4);  // number of iterations
	double dt = pow(10,-4);  // time step for integration (ps)
	double mass = 27 * 1.0364 * pow(10, -4); // mass of Al [eV (ps)^2 / Å^2]

	// initialize variables
	double latparam; 
	int natoms = 4 * ndim * ndim * ndim;  // total number of atoms
	double pos[natoms][3];  // positions of atoms
	//double forces[natoms][3];  // forces on atoms
	double acc[natoms][3];  // accelerations for atoms
	double vel[natoms][3];  // velocities for atoms
	double current_time;
	double totEn;
	gsl_rng *my_rng = initialize_rng();  // random number generator
	int i, j, t;  // iterator variables

	// allocate for variables to save over time
	double *potEn = malloc((tsteps + 1) * sizeof(double));
	double *kinEn = malloc((tsteps + 1) * sizeof(double));
	double (*pos1)[6] = malloc(sizeof(double[tsteps+1][6]));

	// to figure out equilibrium spacing
	double latparamMin = 3.0;
	double latparamMax = 4.0;
	double sweepStep = 0.01;
	int nSweep = (int)((latparamMax - latparamMin)/sweepStep);
	printf("Checking from %.2f to %.2f in steps of %.4f.\n", 
		latparamMin, latparamMax, sweepStep);
	double latParamList[nSweep];
	double energyLatParam[nSweep];
	double minEnergy = pow(10,20); // something large enough

	// initialize file variables
	FILE *file_pos;
	FILE *file_latEn;
	FILE *file_energy;
	FILE *file_pos1;
	printf("Let's fucking go!\n");

	// inititalize positions without displacements for plot
	init_fcc(pos,ndim,latparamMin);
	file_pos = fopen("pos.dat","w");
	for (i = 0; i < natoms; i++) {
		fprintf(file_pos, "%.4f\t%.4f\t%.4f\n", pos[i][0], pos[i][1], pos[i][2]);
	}
	fclose(file_pos);

	// sweep though values of latparam to find equilibrium spacing (min pot energy)
	file_latEn = fopen("latticeparameter-vs-energy.dat","w");
	for (i = 0; i < nSweep; i++) {
		latParamList[i] = latparamMin + i*sweepStep;
		init_fcc(pos, ndim, latParamList[i]);
		energyLatParam[i] = get_potEn_AL(pos, latParamList[i], natoms);
		fprintf(file_latEn, "%.10f\t%.10f\n", latParamList[i], energyLatParam[i]);
		if (energyLatParam[i] < minEnergy) {
			minEnergy = energyLatParam[i];
			latparam = latParamList[i];
		}
	}
	fclose(file_latEn);
	printf("Minimum potential energy is %e at %.4f Ångström.\n", 
		minEnergy, latparam);

	latparam = 4.05;

	// initialize positions with displacements, uniform dist +/- 0.05 lattice spacing
	init_fcc(pos,ndim,latparamMin);
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			pos[i][j] += (gsl_rng_uniform(my_rng) / 10.0 - 0.05) * latparam;
		//printf("%.2f\n", (gsl_rng_uniform(my_rng) / 5.0 - 0.1)); // for checking
		}
		// printf("%.4f\t%.4f\t%.4f\n", pos[i][0], pos[i][1], pos[i][2]);
	}

	// initial accelerations
	get_acc_AL(acc, pos, latparam, mass, natoms);
	//get_forces_AL(acc, pos, latparam, natoms);

	// initial energies
	potEn[0] = get_potEn_AL(pos, latparam, natoms);
	kinEn[0] = 0.0;  // since we initialized with zero velocity

	// initial velocities (all zero)
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			vel[i][j] = 0.0;
		}
	}


	// timesteps according to velocity Verlet algorithm
	for (t = 1; t < tsteps + 1; t++) {
		// v(t+dt/2)
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				vel[i][j] += dt * 0.5 * acc[i][j];
			}
		}
		
		// q(t+dt)
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				pos[i][j] += dt * vel[i][j];
			}
		}
		
		// a(t+dt)
		get_acc_AL(acc, pos, latparam, mass, natoms);
		//get_forces_AL(acc, pos, latparam, natoms);
		
		// v(t+dt)
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				vel[i][j] += dt * 0.5 * acc[i][j];
			}
		}

		// save energies
		potEn[t] = get_potEn_AL(pos, latparam, natoms);
		kinEn[t] = get_kinEn_AL(vel, mass, natoms);

		// save position for first atom
		pos1[t][0] = pos[0][0];
		pos1[t][1] = pos[0][1];
		pos1[t][2] = pos[0][2];
		pos1[t][3] = pos[1][0];
		pos1[t][4] = pos[1][1];
		pos1[t][5] = pos[1][2];
	}

	// write energies to file
	file_energy = fopen("energy.dat","w");
	for (i = 0; i < tsteps + 1; i++) {
		current_time = i * dt;
		totEn = potEn[i] + kinEn[i];
		fprintf(file_energy, "%.4f\t%.4f\t%.4f\t%.4f\n",
			current_time, potEn[i], kinEn[i], totEn);
	}
	fclose(file_energy);

	// write position of first atom to file
	file_pos1 = fopen("pos1.dat","w");
	for (i = 0; i < tsteps + 1; i++) {
		current_time = i * dt;
		fprintf(file_pos1, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",
			current_time, pos1[i][0], pos1[i][1], pos1[i][2],
			pos1[i][3], pos1[i][4], pos1[i][5]);
	}
	fclose(file_pos1);

	// free allocated memory
	free(potEn); potEn = NULL;
	free(kinEn); kinEn = NULL;
	free(pos1); pos1 = NULL;

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
