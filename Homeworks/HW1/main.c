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
	int tsteps = pow(10,5);  // number of iterations
	double dt = pow(10,-2);  // time step for integration (ps)
	int startEquilibration = pow(10,2);  // when to start equilibrating
	int endEquilibration = pow(10,5);  // for how long to equilibrate
	int startAveraging = pow(10,3);  // when to start taking time averages
	double T_eq = 500.0 + 273.15;  // Equilibrium temperature in Kelvin (corresponds to 500 deg C)
	double mass = 27 * 1.0364 * pow(10, -4); // mass of Al [eV (ps)^2 / Å^2]
	double boltz = 8.617 * pow(10,-5);  // Boltzmann constant [eV / K]
	double P_eq = 101.325;  // Equilibrium pressure (101.325 kPa at sea level)
	double rand_disp = 0.05;  // initial dispacements percentage of latparam
	double tauT = pow(10,-1);  // decay time constant for temperature scaling 
	double kappaT_tauP_ratio = pow(10,-10);  // only need to choose appropriate ratio for pressure scaling

	// initialize variables
	double latparam; 
	double cellLength;
	int natoms = 4 * ndim * ndim * ndim;  // total number of atoms
	double pos[natoms][3];  // positions of atoms
	double vel[natoms][3];  // velocities for atoms
	double acc[natoms][3];  // accelerations for atoms
	double current_time;
	double totEn;
	double vsquare;
	double alphaT;
	double alphaP;
	double sumT;
	double sumP;
	gsl_rng *my_rng = initialize_rng();  // random number generator
	int i, j, t;  // iterator variables
	if (tsteps < endEquilibration) {endEquilibration = tsteps;}
	if (startEquilibration > tsteps) {startEquilibration = tsteps;}

	// allocate for variables to save over time
	double *potEn = malloc((tsteps + 1) * sizeof(double));
	double *kinEn = malloc((tsteps + 1) * sizeof(double));
	double *temp = malloc((tsteps + 1) * sizeof(double));
	double *pres = malloc((tsteps + 1) * sizeof(double));
	double (*pos1)[6] = malloc(sizeof(double[tsteps+1][6]));

	// to figure out equilibrium spacing
	double latparamMin = 4.032;
	double latparamMax = 4.033;
	double sweepStep = 0.0001;
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
		energyLatParam[i] = get_potEn_AL(pos, latParamList[i] * ndim, natoms);
		fprintf(file_latEn, "%.10f\t%.10f\n", latParamList[i], energyLatParam[i]);
		if (energyLatParam[i] < minEnergy) {
			minEnergy = energyLatParam[i];
			latparam = latParamList[i];
		}
	}
	fclose(file_latEn);
	printf("Minimum potential energy is %e at %.4f Ångström.\n", 
		minEnergy, latparam);

	// initialize positions with displacements, uniform dist +/- 0.05 lattice spacing
	cellLength = latparam * ndim;
	init_fcc(pos,ndim,latparam);
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			pos[i][j] += (gsl_rng_uniform(my_rng) * rand_disp * 2 - rand_disp) * latparam;
		}
	}
	// save initial position for first and second atoms
	pos1[0][0] = pos[0][0];
	pos1[0][1] = pos[0][1];
	pos1[0][2] = pos[0][2];
	pos1[0][3] = pos[1][0];
	pos1[0][4] = pos[1][1];
	pos1[0][5] = pos[1][2];

	// initial accelerations
	get_acc_AL(acc, pos, cellLength, mass, natoms);
	//get_forces_AL(acc, pos, latparam, natoms);

	// initial energies
	potEn[0] = get_potEn_AL(pos, cellLength, natoms);
	kinEn[0] = 0.0;  // since we initialized with zero velocity
	temp[0] = 0.0;  // initial temperature (since velocities all zero)
	pres[0] = 0.0;  // initial pressure (since temp is zero)

	// initial velocities (all zero)
	for (i = 0; i < natoms; i++) {
		for (j = 0; j < 3; j++) {
			vel[i][j] = 0.0;
		}
	}
	sumT = 0.0;
	sumP = 0.0;

	// timesteps according to velocity Verlet algorithm
	for (t = 1; t < tsteps + 1; t++) {
		// print progress
		if (t % (tsteps / 100) == 0) {
			printf("\rProgress: %3d%% ", (int) (100.0 * t / tsteps));
			fflush(stdout);
		}

		// steps 1 and 2 velocity Verlet
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				vel[i][j] += dt * 0.5 * acc[i][j];  // velocities after dt/2
				pos[i][j] += dt * vel[i][j];  // new positions
			}
		}
		
		// step 3 velocity Verlet
		get_acc_AL(acc, pos, cellLength, mass, natoms);
		
		// step 4 velocity Verlet
		vsquare = 0.0;
		for (i = 0; i < natoms; i++) {
			for (j = 0; j < 3; j++) {
				vel[i][j] += dt * 0.5 * acc[i][j];  // new velocities
				vsquare += vel[i][j] * vel[i][j];  // get v^2 for temp
			}
		}

		temp[t] = vsquare * mass / (3.0 * natoms * boltz);  // current temp in K
		pres[t] = (natoms * boltz * temp[t] + get_virial_AL(pos, cellLength, natoms)) /
				pow(cellLength, 3);  // pressure in eV / Å^3
		pres[t] = 1.6022 * pow(10,8) * pres[t];  // pressure in kPa

		// for calculating time averages of temperature and pressure
		if (t > startAveraging) {
			sumT += temp[t];
			sumP += pres[t];
		}

		// equilibration
		if (t > startEquilibration && t < endEquilibration) {

			// calculate alpha scaling factors
			alphaT = 1.0 + dt / tauT * (T_eq - temp[t]) / temp[t];
			alphaP = 1.0 - kappaT_tauP_ratio * dt * (P_eq - pres[t]);

			// let temperature and pressure decay
			for (i = 0; i < natoms; i++) {
				for (j = 0; j < 3; j++) {
					vel[i][j] = sqrt(alphaT) * vel[i][j];  // adjusted velocities
					pos[i][j] = pow(alphaP, 1.0/3.0) * pos[i][j];  // adjusted pressure
					latparam = pow(alphaP, 1.0/3.0) * latparam;  // adjust the latparam accordingly
					cellLength = pow(alphaP, 1.0/3.0) * cellLength;  // adjust the cell length accordingly
				}
			}
		}

		// save energies
		potEn[t] = get_potEn_AL(pos, cellLength, natoms);
		kinEn[t] = get_kinEn_AL(vel, mass, natoms);

		// save position for first and second atoms
		pos1[t][0] = pos[0][0];
		pos1[t][1] = pos[0][1];
		pos1[t][2] = pos[0][2];
		pos1[t][3] = pos[1][0];
		pos1[t][4] = pos[1][1];
		pos1[t][5] = pos[1][2];
	} // end velocity Verlet loop

	printf("\nAveraged over %d time steps.\n", tsteps - startAveraging);
	printf("Time average temperature: %.5f Celsius.\n", sumT / (tsteps - startAveraging) - 273.15);
	printf("Time average pressure:    %.5f kPa.\n", sumP / (tsteps - startAveraging));

	// write energies to file
	file_energy = fopen("energy.dat","w");
	for (i = 0; i < tsteps + 1; i++) {
		current_time = i * dt;
		totEn = potEn[i] + kinEn[i];
		fprintf(file_energy, "%e\t%e\t%e\t%e\t%e\t%e\n",
			current_time, potEn[i], kinEn[i], totEn, temp[i], pres[i]);
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
	free(temp); temp = NULL;
	free(pres); pres = NULL;

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
