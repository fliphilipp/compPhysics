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
	printf("--> Initialize variables.\n");

	// parameters
	int ndim = 4;  // grid size in each dimension  
	int tsteps = 15 * pow(10,3);  // number of iterations
	double dt = pow(10,-2);  // time step for integration (ps)
	int solid = 0;  // 1 for solid, 0 for fluid
	int startEquilibration = pow(10,2);  // when to start equilibrating
	int endEquilibration = 3 * pow(10,3);  // for how long to equilibrate
	int startAveraging = 3 * pow(10,3);  // when to start taking time averages
	int nDistanceMeasures = 100;  // how many times to calculate all distances for radial dist func
	int distanceMeasurementInterval = (int) (1.0 / dt); // measure each picosecond for radial dist func
	int boxesToMeasure = 1;  // measure distances up to twice the cell length for radial dist func
	double endIntermediateT = 5*pow(10,2); // until when we enforce intermediate temperature
	double mass = 27 * 1.0364 * pow(10, -4); // mass of Al [eV (ps)^2 / Å^2]
	double boltz = 8.617 * pow(10,-5);  // Boltzmann constant [eV / K]
	double P_eq = 101.325;  // Equilibrium pressure (101.325 kPa at sea level)
	double rand_disp = 0.05;  // initial dispacements percentage of latparam
	double tauT = 0.1;  // decay time constant for temperature scaling 
	double kappaT_tauP_ratio = pow(10,-10);  // only need to choose appropriate ratio for pressure scaling

	// adjust temperatures for solid vs. fluid case
	double intermediateT;
	double T_eq;
	if (solid == 1) {
		T_eq = 20.0 + 273.15;  // Equilibrium temperature in Kelvin (500 deg C)
		intermediateT = T_eq;  // no need for different intermediate temperature
	} else {
		T_eq = 700.0 + 273.15;  // Equilibrium temperature in Kelvin (700 deg C)
		intermediateT = 1000.0 + 273.15;  // for melting the system (1000 deg C)
	}

	// initialize variables
	double latparam;  // lattice spacing
	double cellLength;  // size of full simulation box
	double cellLengthLargeBox; 
	int natoms = 4 * ndim * ndim * ndim;  // total number of atoms
	double pos[natoms][3];  // positions of atoms
	double vel[natoms][3];  // velocities for atoms
	double acc[natoms][3];  // accelerations for atoms
	int nDistances = natoms * (natoms - 1) / 2;  // number pairwise distances
	int nBoxes = (int) pow(((2 * boxesToMeasure) + 1), 3);  // number of copies we need for radial dist func
	double largeBoxPos[natoms * nBoxes][3];  // to get all neighboring boxes
	int nDistLargeBox = natoms * (nBoxes * natoms - (natoms + 1) / 2); // number distances to large box for radial dist func
	int distancesMeasured;  // total number of distances measured (for info)
	double current_time;  // time in ps for printing to file
	double totEn;  // total energy for writing to file
	double vsquare;  // squared velocities for pressure calculation
	double alphaT;  // scaling factor for temperature equilibration
	double alphaP;  // scaling factor for pressure equilibration
	double sumT;  // sum of temperatures for temperature time average
	double sumP;  // sum of temperatures for pressure time average
	double desiredT;  // the desired temperature to decay to at a certain time
	double tAverage;  // time average of temperature
	double pAverage;  // time average of pressure
	double sumPotEn;  // sum of potential energies for time average
	double averagePotEn;  // time average of potential energy
	double sumFluctuationsPotEnSquare;  // sum of squared potential energy fluctuations for time average
	double averageFluctuationPotEnSquare;  // time average of squared fluctuations in potential energy
	double heatCapacity;  // calculated heat capaciy
	double molarHeatCapacity;  // to be able to compare to literature
	double sumTotEn;
	double averageTotEn;
	gsl_rng *my_rng = initialize_rng();  // random number generator
	double dx, dy, dz;  // for adding neighboring boxes
	int i, j, t, counter;  // iterator variables
	if (tsteps < endEquilibration) {endEquilibration = tsteps;}
	if (startEquilibration > tsteps) {startEquilibration = tsteps;}

	// allocate for variables to save over time
	double *potEn = malloc((tsteps + 1) * sizeof(double));
	double *kinEn = malloc((tsteps + 1) * sizeof(double));
	double *temp = malloc((tsteps + 1) * sizeof(double));
	double *pres = malloc((tsteps + 1) * sizeof(double));
	double (*pos1)[6] = malloc(sizeof(double[tsteps+1][6]));
	double *distLargeBox = malloc(nDistLargeBox * sizeof(double)); // distance to large box
	double *distances = malloc(nDistances * sizeof(double)); // pairwise distances of atoms

	printf("--> Let's fucking go!\n");

	// to figure out equilibrium spacing
	double latparamMin = 4.032;
	double latparamMax = 4.033;
	double sweepStep = 0.0001;
	int nSweep = (int)((latparamMax - latparamMin)/sweepStep);
	printf("--> Checking potential from %.5f to %.5f in steps of %.5f.\n", 
		latparamMin, latparamMax, sweepStep);
	double latParamList[nSweep];
	double energyLatParam[nSweep];
	double minEnergy = pow(10,20); // something large enough

	// initialize file variables
	FILE *file_pos;
	FILE *file_latEn;
	FILE *file_energy;
	FILE *file_pos1;
	FILE *file_largeBox;
	FILE *file_distancesLargeBox;

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
	printf("--> Minimum potential energy is %e eV at %.4f Ångström.\n", 
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
	sumPotEn = 0.0;
	sumTotEn = 0.0;

	// clear file for distances
	file_distancesLargeBox = fopen("distances-largebox.dat","w");
	fclose(file_distancesLargeBox);

	distancesMeasured = 0;
	// timesteps according to velocity Verlet algorithm
	for (t = 1; t < tsteps + 1; t++) {
		// print progress
		if (t % (tsteps / 100) == 0) {
			printf("\r--> Simulation progress: %3d%% ", (int) (100.0 * t / tsteps));
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

		temp[t] = vsquare * mass / (3.0 * natoms * boltz);  // save current temp in K
		pres[t] = (natoms * boltz * temp[t] + get_virial_AL(pos, cellLength, natoms)) /
				pow(cellLength, 3);  // pressure in eV / Å^3
		pres[t] = 1.6022 * pow(10,8) * pres[t];  // save pressure in kPa
		potEn[t] = get_potEn_AL(pos, cellLength, natoms);  // save E_pot
		kinEn[t] = get_kinEn_AL(vel, mass, natoms);  // save E_kin

		// for calculating time averages of temperature and pressure
		if (t > startAveraging) {
			sumT += temp[t];
			sumP += pres[t];
			sumPotEn += potEn[t];
			sumTotEn += potEn[t] + kinEn[t];
		}

		// equilibration
		if (t > startEquilibration && t < endEquilibration) {

			// set intermediate temperature for melting the system 
			if (t > startEquilibration && t < endIntermediateT) {
				desiredT = intermediateT;
			} else {
				desiredT = T_eq;
			}

			// calculate alpha scaling factors
			alphaT = 1.0 + dt / tauT * (desiredT - temp[t]) / temp[t];
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

		if (t > startAveraging && t % distanceMeasurementInterval == 0 && distancesMeasured < nDistanceMeasures) {
			// get positions of all surrounding boxes
			counter = 0;
			add_box(largeBoxPos, pos, 0.0, 0.0, 0.0, natoms, &counter);
			for (dx = - 1.0 * boxesToMeasure * cellLength; dx < 1.0001 * boxesToMeasure * cellLength; dx += cellLength) {
				for (dy = - 1.0 * boxesToMeasure * cellLength; dy < 1.0001 * boxesToMeasure * cellLength; dy += cellLength) {
					for (dz = - 1.0 * boxesToMeasure * cellLength; dz < 1.0001 * boxesToMeasure * cellLength; dz += cellLength) {
						if (dx != 0.0 || dy != 0.0 || dz != 0.0) {
							add_box(largeBoxPos, pos, dx, dy, dz, natoms, &counter);
							//printf("-----> dx = %.2f, dy = %.2f, dz = %.2f, counter = %d\n",dx, dy, dz, counter);
						}
					}
				}
			}

			// get distances to all atoms in the large box
			cellLengthLargeBox = 3.0 * cellLength;
			get_distances_largebox(largeBoxPos, cellLengthLargeBox, natoms, nBoxes);
			distancesMeasured ++;
		}


		// save position for first and second atoms
		pos1[t][0] = pos[0][0];
		pos1[t][1] = pos[0][1];
		pos1[t][2] = pos[0][2];
		pos1[t][3] = pos[1][0];
		pos1[t][4] = pos[1][1];
		pos1[t][5] = pos[1][2];
	} // end velocity Verlet loop

	// print results from time averages for temperature and pressure
	printf("\n--> Averaged over %d time steps.\n", tsteps - startAveraging);
	tAverage = sumT / (tsteps - startAveraging);  // average T in Kelvin
	pAverage = sumP / (tsteps - startAveraging);  // average P in kPa
	printf("--> Time average temperature: %.5f Celsius.\n", tAverage - 273.15);
	printf("--> Time average pressure:    %.5f kPa.\n", pAverage);

	// calculate heat capacity from fluctuations in potential energy
	averagePotEn = sumPotEn / (tsteps - startAveraging);
	sumFluctuationsPotEnSquare = 0.0;
	for (t = startAveraging; t < tsteps; t++) {
		sumFluctuationsPotEnSquare += (averagePotEn - potEn[t]) * (averagePotEn - potEn[t]);
	}
	averageFluctuationPotEnSquare = sumFluctuationsPotEnSquare / (tsteps - startAveraging);
	heatCapacity = (3.0 * natoms * boltz / 2.0) * 1.0 / (1.0 - 2.0 * averageFluctuationPotEnSquare 
		/ (3 * natoms * boltz * boltz * tAverage * tAverage));
	molarHeatCapacity = heatCapacity * 1.6022 * pow(10,-19) * 6.022 * pow(10,23) / 256.0;
	printf("--> Average squared fluctuation in potential energy: %.5f\n", 
		averageFluctuationPotEnSquare);
	printf("--> Heat capacity: %.5f ev/K.\n", heatCapacity);
	printf("--> Molar heat capacity: %.5f J / (K * mol).\n", molarHeatCapacity);
	averageTotEn = sumTotEn / (tsteps - startAveraging);

	printf("--> Cell length: %.5f\n", cellLength);
	file_distancesLargeBox = fopen("distances-largebox.dat","a"); // write cell length for use in matlab
		fprintf(file_distancesLargeBox, "%.10f\n", cellLength);
		fprintf(file_distancesLargeBox, "%d\n", nDistanceMeasures);
	fclose(file_distancesLargeBox);

	/*
	// get distances of simulating box
	get_pairwise_distances(distances, pos, cellLength,natoms);
	file_distances = fopen("distances.dat","w");
	fprintf(file_distances, "%e\n",cellLength);
	for (i = 0; i < nDistances; i++) {
		fprintf(file_distances, "%e\n",distances[i]);
	}
	fclose(file_distances);
	*/

	// write positions of large box to file
	file_largeBox = fopen("largeBox.dat","w");
	for (i = 0; i < 27 * natoms; i++) {
		fprintf(file_largeBox, "%.4f\t%.4f\t%.4f\n",
			largeBoxPos[i][0], largeBoxPos[i][1], largeBoxPos[i][2]);
	}
	fclose(file_largeBox);


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
	free(distLargeBox); distLargeBox = NULL;
	free(distances); distances = NULL;

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
