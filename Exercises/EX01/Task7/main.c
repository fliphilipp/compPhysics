/*
E1_main.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2014-10-20.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#include "fft_func.h"
#define PI 3.141592653589
#define nbr_of_timesteps 65535 /* nbr_of_timesteps+1 = power of 2, for best speed */
#define nbr_of_particles 3 /* The number of particles is 3 */

/* Main program */
int main()
{
	/* Declartion of variables */
	double timestep;
	int i,j;
	double timestep_sq,current_time;
	double m[nbr_of_particles];
	double kappa;

	/* declare file variable */
	FILE *file;
    FILE *file2;
    FILE *file3;
    FILE *file4;

	/* displacement, velocity and acceleration */
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles]; 


	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double *q_1 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_2 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_3 = malloc((nbr_of_timesteps+1) * sizeof (double));
    double *pe = malloc((nbr_of_timesteps+1) * sizeof (double));
    double *ke = malloc((nbr_of_timesteps+1) * sizeof (double));
    double *pspec1 = malloc((nbr_of_timesteps+1) * sizeof (double));
    double *pspec2 = malloc((nbr_of_timesteps+1) * sizeof (double));
    double *pspec3 = malloc((nbr_of_timesteps+1) * sizeof (double));
    double *freq = malloc((nbr_of_timesteps+1) * sizeof (double));

	/* Set variables */
	timestep = pow(10,-4);
	m[0] = 1.6582 * pow(10,-3);  // for C-atom set to 16u =  1.6582 * 10^−3 eV(ps)^2 / Å^2
    m[1] = 1.2437 * pow(10,-3);  // for O-atom set to 12u =  1.2437 * 10^−3 eV(ps)^2 / Å^2
    m[2] = 1.6582 * pow(10,-3);  // for C-atom set to 16u =  1.6582 * 10^−3 eV(ps)^2 / Å^2
	kappa = 99.8627;  // The spring constant for CO2 is κ = 1.6 kN / m = 99.8627 eV / Å
	timestep_sq = timestep * timestep;

	
	/* Initial conditions */
	/* Set initial displacements and velocites */
	q[0] = 0.01 * 1.16;  // some reasonable initial displacement??
	v[0] = 0;

	for (i = 1; i < nbr_of_particles; i++) {
		q[i] = 0;
		v[i] = 0;
	}
	q_1[0] = q[0];
	q_2[0] = q[1];
	q_3[0] = q[2];
    
    // calculate initial kinetic and potential energies
    pe[0] = calc_pe(q, kappa, nbr_of_particles);
    ke[0] = calc_ke(v, nbr_of_particles, m);

	/* Calculate initial accelerations based on initial displacements */
	calc_acc(a, q, m, kappa, nbr_of_particles);

	/* timesteps according to velocity Verlet algorithm */
	for (i = 1; i < nbr_of_timesteps + 1; i++) {
		/* v(t+dt/2) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* q(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    q[j] += timestep * v[j];
		}

		/* a(t+dt) */
		calc_acc(a, q, m, kappa, nbr_of_particles);

		/* v(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* Save the displacement of the three atoms */
		q_1[i] = q[0];
		q_2[i] = q[1];
		q_3[i] = q[2];
        
        // save potential and kinetic energies
        pe[i] = calc_pe(q, kappa, nbr_of_particles);
        ke[i] = calc_ke(v, nbr_of_particles, m);
	}
    
    printf("displacements at end: %e \t %e \t %e.\n", q[0], q[1], q[2]);
    
    /* make FFT (powerspectrum) */
    powerspectrum(q_1, pspec1, nbr_of_timesteps);
    powerspectrum(q_2, pspec2, nbr_of_timesteps);
    powerspectrum(q_3, pspec3, nbr_of_timesteps);
    powerspectrum_shift(pspec1, nbr_of_timesteps);
    powerspectrum_shift(pspec2, nbr_of_timesteps);
    powerspectrum_shift(pspec3, nbr_of_timesteps);
    fft_freq_shift(freq, timestep, nbr_of_timesteps);

	/* Print displacement data to output file */
	file = fopen("disp.dat","w");
	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i] );	
		fprintf(file, "\n");
	}
	fclose(file);
    
    // write potential energy to file
    file2 = fopen("pe.dat","w");
    for (i = 0; i < nbr_of_timesteps + 1; i++) {
        current_time = i * timestep;
        fprintf(file2, "%.4f \t %e", current_time, pe[i]);
        fprintf(file2, "\n");
    }
    fclose(file2);
    
    // write kinetic energy to file
    file3 = fopen("ke.dat","w");
    for (i = 0; i < nbr_of_timesteps + 1; i++) {
        current_time = i * timestep;
        fprintf(file3, "%.4f \t %e", current_time, ke[i]);
        fprintf(file3, "\n");
    }
    fclose(file3);
    
    /* Print powerspectrum data to output file */
    file4 = fopen("pspec.dat","w");
    for (i = 0; i < nbr_of_timesteps + 1; i++) {
        fprintf(file4, "%.4f \t %e \t %e \t %e", freq[i], pspec1[i], pspec2[i], pspec3[i] );
        fprintf(file4, "\n");
    }
    fclose(file4);

	/* Free allocated memory */ 
	free(q_1); q_1 = NULL;
	free(q_2); q_2 = NULL;
	free(q_3); q_3 = NULL;
    free(pe); pe = NULL;
    free(ke); ke = NULL;
    free(pspec1); pspec1 = NULL;
    free(pspec2); pspec2 = NULL;
    free(pspec3); pspec3 = NULL;
    free(freq); freq = NULL;
	return 0;    
}
