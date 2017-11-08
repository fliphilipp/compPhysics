/*
 E2.c
 
 Created by Anders Lindman on 2014-11-04.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "func.h"
#define PI 3.141592653589
#define npart 32  // number of particles

/* Main program */
int main()
{
    // parameters
    int tmax = pow(10,6);     // legth of simulation
    double dt = 0.1;    // time step
    double a0 = 1.0;    // equilibrium separation between the particles
    double m = 1.0;     // mass
    double kappa = 1.0; // spring constant
    double E0 = npart;    // initial energy
    double alpha = 0.1;  // nonlinearity factor
    int modesToPlot = npart; // number of modes to plot
    int saveEach = 100;  // interval in which to save energy modes
    
    double sum;
    double current_time;
    double omega;
    int i,j,k;
    int tsteps = tmax / dt;  // number of iterations
    
    // file variable
    FILE *file_disp;
    FILE *file_energy;
    FILE *file_modes;
    FILE *file_tavg;
    
    // normal coordinates and moments
    double Q[npart];
    double P[npart];
    
    // ordinary coordinates and moments, and acceleration
    double q[npart];
    double p[npart];
    double a[npart];
    
    // allocating memory for displacements and moments: size (tsteps x npart)
    double (*q_save)[npart] = malloc(sizeof(double[tsteps/saveEach+1][npart]));
    double (*p_save)[npart] = malloc(sizeof(double[tsteps/saveEach+1][npart]));
    double *kin_en = malloc((tsteps/saveEach+1) * sizeof (double));
    double *pot_en = malloc((tsteps/saveEach+1) * sizeof (double));
    double (*E_save)[npart] = malloc(sizeof(double[tsteps/saveEach+1][npart]));
    double (*E_tavg)[npart] = malloc(sizeof(double[tsteps/saveEach+1][npart]));
    double tot_en;
    
    printf("Let's fucking go!\n");
    
    // initial condition
    Q[0] = 0;
    P[0] = sqrt(2.0 * E0);
    for (i = 1; i < npart; i++) {
        Q[i] = 0;
        P[i] = 0;
    }
    
    // construct the transformation matrix
    double factor;
    double trans_matrix[npart][npart];
    factor = 1 / ((double) npart + 1);
    for (i=0; i < npart; i++) {
        for (j=0; j < npart; j++) {
            trans_matrix[i][j] = sqrt(2.0 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        }
    }
    
    // from normal Q to q
    for (i = 0; i < npart; i++){
        sum = 0;
        for (j = 0; j < npart; j++){
            sum += Q[j] / sqrt(m) * trans_matrix[i][j];
        }
        q[i] = sum;
    }
    
    // from normal P to p
    for (i = 0; i < npart; i++){
        sum = 0;
        for (j = 0; j < npart; j++){
            sum += P[j] * sqrt(m) * trans_matrix[i][j];
        }
        p[i] = sum;
    }
    
    // initialize acceleration (zero since zero displacements)
    for (i = 0; i < npart; i++) {
        a[i] = 0;
    }
    
    // calculate initial kinetic and potential energies
    pot_en[0] = calc_pe(q, kappa, npart);
    kin_en[0] = calc_ke(p, npart, m);
    for (j = 0; j < npart; j++) {
        omega = 2.0*sqrt(kappa/m)*sin(0.5*(j+1)*PI*factor);
        E_save[0][j] = 1/2.0 * (P[j]*P[j] + omega*omega*Q[j]*Q[j]);
        E_tavg[0][j] = E_save[0][j];
    }
    
    // timesteps according to velocity Verlet algorithm
    for (i = 1; i < tsteps + 1; i++) {
        // p(t+dt/2)
        for (j = 0; j < npart; j++) {
            p[j] += dt * 0.5 * a[j];
        }
        
        // q(t+dt)
        for (j = 0; j < npart; j++) {
            q[j] += dt * p[j];
        }
        
        // a(t+dt)
        calc_acc(a, q, m, kappa, alpha, npart);
        
        // p(t+dt)
        for (j = 0; j < npart; j++) {
            p[j] += dt * 0.5 * a[j];
        }
        
        if (i % saveEach == 0) {
            // Save the displacement of the atoms
            for (j = 0; j < npart; j++) {
                q_save[i/saveEach][j] = q[j];
            }
            
            // save potential and kinetic energies
            pot_en[i/saveEach] = calc_pe(q, kappa, npart);
            kin_en[i/saveEach] = calc_ke(p, npart, m);
        }//end saving displacements and energies
        
        if (i % saveEach == 0) {
            // from ordinary q to normal Q
            for (k = 0; k < npart; k++){
                sum = 0;
                for (j = 0; j < npart; j++){
                    sum += q[j] * sqrt(m) * trans_matrix[k][j];
                }
                Q[k] = sum;
            }
            
            // from ordinary p to normal P
            for (k = 0; k < npart; k++){
                sum = 0;
                for (j = 0; j < npart; j++){
                    sum += p[j] / sqrt(m) * trans_matrix[k][j];
                }
                P[k] = sum;
            }
            
            // save energy modes
            for (j = 0; j < npart; j++) {
                omega = 2.0*sqrt(kappa/m)*sin(0.5*(j+1)*PI*factor);
                E_save[i/saveEach][j] = 1/2.0 * (P[j]*P[j] + omega*omega*Q[j]*Q[j]);
                E_tavg[i/saveEach][j] = E_tavg[i/saveEach-1][j] + E_save[i/saveEach][j]; // just add up for now
            }
        }//end saving energy modes
        
    }//end Verlet
    
    // save energy modes to output file
    file_modes = fopen("modes.dat","w");
    for (i = 0; i < tsteps + 1; i = i + saveEach) {
        current_time = i * dt;
        fprintf(file_modes, "%.4f\t", current_time);
        for (j = 0; j < modesToPlot; j++) {
            fprintf(file_modes, "%.4f\t", E_save[i/saveEach][j]);
        }
        fprintf(file_modes, "\n");
    }
    fclose(file_modes);
    
    // save energy mode time average to output file
    file_tavg = fopen("mode_tavg.dat","w");
    for (i = 0; i < tsteps + 1; i = i + saveEach) {
        current_time = i * dt;
        fprintf(file_tavg, "%.4f\t", current_time);
        for (j = 0; j < modesToPlot; j++) {
            // devide by number of records up to date
            fprintf(file_tavg, "%.35f\t", E_tavg[i/saveEach][j] / ((double) i/saveEach + 1));
        }
        fprintf(file_tavg, "\n");
    }
    fclose(file_tavg);
    

    // save displacements to output file
    file_disp = fopen("disp.dat","w");
    for (i = 0; i < tsteps + 1; i = i + saveEach) {
        current_time = i * dt;
        fprintf(file_disp, "%.4f\t", current_time);
        for (j = 0; j < npart; j++) {
            fprintf(file_disp, "%.4f\t", q_save[i/saveEach][j]);
        }
        fprintf(file_disp, "\n");
    }
    fclose(file_disp);
    
    // save energies to file
    file_energy = fopen("energy.dat","w");
    for (i = 0; i < tsteps + 1; i = i + saveEach) {
        current_time = i * dt;
        tot_en = pot_en[i/saveEach] + kin_en[i/saveEach];
        fprintf(file_energy, "%.4f\t%.4f\t%.4f\t%.4f", current_time, pot_en[i/saveEach], kin_en[i/saveEach], tot_en);
        fprintf(file_energy, "\n");
    }
    fclose(file_energy);
    
    // free allocated memory
    free(q_save); q_save = NULL;
    free(p_save); p_save = NULL;
    free(kin_en); kin_en = NULL;
    free(pot_en); pot_en = NULL;
    free(E_save); E_save = NULL;
    free(E_tavg); E_save = NULL;
    
    return 0;
} // end main
