/*
 E3.c
 This program reads data from the file MC.txt.
 Created by Anders Lindman on 2015-11-12.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include "helper_funcs.h"
#define PI 3.141592653589

gsl_rng* initialize_rng();

// Task 1
// Evaluate the integral with the Monte Carlo method
void one_dim_integral(){

  // Initialize vars
  int x, t, N;
  double error, sum, avg, var_avg, square_sum;
  // random number generator
  gsl_rng *my_rng = initialize_rng();

  // Print expected value to compare
  printf("*******************************************\n");
  printf("************UNIFORM DISTRIBUTION***********\n");
  printf("*******************************************\n");
  printf("The expected value of the integral is approx: %.8F\n", (double) 1/6);
  printf("*******************************************\n");

  // Run for 10, 10^2, 10^3 and 10^4
  for(x = 0; x < 4; x++){
    
    N = pow(10,x+1);
    // Reset according to new length N
    double x_var[N];
    sum = 0;
    square_sum = 0;

    // Integrate!
    for(t = 0; t < N; t++){

      // Get a random sample
      x_var[t] = gsl_rng_uniform(my_rng);

      
      // Calculate the big sum
      sum += x_var[t] * (1 - x_var[t]);

      // Calculate the variance
      square_sum += pow(x_var[t] * (1 - x_var[t]),2);

    }

    // Take the average of the big sum
    avg = sum / N;

    // Calculate the error margin
    var_avg = square_sum / N;
    error = (var_avg - pow(avg,2));
    error = sqrt(error);

    // Results!!
    printf("N = 10^%i | Integral value = %.4F | Error margin = %.4F \n", x+1, avg, error );

  }
  printf("*******************************************\n");
}

// Task 2
// Monte Carlo sampling from sine distribution [0,1]
void sine_integral(){

  // Initialize vars
  int x, t, N;
  double r, error, sum, avg, var_avg, square_sum;
  // Set x array size to max N
  double x_var[4][10000];
  // random number generator
  gsl_rng *my_rng = initialize_rng();

  // Print expected value to compare
  printf("*************SINE DISTRIBUTION*************\n");
  printf("*******************************************\n");
  printf("The expected value of the integral is approx: %.8F\n", (double) 1/6);
  printf("*******************************************\n");

  // Run for 10, 10^2, 10^3 and 10^4
  for(x = 0; x < 4; x++){
    
    N = pow(10,x+1);
    sum = 0;
    square_sum = 0;

    // Integrate!
    for(t = 0; t < N; t++){

      // Get a random sample
      r = gsl_rng_uniform(my_rng);

      // Push the random sample through the sine distribution
      x_var[x][t] = acos(1 - 2*r) / PI;
      
      // Calculate the big sum
      sum += x_var[x][t] * (1 - x_var[x][t]) * 2 / sin(PI * x_var[x][t]) / PI;

      // Calculate the variance
      square_sum += x_var[x][t] * (1 - x_var[x][t]) * 2 \
                    / sin(PI * x_var[x][t]) * x_var[x][t] * (1 - x_var[x][t]) * 2 \
                    / sin(PI * x_var[x][t]) \
                    / PI \
                    / PI;

    }

    // Take the average of the big sum
    avg = sum / N;

    // Calculate the error margin
    var_avg = square_sum / N;
    error = (var_avg - pow(avg,2))/N;
    error = sqrt(error);

    // Results!!
    printf("N = 10^%i | Integral value = %.4F | Error margin = %.4F \n", x+1, avg, error );

  }
  printf("*******************************************\n");

  // Save results in file:
  FILE *out_file;
  out_file = fopen("task2.data", "w");

  for(x = 0; x < 10000; x++){
    fprintf(out_file, "%F \t %F \t %F \t %F \n", x_var[0][x], x_var[1][x], x_var[2][x], x_var[3][x]);
  }

  fclose(out_file);
}

void metropolis(){
  // Initialize vars
  int N = 10000, t, x;
  double p_max = 0.01;
  double x_var[N], p_var, sum, square_sum, trial, r;
  double error, avg, var_avg;
  double delta = 0.55;
  int total = 0;

  // random number generator
  gsl_rng *my_rng = initialize_rng();

  printf("************METROPOLS INTEGRAL*************\n");
  printf("*******************************************\n");
  printf("The expected value of the integral is approx: %.8F\n", (double) 1/6);
  printf("*******************************************\n");

  // Prepare variables before looping
  x_var[0] = 0.5;
  p_var = 2 / sin(PI * x_var[0]) / PI;

  for(t = 1; t < N; t++){
    
    // Let the trials begin!
    do {

      // Get random
      do {
        r = gsl_rng_uniform(my_rng);
      } while (r > 1.0);

      // Calculate next value of x for the trial
      x_var[t] =  x_var[t - 1] + delta * (r - 0.5);

      // Find the probability
      p_var = 2 / sin(PI * x_var[t]) / PI;

      trial = p_var / p_max;

      // Get random
      do {
        r = gsl_rng_uniform(my_rng);
      } while (r > 1.0);

      total++;

    // Test for the Rejection Method
    } while (trial < r);
    
    // Calculate the big sum
    sum += x_var[t] * (1 - x_var[t]) * 2 / sin(PI * x_var[t]) / PI;

    // Calculate the variance
    square_sum += x_var[t] * (1 - x_var[t]) * 2 \
                  / sin(PI * x_var[t]) * x_var[t] * (1 - x_var[t]) * 2 \
                  / sin(PI * x_var[t]) \
                  / PI \
                  / PI;

  }


  // Take the average of the big sum
  avg = sum / N;

  // Calculate the error margin
  var_avg = square_sum / N;
  error = (var_avg - pow(avg,2))/N;
  error = sqrt(error);

  // Results!!
  printf("N = %i | Integral value = %.4F | Error margin = %.4F \n", N, avg, error );
  printf("Percentage of approved trials = %.1f\n", 100*(N/total) );
  printf("%i\n", total);

  // Save results in file:
  FILE *out_file;
  out_file = fopen("task3.data", "w");

  for(x = 0; x < N; x++){
    fprintf(out_file, "%F \n", x_var[x]);
  }

  fclose(out_file);
}

int main()
{
  int i, nbr_of_lines;
  FILE *in_file;
  
  nbr_of_lines = 1e6; /* The number of lines in MC.txt. */
  double *data = malloc((nbr_of_lines) * sizeof (double));
  
  /* Read data from file. */
  in_file = fopen("MC.txt","r");
  for (i=0; i<nbr_of_lines; i++) {
    fscanf(in_file,"%lf",&data[i]);
  }
  fclose(in_file);


  // Run Task 1
  one_dim_integral();

  // Run Task 2
  sine_integral();

  // Run Task 3
  metropolis();
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
