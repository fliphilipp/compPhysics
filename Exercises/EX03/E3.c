/*
 E3.c
 This program reads data from the file MC.txt.
 Created by Anders Lindman on 2015-11-12.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include "helper_funcs.h"
#define PI 3.141592653589

gsl_rng* initialize_rng();

// Task 1
// Evaluate the integral with the Monte Carlo method
void one_dim_integral(){

  // Initialize vars
  int x, t, N;
  double error, sum, avg, var_avg, square_sum, new;
  double x_var[4][10000];
  // random number generator
  gsl_rng *my_rng = initialize_rng();

  // Print expected value to compare
  printf("*******************************************\n");
  printf("*******************TASK 1******************\n");
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
      x_var[x][t] = gsl_rng_uniform(my_rng);

      
      // Calculate the big sum
      new = x_var[x][t] * (1 - x_var[x][t]);
      sum += new;

      // Calculate the variance
      square_sum += pow(new,2);

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

  // Save results in file:
  FILE *out_file;
  out_file = fopen("task1.data", "w");

  for(x = 0; x < 10000; x++){
    fprintf(out_file, "%F \t %F \t %F \t %F \n", x_var[0][x], x_var[1][x], x_var[2][x], x_var[3][x]);
  }

  fclose(out_file);

}

// Task 2
// Monte Carlo sampling from sine distribution [0,1]
void sine_integral(){

  // Initialize vars
  int x, t, N;
  double r, error, sum, avg, var_avg, square_sum, new;
  // Set x array size to max N
  double x_var[4][10000];
  // random number generator
  gsl_rng *my_rng = initialize_rng();

  // Print expected value to compare
  printf("*******************TASK 2******************\n");
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
      x_var[x][t] = asin(r) / PI;
      
      // Calculate the big sum
      new = x_var[x][t] * (1 - x_var[x][t]);
      sum += new;

      // Calculate the variance
      square_sum += pow(new,2);

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
  int N = 10000, t=0, x;
  double p_m = 0.5; // Initialize probability of state null to 0.5
  double x_var[N], y_var[N], z_var[N]; 
  double p_n[N], sum, square_sum = 0, trial, r;
  double error, avg, var_avg, new;
  double delta = -0.9;
  int total = 0;

  // random number generator
  gsl_rng *my_rng = initialize_rng();

  printf("*******************TASK 3******************\n");
  printf("*******************************************\n");

  // Prepare variables before looping
  x_var[0] = 0;
  y_var[0] = 0;
  z_var[0] = 0;
  p_n[0] = sqrt(-x_var[t] * x_var[t] - log(pow(PI,-3/2) * y_var[t]) \
          - z_var[t] * z_var[t]);
  sum = 0;
  square_sum = 0;

  for(t = 1; t < N; t++){
    
    // Let the trials begin!

    // Get random
    do {
      r = gsl_rng_uniform(my_rng);
    } while (r > 1.0); // Just to be safe

    // Calculate next value of x for the trial
    x_var[t] =  x_var[t - 1] + delta * (r - 0.5);
    y_var[t] =  y_var[t - 1] + delta * (r - 0.5);
    z_var[t] =  z_var[t - 1] + delta * (r - 0.5);

    // Transform through the inverse of the weight function!
    p_n[t] = sqrt( - pow(x_var[t],2) - log(pow(PI,-3/2) * y_var[t]) \
          - pow(z_var[t],2));
    // Ratio of current state prob over previous state prob
    trial = p_n[t] / p_m;

    // Get new random
    do {
      r = gsl_rng_uniform(my_rng);
    } while (r > 1.0); // Just to be safe


    // Test for the Rejection Method
    if(trial < r){
      // Not accepted, keep previous value
      x_var[t] = x_var[t - 1];
      y_var[t] = y_var[t - 1];
      z_var[t] = z_var[t - 1];
      total++;
    } else {
      // Accepted
      // Update the probability of the previous state
      // to the one of the current state
      p_m = p_n[t];
    }
    
    // Calculate the big sum
    new = pow(PI,-3/2) * pow(x_var[t],2) + pow(x_var[t],2) * pow(y_var[t],2) \
          + pow(x_var[t],2) * pow(y_var[t],2) * pow(z_var[t],2) \
          * exp( - (pow(x_var[t],2) + pow(y_var[t],2) + pow(z_var[t],2)) );

    sum += new;

    // Calculate the variance
    square_sum += pow(new,2);

  }


  // Take the average of the big sum
  avg = sum / N;

  // Calculate the error margin
  var_avg = square_sum / N;
  error = (var_avg - pow(avg,2))/N;
  error = sqrt(error);

  // Results!!
  printf("N = %i | Integral value = %.4F | Error margin = %.4F \n", N, avg, error );
  printf("*******************************************\n");
  printf("Percentage of approved trials = %.1f\n", (double) 100*(((double) N - (double) total) / (double) N) );
  printf("*******************************************\n");
  
  // Save results in file:
  FILE *out_file;
  out_file = fopen("task3.data", "w");

  for(x = 0; x < N; x++){
    fprintf(out_file, "%F \t %F \t %F \t %F\n", x_var[x], y_var[x], z_var[x], p_n[x]);
  }

  fclose(out_file);
}

// Task 4a
void correlation_function(double *array, int array_length){
  printf("******************TASK 4a******************\n");
  printf("*******************************************\n");

  // Initialize vars
  int steps = 300; // Steps for the correlation function
  int decay_time;
  double mean_fi = 0, square_mean_fi = 0;
  double correlation_func[steps], fik[steps];
  double stat_ineff = 0, relaxation_time = 0, total_error;

  // Array of zeros
  for (int k = 0; k < steps; ++k){
    fik[k] = 0.0;
  }

  // Calculate <f_i> and <f_i>^2
  for (int i = 0; i < array_length - steps; ++i){
    mean_fi += array[i] / (array_length - steps);
    square_mean_fi += pow(array[i],2) / (array_length - steps);
  }

  // Calculate <f_(i+k)*f_i> 
  for (int i = 0; i < array_length - steps; ++i){
    for (int k = 0; k < steps; ++k){
      fik[k] += (array[i + k] * array[i]) / (array_length - steps);
    } 
  }

  // Calculate correlation function
  for (int k = 0; k < steps; ++k){
    correlation_func[k] = ( fik[k] - pow(mean_fi,2) ) \
                          / ( square_mean_fi - pow(mean_fi,2)); 
  }

  // Calculate the time when the function will have decayed by approximately 10%
  // aka statistical inefficiency
  decay_time = 0;
  while(correlation_func[decay_time] >= exp(-2)){
    relaxation_time += correlation_func[decay_time];
    decay_time++;
  }

  stat_ineff = relaxation_time * 2;

  printf("Function reaches 10%% of it's initial value %i times!\n", decay_time);
  total_error = sqrt((square_mean_fi - mean_fi*mean_fi)/steps*stat_ineff);
  printf("Statistical inefficiency: %F \n", stat_ineff);
  printf("Result: %.8e | Error margin: %.10e \n", mean_fi, total_error);

  printf("*******************************************\n");
}

// Task 4b
void block_averaging(double *array, int array_length){
  printf("******************TASK 4b******************\n");
  printf("*******************************************\n");



  printf("*******************************************\n");
}

int main()
{
  int i, nbr_of_lines;
  FILE *in_file;
  
  nbr_of_lines = 1e6; /* The number of lines in MC.txt. */
  double *data = malloc((nbr_of_lines) * sizeof (double));
  double temp;
  
  /* Read data from file. */
  in_file = fopen("MC.txt","r");
  for (i=0; i<nbr_of_lines; i++) {
    temp = fscanf(in_file,"%lf",&data[i]);
  }
  fclose(in_file);


  // Run Task 1
  one_dim_integral();

  // Run Task 2
  sine_integral();

  // Run Task 3
  metropolis();

  // Run Task 4a
  correlation_function(data, nbr_of_lines);

  // Run Task 4b
  block_averaging(data, nbr_of_lines);
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
