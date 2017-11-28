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
  printf("The expected value of the integral is approx: %.4F\n", (double) 1/6);
  printf("*******************************************\n");

  // Run for 10, 10^2, 10^3 and 10^4
  for(x = 0; x < 4; x++){
    
    N = pow(10, x+1);
    sum = 0.0;
    square_sum = 0.0;

    // Integrate!
    for(t = 0; t < N; t++){

      // Get a random sample
      x_var[x][t] = gsl_rng_uniform(my_rng);

      
      // Calculate the big sum
      new = x_var[x][t] * (1.0 - x_var[x][t]);
      sum += new;

      // Calculate the variance
      square_sum += pow(new, 2.0);

    }

    error = sqrt(square_sum / N - pow(sum / N, 2.0)) / sqrt(N);

    // Results!!
    printf("N = 10^%i | Integral value = %.6F +/- %.6F \n", x+1, sum / N, error );

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
  double r, error, sum, avg, var_avg, square_sum, new, fx, px;
  // Set x array size to max N
  double x_var[4][10000];
  // random number generator
  gsl_rng *my_rng = initialize_rng();

  // Print expected value to compare
  printf("*******************TASK 2******************\n");
  printf("*******************************************\n");
  printf("The expected value of the integral is approx: %.4F\n", (double) 1/6);
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

      // Push the random sample through the reverse of the integral of the weight function 
      // Integral of sin(PI * r) => -cos(PI * r) / PI
      // Inverse => acos(- PI * r) / PI
      x_var[x][t] = acos(- r) / PI;
      x_var[x][t] = acos(1.0 - 2.0*r) / PI;
      
      // Calculate the big sum
      fx = x_var[x][t] * (1.0 - x_var[x][t]);
      px = PI * sin(x_var[x][t] * PI) / 2.0;
      new = fx / px;
      sum += new;

      // Calculate the variance
      square_sum += pow(new,2);

    }

    error = sqrt(square_sum / N - pow(sum / N, 2.0)) / sqrt(N);

    // Results!!
    printf("N = 10^%i | Integral value = %.6F +/- %.6F \n", x+1, sum / N, error);

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
  int N = pow(10,8), t=0, i;

  // Initialize probability of state null to a low value
  double p_m = 0.01;

  // Allocate memory for big arrays
  double *x_var = malloc((N) * sizeof(double));
  double *y_var = malloc((N) * sizeof(double));
  double *z_var = malloc((N) * sizeof(double));
  double *p_n = malloc((N) * sizeof(double));

  double sum, square_sum = 0.0, trial, r, ft, px, x, y, z;
  double error, avg, var_avg, new;
  double delta = 2.0;
  int total = 0;

  // random number generator
  gsl_rng *my_rng = initialize_rng();

  printf("*******************TASK 3******************\n");
  printf("*******************************************\n");

  // Prepare variables before looping
  x_var[0] = 0.0;
  y_var[0] = 0.0;
  z_var[0] = 0.0;
  x = 0.0;
  y = 0.0;
  z = 0.0;
  p_n[0] = pow(PI,-3.0/2.0) * exp(- (x*x + y*y + z*z));
  sum = 0.0;
  square_sum = 0.0;

  for(t = 1; t < N; t++){
    
    // Let the trials begin!

    // Get random
    r = gsl_rng_uniform(my_rng);
    // Calculate next value of x for the trial
    x_var[t] =  x_var[t - 1] + delta * (r - 0.5);

    // Get random
    r = gsl_rng_uniform(my_rng);
    // Calculate next value of y for the trial
    y_var[t] =  y_var[t - 1] + delta * (r - 0.5);

    // Get random
    r = gsl_rng_uniform(my_rng);
    // Calculate next value of y for the trial
    z_var[t] =  z_var[t - 1] + delta * (r - 0.5);

    // Transform through the weight function!
    p_n[t] = pow(PI,-3.0/2.0) * exp(- (pow(x_var[t],2) + pow(y_var[t],2) + pow(z_var[t],2)) );
    // Ratio of current state prob over previous state prob
    trial = p_n[t] / p_m;

    // Get new random
    r = gsl_rng_uniform(my_rng);


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

    x = x_var[t];
    y = y_var[t];
    z = z_var[t];
    
    // Calculate the big sum
    ft = pow(PI,-3.0/2.0) * (x*x + x*x * y*y +  x*x * y*y * z*z) * exp( -(x*x + y*y + z*z));
    new = ft / p_m;

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
  printf("N = %i | Integral value = %.6F | Error margin = %.6F \n", N, avg, error );
  printf("*******************************************\n");
  printf("Percentage of approved trials = %.1f\n", (double) 100*(((double) N - (double) total) / (double) N) );
  printf("*******************************************\n");
  
  // Save results in file:
  FILE *out_file;
  out_file = fopen("task3.data", "w");

  for(i = 0; i < 10000; i++){
    fprintf(out_file, "%F \t %F \t %F \t %F\n", x_var[i], y_var[i], z_var[i], p_n[i]);
  }

  fclose(out_file);

  // free memory
  free(x_var); free(y_var); free(z_var); free(p_n);
  x_var = NULL; y_var = NULL; z_var = NULL; p_n = NULL;
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
  // Also calculate the sum of the function values to that point aka relaxation time
  decay_time = 0;
  while(correlation_func[decay_time] >= exp(-2)){
    relaxation_time += correlation_func[decay_time];
    decay_time++;
  }

  // Statistical inefficiency approx is twice the relaxation time 
  stat_ineff = (double) relaxation_time * 2;

  printf("Function reaches 10%% of its initial value at time t = %i\n", decay_time);
  printf("Statistical inefficiency s = %F \n", stat_ineff);
  total_error = sqrt((square_mean_fi - mean_fi*mean_fi)/steps*stat_ineff);
  printf("Result: %.4f | Error margin: %.4e \n", mean_fi, total_error);

  printf("*******************************************\n");

  // Save results in file:
  FILE *out_file;
  out_file = fopen("task4a.data", "w");

  for(int x = 0; x < steps; x++){
    fprintf(out_file, "%F\n", correlation_func[x]);
  }

  fclose(out_file);

}

// Task 4b
void block_averaging(double *array, int array_length){
  printf("******************TASK 4b******************\n");
  printf("*******************************************\n");

  // Declaration and initiation of variables
  int i, j;
  int min_block_size = 5;
  int max_block_size = 1000;
  int block_size = max_block_size;
  int nbr_blocks = array_length/block_size;
  double block_means[nbr_blocks];
  double square_block_means[nbr_blocks];
  double total_error, var_F;
  double mean = 0.0, square_mean = 0.0, var_f, s[max_block_size];
  double new, mean_var_F = 0.0, square_mean_var_F = 0.0;


  mean = 0.0;
  square_mean = 0.0;
  // Determine the variance for the whole array
  for(i = 0; i < array_length; i++){
    mean += array[i] / array_length;
    square_mean += pow(array[i],2) / array_length;
  }

  var_f = square_mean - pow(mean,2);

  // Calculate statistical inefficiency for all values of block sizes
  for (int block_size = min_block_size; block_size < max_block_size + 1; ++block_size){

    nbr_blocks = (int) array_length/block_size;
    double block_means[nbr_blocks];
    double square_block_means[nbr_blocks];
    mean_var_F = 0.0;
    square_mean_var_F = 0.0;
    var_F = 0.0;
    
    // Determine the average in each block
    for(i = 0; i < nbr_blocks; i++){
      
      block_means[i] = 0.0;
      square_block_means[i] = 0.0;

      for(j = 0; j < block_size; j++){
        block_means[i] += array[( i * block_size + j )] / block_size;
        square_block_means[i] += pow(array[( i * block_size + j )],2) / block_size;
      }
      
      var_F += pow(block_means[i] - mean,2) / nbr_blocks;
    }


    s[block_size] = block_size * var_F / var_f;
  }

  printf("Statistical inefficiency s = %f\n", s[max_block_size]);

  total_error = pow(var_F,2) / sqrt(max_block_size);
  printf("Variance from original mean: %.4e \n", var_F);
  printf("*******************************************\n");

  // Save results in file:
  FILE *out_file;
  out_file = fopen("task4b.data", "w");

  for(int x = 0; x < block_size; x++){
    fprintf(out_file, "%F\n", s[x]);
  }

  fclose(out_file);
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
