/*
 * For all of these functions, the first column of the data_matrix is the 
 * dependent variable, and the remainder are the independent.
 * The last argument to the maximum_likelihood function is 1=verbose; 0=not.
 */

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include "linear_algebra.h"


double probit_likelihood(const gsl_vector *beta, void *d);
void d_probit_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient);
void probit_fdf( const gsl_vector *beta, void *d, double *f, gsl_vector *df);
//usage:
//gsl_vector * probit_parameter = gsl_vector_alloc(data_size-1);
//maximum_likelihood_w_d(data_matrix, &probit_parameter, data_size-1, probit_likelihood, d_probit_likelihood, probit_fdf, 0);

double yule_likelihood(const gsl_vector *beta, void *d);
//usage:
//gsl_vector * yule_parameter = gsl_vector_alloc(1);
//maximum_likelihood(data_matrix, &yule_parameter, 1, yule_likelihood, 0);



double zipf_likelihood(const gsl_vector *beta, void *d);
void d_zipf_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient);
void zipf_fdf(const gsl_vector *beta, void *d, double *f, gsl_vector *df);

//usage:
//gsl_vector * zipf_parameter = gsl_vector_alloc(1);
//maximum_likelihood_w_d(data_matrix, &zipf_parameter, 1, zipf_likelihood, d_zipf_likelihood, zipf_fdf, 0);


#define MAX_ITERATIONS 		500000
#define MAX_ITERATIONS_w_d	500000


//The maximum likelihood functions themselves. Call them using one of the likelihood fns as above.
void	maximum_likelihood_w_d(void * data, gsl_vector **betas, int betasize,
					double (* likelihood)(const gsl_vector *beta, void *d),
					void (* d_likelihood)(const gsl_vector *beta, void *d, gsl_vector *df), 
					void (* fdf)(const gsl_vector *beta, void *d, double *f, gsl_vector *df), int verbose);

double	maximum_likelihood(void * data, gsl_vector **betas, int betasize,
					double (* likelihood)(const gsl_vector *beta, void *d), int verbose);
//Feed in data, the parameters (to be output), the # of parameters,
//and a pointer to the likelihood fn itself (which takes beta first,
//then data).  You'll get the most likely betas back out.
