#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include "apop_linear_algebra.h"
#include "apop_conversions.h"

#define MAX_ITERATIONS 		500
#define MAX_ITERATIONS_w_d	500

double mle_probit(gsl_matrix *data, gsl_vector **beta, double *starting_pt, double step_size, int verbose);

double mle_waring(gsl_matrix *data, gsl_vector **beta, double *starting_pt, double step_size, int verbose);
double mle_yule(gsl_matrix *data, gsl_vector **beta, double *starting_pt, double step_size, int verbose);
double mle_zipf(gsl_matrix *data, gsl_vector **beta, double *starting_pt, double step_size, int verbose);

/*
For the Probit, the first column of the data matrix is the dependent
variable, and the remaining variables are the independent. This means
that the beta which will be output will be of size (data->size2 - 1).

For the Waring, Yule, and Zipf estimates, each row of the data matrix
is one observation. The first column is the number of elements with one
link, the second is the number of elements with two links, et cetera.

Beta needs to be declared but should not be allocated; these functions
do the appropriate allocation for you. When the function returns, Beta
will have the most likely estimate.

starting_pt is a vector of the appropriate size which indicates your
best initial guess for beta. if starting_pt=NULL, then (0,0,...0) will
be assumed.

step_size is the scale of the initial steps the maximization algorithm
will take. Currently, it is a scalar, so every dimension will have the
same step_size.

verbose is zero or one depending on whether you want to see the
maximizer's iterations.

For each function, the return value is the log likelihood (i.e.,
not the MLE, which is returned in beta, but the probability that the
MLE is true).

Sample usage:

gsl_vector * 	waring_parameter; 	//do not allocate.
double 		starting_pt[2] = {3, 0};
double		likelihood;
likelihood	= mle_waring(data, &zipf_parameter, starting_pt, .01, 0);
printf("Your most likely waring parameter is %g, with likelihood %g", 
				gsl_vector_get(waring_parameter, 0), likelihood);
gsl_vector_free waring_parameter; 	//Don't forget to clean up when you're done.

*/
