//stats.h			  	Copyright 2005 by Ben Klemens. Licensed under the GNU GPL.
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h>
#include "linear_algebra.h"


	//The following are just convenient hooks to gsl vector functions.
	//var_m lets you input a mean if you've already calculated it, saving
	//some repetition.
inline double apop_mean(gsl_vector *in);
inline double apop_var(gsl_vector *in);
inline double apop_var_m(gsl_vector *in, double mean);
inline double apop_covar(gsl_vector *ina, gsl_vector *inb);
inline double apop_cov(gsl_vector *ina, gsl_vector *inb);
inline double apop_correlation(gsl_vector *ina, gsl_vector *inb);
inline double apop_kurtosis(gsl_vector *in);
inline double apop_kurt(gsl_vector *in);


void apop_normalize_vector(gsl_vector *in, gsl_vector **out, int in_place, int normalization_type);
/*
	normalization_type:
		1: out will have mean = 0, std deviation = 1 [Uses sample variance (stuff/N-1, not stuff/N).]
		2: out will have min = 0, max = 1;
	in_place:
		0: out will be allocated and filled, in will be	unchanged.
		1: in will be normalized in place. 
	call: 
		gsl_vector  *unnormed, *normed;
		[allocate and fill unnormed with data; do not allocate normed.]
		normalize_vector(unnormed, &normed, 1, 1);
		[or:]
		normalize_vector(unnormed, NULL, 0, 1);
*/

void apop_normalize_matrix(gsl_matrix *data);
	//Regression methods often require that the mean of each colum of the data matrix have mean zero.

inline double apop_test_chi_squared_var_not_zero(gsl_vector *in);
	//As described: give it a vector, and it'll tell you the confidence 
	//with which you can say that the vector is not zero.

inline double apop_double_abs(double a);
	//This has to exist somewhere...

double apop_random_beta(double m, double v, gsl_rng *r) ;
	/*Give me mean m and variance v, and I'll give you
	 * n draws from the appropriate beta dist.
	 * remember: 0<m<1, and v is tiny (<<1/12). You get NaNs if no
	 * appropriate distribution exists.*/

double apop_multivariate_normal_prob(gsl_vector *x, gsl_vector* mu, gsl_matrix* sigma, int first_use);
	//Evaluate a multivariate normal(mu, sigma) at the point x.
//The equation:
//	exp(-1/2 (X-mu)' sigma^-1 (x-mu))
//	--------------------------
//	sqrt((2 Pi)^n det(sigma))
//
//The inverse and determinant are expensive, so keep them around where possible: on the first call, set 
//first_use to 1, then feed in as many new values of X as you want.

double apop_random_double(double min, double max, gsl_rng *r);

//produce a 101-element vector of percentiles.
double * apop_percentiles(gsl_vector *data, char rounding);

long double apop_matrix_sum(gsl_matrix *m);
double apop_matrix_mean(gsl_matrix *data);
double apop_matrix_var_m(gsl_matrix *data, double mean);
void apop_matrix_mean_and_var(gsl_matrix *data, double *mean, double *var);
