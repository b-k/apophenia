/** \file apop_normal.c

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/

//The default list. Probably don't need them all.
#include "types.h"
#include "bootstrap.h"
#include "regression.h"
#include "conversions.h"
#include "likelihoods.h"
#include <apophenia/model.h>
#include <apophenia/stats.h>
#include <apophenia/linear_algebra.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <assert.h>
static double normal_log_likelihood(const gsl_vector *beta, apop_data *d);


//////////////////
//The Normal (gaussian) distribution
//////////////////

/** The normal estimate
\todo Get off my ass and check the closed-form var-covar matrix, instead of using the inverse hessian. */
static apop_estimate * normal_estimate(apop_data * data, void *parameters){
apop_estimate 	*est	    = apop_estimate_alloc(data,apop_normal, parameters);
double		mean, var;
	apop_matrix_mean_and_var(data->matrix, &mean, &var);	
	gsl_vector_set(est->parameters->vector, 0, mean);
	gsl_vector_set(est->parameters->vector, 1, var);
	if (est->estimation_params.uses.log_likelihood)
		est->log_likelihood	= normal_log_likelihood(est->parameters->vector, data->matrix);
	if (est->estimation_params.uses.covariance)
		apop_numerical_var_covar_matrix(apop_normal, est, data->matrix);
	return est;
}


static double beta_1_greater_than_x_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
double  limit       = 0,
        tolerance   = 1e-1;
double  mu          = gsl_vector_get(beta, 1);
    if (mu > limit) 
        return 0;
    //else:
    gsl_vector_memcpy(returned_beta, beta);
    gsl_vector_set(returned_beta, 1, limit + tolerance);
    return limit - mu;    
}

/* The log likelihood function for the Normal.

The log likelihood function and dlog likelihood don't care about your
rows of data; if you have an 8 x 7 data set, it will give you the log
likelihood of those 56 observations given the mean and variance (i.e.,
\f$\sigma^2\f$, not std deviation=\f$\sigma\f$) you provide.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-(x-\mu)^2 / 2\sigma^2)\f$
\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\param beta	beta[0]=the mean; beta[1]=the variance
\param d	the set of data points; see notes.
*/
static double normal_log_likelihood(const gsl_vector *beta, apop_data *d){
double 		mu	= gsl_vector_get(beta,0),
	    	ss	= 2 * gsl_vector_get(beta,1),
	    	ll	= 0;
size_t		i,j;
gsl_matrix	*data	= d->matrix;
	for (i=0;i< data->size1; i++)
		for (j=0;j< data->size2; j++)
			ll	-= gsl_pow_2(gsl_matrix_get(data, i, j) - mu) / ss;
			//Or use the stock gsl function (but then just return -ll).
			//ll	+= log(gsl_ran_gaussian_pdf((gsl_matrix_get(data, i, j) - mu), gsl_vector_get(beta,1)));
	ll 	-= data->size1 * data->size2 * log(ss * M_PI)/2.0; //second term is positive because it subtracts -log.
	return ll;
}

/** Gradient of the log likelihood function

To tell you the truth, I have no idea when anybody would need this, but it's here for completeness. 
\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$
\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$
\todo Add constraint that \f$\sigma^2>0\f$.
 */
static void normal_dlog_likelihood(const gsl_vector *beta, apop_data *d, gsl_vector *gradient){
double 		mu	    = gsl_vector_get(beta,0),
	    	ss	    = gsl_vector_get(beta,1),
	    	dll	    = 0,
	    	sll	    = 0,
	    	x;
int	    	i,j;
gsl_matrix	*data	= d->matrix;
	for (i=0;i< data->size1; i++)
		for (j=0;j< data->size2; j++){
			x	 = gsl_matrix_get(data, i, j);
			dll	+= (x - mu);
			sll	+= gsl_pow_2(x - mu);
		}
	gsl_vector_set(gradient, 0, dll/ss);
	gsl_vector_set(gradient, 1, sll/(2*gsl_pow_2(ss))+ data->size1 * data->size2 * 0.5/ss);
}

/** An apophenia wrapper for the GSL's Normal RNG.

Two differences: this one asks explicitly for a mean, and the GSL
assumes zero and makes you add the mean yourself; Apophenia tends to
prefer the variance (\f$\sigma^2\f$) wherever possible, while the GSL
uses the standard deviation here (\f$\sigma\f$)

\param r	a gsl_rng already allocated
\param a	the mean and the variance
 */
static double normal_rng(gsl_rng *r, double *a){
	//return gsl_ran_gaussian(r, sqrt(a[1])) + a[0];
	return gsl_ran_gaussian(r, sqrt(a[1])) + a[0];
}

/** You know it, it's your attractor in the limit, it's the Gaussian distribution.

Where possible, Apophenia tries to report the variance, \f$\sigma^2\f$,
not the standard deviation, \f$\sigma\f$. So the first parameter for
this model is the mean, and the second is the variance.

The log likelihood function and dlog likelihood don't care about your
rows of data; if you have an 8 x 7 data set, it will give you the log
likelihood of those 56 observations given the mean and variance you provide.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2)\f$

\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$

\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$
\ingroup models
*/
apop_model apop_normal = {"Normal", 2, 
 normal_estimate, normal_log_likelihood, normal_dlog_likelihood, NULL,beta_1_greater_than_x_constraint, normal_rng};

/** This is a synonym for \ref apop_normal, q.v.
\ingroup models
*/
apop_model apop_gaussian = {"Gaussian", 2, 
 normal_estimate, normal_log_likelihood, normal_dlog_likelihood, NULL,beta_1_greater_than_x_constraint, normal_rng};
