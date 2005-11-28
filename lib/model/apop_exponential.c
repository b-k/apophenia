/** \file apop_exponential.c

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/

//The default list. Probably don't need them all.
#include "name.h"
#include "stats.h"
#include "model.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <stdio.h>
#include <assert.h>

static double exponential_log_likelihood(const gsl_vector *beta, void *d);

static double keep_away(double value, double limit,  double base){
	return (50000+fabs(value - limit)) * base;
}

static apop_estimate * exponential_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
	apop_inventory_filter(uses, apop_exponential.inventory_filter);
apop_estimate 	*est	= apop_estimate_alloc(data->size1,1,NULL, *uses);
	gsl_vector_set(est->parameters, 0, apop_matrix_mean(data));
	if (est->uses.log_likelihood)
		est->log_likelihood	= exponential_log_likelihood(est->parameters, data);
	if (est->uses.covariance)
		apop_numerical_var_covar_matrix(apop_exponential, est, data);
	return est;
}

static double beta_greater_than_x_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
double  limit       = 0,
        tolerance   = 1e-1;
double  mu          = gsl_vector_get(beta, 0);
    if (mu > limit) 
        return 0;
    //else:
    gsl_vector_memcpy(returned_beta, beta);
    gsl_vector_set(returned_beta, 0, limit + tolerance);
    return limit - mu;    
}

/** The exponential distribution. A one-parameter likelihood fn.
\f$Z(\mu,k) 		= \sum_k 1/\mu e^{-k/\mu} 			\f$ <br>
\f$ln Z(\mu,k) 		= \sum_k -\ln(\mu) - k/\mu			\f$ <br>
\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>

Some folks write the function as:
\f$Z(C,k) dx = \ln C C^{-k}. \f$
If you prefer this form, just convert your parameter via \f$\mu = {1\over
\ln C}\f$ (and convert back from the parameters this function gives you
via \f$C=\exp(1/\mu)\f$.

\ingroup likelihood_fns
\todo Set up an exponential object which makes use of the GSL.
\todo Check that the borderline work here is correct.
*/
static double exponential_log_likelihood(const gsl_vector *beta, void *d){
gsl_matrix	*data	= d;
double		mu		= gsl_vector_get(beta, 0),
		llikelihood;
	llikelihood	 = -apop_matrix_sum(data);
	llikelihood	/= mu;
	llikelihood	-= data->size1 * data->size2 * log(mu);
	return llikelihood;
}

/*
static double exponential_log_likelihood(const gsl_vector *beta, void *d){
double		mu		= gsl_vector_get(beta, 0),
		llikelihood;
static double	ka		= 0;
gsl_matrix	*data		= d;
	if (mu <= 0) {		//run away
		if (ka ==0){
			gsl_vector *	b_ka	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka, 0, GSL_DBL_EPSILON);
		 	ka	= exponential_log_likelihood(b_ka, d);
			gsl_vector_free (b_ka);
		}
		return keep_away(mu, 0, ka);
	}//else:
	llikelihood	 = -apop_matrix_sum(data);
	llikelihood	/= mu;
	llikelihood	-= data->size1 * data->size2 * log(mu);
	return llikelihood;
}
*/


/** The exponential distribution. A one-parameter likelihood fn.

\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>
\todo Check that the borderline work here is correct too.
*/
static void exponential_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double		mu		= gsl_vector_get(beta, 0);
static double	dka		= 0;
gsl_matrix	*data		= d;
double 		d_likelihood;
	if (mu <= 0) {		//keep away
		if (dka ==0){
			gsl_vector 	*b_ka	= gsl_vector_alloc(1);
			gsl_vector 	*b_kg	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, GSL_DBL_EPSILON);
		 	exponential_dlog_likelihood(b_ka , d, b_kg);
			dka	= gsl_vector_get(b_kg, 0);
			gsl_vector_free (b_ka);
			gsl_vector_free (b_kg);
		}
		gsl_vector_set(gradient,0, -keep_away(mu, 1, dka));
		return;
	}			//else:
	d_likelihood	 = apop_matrix_sum(data);
	d_likelihood	/= gsl_pow_2(mu);
	d_likelihood	-= data->size1 * data->size2 /mu;
	gsl_vector_set(gradient,0, d_likelihood);
}

/* Just a wrapper for gsl_ran_exponential.

   cut & pasted from the GSL documentation:
\f$p(x) dx = {1 \over \mu} \exp(-x/\mu) dx \f$

See the notes for \ref apop_exponential on a popular alternate form.
*/
static double exponential_rng(gsl_rng* r, double * a){
	//This fn exists because the GSL requires a double, 
	//while the apop_model structure requires a double*. 
	return gsl_ran_exponential(r, *a);
}


/** The exponential distribution. A one-parameter likelihood fn.

Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.


\f$Z(\mu,k) 	= \sum_k 1/\mu e^{-k/\mu} 			\f$ <br>
\f$ln Z(\mu,k) 	= \sum_k -\ln(\mu) - k/\mu			\f$ <br>
\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>

Some write the function as:
\f$Z(C,k) dx = \ln C C^{-k}. \f$
If you prefer this form, just convert your parameter via \f$\mu = {1\over
\ln C}\f$ (and convert back from the parameters this function gives you
via \f$C=\exp(1/\mu)\f$.

\ingroup models
\todo Check that the borderline work here is correct.
\todo Write a second object for the plain old not-network data Exponential.
*/
apop_model apop_exponential = {"Exponential", 1,
{	1,	//parameters
	1,	//covariance
	1,	//confidence
	0,	//predicted
	0,	//residuals
	1,	//log_likelihood
	1	//names;
},
	 exponential_estimate, exponential_log_likelihood, exponential_dlog_likelihood, NULL, beta_greater_than_x_constraint, exponential_rng};


