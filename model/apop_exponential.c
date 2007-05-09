/** \file apop_exponential.c

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/

//The default list. Probably don't need them all.
#include "types.h"
#include "stats.h"
#include "model.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

static double exponential_log_likelihood(const apop_data *beta, apop_data *d, apop_params *p);

static apop_params * exponential_estimate(apop_data * data,  apop_params *parameters){
apop_params 	*est	    = apop_params_alloc(data, apop_exponential, parameters, NULL);
	gsl_vector_set(est->parameters->vector, 0, apop_matrix_mean(data->matrix));
    est->log_likelihood	= exponential_log_likelihood(est->parameters, data, NULL);
	//if (est->ep.uses.covariance)
		//apop_numerical_var_covar_matrix(apop_exponential, est, data->matrix);
	return est;
}

static double beta_greater_than_x_constraint(const apop_data *beta, apop_data *returned_beta, apop_params *v){
    //constraint is 0 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint) constraint = apop_data_calloc(1,1,1);
    apop_data_set(constraint, 0, 0, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
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
static double exponential_log_likelihood(const apop_data *beta, apop_data *d, apop_params *p){
  gsl_matrix	*data	= d->matrix;
  double		mu		= gsl_vector_get(beta->vector, 0),
		llikelihood;
	llikelihood	 = -apop_matrix_sum(data)/ mu;
	llikelihood	-= data->size1 * data->size2 * log(mu);
	return llikelihood;
}

static double exp_p(const apop_data *beta, apop_data *d, apop_params *p){
    return exp(exponential_log_likelihood(beta, d, p));
}

/** The exponential distribution. A one-parameter likelihood fn.

\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>
\todo Check that the borderline work here is correct too.
*/
static void exponential_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_params *p){
  double		mu	    = gsl_vector_get(beta->vector, 0);
  gsl_matrix	*data	= d->matrix;
  double 		d_likelihood;
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
static void exponential_rng(double *out, gsl_rng* r, apop_params *p){
	*out = gsl_ran_exponential(r, p->parameters->vector->data[0]);
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
apop_model apop_exponential = {"Exponential distribution", 1,0,0,
	 .estimate = exponential_estimate, .p = exp_p, .log_likelihood = exponential_log_likelihood, 
     .score = exponential_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
     .draw = exponential_rng};


