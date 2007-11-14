/** \file apop_exponential.c

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

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

static double exponential_log_likelihood(const apop_data *d, apop_model *p);




/* Let k be the rank, and x_k be the number of elements at that rank;
 then the mean rank (and therefore the most likely estimate for the
 exponential parameter) is sum(k * x_k)/sum(x) */
static apop_model * rank_exponential_estimate(apop_data * data, apop_model *parameters){
  double          colsum,
                  numerator   = 0,
                  grand_total = 0;
  apop_model 	*est= parameters ? parameters : apop_model_copy(apop_exponential);
  apop_model_clear(data, est);
  int             i;
    for(i=0; i< data->matrix->size2; i++){
        APOP_MATRIX_COL(data->matrix, i, v);
        colsum       = apop_sum(v);
        numerator   += colsum * (i+1);
        grand_total += colsum;
    }
	gsl_vector_set(est->parameters->vector, 0, numerator/grand_total);
    est->llikelihood	= apop_exponential.log_likelihood(data, parameters);
	/*if (est->uses.covariance)
		apop_numerical_covariance_matrix(apop_exponential_rank, est, data);*/
	return est;
}

static double rank_exponential_log_likelihood(const apop_data *d, apop_model *params){
  if (!params->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  double          b		    = gsl_vector_get(params->parameters->vector, 0),
		          p,
		          llikelihood = 0,
		          ln_b		= log(b);
  gsl_matrix	  *data		= d->matrix;
  int 		      k;
  gsl_vector      v;
	for (k=0; k< data->size2; k++){
		p	            = -ln_b - (k+1)/b;
        v               = gsl_matrix_column(data, k).vector;
		llikelihood    += apop_sum(&v) * p; 
	}
	return llikelihood;
}

static void rank_exponential_dlog_likelihood(const apop_data *d, gsl_vector *gradient, apop_model *params){
  if (!params->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  double		    bb		        = gsl_vector_get(params->parameters->vector, 0);
  int 		    k;
  gsl_matrix	    *data		    = d->matrix;
  double 		    d_likelihood 	= 0,
		          one_over_ln_b	    = 1/log(bb),
		          p;
  gsl_vector      v;
	for (k=0; k< data->size2; k++) {
		p	            = (one_over_ln_b -(k+1)) /bb;
        v               = gsl_matrix_column(data, k).vector;
		d_likelihood   += apop_sum(&v) * p; 
    }
	gsl_vector_set(gradient,0, d_likelihood);
}






static apop_model * exponential_estimate(apop_data * data,  apop_model *m){
  if (m->model_settings && (!strcmp((char *)m->model_settings, "r") || !strcmp((char *)m->model_settings, "R")))
      return rank_exponential_estimate(data, m);
  apop_model 	*est= apop_model_copy(apop_exponential);
  apop_model_clear(data, est);
	gsl_vector_set(est->parameters->vector, 0, apop_matrix_mean(data->matrix));
    est->llikelihood	= exponential_log_likelihood(data, est);
	//if (est->ep.uses.covariance)
		//apop_numerical_var_covar_matrix(apop_exponential, est, data->matrix);
	return est;
}

static double beta_greater_than_x_constraint(const apop_data *data, apop_model *v){
    //constraint is 0 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint){ 
        constraint = apop_data_calloc(1,1,1);
        apop_data_set(constraint, 0, 0, 1);
        }
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
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
static double exponential_log_likelihood(const apop_data *d, apop_model *p){
  if (!p->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  if (p->model_settings && (!strcmp((char *)p->model_settings, "r") || !strcmp((char *)p->model_settings, "R")))
      return rank_exponential_log_likelihood(d, p);
  gsl_matrix	*data	= d->matrix;
  double		mu		= gsl_vector_get(p->parameters->vector, 0),
		llikelihood;
	llikelihood	 = -apop_matrix_sum(data)/ mu;
	llikelihood	-= data->size1 * data->size2 * log(mu);
	return llikelihood;
}

static double exp_p(const apop_data *d, apop_model *p){
    return exp(exponential_log_likelihood(d, p));
}

/** The exponential distribution. A one-parameter likelihood fn.

\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>
\todo Check that the borderline work here is correct too.
*/
static void exponential_dlog_likelihood(const apop_data *d, gsl_vector *gradient, apop_model *p){
  if (!p->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  if (p->model_settings && (!strcmp((char *)p->model_settings, "r") || !strcmp((char *)p->model_settings, "R")))
      return rank_exponential_dlog_likelihood(d, gradient, p);
  double		mu	    = gsl_vector_get(p->parameters->vector, 0);
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
static void exponential_rng(double *out, gsl_rng* r, apop_model *p){
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
