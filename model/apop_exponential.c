/** \file 
        The Exponential distribution.
        */
/* Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "settings.h"
#include "internal.h"
#include "likelihoods.h"

static double beta_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1
    return apop_linear_constraint(v->parameters->vector, .margin = 1e-3);
}

static double exponential_log_likelihood(apop_data *d, apop_model *p){
  Nullcheck_m(p); Nullcheck_p(p);
  gsl_matrix	*data	= d->matrix;
  double		mu		= gsl_vector_get(p->parameters->vector, 0),
		llikelihood;
	llikelihood	 = -apop_matrix_sum(data)/ mu;
	llikelihood	-= data->size1 * data->size2 * log(mu);
	return llikelihood;
}

static apop_model * exponential_estimate(apop_data * data,  apop_model *est){
	gsl_vector_set(est->parameters->vector, 0, apop_matrix_mean(data->matrix));
    apop_data_add_named_elmt(est->info, "log likelihood", exponential_log_likelihood(data, est));
	return est;
}

static double expo_cdf(apop_data *d, apop_model *params){
  Nullcheck_m(params) Nullcheck_p(params) Nullcheck_d(d) 
  Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double lambda = gsl_vector_get(params->parameters->vector, 0);
    return gsl_cdf_exponential_P(val, lambda);
}

/* The exponential distribution. A one-parameter likelihood fn.

\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>
*/
static void exponential_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  Nullcheck_m(p); Nullcheck_p(p);
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

apop_model apop_exponential = {"Exponential distribution", 1,0,0,.dsize=1,
	 .estimate = exponential_estimate, .log_likelihood = exponential_log_likelihood, 
     .score = exponential_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
     .draw = exponential_rng, .cdf = expo_cdf};
