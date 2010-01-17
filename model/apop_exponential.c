/** \file apop_exponential.c 
        The Exponential distribution.
        */
/* Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "settings.h"
#include "likelihoods.h"

/** If this settings group is present, models that can take rank data
  will read the input data as such.  Allocation is thus very simple, e.g.
  \code
  Apop_model_group_add(your_model, apop_rank);
  \endcode */
apop_rank_settings *apop_rank_settings_init(apop_rank_settings in){
    apop_rank_settings *out = malloc(sizeof(apop_rank_settings));
    out->rank_data = 'r';
    return out;
}

apop_rank_settings *apop_rank_settings_alloc(void *ignoreme){
    return apop_rank_settings_init((apop_rank_settings){}); }

void apop_rank_settings_free(apop_rank_settings *in){free(in);}

void *apop_rank_settings_copy(apop_rank_settings *in){
    apop_rank_settings *out = apop_rank_settings_alloc(NULL);
    out->rank_data = in->rank_data;
    return out;
}

static double exponential_log_likelihood(apop_data *d, apop_model *p);
static double rank_exponential_log_likelihood(const apop_data *d, apop_model *params);

/* Let k be the rank, and x_k be the number of elements at that rank;
 then the mean rank (and therefore the most likely estimate for the
 exponential parameter) is sum(k * x_k)/sum(x) */
static apop_model * rank_exponential_estimate(apop_data * data, apop_model *est){
  double          colsum,
                  numerator   = 0,
                  grand_total = 0;
    for(size_t i=0; i< data->matrix->size2; i++){
        APOP_MATRIX_COL(data->matrix, i, v);
        colsum       = apop_sum(v);
        numerator   += colsum * (i+1);
        grand_total += colsum;
    }
	gsl_vector_set(est->parameters->vector, 0, numerator/grand_total);
    est->llikelihood	= rank_exponential_log_likelihood(data, est);
	return est;
}

static double rank_exponential_log_likelihood(const apop_data *d, apop_model *params){
  apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  double          b		    = gsl_vector_get(params->parameters->vector, 0),
		          p,
		          llikelihood = 0,
		          ln_b		= log(b);
  gsl_matrix	  *data		= d->matrix;
  gsl_vector      v;
	for (size_t k=0; k< data->size2; k++){
		p	            = -ln_b - (k+1)/b;
        v               = gsl_matrix_column(data, k).vector;
		llikelihood    += apop_sum(&v) * p; 
	}
	return llikelihood;
}

static void rank_exponential_dlog_likelihood(const apop_data *d, gsl_vector *gradient, apop_model *params){
  apop_assert_void(params->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
  double		    bb		        = gsl_vector_get(params->parameters->vector, 0);
  gsl_matrix	    *data		    = d->matrix;
  double 		    d_likelihood 	= 0,
		          one_over_ln_b	    = 1/log(bb),
		          p;
  gsl_vector      v;
	for (size_t k=0; k< data->size2; k++) {
		p	            = (one_over_ln_b -(k+1)) /bb;
        v               = gsl_matrix_column(data, k).vector;
		d_likelihood   += apop_sum(&v) * p; 
    }
	gsl_vector_set(gradient,0, d_likelihood);
}

static apop_model * exponential_estimate(apop_data * data,  apop_model *est){
    if (apop_settings_get_group(est, "apop_rank"))
      return rank_exponential_estimate(data, est);
	gsl_vector_set(est->parameters->vector, 0, apop_matrix_mean(data->matrix));
    est->llikelihood	= exponential_log_likelihood(data, est);
	return est;
}

static double beta_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1
    return apop_linear_constraint(v->parameters->vector, .margin = 1e-3);
}

static double exponential_log_likelihood(apop_data *d, apop_model *p){
  apop_assert(p->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(p, "apop_rank"))
      return rank_exponential_log_likelihood(d, p);
  gsl_matrix	*data	= d->matrix;
  double		mu		= gsl_vector_get(p->parameters->vector, 0),
		llikelihood;
	llikelihood	 = -apop_matrix_sum(data)/ mu;
	llikelihood	-= data->size1 * data->size2 * log(mu);
	return llikelihood;
}

/* The exponential distribution. A one-parameter likelihood fn.

\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>
*/
static void exponential_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  apop_assert_void(p->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(p, "apop_rank"))
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

apop_model apop_exponential = {"Exponential distribution", 1,0,0,
	 .estimate = exponential_estimate, .log_likelihood = exponential_log_likelihood, 
     .score = exponential_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
     .draw = exponential_rng};
