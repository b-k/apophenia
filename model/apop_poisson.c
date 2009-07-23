/** \file apop_poisson.c

  The poisson distribution.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

//The default list. Probably don't need them all.
#include "model.h"
#include "stats.h"
#include "types.h"
#include "mapply.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>

static double poisson_log_likelihood(apop_data *d, apop_model *);

static apop_model * poisson_estimate(apop_data * data,  apop_model *parameters){
  apop_model 	*est= apop_model_copy(parameters ? *parameters : apop_poisson);
  double		mean    = apop_matrix_mean(data->matrix);
    if (!est->parameters) 
        est->parameters   = apop_data_alloc(1,0,0);
	gsl_vector_set(est->parameters->vector, 0, mean);
    est->llikelihood	= poisson_log_likelihood(data, est);
    //to prevent an infinite loop, the jackknife needs to be flagged to
    //not run itself. We free-ride of the apop_ls_settings struct to signal.
    apop_ls_settings *dummy = apop_settings_get_group(est, "apop_ls");
    if (!dummy || dummy->want_cov){
        if (!dummy)
            Apop_settings_add_group(est, apop_ls, data);
        Apop_settings_add(est, apop_ls, want_cov, 0);
        est->covariance = apop_jackknife_cov(data, *est);
    }
	return est;
}

static double beta_zero_greater_than_x_constraint(apop_data *returned_beta, apop_model *v){
    //constraint is 0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_calloc(1,1,1);
        apop_data_set(constraint, 0, 0, 1);
    }
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
}

static double ln_l;

static double apply_me(gsl_vector *v){
  int       k;
  double    x,
            llikelihood = 0;
    for (k=0; k< v->size; k++){
        x	= gsl_vector_get(v, k);
        if (x!=0)
            llikelihood    += ln_l *x - gsl_sf_lngamma(x+1);
    }
    return llikelihood;
}

static double poisson_log_likelihood(apop_data *d, apop_model * p){
  apop_assert(p->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  double        lambda      = gsl_vector_get(p->parameters->vector, 0);
    ln_l 	= log(lambda);
  gsl_vector *  v           = apop_matrix_map(d->matrix, apply_me);
  double        llikelihood = apop_vector_sum(v);
    gsl_vector_free(v);
    return llikelihood - d->matrix->size1*d->matrix->size2*lambda;
}

static void poisson_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  apop_assert_void(p->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
  double       	lambda  = gsl_vector_get(p->parameters->vector, 0);
  gsl_matrix      *data	= d->matrix;
  float           d_a;
    d_a  = apop_matrix_sum(data)/lambda;
    d_a -= data->size1*data->size2;
    gsl_vector_set(gradient,0, d_a);
}

/* Just a wrapper for gsl_ran_poisson.

   cut & pasted from the GSL documentation:
\f$       
p(k) = {\mu^k \over k!} \exp(-\mu), \f$

where \f$k\geq 0\f$.
*/
static void poisson_rng(double *out, gsl_rng* r, apop_model *p){
    *out = gsl_ran_poisson(r, *p->parameters->vector->data);
}


/** The poisson distribution

Location of data in the grid is not relevant; send it a 1 x N, N x 1, or N x M and it will all be the same.

\f$p(k) = {\mu^k \over k!} \exp(-\mu), \f$

\ingroup models
*/
apop_model apop_poisson = {"poisson", 1, 0,0, 
     .estimate = poisson_estimate, .log_likelihood = poisson_log_likelihood, 
     .score = poisson_dlog_likelihood, .constraint = beta_zero_greater_than_x_constraint, 
     .draw = poisson_rng};
