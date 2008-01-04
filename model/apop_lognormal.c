/** \file apop_normal.c

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

//The default list. Probably don't need them all.
#include "types.h"
#include "mapply.h"
#include "conversions.h"
#include "likelihoods.h"
#include <apophenia/model.h>
#include <apophenia/stats.h>
#include <apophenia/linear_algebra.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <assert.h>
static double lognormal_log_likelihood(apop_data *d, apop_model *params);



/** The normal estimate */
static apop_model * lognormal_estimate(apop_data * data, apop_model *parameters){
  apop_model 	*est;
  apop_ls_settings *p;
  double   mean    = 0,
           var     = 0; 
    if (!parameters) {
        p   = apop_ls_settings_alloc(data, apop_normal);
        est = p->model;
    } else {
        p   = parameters->model_settings;
        est	= apop_model_copy(*parameters);
    }
    apop_matrix_mean_and_var(data->matrix, &mean, &var);
    if (!est->parameters)
        est->parameters = apop_data_alloc(2, 0, 0);
    double sigsq   = log(1+ var/gsl_pow_2(mean));
	gsl_vector_set(est->parameters->vector, 0, log(mean)- sigsq/2);
	gsl_vector_set(est->parameters->vector, 1, sqrt(sigsq));
    est->llikelihood	= lognormal_log_likelihood(data, est);
	return est;
}

static double beta_1_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint) {
        constraint = apop_data_calloc(1,1,2);
        apop_data_set(constraint, 0, 1, 1);
    }
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-3);
}

static double   mu, sd;

static double apply_me2a(gsl_vector *v){
  int           i;
  long double    ll  = 0;
    for(i=0; i< v->size; i++)
	    ll	+= log(gsl_vector_get(v, i));
    return ll;
}

static double apply_me2b(gsl_vector *v){
  int           i;
  long double    ll  = 0;
    for(i=0; i< v->size; i++)
	    ll	+= gsl_pow_2(log(gsl_vector_get(v, i)) - mu);
    return ll;
}

/* The log likelihood function for the lognormal distribution.

\$f = exp(-(ln(x)-\mu)^2/(2\sigma^2))/ (x\sigma\sqrt{2\pi})\$
\$ln f = -(ln(x)-\mu)^2/(2\sigma^2) - ln(x) - ln(\sigma\sqrt{2\pi})\$

\param beta	beta[0]=the mean; beta[1]=the variance
\param d	the set of data points; see notes.
*/
static double lognormal_log_likelihood(apop_data *d, apop_model *params){
  apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    mu	        = gsl_vector_get(params->parameters->vector,0);
    sd          = gsl_vector_get(params->parameters->vector,1);
    gsl_vector *  v       = apop_matrix_map(d->matrix, apply_me2b);//sum of (ln(x)-mu)^2
    long double   ll      = -apop_vector_sum(v)/(2*gsl_pow_2(sd));
    gsl_vector *  v1      = apop_matrix_map(d->matrix, apply_me2a);//sum of ln(x)
                ll       -= apop_vector_sum(v1);
                ll       -= d->matrix->size1*d->matrix->size2*(M_LNPI+M_LN2+log(sd));
    gsl_vector_free(v);
    gsl_vector_free(v1);
	return ll;
}

static double lognormal_p(apop_data *d, apop_model *params){
  apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  int   i, j;
  long double ll    = 1;
    mu	= gsl_vector_get(params->parameters->vector,0);
    sd  = gsl_pow_2(gsl_vector_get(params->parameters->vector,1));//really, var.
    for (i=0; i< d->matrix->size1; i++)
        for (j=0; j< d->matrix->size1; j++){
            double x = gsl_matrix_get(d->matrix, i, j);
            ll      *=  exp(-gsl_pow_2(log(x)-mu)/(2*sd))/(x*sd * M_LNPI+M_LN2);
        }
	return ll;
}

/* This is copied from the Normal. The first one who needs it gets to
 * fix it. 
static void lognormal_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *params){    
              mu      = gsl_vector_get(params->parameters->vector,0);
  double      sd      = gsl_vector_get(params->parameters->vector,1),
              dll     = 0,
              sll     = 0,
              x;
  int         i,j;
  gsl_matrix  *data   = d->matrix;
    for (i=0;i< data->size1; i++)
        for (j=0;j< data->size2; j++){
            x    = gsl_matrix_get(data, i, j);
            dll += (x - mu);
            sll += gsl_pow_2(x - mu);
        }
    gsl_vector_set(gradient, 0, dll/gsl_pow_2(sd));
    //gsl_vector_set(gradient, 1, sll/(2*gsl_pow_2(ss))- data->size1 * data->size2 * 0.5/ss);
    gsl_vector_set(gradient, 1, sll/gsl_pow_3(sd)- data->size1 * data->size2 /sd);
}
*/


/** An Apophenia wrapper for the GSL's Normal RNG, logged.

\param r	a gsl_rng already allocated
\param *out	To where I will write the drawn number
\param *p   A pointer to the model.
 */
static void lognormal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = exp(gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0]);
}

/** You know it, it's your attractor in the limit, it's the Gaussian distribution.

  As is custom, the first parameter is the mean, the second is the standard deviation (i.e., the square root of the variance).

The log likelihood function and dlog likelihood don't care about your
rows of data; if you have an 8 x 7 data set, it will give you the log
likelihood of those 56 observations given the mean and variance you provide.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2)\f$

\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$

\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$
\ingroup models
*/
apop_model apop_lognormal = {"Lognormal distribution", 2, 0, 0,
 .estimate = lognormal_estimate, .p = lognormal_p, .log_likelihood = lognormal_log_likelihood, /*.score = lognormal_dlog_likelihood,*/ 
 .constraint = beta_1_greater_than_x_constraint, .draw = lognormal_rng};
