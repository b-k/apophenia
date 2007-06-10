/** \file apop_normal.c

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

//The default list. Probably don't need them all.
#include "types.h"
#include "mapply.h"
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
static double normal_log_likelihood(const apop_data *beta, apop_data *d, apop_model *params);


//////////////////
//The Normal (gaussian) distribution
//////////////////

/** The normal estimate */
static apop_model * normal_estimate(apop_data * data, apop_model *parameters){
  double		mean, var;
  apop_model 	*est	    = apop_model_copy(*parameters);
  apop_OLS_params *p;
    if (!parameters) {
        p = apop_OLS_params_alloc(data, apop_normal);
        est   = p->model;
    } else p = parameters->model_params;
	apop_matrix_mean_and_var(data->matrix, &mean, &var);	
    if (!est->parameters)
        est->parameters = apop_data_alloc(2, 0, 0);
	gsl_vector_set(est->parameters->vector, 0, mean);
	gsl_vector_set(est->parameters->vector, 1, sqrt(var));
    est->llikelihood	= normal_log_likelihood(est->parameters, data, NULL);
	if (!p || p->want_cov){
        est->covariance   = apop_data_calloc(0, 2, 2);
        int ct  = data->matrix->size1 * data->matrix->size2;
        apop_data_set(est->covariance, 0, 0, mean/ct);
        apop_data_set(est->covariance, 1, 1, 2*gsl_pow_2(var)/(ct-1));
    }
	return est;
}

static double beta_1_greater_than_x_constraint(const apop_data *beta, apop_data *returned_beta, apop_model *v){
    //constraint is 0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint) {
        constraint = apop_data_calloc(1,1,2);
        apop_data_set(constraint, 0, 1, 1);
    }
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}

static double   mu, sd;

static double apply_me(gsl_vector *v){
  int           i;
  long double    ll  = 1;
    for(i=0; i< v->size; i++)
	    ll	*= gsl_ran_gaussian_pdf((gsl_vector_get(v, i) - mu), sd);
    return ll;
}

//This just takes the sum of (x-mu)^2. Using gsl_ran_gaussian_pdf
//would be to calculate log(exp((x-mu)^2)) == slow.
static double apply_me2(gsl_vector *v){
  int           i;
  long double    ll  = 0;
    for(i=0; i< v->size; i++)
	    ll	+= gsl_pow_2(gsl_vector_get(v, i) - mu);
    return ll;
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
static double normal_log_likelihood(const apop_data *beta, apop_data *d, apop_model *params){
    mu	        = gsl_vector_get(beta->vector,0);
    sd          = gsl_vector_get(beta->vector,1);
  gsl_vector *  v       = apop_matrix_map(d->matrix, apply_me2);//sum of (x-mu)^2
  long double   ll      = -apop_vector_sum(v)/(2*gsl_pow_2(sd)) - d->matrix->size1*d->matrix->size2*(M_LNPI+M_LN2+log(sd));
    gsl_vector_free(v);
	return ll;
}

static double normal_p(const apop_data *beta, apop_data *d, apop_model *params){
    mu	        = gsl_vector_get(beta->vector,0);
    sd          = gsl_vector_get(beta->vector,1);
  gsl_vector *  v       = apop_matrix_map(d->matrix, apply_me);
  int           i;
  long double   ll      = 1;
    for(i=0; i< v->size; i++)
        ll  *= gsl_vector_get(v, i);
    gsl_vector_free(v);
	return ll;
}

/** Gradient of the log likelihood function

To tell you the truth, I have no idea when anybody would need this, but it's here for completeness.
\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$
\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$
\f$d\ln N(\mu,\sigma)/d\sigma = ((x-\mu)^2 / \sigma^3) - 1/\sigma \f$
 */
static void normal_dlog_likelihood(const apop_data *beta, apop_data *d, 
                                    gsl_vector *gradient, apop_model *params){    
              mu      = gsl_vector_get(beta->vector,0);
  double      sd      = gsl_vector_get(beta->vector,1),
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


/** An apophenia wrapper for the GSL's Normal RNG.

Two differences: this one asks explicitly for a mean, and the GSL
assumes zero and makes you add the mean yourself; Apophenia tends to
prefer the variance (\f$\sigma^2\f$) wherever possible, while the GSL
uses the standard deviation here (\f$\sigma\f$)

\param r	a gsl_rng already allocated
\param a	the mean and the variance
 */
static void normal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0];
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
apop_model apop_normal = {"Normal", 2, 0, 0,
 .estimate = normal_estimate, .p = normal_p, .log_likelihood = normal_log_likelihood, .score = normal_dlog_likelihood, 
 .constraint = beta_1_greater_than_x_constraint, .draw = normal_rng};

/** This is a synonym for \ref apop_normal, q.v.
\ingroup models
*/
apop_model apop_gaussian = {"Gaussian", 2,0,0,
 .estimate = normal_estimate, .p = normal_p, .log_likelihood = normal_log_likelihood, 
 .score = normal_dlog_likelihood, .constraint = beta_1_greater_than_x_constraint, .draw = normal_rng};
