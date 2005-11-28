/** \file apop_gamma.c

  The gamma distribution.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/

#include "model.h"


//The default list. Probably don't need them all.
#include "name.h"
#include "bootstrap.h"
#include "regression.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>

static double keep_away(double value, double limit,  double base){
    return (50000+fabs(value - limit)) * base;
}


static apop_estimate * gamma_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
    apop_inventory_filter(uses, apop_gamma.inventory_filter);
    return apop_maximum_likelihood(data, uses, apop_gamma, *(apop_estimation_params *)parameters);
}

static double beta_zero_and_one_greater_than_x_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
double  limit0      = 0,
        limit1      = 0,
        tolerance   = 1e-1;
double  beta0       = gsl_vector_get(beta, 0),
        beta1       = gsl_vector_get(beta, 1);
    if (beta0 > limit0 && beta1 > limit1) 
        return 0;
    //else:
    gsl_vector_memcpy(returned_beta, beta);
    gsl_vector_set(returned_beta, 0, GSL_MAX(limit0 + tolerance, beta0));
    gsl_vector_set(returned_beta, 1, GSL_MAX(limit1 + tolerance, beta1));
    return GSL_MAX(limit0 - beta0, 0) + GSL_MAX(limit1 - beta1, 0);    
}

static double gamma_log_likelihood(const gsl_vector *beta, void *d){
float         a    = gsl_vector_get(beta, 0),
              b    = gsl_vector_get(beta, 1);
int           i, k;
gsl_matrix    *data        = d;
float         llikelihood  = 0,
        ln_ga 	= gsl_sf_lngamma(a),
        ln_b	= log(b),
        a_ln_b	= a * ln_b,
        x;
    for (i=0; i< data->size1; i++)
        for (k=0; k< data->size2; k++){
            x	= gsl_matrix_get(data, i, k);
            if (x!=0)
                llikelihood    += -ln_ga - a_ln_b + (a-1) * log(x) - x/b;
        }
    return llikelihood;
}

/** The derivative of the Gamma distribution, for use in likelihood
 * minimization. You'll probably never need to call this directly.*/
static void gamma_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
float       	a    	= gsl_vector_get(beta, 0),
        		b    	= gsl_vector_get(beta, 1);
    //if (a <= 0 || b <= 0 || gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;    
                        //a sign to the minimizer to look elsewhere.
int             i, k;
gsl_matrix      *data	= d;
float           d_a     = 0,
        d_b		= 0,
        psi_a	= gsl_sf_psi(a),
        ln_b    = log(b),
        x;
    for (i=0; i< data->size1; i++)
        for (k=0; k< data->size2; k++){
            x	= gsl_matrix_get(data, i, k);
            if (x!=0){
                d_a    += -psi_a - ln_b + log(x);
                d_b    += -a/b - x;
            }
        }
    gsl_vector_set(gradient,0, d_a);
    gsl_vector_set(gradient,1, d_b);
}

/* Just a wrapper for gsl_ran_gamma.

   cut & pasted from the GSL documentation:
\f$          p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx \f$

See the notes for \ref apop_exponential on a popular alternate form.
*/
static double gamma_rng(gsl_rng* r, double * a){
    //This fn exists because the GSL requires a double, 
    //while the apop_model structure requires a double*. 
    return gsl_ran_gamma(r, a[0], a[1]);
}


/** The Gamma distribution

Location of data in the grid is not relevant; send it a 1 x N, N x 1, or N x M and it will all be the same.

apop_gamma.estimate() is an MLE, so feed it appropriate \ref apop_estimation_params.
  
If you have frequency or ranking data, you probably mean to be using \ref apop_gamma_rank.

\f$G(x, a, b)     = 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da    =  -\psi(a) - ln b + ln(x) \f$    (also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db    =  -a/b - x \f$
\ingroup models
*/
apop_model apop_gamma = {"Gamma", 2,
{
    1,    //parameters
    1,    //covariance
    1,    //confidence
    0,    //predicted
    0,    //residuals
    1,    //log_likelihood
    1    //names;
},
     gamma_estimate, gamma_log_likelihood, gamma_dlog_likelihood, NULL,  beta_zero_and_one_greater_than_x_constraint, gamma_rng};
