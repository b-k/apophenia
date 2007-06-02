/** \file apop_gamma.c

  The gamma distribution.

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"


//The default list. Probably don't need them all.
#include "types.h"
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


static apop_model * gamma_estimate(apop_data * data, apop_model *parameters){
    return apop_maximum_likelihood(data, *parameters);
}

static double beta_zero_and_one_greater_than_x_constraint(const apop_data *beta, apop_data *returned_beta, apop_model *v){
    //constraint is 0 < beta_1 and 0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint) constraint = apop_data_calloc(2,2,2);
    apop_data_set(constraint, 0, 0, 1);
    apop_data_set(constraint, 1, 1, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}

static double gamma_log_likelihood(const apop_data *beta, apop_data *d, apop_model *p){
float         a    = gsl_vector_get(beta->vector, 0),
              b    = gsl_vector_get(beta->vector, 1);
int           i, k;
gsl_matrix    *data        = d->matrix;
double        llikelihood  = 0,
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

static double gamma_p(const apop_data *beta, apop_data *d, apop_model *p){
    return exp(gamma_log_likelihood(beta, d, p));
}

/** The derivative of the Gamma distribution, for use in likelihood
 * minimization. You'll probably never need to call this directly.*/
static void gamma_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_model *p){
float       	a    	= gsl_vector_get(beta->vector, 0),
        		b    	= gsl_vector_get(beta->vector, 1);
    //if (a <= 0 || b <= 0 || gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;    
                        //a sign to the minimizer to look elsewhere.
int             i, k;
gsl_matrix      *data	= d->matrix;
double          d_a     = 0,
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
    gsl_vector_set(gradient, 0, d_a);
    gsl_vector_set(gradient, 1, d_b);
}

/* Just a wrapper for gsl_ran_gamma.

   cut & pasted from the GSL documentation:
\f$          p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx \f$

See the notes for \ref apop_exponential on a popular alternate form.
*/
static void gamma_rng( double *out, gsl_rng* r, apop_model *p){
    *out    = gsl_ran_gamma(r, gsl_vector_get(p->parameters->vector, 0), gsl_vector_get(p->parameters->vector, 1));
}


/** The Gamma distribution

Location of data in the grid is not relevant; send it a 1 x N, N x 1, or N x M and it will all be the same.

apop_gamma.estimate() is an MLE, so feed it appropriate \ref apop_params.
  
If you have frequency or ranking data, you probably mean to be using \ref apop_gamma_rank.

\f$G(x, a, b)     = 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da    =  -\psi(a) - ln b + ln(x) \f$    (also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db    =  -a/b - x \f$
\ingroup models
*/
apop_model apop_gamma = {"Gamma distribution", 2,0,0,
     .estimate = gamma_estimate, .p = gamma_p, .log_likelihood = gamma_log_likelihood, 
     .score = gamma_dlog_likelihood, .constraint = beta_zero_and_one_greater_than_x_constraint, 
     .draw = gamma_rng};
