/** \file apop_gamma_rank.c


Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/


//The default list. Probably don't need them all.
#include "types.h"
#include "stats.h"
#include "model.h"
#include "bootstrap.h"
#include "regression.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <stdio.h>
#include <assert.h>

static apop_estimate * gamma_rank_estimate(apop_data * data,  void *parameters){
    return apop_maximum_likelihood(data, apop_gamma_rank, parameters);
}

static double beta_zero_and_one_greater_than_x_constraint(gsl_vector *beta, void * d, gsl_vector *returned_beta){
double  limit0      = 0,
        limit1      = 0,
        tolerance   = 1e-1,
        beta0       = gsl_vector_get(beta, 0),
        beta1       = gsl_vector_get(beta, 1);
    if (beta0 > limit0 && beta1 > limit1) 
        return 0;
    //else:
    gsl_vector_memcpy(returned_beta, beta);
    gsl_vector_set(returned_beta, 0, GSL_MAX(limit0 + tolerance, beta0));
    gsl_vector_set(returned_beta, 1, GSL_MAX(limit1 + tolerance, beta1));
    return GSL_MAX(limit0 - beta0, 0) + GSL_MAX(limit1 - beta1, 0);    
}

static double gamma_rank_log_likelihood(const gsl_vector *beta, void *d){
float           a           = gsl_vector_get(beta, 0),
                b           = gsl_vector_get(beta, 1);
    //assert (a>0 && b>0);
    if (a<=0 || b<=0) printf("assertion failed.\n");
    if (gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;    
int             k;
gsl_matrix      *data       = d;
double          llikelihood = 0,
                ln_ga       = gsl_sf_lngamma(a),
                a_ln_b      = a * log(b),
                ln_k;
gsl_vector      v;
    for (k=0; k< data->size2; k++){
            ln_k            = log(k+1);
            v               = gsl_matrix_column(data, k).vector;
            llikelihood    += apop_sum(&v) * (-ln_ga - a_ln_b + (a-1) * ln_k  - (k+1)/b);
        }
    return llikelihood;
}

/** The derivative of the Gamma distribution, for use in likelihood
 * minimization. You'll probably never need to call this directly.*/
static void gamma_rank_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double          a       = gsl_vector_get(beta, 0),
                b       = gsl_vector_get(beta, 1);
int             k;
gsl_matrix     *data    = d;
double          d_a     = 0,
                d_b     = 0,
                psi_a   = gsl_sf_psi(a),
                ln_b    = log(b),
                x, ln_k;
gsl_vector_view v;
    for (k=0; k< data->size2; k++){
        ln_k    = log(k +1);
        v       = gsl_matrix_column(data, k);
        x       = apop_sum(&(v.vector));
        d_a    += x * (-psi_a - ln_b + ln_k);
        d_b    += x * (a/b - (k+1)/gsl_pow_2(b));
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
    return gsl_ran_gamma(r, a[0], a[1]);
}


/** The Gamma distribution for frequency ranks

  The usual \ref apop_gamma distribution assumes the data is just a
  list of numbers to which \f$\alpha\f$ and \f$\beta\f$ will be fit,
  that hapens to be in grid format.

apop_gamma_rank.estimate() is an MLE, so feed it appropriate \ref apop_estimation_params.

Here, we assume that the data is ranking frequencies: data[7][0] is
the number of times the first-ranked item appears in data set number
seven, data[7][1] is the number of times the second-ranked item appears,
et cetera. Below, \f$x_k\f$ is the count (or percentage) of the \f$k\f$th ranked item.
Also, \f$d ln \gamma(k) \equiv \psi(k)\f$.

\f$G(x, a, b)       = \prod_{k=1}^n \left(1/(\Gamma(a) b^a)  k^{a-1} e^{-k/b}\right)^{x_k}\f$

\f$ln G(x, a, b)    = \sum_{k=1}^n x_k (-ln \Gamma(a) - a ln b + (a-1)ln(k) + -k/b)\f$

\f$d ln G/ da       = \sum_{k=1}^n x_k ( -\psi(a) - ln b + ln(k)) \f$

\f$d ln G/ db       = \sum_{k=1}^n x_k (-a/b - k) \f$

\ingroup models
*/
apop_model apop_gamma_rank = {"Gamma, rank data", 2, 
{
    1,    //parameters
    1,    //covariance
    1,    //confidence
	0,	//dependent
	0,	//predicted
    1,    //log_likelihood
    0    //names;
},
    gamma_rank_estimate, gamma_rank_log_likelihood, gamma_rank_dlog_likelihood, NULL, 
    beta_zero_and_one_greater_than_x_constraint,  gamma_rng};
