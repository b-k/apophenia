/** \file apop_exponential_rank.c

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"

//The default list. Probably don't need them all.
#include "types.h"
#include "bootstrap.h"
#include "regression.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <apophenia/stats.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>

/* Let k be the rank, and x_k be the number of elements at that rank;
 then the mean rank (and therefore the most likely estimate for the
 exponential parameter) is sum(k * x_k)/sum(x) */
static apop_model * exponential_rank_estimate(apop_data * data, apop_model *parameters){
  double          colsum,
                  numerator   = 0,
                  grand_total = 0;
  apop_model 	*est= parameters ? parameters : apop_model_copy(apop_exponential_rank);
  apop_model_clear(data, est);
  int             i;
    for(i=0; i< data->matrix->size2; i++){
        APOP_MATRIX_COL(data->matrix, i, v);
        colsum       = apop_sum(v);
        numerator   += colsum * (i+1);
        grand_total += colsum;
    }
	gsl_vector_set(est->parameters->vector, 0, numerator/grand_total);
    est->llikelihood	= apop_exponential_rank.log_likelihood(est->parameters, data, parameters);
	/*if (est->uses.covariance)
		apop_numerical_covariance_matrix(apop_exponential_rank, est, data);*/
	return est;
}

static double beta_greater_than_x_constraint(const apop_data *beta, apop_data *returned_beta, apop_model *v){
    //constraint is 0 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint) constraint = apop_data_calloc(1,1,1);
    apop_data_set(constraint, 0, 0, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}


static double rank_exponential_log_likelihood(const apop_data *beta, apop_data *d, apop_model *params){
  double		    b		    = gsl_vector_get(beta->vector, 0),
		          p,
		          llikelihood = 0,
		          ln_b		= log(b);
  gsl_matrix	    *data		= d->matrix;
  int 		    k;
  gsl_vector      v;
	for (k=0; k< data->size2; k++){
		p	            = -ln_b - (k+1)/b;
        v               = gsl_matrix_column(data, k).vector;
		llikelihood    += apop_sum(&v) * p; 
	}
	return llikelihood;
}

static double rank_exp_p(const apop_data *beta, apop_data *d, apop_model *p){
    return exp(rank_exponential_log_likelihood(beta, d, p));
}

// The exponential distribution. A one-parameter likelihood fn.
//
//\todo Check that the borderline work here is correct too.
static void rank_exponential_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_model *params){
  double		    bb		        = gsl_vector_get(beta->vector, 0);
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


/** The exponential distribution. A one-parameter likelihood fn.

Right now, it is keyed toward network analysis, meaning that the data
structure requires that the first column be the percentage of observations
which link to the most popular, the second column the percentage of
observations which link to the second-most popular, et cetera.


\f$Z(\mu,k) 	= 1/\mu e^{-k/\mu} 			\f$ <br>
\f$ln Z(\mu,k) 	= -\ln(\mu) - k/\mu			\f$ <br>

Some folks write the function as:
\f$Z(C,k) dx = \ln C C^{-k}. \f$
If you prefer this form, just convert your parameter via \f$\mu = {1\over
\ln C}\f$ (and convert back from the parameters this function gives you
via \f$C=\exp(1/\mu)\f$.

apop_exponential.estimate() is an MLE, so feed it appropriate \ref apop_params.

\ingroup models
\todo Check that the borderline work here is correct.
\todo Write a second object for the plain old not-network data Exponential.
*/
apop_model apop_exponential_rank = {"Exponential, rank data",  1,0,0,
	 .estimate = exponential_rank_estimate, .p = rank_exp_p, .log_likelihood = rank_exponential_log_likelihood, 
     .score = rank_exponential_dlog_likelihood, .constraint = beta_greater_than_x_constraint};
