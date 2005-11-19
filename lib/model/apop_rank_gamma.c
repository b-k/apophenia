/** \file apop_rank_gamma.c


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


double apop_gamma_rank_log_likelihood(const gsl_vector *beta, void *d){
float		a	= gsl_vector_get(beta, 0),
		b	= gsl_vector_get(beta, 1);
	if (a <= 0 || b <= 0 || gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;	
static double		ka = 0;
	if (b <0 || a < 0) {
		if (ka==0){
		gsl_vector *	b_ka	= gsl_vector_alloc(2);
			gsl_vector_set(b_ka, 0, GSL_MAX(b, 0) + 1e-6);
			gsl_vector_set(b_ka, 0, GSL_MAX(a, 0) + 1e-6);
	 		ka	= apop_gamma_log_likelihood(b_ka, d);
			gsl_vector_free (b_ka);
		}
		if (b<=0) 	return keep_away(b, 0, ka);
		else 		return keep_away(a, 0, ka);
	}			//else:
int 		i, k;
gsl_matrix	*data		= d;
float 		llikelihood 	= 0,
		ln_ga		= gsl_sf_lngamma(a),
		ln_b		= log(b),
		a_ln_b		= a * ln_b,
		ln_k;
	for (k=0; k< data->size2; k++){
		ln_k	= log(k+1);
		for (i=0; i< data->size1; i++)
			llikelihood	+= gsl_matrix_get(data, i, k) * (-ln_ga - a_ln_b + (a-1) * ln_k  - (k+1)/b);
	}
	return llikelihood;
}

/** The derivative of the Gamma distribution, for use in likelihood
 * minimization. You'll probably never need to call this directly.*/
void apop_gamma_rank_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
float		a	= gsl_vector_get(beta, 0),
		b	= gsl_vector_get(beta, 1);
	//if (a <= 0 || b <= 0 || gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;	
						//a sign to the minimizer to look elsewhere.
int 		i, k;
gsl_matrix	*data	= d;
float 		d_a 	= 0,
		d_b	= 0,
		psi_a	= gsl_sf_psi(a),
		ln_b	= log(b),
		x, ln_k;
	for (k=0; k< data->size2; k++){
		ln_k	= log(k +1);
		for (i=0; i< data->size1; i++){
			x	= gsl_matrix_get(data, i, k);
			d_a	+= x * (-psi_a - ln_b + ln_k);
			d_b	+= x *(-a/b - (k+1));
		}
		}
	gsl_vector_set(gradient,0, d_a);
	gsl_vector_set(gradient,1, d_b);
}


/** The Gamma distribution for frequency ranks

  The usual \ref apop_gamma distribution assumes the data is just a
  list of numbers to which \f$\alpha\f$ and \f$\beta\f$ will be fit,
  that hapens to be in grid format.

  Here, we assume that the data is ranking frequencies: data[7][0] is the number of times the first-ranked item appears in data set number seven, data[7][1] is the number of times the second-ranked item appears, et cetera. 

\ingroup likelihood_fns
*/
apop_model apop_gamma_rank = {"Gamma", 2, apop_gamma_rank_log_likelihood, apop_gamma_rank_dlog_likelihood, NULL, 0, NULL, NULL};
//apop_model apop_gamma_rank = {"Gamma", 2, apop_gamma_rank_log_likelihood, NULL, NULL, 0, NULL, NULL};

