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


/** The Gamma distribution

\f$G(x, a, b) 	= 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x)  -x/b\f$

\f$d ln G/ da	=  -\psi(a) - ln b + ln(x) \f$	(also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db	=  -a/b - x \f$

If you have frequency or ranking data, you probably mean to be using \ref apop_gamma_rank.
\ingroup likelihood_fns
*/
double apop_gamma_log_likelihood(const gsl_vector *beta, void *d){
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
		x;
	for (i=0; i< data->size1; i++)
		for (k=0; k< data->size2; k++){
			x		 = gsl_matrix_get(data, i, k);
			if (x!=0)
				llikelihood	+= -ln_ga - a_ln_b + (a-1) * log(x) - x/b;
		}
	return llikelihood;
}

/** The derivative of the Gamma distribution, for use in likelihood
 * minimization. You'll probably never need to call this directly.*/
void apop_gamma_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
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
		x;
	for (i=0; i< data->size1; i++)
		for (k=0; k< data->size2; k++){
			x		 = gsl_matrix_get(data, i, k);
			if (x!=0){
				d_a	+= -psi_a - ln_b + log(x);
				d_b	+= -a/b - x;
			}
		}
	gsl_vector_set(gradient,0, d_a);
	gsl_vector_set(gradient,1, d_b);
}


/** The Gamma distribution

The data set needs to be in rank-form. The first column is the frequency of the most common item, the second is the frequency of the second most common item, &c.

\f$G(x, a, b) 	= 1/(\Gamma(a) b^a)  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da	=  -\psi(a) - ln b + ln(x) \f$	(also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db	=  -a/b - x \f$
\ingroup likelihood_fns
*/
apop_model apop_gamma = {"Gamma", 2, apop_gamma_log_likelihood, apop_gamma_dlog_likelihood, NULL, NULL};
