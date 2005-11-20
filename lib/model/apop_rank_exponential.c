/** \file apop_rank_exponential.c

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

static apop_estimate * exponential_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
	apop_inventory_filter(uses, apop_exponential.inventory_filter);
	return apop_maximum_likelihood(data, uses, apop_exponential, *(apop_estimation_params *)parameters);
}

/** The exponential distribution. A one-parameter likelihood fn.
\f$Z(\mu,k) 	= 1/\mu e^{-k/\mu} 			\f$ <br>
\f$ln Z(\mu,k) 	= -\ln(\mu) - k/\mu			\f$ <br>

Some folks write the function as:
\f$Z(C,k) dx = \ln C C^{-k}. \f$
If you prefer this form, just convert your parameter via \f$\mu = {1\over
\ln C}\f$ (and convert back from the parameters this function gives you
via \f$C=\exp(1/\mu)\f$.

\ingroup likelihood_fns
\todo Set up an exponential object which makes use of the GSL.
\todo Check that the borderline work here is correct.
*/
static double rank_exponential_log_likelihood(const gsl_vector *beta, void *d){
double		bb		= gsl_vector_get(beta, 0),
		p,
		llikelihood 	= 0,
		ln_c		= log(bb);
		//ln_ln_c	= log(ln_c);
static double	ka		= 0;
gsl_matrix	*data		= d;
int 		i, k;
	if (bb <= 0) {		//run away
		if (ka ==0){
			gsl_vector *	b_ka	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka, 0, GSL_DBL_EPSILON);
		 	ka	= rank_exponential_log_likelihood(b_ka, d);
			gsl_vector_free (b_ka);
		}
		return keep_away(bb, 0, ka);
	}			//else:
	for (k=0; k< data->size2; k++){
		//p	= ln_ln_c - ln_c * (k+1);
		p	= -ln_c - (k)/bb;
		for (i=0; i< data->size1; i++)
			llikelihood	+=  gsl_matrix_get(data, i, k) * p;
	}
	return llikelihood;
}


/** The exponential distribution. A one-parameter likelihood fn.

\todo Check that the borderline work here is correct too.
*/
static void rank_exponential_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double		bb		= gsl_vector_get(beta, 0);
int 		i, k;
static double	dka		= 0;
gsl_matrix	*data		= d;
double 		d_likelihood 	= 0,
		one_over_ln_c	= 1/log(bb),
		p;
	if (bb <= 0) {		//keep away
		if (dka ==0){
			gsl_vector 	*b_ka	= gsl_vector_alloc(1);
			gsl_vector 	*b_kg	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, GSL_DBL_EPSILON);
		 	rank_exponential_dlog_likelihood(b_ka , d, b_kg);
			dka	= gsl_vector_get(b_kg, 0);
			gsl_vector_free (b_ka);
			gsl_vector_free (b_kg);
		}
		gsl_vector_set(gradient,0, -keep_away(bb, 1, dka));
		return;
	}			//else:
	for (k=0; k< data->size2; k++) {
		p	= (one_over_ln_c -(k+1)) /bb;
		for (i=0; i< data->size1; i++)			
			d_likelihood	+= gsl_matrix_get(data, i, k) * p; 
	}
	gsl_vector_set(gradient,0, d_likelihood);
}

/* Just a wrapper for gsl_ran_exponential.

   cut & pasted from the GSL documentation:
\f$p(x) dx = {1 \over \mu} \exp(-x/\mu) dx \f$

See the notes for \ref apop_exponential_rng on a popular alternate form.
*/
static double rank_exponential_rng(gsl_rng* r, double * a){
	//This fn exists because the GSL requires a double, 
	//while the apop_model structure requires a double*. 
	return gsl_ran_exponential(r, *a);
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

apop_exponential.estimate() is an MLE, so feed it appropriate \ref apop_estimation_params.

\ingroup models
\todo Check that the borderline work here is correct.
\todo Write a second object for the plain old not-network data Exponential.
*/
apop_model apop_rank_exponential = {"Exponential", 1,
	{
	1,	//parameters
	1,	//covariance
	1,	//confidence
	0,	//predicted
	0,	//residuals
	1,	//log_likelihood
	1	//names;
},
	 exponential_estimate, rank_exponential_log_likelihood, rank_exponential_dlog_likelihood, NULL,  rank_exponential_rng};


