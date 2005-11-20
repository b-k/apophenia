/** \file apop_yule.c

  The Yule distribution. A special case of the Waring.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/



//The default list. Probably don't need them all.
#include "name.h"
#include "output.h"
#include "bootstrap.h"
#include "regression.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>

/** This function is used to keep the minimizer away from bounds.

If you just return GSL_POSINF at the bounds, it's not necessarily smart
enough to get it.  This helps the minimzer along by providing a (almost)
continuous, steep line which steers the minimizer back to the covered
range. 
\todo Replace this with apop_constraints.
*/
static double keep_away(double value, double limit,  double base){
	return (50000+fabs(value - limit)) * base;
}

static apop_estimate * yule_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
	return apop_maximum_likelihood(data, uses, apop_yule, *(apop_estimation_params *)parameters);
}

static double yule_log_likelihood(const gsl_vector *beta, void *d){
float		bb	= gsl_vector_get(beta, 0);
static double	ka	= 0;
	if (bb < 1) {		//run away
		if (ka ==0){
			gsl_vector *	b_ka	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka, 0, 1.00001);
		 	ka	= yule_log_likelihood(b_ka, d);
			gsl_vector_free (b_ka);
		}
		return keep_away(bb, 1, ka);
	}			//else:
int 		i, k;
gsl_matrix 	*data		= d;
float 		ln_k, ln_bb_k,
		likelihood 	= 0,
		ln_bb		= gsl_sf_lngamma(bb),
		ln_bb_less_1	= log(bb-1);
	for (k=0; k< data->size2; k++)	{
		if (k>=1) 	ln_k	= gsl_sf_lngamma(k+1);
		else		ln_k	= 0;
		ln_bb_k		= gsl_sf_lngamma(k+1+bb);
		for (i=0; i< data->size1; i++){
			likelihood	+= gsl_matrix_get(data,i,k) * (ln_bb_less_1 + ln_k + ln_bb - ln_bb_k);
		}
	}
	return likelihood;
}

static void yule_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
	//Psi is the derivative of the log gamma function.
float		bb		= gsl_vector_get(beta, 0);
static double	dka		= 0;
int 		i, k;
gsl_matrix	*data		= d;
double		bb_minus_one_inv= 1/(bb-1),
		psi_bb		= gsl_sf_psi(bb),
		psi_bb_k,
		p,
		d_bb		= 0;
	if (bb < 1) {		//keep away
		if (dka ==0){
			gsl_vector 	*b_ka	= gsl_vector_alloc(1);
			gsl_vector 	*b_kg	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, 1+GSL_DBL_EPSILON);
		 	yule_dlog_likelihood(b_ka , d, b_kg);
			dka	= gsl_vector_get(b_kg, 0);
			gsl_vector_free (b_ka);
			gsl_vector_free (b_kg);
		}
		gsl_vector_set(gradient,0, keep_away(bb, 1, dka));
		return;
	}			//else:
	for (k=0; k< data->size2; k++){
		psi_bb_k= gsl_sf_psi(k +1 + bb);
		p	= bb_minus_one_inv + psi_bb - psi_bb_k;
		for (i=0; i< data->size1; i++){
			d_bb		+= gsl_matrix_get(data, i, k) * p;
		}
	}
	gsl_vector_set(gradient, 0, d_bb);
}


/** Draw from a Yule distribution with parameter a

Call this fn using \ref apop_yule.rng().

\param	a	The parameter.
\param	r	A gsl_rng that you've already set up.

For example:
\code
gsl_rng *       r;
gsl_rng_env_setup();
r=gsl_rng_alloc(gsl_rng_taus);	//for example. 
apop_yule_rng(r, 1.4);
\endcode

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 553.  */
static double yule_rng(gsl_rng * r, double* a){
double 		e1, e2;
int		x;
	e1	= gsl_ran_exponential(r, 1);
	e2	= gsl_ran_exponential(r, 1);
	x	= - e1  / log(1 - exp(-e2 / (*a -1)));
	return  x + 1;	//we rounded down to floor, but want ceil.
}



/** The Yule distribution

The special case of Waring where \f$ \alpha = 0.	\f$<br>

The data set needs to be in rank-form. The first column is the frequency
of the most common item, the second is the frequency of the second most
common item, &c.

apop_yule.estimate() is an MLE, so feed it appropriate \ref apop_estimation_params.

\f$ Y(x, b) 	= (b-1) \gamma(b) \gamma(k) / \gamma(k+b)			\f$

\f$ \ln Y(x, b)	= \ln(b-1) + ln\gamma(b) + \ln\gamma(k) - \ln\gamma(k+b)	\f$

\f$ d\ln Y/db	= 1/(b-1)  + \psi(b) - \psi(k+b)				\f$
\ingroup models
\todo I'm pretty sure Wikipedia's specification of the Yule is wrong; I should check and fix when I can check references.
*/
apop_model apop_yule = {"Yule", 1, yule_estimate, yule_log_likelihood, yule_dlog_likelihood, NULL, yule_rng};
