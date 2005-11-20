/** \file apop_zipf.c

  The Zipf distribution.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/


//The default list. Probably don't need them all.
#include "name.h"
#include "model.h"
#include "output.h"
#include "regression.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <assert.h>



static apop_estimate * apop_zipf_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
	return apop_maximum_likelihood(data, uses, apop_zipf, *(apop_estimation_params *)parameters);
}



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


///////////////////////
//The Zipf distribution
///////////////////////
#include <gsl/gsl_sf_zeta.h>

/* The Zipf distribution.

\f$Z(a)		= {1\over \zeta(a) * i^a}		\f$<br>

 \todo link this fn in with the object 
static double apop_zipf_likelihood(double a, int i){
double		z	= 1/(gsl_sf_zeta(a) * pow(i, a));
	return z;
}
*/

static double apop_zipf_log_likelihood(const gsl_vector *beta, void *d){
static double	ka	= 0;
gsl_matrix	*data	= d;
double		like	= 0, 
		bb	= gsl_vector_get(beta, 0),
		z;
int 		i, j;
	if (bb <= 1) {		//run away
		if (ka ==0){
			gsl_vector *	b_ka	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka, 0, 1+GSL_DBL_EPSILON);
		 	ka	= apop_zipf_log_likelihood(b_ka, d);
			gsl_vector_free (b_ka);
		}
		return keep_away(bb, 1, ka);
	}			//else:
	for(j=0; j< data->size2; j++){
		z	 = -log(gsl_sf_zeta(bb)) - bb * log(j+1);
		for(i=0; i< data->size1; i++){
			like	+= gsl_matrix_get(data,i,j) * z;
		}
	}
	return like;
}	


/** Dlog likelihood for the zipf distribution.

This fn has a bug, and I can't find it right now. Feel free to look over it.  In the mean time, it's not linked in for the apop_zipf object.

\todo Fix this. */
static void apop_zipf_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double		a	= gsl_vector_get(beta, 0);
static double	ka	= 0;
gsl_matrix	*data	= d;
int 		i, j;
double		dlike	= 0, 
		colsum, dz;
	if (a <= 1) {		//keep away
		if (ka ==0){
			gsl_vector 	*b_ka	= gsl_vector_alloc(1);
			gsl_vector 	*b_kg	= gsl_vector_alloc(1);
			gsl_vector_set(b_ka,0, 1+GSL_DBL_EPSILON);
		 	apop_zipf_dlog_likelihood(b_ka , d, b_kg);
			ka	= gsl_vector_get(b_kg, 0);
			gsl_vector_free (b_ka);
			gsl_vector_free (b_kg);
		}
		gsl_vector_set(gradient,0, -keep_away(a, 1, ka));
		return;
	}			//else:
	for(j=0; j< data->size2; j++){
		dz	 = a*gsl_sf_zeta(a+1)/gsl_sf_zeta(a) - log(j+1);
		colsum	= 0;
		for(i=0; i< data->size1; i++){
			colsum	+= gsl_matrix_get(data,i,j);
		}
		dlike	+= colsum * dz;
	}
	gsl_vector_set(gradient,0,dlike);
}	


/** Draw from a Zipf distribution with parameter \f$ a \f$

Call this fn using \ref apop_zipf.rng().

Returns a ranking: If the population were Zipf distributed, you're most
likely to get the 1st most common item, so this produces a lot of ones,
a great deal of twos, and so on.


For example:
\code
gsl_rng *       r;
gsl_rng_env_setup();
r=gsl_rng_alloc(gsl_rng_taus);	//for example. 
apop_zipf.rng(r, 1.4);
\endcode

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 551.  */
static double apop_zipf_rng(gsl_rng* r, double * a){
if (*a  <= 1)	
	{printf("apop_zipf.rng: Zipf needs a parameter >=1. Returning 0.\n"); return 0;};
int		x;
double		u, v, t, 
		b 	= pow(2, *a-1), 
		ainv	= -(1.0/(*a-1));
	do {
		u	= gsl_rng_uniform(r);
		v	= gsl_rng_uniform(r);
		x	= pow(u, ainv);
		t	= pow((1.0 + 1.0/x), (*a-1));
	} while (v * x * (t-1.0)/(b-1) > t/b);
	return x;
}




/** The Zipf distribution.
Wikipedia has notes on the <a href="http://en.wikipedia.org/wiki/Zipf_distribution">Zipf distribution</a>. 

The data set needs to be in rank-form. The first column is the frequency of the most common item, the second is the frequency of the second most common item, &c.

apop_zipf.estimate() is an MLE, so feed it appropriate \ref apop_estimation_params.

\f$Z(a)		= {1\over \zeta(a) * i^a}		\f$

\f$lnZ(a)	= -(\log(\zeta(a)) + a \log(i))	\f$

\f$dlnZ(a)/da	= -{\zeta(a)\over a \log(\zeta(a-1))} -  \log(i)		\f$
\ingroup likelihood_fns
*/
apop_model apop_zipf = {"Zipf", 1, apop_zipf_estimate, apop_zipf_log_likelihood, apop_zipf_dlog_likelihood, NULL, apop_zipf_rng};
//apop_model apop_zipf = {"Zipf", 1, apop_zipf_log_likelihood, NULL, NULL, 0, NULL, apop_zipf_rng};
