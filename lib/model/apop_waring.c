
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


/////////////////////////
//The Waring distribution
/////////////////////////
/** The Waring distribution

\f$W(x, b,a) 	= (b-1) \gamma(b+a) \gamma(k+a) / [\gamma(a+1) \gamma(k+a+b)]\f$

\f$\ln W(x, b, a) = \ln(b-1) + \ln\gamma(b+a) + \ln\gamma(k+a) - \ln\gamma(a+1) - \ln\gamma(k+a+b)\f$

\f$dlnW/db	= 1/(b-1)  + \psi(b+a) - \psi(k+a+b)\f$

\f$dlnW/da	= \psi(b+a) + \psi(k+a) - \psi(a+1) - \psi(k+a+b)\f$

\ingroup likelihood_fns
*/
static apop_estimate * waring_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
	apop_inventory_filter(uses, apop_waring.inventory_filter);
	return apop_maximum_likelihood(data, uses, apop_waring, *(apop_estimation_params *)parameters);
}

static double waring_log_likelihood(const gsl_vector *beta, void *d){
float		bb	= gsl_vector_get(beta, 0),
		a	= gsl_vector_get(beta, 1);
double		ka;		//recalculated every time.
	if (bb <1 || a < 0) {
		gsl_vector *	b_ka	= gsl_vector_alloc(2);
		gsl_vector_set(b_ka, 0, GSL_MAX(bb, 1) +  1e20*GSL_DBL_EPSILON);
		gsl_vector_set(b_ka, 1, GSL_MAX(a, 0) + 1e20*GSL_DBL_EPSILON);
	 	ka	= waring_log_likelihood(b_ka, d);
		gsl_vector_free (b_ka);
		if (bb<=1) 	return keep_away(bb, 1, ka);
		else 		return keep_away(bb, 0, ka);
	}			//else:
int 		i, k;
gsl_matrix*	data		= d;
double 		ln_a_k, ln_bb_a_k, p,
		likelihood 	= 0,
		ln_bb_a		= gsl_sf_lngamma(bb + a),
		ln_a_mas_1	= gsl_sf_lngamma(a + 1),
		ln_bb_less_1	= log(bb - 1);
	for (k=0; k< data->size2; k++){	//more efficient to go column-by-column
		ln_bb_a_k	 = gsl_sf_lngamma(k +1 + a + bb);
		ln_a_k		 = gsl_sf_lngamma(k +1 + a);
		p		 = ln_bb_less_1 + ln_a_k + ln_bb_a - ln_a_mas_1 - ln_bb_a_k;
		for (i=0; i< data->size1; i++){
			likelihood	+= gsl_matrix_get(data, i, k) *p;
		}
	}
	return likelihood;
}

/** The derivative of the Waring distribution, for use in likelihood
 minimization. You'll probably never need to call this directy.*/
static void waring_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
	//Psi is the derivative of the log gamma function.
float		bb		= gsl_vector_get(beta, 0),
		a		= gsl_vector_get(beta, 1);
int 		i, k;
gsl_matrix	*data		= d;
double		bb_minus_one_inv= 1/(bb-1),
		psi_a_bb	= gsl_sf_psi(bb + a),
		psi_a_mas_one	= gsl_sf_psi(a+1),
		psi_a_k,
		psi_bb_a_k,
		d_bb		= 0,
		d_a		= 0;
	for (k=0; k< data->size2; k++){	//more efficient to go column-by-column
		psi_bb_a_k	 = gsl_sf_psi(k +1 + a + bb);
		psi_a_k		 = gsl_sf_psi(k +1 + a);
		for (i=0; i< data->size1; i++){
			d_bb		+= gsl_matrix_get(data, i, k) *(bb_minus_one_inv + psi_a_bb - psi_bb_a_k);
			d_a		+= gsl_matrix_get(data, i, k) *(psi_a_bb + psi_a_k - psi_a_mas_one - psi_bb_a_k);
		}
	}
	gsl_vector_set(gradient, 0, d_bb);
	gsl_vector_set(gradient, 1, d_a);
}


/** RNG from a Generalized Hypergeometric type B3.

 Devroye uses this as the base for many of his
 distribution-generators, e.g., \ref apop_waring_rng. 
*/ 
double apop_GHgB3_rng(gsl_rng * r, double* a){
if ((a[0]<=0) || (a[1] <= 0) || (a[2] <=0)){
	printf("apop_GHgB3_rng took a zero parameter; bad.\n");
	return 0;
	}
double		aa	= gsl_ran_gamma(r, a[0], 1),
		b	= gsl_ran_gamma(r, a[1], 1),
		c	= gsl_ran_gamma(r, a[2], 1);
int		p;
	p	= gsl_ran_poisson(r, aa*b/c);
	return p;
}

/** Give me parameters, and I'll draw a ranking from the appropriate
Waring distribution. [I.e., if I randomly draw from a Waring-distributed
population, return the ranking of the item I just drew.]

Page seven of:
L. Devroye, <a href="http://cgm.cs.mcgill.ca/~luc/digammapaper.ps">Random
variate generation for the digamma and trigamma distributions</a>, Journal
of Statistical Computation and Simulation, vol. 43, pp. 197-216, 1992.
*/
static double waring_rng(gsl_rng *r, double *a){
//The key to covnert from Devroye's GHgB3 notation to what I
//consider to be the standard Waring notation in \ref apop_waring:
// a = a + 1
// b = 1 
// c = b - 1 
// n = k - 1 , so if it returns 0, that's first rank.
// OK, I hope that clears everything up.
double		x, u,
		params[]	={a[0]+1, 1, a[1]-1};
	do{
		x	= 1+ apop_GHgB3_rng(r, params);
		u	= gsl_rng_uniform(r);
	} while (u >= (x + a[0])/(GSL_MAX(a[0]+1,1)*x));
	return x;
}

/** The Waring distribution
The data set needs to be in rank-form. The first column is the frequency of the most common item, the second is the frequency of the second most common item, &c.

apop_waring.estimate() is an MLE, so feed it appropriate \ref apop_estimation_params.

\f$W(x,k, b,a) 	= (b-1) \gamma(b+a) \gamma(k+a) / [\gamma(a+1) \gamma(k+a+b)]\f$

\f$\ln W(x,k, b, a) = ln(b-1) + lng(b+a) + lng(k+a) - lng(a+1) - lng(k+a+b)\f$

\f$dlnW/db	= 1/(b-1)  + \psi(b+a) - \psi(k+a+b)\f$

\f$dlnW/da	= \psi(b+a) + \psi(k+a) - \psi(a+1) - \psi(k+a+b)\f$
\ingroup models
*/
//apop_model apop_waring = {"Waring", 2, apop_waring_log_likelihood, NULL, NULL, 0, NULL, apop_waring_rng};
apop_model apop_waring = {"Waring", 2, 
 {
	1,	//parameters
	1,	//covariance
	1,	//confidence
	0,	//predicted
	0,	//residuals
	1,	//log_likelihood
	1	//names;
},
	waring_estimate, waring_log_likelihood, waring_dlog_likelihood, NULL,  waring_rng};
