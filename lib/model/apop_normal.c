/** \file apop_normal.c

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
static double normal_log_likelihood(const gsl_vector *beta, void *d);

static double keep_away(double value, double limit,  double base){
	return (50000+fabs(value - limit)) * base;
}


//////////////////
//The Normal (gaussian) distribution
//////////////////

static apop_estimate * normal_estimate(gsl_matrix * data, apop_inventory *uses, void *parameters){
	apop_inventory_filter(uses, apop_normal.inventory_filter);
apop_estimate *est	= apop_estimate_alloc(0,2,NULL, *uses);
double		avg	= 0,
		avg2	= 0;
int		i,j, cnt= 0;

double          x, ratio;
	for(i=0; i < data->size1; i++)
		for(j=0; j < data->size1; j++){
    			x 	= gsl_matrix_get(data, i,j);
    			ratio	= cnt/(cnt+1.0);
    			cnt	++;
    			avg	*= ratio;
    			avg2	*= ratio;
    			avg	+= x/(cnt +0.0);
    			avg2 	+= gsl_pow_2(x)/(cnt +0.0);
  		}
	gsl_vector_set(est->parameters,0, avg);
	gsl_vector_set(est->parameters,1, (avg2 - gsl_pow_2(avg))); //E[x^2] - E^2[x]
	if (est->uses.log_likelihood)
		est->log_likelihood	= normal_log_likelihood(est->parameters, data);
	return est;
}


/** The log likelihood function for the Normal.

The log likelihood function and dlog likelihood don't care about your
rows of data; if you have an 8 x 7 data set, it will give you the log
likelihood of those 56 observations given the mean and variance (i.e.,
\f$\sigma^2\f$, not std deviation=\f$\sigma\f$) you provide.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-(x-\mu)^2 / 2\sigma^2)\f$
\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\param beta	beta[0]=the mean; beta[1]=the variance
\param d	the set of data points; see notes.
*/
static double normal_log_likelihood(const gsl_vector *beta, void *d){
double 		mu	= gsl_vector_get(beta,0),
		ss	= 2 * gsl_vector_get(beta,1),
		ll	= 0;
size_t		i,j;
gsl_matrix	*data	= d;
static double	ka	= 0;
	if (ss <0) {
		gsl_vector *	b_ka	= gsl_vector_alloc(2);
		gsl_vector_set(b_ka, 0, mu);
		gsl_vector_set(b_ka, 1, GSL_DBL_EPSILON);
	 	ka	= normal_log_likelihood(b_ka, d);
		gsl_vector_free (b_ka);
		//return keep_away(ss/2.0, 0, ka);
		return keep_away(gsl_vector_get(beta,1), 0, ka);
	}			//else:

	for (i=0;i< data->size1; i++)
		for (j=0;j< data->size2; j++)
			ll	-= gsl_pow_2(gsl_matrix_get(data, i, j) - mu) / ss;
			//Or use the stock gsl function (but then just return -ll).
			//ll	+= log(gsl_ran_gaussian_pdf((gsl_matrix_get(data, i, j) - mu), gsl_vector_get(beta,1)));
	ll 	-= data->size1 * data->size2 * log(ss * M_PI)/2.0; //second term is positive because it subtracts -log.
	return ll;
}

/** Gradient of the log likelihood function

To tell you the truth, I have no idea when anybody would need this, but it's here for completeness. 
\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$
\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$
\todo Add constraint that \f$\sigma^2>0\f$.
 */
static void normal_dlog_likelihood(const gsl_vector *beta, void *d, gsl_vector *gradient){
double 		mu	= gsl_vector_get(beta,0),
		ss	= gsl_vector_get(beta,1),
		dll	= 0,
		sll	= 0,
		x;
int		i,j;
gsl_matrix	*data	= d;
	for (i=0;i< data->size1; i++)
		for (j=0;j< data->size2; j++){
			x	 = gsl_matrix_get(data, i, j);
			dll	+= (x - mu);
			sll	+= gsl_pow_2(x - mu);
		}
	gsl_vector_set(gradient, 0, dll/ss);
	gsl_vector_set(gradient, 1, sll/(2*gsl_pow_2(ss))+ data->size1 * data->size2 * 0.5/ss);
}

/** An apophenia wrapper for the GSL's Normal RNG.

Two differences: this one asks explicitly for a mean, and the GSL
assumes zero and makes you add the mean yourself; Apophenia tends to
prefer the variance (\f$\sigma^2\f$) wherever possible, while the GSL
uses the standard deviation here (\f$\sigma\f$)

\param r	a gsl_rng already allocated
\param a	the mean and the variance
 */
static double normal_rng(gsl_rng *r, double *a){
	//return gsl_ran_gaussian(r, sqrt(a[1])) + a[0];
	return gsl_ran_gaussian(r, sqrt(a[1])) + a[0];
}

/** You know it, it's your attractor in the limit, it's the Gaussian distribution.

Generally, fitting a Normal distribution via maximum likelihood is silly
(\ref apop_mean and \ref apop_var will find the parameters with maximum
likelihood just as quickly), but there exist reasons for wanting this
as an \ref apop_model object, so here you are.

The log likelihood function and dlog likelihood don't care about your
rows of data; if you have an 8 x 7 data set, it will give you the log
likelihood of those 56 observations given the mean and variance you provide.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2)\f$

\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$

\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$
\ingroup models
*/
apop_model apop_normal = {"Normal", 2, 
{
	1,	//parameters
	1,	//covariance
	1,	//confidence
	0,	//predicted
	0,	//residuals
	1,	//log_likelihood
	1	//names;
}, normal_estimate,
	 normal_log_likelihood, normal_dlog_likelihood, NULL, normal_rng};
//apop_model apop_normal = {"Normal", 2, apop_normal_log_likelihood, NULL, NULL, 0, NULL, apop_normal_rng};

/** This is a synonym for \ref apop_normal, q.v.
\ingroup models
*/
apop_model apop_gaussian = {"Gaussian", 2, 
{
	1,	//parameters
	1,	//covariance
	1,	//confidence
	0,	//predicted
	0,	//residuals
	1,	//log_likelihood
	1	//names;
}, normal_estimate,
	 normal_log_likelihood, normal_dlog_likelihood, NULL, normal_rng};
