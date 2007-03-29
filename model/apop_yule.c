/** \file apop_yule.c

  The Yule distribution. A special case of the Waring.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/



//The default list. Probably don't need them all.
#include "model.h"
#include "types.h"
#include "output.h"
#include "mapply.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>

static apop_estimate * yule_estimate(apop_data * data, void *parameters){
	return apop_maximum_likelihood(data, apop_yule, parameters);
}

static double beta_greater_than_x_constraint(const apop_data *beta, void * d, apop_data *returned_beta, void *v){
    //constraint is 1 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint)constraint= apop_data_calloc(1,1,1);
    apop_data_set(constraint, 0, 0, 1);
    apop_data_set(constraint, 0, -1, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}

static long double  bb;

static double  apply_me(gsl_vector *v){
  long double   likelihood =0, 
                ln_k, ln_bb_k;
  int 	        k;
  double        pt;
    for (k=0; k< v->size; k++)	{
        pt            = gsl_vector_get(v, k);
		if (pt>=1) 	ln_k	= gsl_sf_lngamma(pt);
		else		ln_k	= 0;
        ln_bb_k		  = gsl_sf_lngamma(pt+bb);
        likelihood   +=  ln_k - ln_bb_k;
    }
    return likelihood;
}

static double  dapply_me(gsl_vector *v){
  int           k;
  double        pt; 
  long double   d   = 0;
    for (k=0; k< v->size; k++)	{
        pt   = gsl_vector_get(v, k);
        d	-= gsl_sf_psi(pt + bb);
    }
    return d;
}

static double yule_log_likelihood(const apop_data *beta, apop_data *d, void *p){
    bb	            = gsl_vector_get(beta->vector, 0);
  long double   ln_bb		    = gsl_sf_lngamma(bb),
                ln_bb_less_1    = log(bb-1);
  gsl_vector *  v               = apop_matrix_map(d->matrix, apply_me);
  double        likelihood      = apop_vector_sum(v);
    gsl_vector_free(v);
	return likelihood + (ln_bb_less_1 + ln_bb) * d->matrix->size1 * d->matrix->size2;
}

static void yule_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, void *p){
	//Psi is the derivative of the log gamma function.
    bb		= gsl_vector_get(beta->vector, 0);
double		    bb_minus_one_inv= 1/(bb-1),
		        psi_bb	        = gsl_sf_psi(bb);
  gsl_vector *  v               = apop_matrix_map(d->matrix, dapply_me);
  double        d_bb            = apop_vector_sum(v);
    gsl_vector_free(v);
    d_bb    += (bb_minus_one_inv + psi_bb) * d->matrix->size1 * d->matrix->size2;
	gsl_vector_set(gradient, 0, d_bb);
}

static double yule_p(const apop_data *beta, apop_data *d, void *p){
    return exp(yule_log_likelihood(beta, d, p));
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
static void yule_rng( double *out, apop_data * a, apop_ep *p, gsl_rng * r){
double 	e1, e2;
int		x;
	e1	= gsl_ran_exponential(r, 1);
	e2	= gsl_ran_exponential(r, 1);
	x	= GSL_MAX((int) (- e1  / log(1 - exp(-e2 / ((a->vector->data)[0] -1)))), 0);
	*out =  x + 1;	//we rounded down to floor, but want ceil.
}



/** The Yule distribution

The special case of Waring where \f$ \alpha = 0.	\f$<br>

The data set needs to be in rank-form. The first column is the frequency
of the most common item, the second is the frequency of the second most
common item, &c.

apop_yule.estimate() is an MLE, so feed it appropriate \ref apop_ep.

\f$ Y(x, b) 	= (b-1) \gamma(b) \gamma(k) / \gamma(k+b)			\f$

\f$ \ln Y(x, b)	= \ln(b-1) + ln\gamma(b) + \ln\gamma(k) - \ln\gamma(k+b)	\f$

\f$ d\ln Y/db	= 1/(b-1)  + \psi(b) - \psi(k+b)				\f$
\ingroup models
\todo I'm pretty sure Wikipedia's specification of the Yule is wrong; I should check and fix when I can check references.
*/
apop_model apop_yule = {"Yule", 1, 0,0, 
	yule_estimate, yule_p, yule_log_likelihood, yule_dlog_likelihood, beta_greater_than_x_constraint, yule_rng};
