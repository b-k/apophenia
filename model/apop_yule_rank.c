/** \file apop_yule_rank.c

  The Yule distribution. A special case of the Waring.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/



//The default list. Probably don't need them all.
#include "types.h"
#include "output.h"
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

static apop_model * yule_estimate(apop_data * data, apop_model *parameters){
	return apop_maximum_likelihood(data, *parameters);
}

static double beta_greater_than_x_constraint(const apop_data *beta, apop_data *returned_beta, apop_model *v){
    //constraint is 1 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint)constraint= apop_data_calloc(1,1,1);
    apop_data_set(constraint, 0, 0, 1);
    apop_data_set(constraint, 0, -1, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}

static double yule_log_likelihood(const apop_data *beta, apop_data *d, apop_model *p){
  float         bb	            = gsl_vector_get(beta->vector, 0);
  int 		    k;
  float 	    ln_k, ln_bb_k,
		        likelihood 	    = 0,
		        ln_bb		    = gsl_sf_lngamma(bb),
		        ln_bb_less_1    = log(bb-1);
	for (k=0; k< d->matrix->size2; k++)	{
		if (k>=1) 	ln_k	= gsl_sf_lngamma(k+1);
		else		ln_k	= 0;
		ln_bb_k		  = gsl_sf_lngamma(k+1+bb);
        APOP_COL(d, k, v);
		likelihood   += apop_sum(v) *  (ln_k - ln_bb_k);
	}
	likelihood   += (ln_bb_less_1 + ln_bb) * d->matrix->size1 * d->matrix->size2;
	return likelihood;
}

static void yule_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_model *p){
	//Psi is the derivative of the log gamma function.
  float         bb		        = gsl_vector_get(beta->vector, 0);
  gsl_matrix    *data	        = d->matrix;
  int 		    k;
  double	    bb_minus_one_inv= 1/(bb-1),
		        psi_bb	        = gsl_sf_psi(bb),
		        psi_bb_k,
		        d_bb		= 0;
	for (k=0; k< data->size2; k++){
		psi_bb_k= gsl_sf_psi(k +1 + bb);
        APOP_COL(d, k, v);
		d_bb		-= apop_sum(v) * psi_bb_k;
	}
	d_bb   += (bb_minus_one_inv + psi_bb) * data->size1 * data->size2;
	gsl_vector_set(gradient, 0, d_bb);
}


static double yule_rank_p(const apop_data *beta, apop_data *d, apop_model *p){
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
static void yule_rng( double *out, gsl_rng * r, apop_model *a){
double 	e1, e2;
int		x;
	e1	= gsl_ran_exponential(r, 1);
	e2	= gsl_ran_exponential(r, 1);
	x	= GSL_MAX((int) (- e1  / log(1 - exp(-e2 / (*a->parameters->vector->data -1)))), 0);
	*out =  x + 1;	//we rounded down to floor, but want ceil.
}



/** The Yule distribution

The special case of Waring where \f$ \alpha = 0.	\f$<br>

The data set needs to be in rank-form. The first column is the frequency
of the most common item, the second is the frequency of the second most
common item, &c.

apop_yule.estimate() is an MLE, so feed it appropriate \ref apop_params.

\f$ Y(x, b) 	= (b-1) \gamma(b) \gamma(k) / \gamma(k+b)			\f$

\f$ \ln Y(x, b)	= \ln(b-1) + ln\gamma(b) + \ln\gamma(k) - \ln\gamma(k+b)	\f$

\f$ d\ln Y/db	= 1/(b-1)  + \psi(b) - \psi(k+b)				\f$
\ingroup models
\todo I'm pretty sure Wikipedia's specification of the Yule is wrong; I should check and fix when I can check references.
*/
apop_model apop_yule_rank = {"Yule, rank data", 1,0,0, 
	.estimate = yule_estimate, .p = yule_rank_p, .log_likelihood = yule_log_likelihood, 
    .score = yule_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
    .draw = yule_rng};
