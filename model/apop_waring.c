/** \file apop_waring.c

  The Waring distribution. 

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */


//The default list. Probably don't need them all.
#include "types.h"
#include "stats.h"
#include "model.h"
#include "mapply.h"
#include "settings.h"
#include "conversions.h"
#include "likelihoods.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

//First the rank versions
static double waring_log_likelihood_rank(const apop_data *d, apop_model *m){
  float		      bb	= gsl_vector_get(m->parameters->vector, 0),
    		      a	    = gsl_vector_get(m->parameters->vector, 1);
  int 		      k;
  gsl_matrix      *data	= d->matrix;
  double 		  ln_a_k, ln_bb_a_k,
		          likelihood 	= 0,
		          ln_bb_a		= gsl_sf_lngamma(bb + a),
		          ln_a_mas_1	= gsl_sf_lngamma(a + 1),
		          ln_bb_less_1= log(bb - 1);
	for (k=0; k< data->size2; k++){	//more efficient to go column-by-column
		ln_bb_a_k	 = gsl_sf_lngamma(k +1 + a + bb);
		ln_a_k		 = gsl_sf_lngamma(k +1 + a);
        APOP_COL(d,k, v);
		likelihood   += apop_sum(v) * (ln_a_k - ln_bb_a_k);
	}
    likelihood   +=  (ln_bb_less_1 + ln_bb_a - ln_a_mas_1) * d->matrix->size1 * d->matrix->size2;
	return likelihood;
}

static void waring_dlog_likelihood_rank(const apop_data *d, gsl_vector *gradient, apop_model *m){
  float		      bb		    = gsl_vector_get(m->parameters->vector, 0),
	    	      a		        = gsl_vector_get(m->parameters->vector, 1);
  int 		      k;
  gsl_matrix	  *data		    = d->matrix;
  double		  bb_minus_one_inv= 1/(bb-1),
    		      psi_a_bb	        = gsl_sf_psi(bb + a),
		          psi_a_mas_one	    = gsl_sf_psi(a+1),
		          psi_a_k,
		          psi_bb_a_k,
		          d_bb		        = 0,
		          d_a		            = 0;
	for (k=0; k< data->size2; k++){	//more efficient to go column-by-column
		psi_bb_a_k	 = gsl_sf_psi(k +1 + a + bb);
		psi_a_k		 = gsl_sf_psi(k +1 + a);
        APOP_COL(d, k, v);
		d_bb	    += apop_sum(v) * -psi_bb_a_k;
		d_a		    += apop_sum(v) * (psi_a_k - psi_bb_a_k);
	}
    d_bb	    += (bb_minus_one_inv + psi_a_bb) * d->matrix->size1 * d->matrix->size2;
    d_a		    += (psi_a_bb - psi_a_mas_one) * d->matrix->size1 * d->matrix->size2;
	gsl_vector_set(gradient, 0, d_bb);
	gsl_vector_set(gradient, 1, d_a);
}


static double beta_zero_and_one_greater_than_x_constraint(apop_data *returned_beta, apop_model *m){
    //constraint is 1 < beta_1 and  0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_calloc(2,2,2);
        apop_data_fill(constraint, 1., 1., 0.,
                                   0., 0., 1.);
    }
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-3);
}


double a, bb;
static double apply_me(gsl_vector *data){
  double    likelihood  = 0,
            val, ln_bb_a_k, ln_a_k;
    for (size_t i=0; i< data->size; i++){
        val          = gsl_vector_get(data, i);
        ln_bb_a_k	 = gsl_sf_lngamma(val + a + bb);
        ln_a_k		 = gsl_sf_lngamma(val + a);
        likelihood	+=  ln_a_k - ln_bb_a_k;
    }
    return likelihood;
}

static double waring_log_likelihood(apop_data *d, apop_model *m){
  apop_assert(m->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(m, "apop_rank"))
      return waring_log_likelihood_rank(d, m);
  bb	= gsl_vector_get(m->parameters->vector, 0),
  a	    = gsl_vector_get(m->parameters->vector, 1);
  double 		likelihood 	= 0,
		        ln_bb_a		= gsl_sf_lngamma(bb + a),
                ln_a_mas_1	= gsl_sf_lngamma(a + 1),
                ln_bb_less_1= log(bb - 1);
  gsl_vector * v = apop_matrix_map(d->matrix, apply_me);
    likelihood   = apop_vector_sum(v);
    gsl_vector_free(v);
	likelihood	+= (ln_bb_less_1 + ln_bb_a - ln_a_mas_1)* d->matrix->size1 * d->matrix->size2;
	return likelihood;
}

static void waring_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
	//Psi is the derivative of the log gamma function.
  apop_assert_void(m->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(m, "apop_rank"))
      return waring_dlog_likelihood_rank(d, gradient, m);
  bb		        = gsl_vector_get(m->parameters->vector, 0);
  a		        = gsl_vector_get(m->parameters->vector, 1);
  int 		    i, k;
  gsl_matrix	*data		    = d->matrix;
  double		bb_minus_one_inv= 1/(bb-1), val,
		        psi_a_bb	        = gsl_sf_psi(bb + a),
		        psi_a_mas_one	    = gsl_sf_psi(a+1),
		        psi_a_k,
		        psi_bb_a_k,
		        d_bb		        = 0,
		        d_a		            = 0;
	for (i=0; i< data->size1; i++){
	    for (k=0; k< data->size2; k++){	
            val          = gsl_matrix_get(data, i, k);
		    psi_bb_a_k	 = gsl_sf_psi(val + a + bb);
		    psi_a_k		 = gsl_sf_psi(val + a);
			d_bb	    -= psi_bb_a_k;
			d_a		    += (psi_a_k - psi_bb_a_k);
		}
	}
	d_bb	+= (bb_minus_one_inv + psi_a_bb) * data->size1 * data->size2;
	d_a	    += (psi_a_bb- psi_a_mas_one) * data->size1 * data->size2;
	gsl_vector_set(gradient, 0, d_bb);
	gsl_vector_set(gradient, 1, d_a);
}

/** Give me parameters, and I'll draw a ranking from the appropriate
Waring distribution. [I.e., if I randomly draw from a Waring-distributed
population, return the ranking of the item I just drew.]

Page seven of:
L. Devroye, <a href="http://cgm.cs.mcgill.ca/~luc/digammapaper.ps">Random
variate generation for the digamma and trigamma distributions</a>, Journal
of Statistical Computation and Simulation, vol. 43, pp. 197-216, 1992.
*/
static void waring_rng(double *out, gsl_rng *r, apop_model *eps){
//The key to covnert from Devroye's GHgB3 notation to what I
//consider to be the standard Waring notation in \ref apop_waring:
// a = a + 1
// b = 1 
// c = b - 1 
// n = k - 1 , so if it returns 0, that's first rank.
// OK, I hope that clears everything up.
  double		x, u,
                b   = gsl_vector_get(eps->parameters->vector, 0),
                a   = gsl_vector_get(eps->parameters->vector, 1),
		params[]	={a+1, 1, b-1};
	do{
		x	= apop_rng_GHgB3(r, params);
		u	= gsl_rng_uniform(r);
	} while (u >= (x + a)/(GSL_MAX(a+1,1)*x));
	*out = x+1;
}

/** The Waring distribution
Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.

apop_waring.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.

\f$W(x,k, b,a) 	= (b-1) \gamma(b+a) \gamma(k+a) / [\gamma(a+1) \gamma(k+a+b)]\f$

\f$\ln W(x, b, a) = \ln(b-1) + \ln\gamma(b+a) + \ln\gamma(k+a) - \ln\gamma(a+1) - \ln\gamma(k+a+b)\f$

\f$dlnW/db	= 1/(b-1)  + \psi(b+a) - \psi(k+a+b)\f$

\f$dlnW/da	= \psi(b+a) + \psi(k+a) - \psi(a+1) - \psi(k+a+b)\f$

To specify that you have frequency or ranking data, use 
\code
Apop_settings_add_group(your_model, apop_rank, NULL);
\endcode

\ingroup models
\todo This function needs better testing.
*/
apop_model apop_waring = {"Waring", 2,0,0, 
	 .log_likelihood =  waring_log_likelihood, .score = waring_dlog_likelihood, 
     .constraint =  beta_zero_and_one_greater_than_x_constraint,  .draw = waring_rng};
//estimate via the default MLE 
