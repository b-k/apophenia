/** \file apop_waring.c

  The Waring distribution.  */
/* Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "mapply.h"
#include "likelihoods.h"

typedef struct {
double a, bb;
} ab_type;

static double waring_apply(double in, void *param, int k){
    ab_type *ab = param;

		double ln_bb_a_k	 = gsl_sf_lngamma(k +1 + ab->a + ab->bb);
		double ln_a_k		 = gsl_sf_lngamma(k +1 + ab->a);
		return in * (ln_a_k - ln_bb_a_k);
}

//First the rank versions
static double waring_log_likelihood_rank(const apop_data *d, apop_model *m){
  ab_type ab;
  ab.bb	= gsl_vector_get(m->parameters->vector, 0),
  ab.a	= gsl_vector_get(m->parameters->vector, 1);
  double 		  likelihood,
		          ln_bb_a		= gsl_sf_lngamma(ab.bb + ab.a),
		          ln_a_mas_1	= gsl_sf_lngamma(ab.a + 1),
		          ln_bb_less_1= log(ab.bb - 1);
        likelihood = apop_map_sum((apop_data*)d, .fn_dpi = waring_apply, .part='c');
    likelihood   +=  (ln_bb_less_1 + ln_bb_a - ln_a_mas_1) * d->matrix->size1 * d->matrix->size2;
	return likelihood;
}

static void waring_dlog_likelihood_rank(const apop_data *d, gsl_vector *gradient, apop_model *m){
  double	      bb		    = gsl_vector_get(m->parameters->vector, 0),
	    	      a		        = gsl_vector_get(m->parameters->vector, 1);
  gsl_matrix	  *data		    = d->matrix;
  double		  bb_minus_one_inv= 1/(bb-1),
    		      psi_a_bb	        = gsl_sf_psi(bb + a),
		          psi_a_mas_one	    = gsl_sf_psi(a+1),
		          psi_a_k,
		          psi_bb_a_k,
		          d_bb		        = 0,
		          d_a		            = 0;
	for (size_t k=0; k< data->size2; k++){	//more efficient to go column-by-column
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
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-4);
}

static double apply_me(double val, void *in){
    ab_type *ab = in;
        double ln_a_k		 = gsl_sf_lngamma(val + ab->a);
        double ln_bb_a_k	 = gsl_sf_lngamma(val + ab->a + ab->bb);
        return  ln_a_k - ln_bb_a_k;
}

static double waring_log_likelihood(apop_data *d, apop_model *m){
  apop_assert(m->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(m, "apop_rank"))
      return waring_log_likelihood_rank(d, m);
  ab_type abstruct;
  abstruct.bb	= gsl_vector_get(m->parameters->vector, 0),
  abstruct.a	    = gsl_vector_get(m->parameters->vector, 1);
  double 		likelihood 	= 0,
		        ln_bb_a		= gsl_sf_lngamma(abstruct.bb + abstruct.a),
                ln_a_mas_1	= gsl_sf_lngamma(abstruct.a + 1),
                ln_bb_less_1= log(abstruct.bb - 1);
    likelihood= apop_map_sum(d, .fn_dp = apply_me, .param=&abstruct, .part='m');
	likelihood	+= (ln_bb_less_1 + ln_bb_a - ln_a_mas_1)* d->matrix->size1 * d->matrix->size2;
	return likelihood;
}

static void waring_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
	//Psi is the derivative of the log gamma function.
  apop_assert_void(m->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
    if (apop_settings_get_group(m, "apop_rank"))
      return waring_dlog_likelihood_rank(d, gradient, m);
  double bb		        = gsl_vector_get(m->parameters->vector, 0);
  double a		        = gsl_vector_get(m->parameters->vector, 1);
  gsl_matrix	*data		    = d->matrix;
  double		bb_minus_one_inv= 1/(bb-1), val,
		        psi_a_bb	        = gsl_sf_psi(bb + a),
		        psi_a_mas_one	    = gsl_sf_psi(a+1),
		        psi_a_k,
		        psi_bb_a_k,
		        d_bb		        = 0,
		        d_a		            = 0;
	for (size_t i=0; i< data->size1; i++){
	    for (size_t k=0; k< data->size2; k++){	
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
                a   = gsl_vector_get(eps->parameters->vector, 1);
  double		params[]	={a+1, 1, b-1};
/*	do{
		x	= apop_rng_GHgB3(r, params);
		u	= gsl_rng_uniform(r);
	} while (u >= (x + a)/(GSL_MAX(a+1,1)*x));
	*out = x+1;*/
	do{
		x	= apop_rng_GHgB3(r, params)+1;
		u	= gsl_rng_uniform(r);
	} while (u >= (x + a)/(GSL_MAX(a+1,1)*x));
	*out = x;
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

\hideinitializer
\ingroup models
\todo This function needs better testing.
*/
apop_model apop_waring = {"Waring", 2,0,0, 
	 .log_likelihood =  waring_log_likelihood, .score = waring_dlog_likelihood, 
     .constraint =  beta_zero_and_one_greater_than_x_constraint,  .draw = waring_rng};
//estimate via the default MLE 
