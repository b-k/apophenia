/** \file apop_waring.c

  The Waring distribution.  */
/* Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

/** RNG from a Generalized Hypergeometric type B3.

 Devroye uses this as the base for many of his
 distribution-generators, including the Waring.
*/  //Header in stats.h
double apop_rng_GHgB3(gsl_rng * r, double* a){
if ((a[0]<=0) || (a[1] <= 0) || (a[2] <=0)){
	printf("apop_GHgB3_rng took a zero parameter; bad.\n");
	return 0;
	}
double		aa	= gsl_ran_gamma(r, a[0], 1),
		b	= gsl_ran_gamma(r, a[1], 1),
		c	= gsl_ran_gamma(r, a[2], 1);
int		p	= gsl_ran_poisson(r, aa*b/c);
	return p;
}

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
  Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
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
  Nullcheck_v(d); Nullcheck_mv(m); Nullcheck_pv(m);
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
  Get_vmsizes(d) //tsize
  Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
    if (apop_settings_get_group(m, apop_rank))
      return waring_log_likelihood_rank(d, m);
  ab_type abstruct;
  abstruct.bb	= gsl_vector_get(m->parameters->vector, 0),
  abstruct.a	    = gsl_vector_get(m->parameters->vector, 1);
  double 		likelihood 	= 0,
		        ln_bb_a		= gsl_sf_lngamma(abstruct.bb + abstruct.a),
                ln_a_mas_1	= gsl_sf_lngamma(abstruct.a + 1),
                ln_bb_less_1= log(abstruct.bb - 1);
    likelihood= apop_map_sum(d, .fn_dp = apply_me, .param=&abstruct);
	likelihood	+= (ln_bb_less_1 + ln_bb_a - ln_a_mas_1)* tsize;
	return likelihood;
}

static void waring_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
	//Psi is the derivative of the log gamma function.
  Get_vmsizes(d) //tsize
  Nullcheck_v(d); Nullcheck_mv(m); Nullcheck_pv(m);
  int min = vsize ? -1 : 0;
  int max = msize2 ? msize2 : 0;
    if (apop_settings_get_group(m, apop_rank))
      return waring_dlog_likelihood_rank(d, gradient, m);
  double bb		        = gsl_vector_get(m->parameters->vector, 0);
  double a		        = gsl_vector_get(m->parameters->vector, 1);
  double		bb_minus_one_inv= 1/(bb-1), val,
		        psi_a_bb	        = gsl_sf_psi(bb + a),
		        psi_a_mas_one	    = gsl_sf_psi(a+1),
		        psi_a_k,
		        psi_bb_a_k,
		        d_bb		        = 0,
		        d_a		            = 0;
	for (size_t i=0; i< GSL_MAX(msize1, vsize); i++){
	    for (int k=min; k< max; k++){	
            val          = apop_data_get(d, i, k);
		    psi_bb_a_k	 = gsl_sf_psi(val + a + bb);
		    psi_a_k		 = gsl_sf_psi(val + a);
			d_bb	    -= psi_bb_a_k;
			d_a		    += (psi_a_k - psi_bb_a_k);
		}
	}
	d_bb	+= (bb_minus_one_inv + psi_a_bb) *tsize;
	d_a	    += (psi_a_bb- psi_a_mas_one) * tsize;
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
//The key to convert from Devroye's GHgB3 notation to what I
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

apop_model apop_waring = {"Waring", 2,0,0, .dsize=1,
	 .log_likelihood =  waring_log_likelihood, .score = waring_dlog_likelihood, 
     .constraint =  beta_zero_and_one_greater_than_x_constraint,  .draw = waring_rng};
//estimate via the default MLE 
