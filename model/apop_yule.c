/* The Yule distribution. A special case of the Waring.

Copyright (c) 2005--2007, 2009, 2011 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2. 

\amodel apop_yule
The special case of the \ref apop_waring "Waring" where \f$ \alpha = 0.	\f$<br>

\f$ Y(x, b) 	= (b-1) \gamma(b) \gamma(k) / \gamma(k+b)			\f$

\f$ \ln Y(x, b)	= \ln(b-1) + ln\gamma(b) + \ln\gamma(k) - \ln\gamma(k+b)	\f$

\f$ d\ln Y/db	= 1/(b-1)  + \psi(b) - \psi(k+b)				\f$

apop_yule.estimate() is an MLE, so feed it appropriate \ref apop_mle_settings.

\adoc    Input_format     
Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.

See also \ref apop_data_rank_compress for means of dealing with one more input data format.
\adoc    Parameter_format  One element at the top of the parameter set's vector.
\adoc    settings   MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings    */

#include "apop_internal.h"

static long double yule_constraint(apop_data *returned_beta, apop_model *m){
  Nullcheck_mp(m, GSL_NAN);
    //constraint is 1 < beta_1
  Staticdef(apop_data *, constraint, apop_data_falloc((1,1,1), 1, 1));
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-4);
}

static double apply_me(double pt, void *bb){
    double ln_k = (pt>=1) 
                   ? gsl_sf_lngamma(pt)
                   : 0;
    double ln_bb_k = gsl_sf_lngamma(pt+*(double*)bb);
    return ln_k - ln_bb_k;
}

static double dapply_me(double pt, void *bb){ return -gsl_sf_psi(pt+*(double*)bb); }

static long double yule_log_likelihood(apop_data *d, apop_model *m){
  Nullcheck_mpd(d, m, GSL_NAN);
  Get_vmsizes(d) //tsize
    double bb = gsl_vector_get(m->parameters->vector, 0);
    long double ln_bb        = gsl_sf_lngamma(bb),
                ln_bb_less_1 = log(bb-1);
    double      likelihood   = apop_map_sum(d, .fn_dp = apply_me,.param= &bb);
	return likelihood + (ln_bb_less_1 + ln_bb) * tsize;
}

static void yule_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
  Nullcheck_mpd(d, m, );
  Get_vmsizes(d) //tsize
	//Psi is the derivative of the log gamma function.
    double bb  = gsl_vector_get(m->parameters->vector, 0);
    long double bb_minus_one_inv= 1/(bb-1),
		        psi_bb	        = gsl_sf_psi(bb);
    double d_bb  = apop_map_sum(d, .fn_dp=dapply_me, .param=&bb);
    d_bb += (bb_minus_one_inv + psi_bb) * tsize;
	gsl_vector_set(gradient, 0, d_bb);
}

/* \adoc RNG Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 553.  */
static int yule_rng( double *out, gsl_rng * r, apop_model *a){
	double e1 = gsl_ran_exponential(r, 1);
	double e2 = gsl_ran_exponential(r, 1);
	int x = GSL_MAX((int) (- e1  / log(1 - exp(-e2 / (*a->parameters->vector->data -1)))), 0);
	*out =  x + 1;	//we rounded down to floor, but want ceil.
    return 0;
}

static void yule_prep(apop_data *data, apop_model *params){
    apop_score_vtable_add(yule_dlog_likelihood, apop_yule);
    apop_model_clear(data, params);
}

apop_model *apop_yule = &(apop_model){"Yule distribution", 1,0,0, .dsize=1, .log_likelihood = yule_log_likelihood, 
    .prep = yule_prep, .constraint = yule_constraint, .draw = yule_rng};
