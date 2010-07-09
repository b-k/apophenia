/** \file apop_yule.c

  The Yule distribution. A special case of the Waring.*/ 
/*Copyright (c) 2005--2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "types.h"
#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

static double beta_greater_than_x_constraint(apop_data *returned_beta, apop_model *m){
  Nullcheck_m(m); Nullcheck_p(m);
    //constraint is 1 < beta_1
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_calloc(1,1,1);
        apop_data_fill(constraint, 1, 1, 0);
        }
    return apop_linear_constraint(m->parameters->vector, constraint, 1e-4);
}

static double  apply_me(double pt, void *bb){
    double ln_k = (pt>=1) 
                   ? gsl_sf_lngamma(pt)
                   : 0;
    double ln_bb_k = gsl_sf_lngamma(pt+*(double*)bb);
    return ln_k - ln_bb_k;
}

static double  dapply_me(double pt, void *bb){ return -gsl_sf_psi(pt+*(double*)bb); }

static double yule_log_likelihood(apop_data *d, apop_model *m){
  Get_vmsizes(d) //tsize
  Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
    double bb = gsl_vector_get(m->parameters->vector, 0);
    long double ln_bb        = gsl_sf_lngamma(bb),
                ln_bb_less_1 = log(bb-1);
    double      likelihood   = apop_map_sum(d, .fn_dp = apply_me,.param= &bb);
	return likelihood + (ln_bb_less_1 + ln_bb) * tsize;
}

static void yule_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
  Get_vmsizes(d) //tsize
  Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
	//Psi is the derivative of the log gamma function.
    double bb  = gsl_vector_get(m->parameters->vector, 0);
    long double bb_minus_one_inv= 1/(bb-1),
		        psi_bb	        = gsl_sf_psi(bb);
    double d_bb  = apop_map_sum(d, .fn_dp=dapply_me, .param=&bb);
    d_bb += (bb_minus_one_inv + psi_bb) * tsize;
	gsl_vector_set(gradient, 0, d_bb);
}

/** Draw from a Yule distribution with parameter a

Call this fn using <tt> apop_draw(*out, r, apop_yule)</tt>.

\param	a	The parameter.
\param	r	A gsl_rng that you've already set up.

Cribbed from <a href="http://cgm.cs.mcgill.ca/~luc/mbookindex.html>Devroye (1986)</a>, p 553.  */
static void yule_rng( double *out, gsl_rng * r, apop_model *a){
double 	e1, e2;
int		x;
	e1	= gsl_ran_exponential(r, 1);
	e2	= gsl_ran_exponential(r, 1);
	x	= GSL_MAX((int) (- e1  / log(1 - exp(-e2 / (*a->parameters->vector->data -1)))), 0);
	*out =  x + 1;	//we rounded down to floor, but want ceil.
}

apop_model apop_yule = {"Yule", 1,0,0, .dsize=1, .log_likelihood = yule_log_likelihood, 
    .score = yule_dlog_likelihood, .constraint = beta_greater_than_x_constraint, 
    .draw = yule_rng};
//estimate via the default MLE method
