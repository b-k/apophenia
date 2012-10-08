/* The gamma distribution.
Copyright (c) 2005--2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_gamma The Gamma distribution

\f$G(x, a, b)     = {1\over (\Gamma(a) b^a)}  x^{a-1} e^{-x/b}\f$

\f$ln G(x, a, b)= -ln \Gamma(a) - a ln b + (a-1)ln(x) + -x/b\f$

\f$d ln G/ da    =  -\psi(a) - ln b + ln(x) \f$    (also, \f$d ln \gamma = \psi\f$)

\f$d ln G/ db    =  -a/b + x/(b^2) \f$

\adoc    Input_format     
Location of data in the grid is not relevant; send it a 1 x N, N x 1, or N x M and it will all be the same.     

\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.
\adoc    Parameter_format   First two elements of the vector.   
\adoc    settings    MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings  
  */

#include "apop_internal.h"

static double gamma_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1 and 0 < beta_2
    return apop_linear_constraint(v->parameters->vector, .margin= 1e-5);
}

typedef struct {double a, b, ln_ga_plus_a_ln_b;} abstruct;

static double apply_for_gamma(double x, void *abin) { 
    abstruct *ab = abin;
    return x ? ((ab->a-1)*log(x) - x/ab->b - ab->ln_ga_plus_a_ln_b) : 0; 
}

static double gamma_log_likelihood(apop_data *d, apop_model *p){
    Nullcheck_mpd(d, p, GSL_NAN) 
    Get_vmsizes(d)
    abstruct ab = {.a = gsl_vector_get(p->parameters->vector, 0),
                   .b = gsl_vector_get(p->parameters->vector, 1) };
    double llikelihood  = 0,
        ln_ga  = gsl_sf_lngamma(ab.a),
        ln_b   = log(ab.b),
        a_ln_b = ab.a * ln_b;
    ab.ln_ga_plus_a_ln_b = ln_ga + a_ln_b;
    llikelihood = apop_map_sum(d, .fn_dp = apply_for_gamma, .param = &ab);
    return llikelihood;
}

static double a_callback(double x, void *ab){ return log(x)- *(double*)ab; }

static double b_callback(double x, void *abv){ 
    double *ab = abv;
    return x/gsl_pow_2(ab[0]) - ab[1]; 
}

static void gamma_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
    Nullcheck_mp(p, ) 
    double  a = gsl_vector_get(p->parameters->vector, 0),
        	b = gsl_vector_get(p->parameters->vector, 1);
    double psi_a_ln_b  = gsl_sf_psi(a) + log(b);
    double b_and_ab[2] = {b, a/b};
    gsl_vector_set(gradient, 0, apop_map_sum(d, .fn_dp = a_callback, .param=&psi_a_ln_b));
    gsl_vector_set(gradient, 1, apop_map_sum(d, .fn_dp = b_callback, .param=&b_and_ab));
}

/* \adoc RNG Just a wrapper for \c gsl_ran_gamma.

See the notes for \ref apop_exponential on a popular alternate form.  */
static void gamma_rng( double *out, gsl_rng* r, apop_model *p){
    *out    = gsl_ran_gamma(r, gsl_vector_get(p->parameters->vector, 0), gsl_vector_get(p->parameters->vector, 1));
}

static double gamma_cdf(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN)
    Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double alpha = gsl_vector_get(params->parameters->vector, 0);
    double beta = gsl_vector_get(params->parameters->vector, 1);
    return gsl_cdf_gamma_P(val, alpha, beta);
}

apop_model apop_gamma = {"Gamma distribution", 2,0,0, .dsize=1, 
      .log_likelihood = gamma_log_likelihood, .score = gamma_dlog_likelihood, 
      .constraint = gamma_constraint, .cdf = gamma_cdf, .draw = gamma_rng};
