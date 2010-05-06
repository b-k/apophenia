/** \file apop_gamma.c

  The gamma distribution.*/
/*Copyright (c) 2005--2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "types.h"
#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

static double gamma_rank_log_likelihood(apop_data *d, apop_model *p){
  Nullcheck_p(p) 
  float           a           = gsl_vector_get(p->parameters->vector, 0),
                  b           = gsl_vector_get(p->parameters->vector, 1);
    apop_assert(a>0 && b>0, 0, 0,'c', "The Gamma's log likelihood needs positive params, and you gave me %g and %g. Returning zero.", a, b);
    if (gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;    
  gsl_matrix      *data       = d->matrix;
  double          llikelihood = 0,
                  ln_ga       = gsl_sf_lngamma(a),
                  a_ln_b      = a * log(b),
                  ln_k;
  gsl_vector      v;
    for (size_t k=0; k< data->size2; k++){
            ln_k            = log(k+1);
            v               = gsl_matrix_column(data, k).vector;
            llikelihood    += apop_sum(&v) * (-ln_ga - a_ln_b + (a-1) * ln_k  - (k+1)/b);
        }
    return llikelihood;
}

static void gamma_rank_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  Nullcheck_p(p) 
  double          a       = gsl_vector_get(p->parameters->vector, 0),
                  b       = gsl_vector_get(p->parameters->vector, 1);
  gsl_matrix     *data    = d->matrix;
  double          d_a     = 0,
                  d_b     = 0,
                  psi_a   = gsl_sf_psi(a),
                  ln_b    = log(b),
                  x, ln_k;
    for (size_t k=0; k< data->size2; k++){
        ln_k    = log(k +1);
        Apop_col(d, k, v);
        x       = apop_sum(v);
        d_a    += x * (-psi_a - ln_b + ln_k);
        d_b    += x * (a/b - (k+1)/gsl_pow_2(b));
    }
    gsl_vector_set(gradient,0, d_a);
    gsl_vector_set(gradient,1, d_b);
}

static double beta_zero_and_one_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1 and 0 < beta_2
    return apop_linear_constraint(v->parameters->vector, .margin= 1e-5);
}

typedef struct {double a, b, ln_ga_plus_a_ln_b;} abstruct;

static double apply_for_gamma(double x, void *abin) { 
    abstruct *ab = abin;
     return x ? ((ab->a-1)*log(x) - x/ab->b - ab->ln_ga_plus_a_ln_b) : 0; 
}

static double gamma_log_likelihood(apop_data *d, apop_model *p){
  Nullcheck_m(p) 
  Nullcheck_p(p) 
  Get_vmsizes(d)
    if (apop_settings_get_group(p, apop_rank))
        return gamma_rank_log_likelihood(d, p);
  abstruct ab = {
      .a    = gsl_vector_get(p->parameters->vector, 0),
      .b    = gsl_vector_get(p->parameters->vector, 1)
  };
  double        llikelihood  = 0,
        ln_ga 	= gsl_sf_lngamma(ab.a),
        ln_b	= log(ab.b),
        a_ln_b	= ab.a * ln_b;
    ab.ln_ga_plus_a_ln_b = ln_ga + a_ln_b;
    llikelihood    = apop_map_sum(d, .fn_dp = apply_for_gamma, .param = &ab);
    return llikelihood;
}

//static double a_callback(double x, void *ab){ return x ? log(x)- *(double*)ab : 0; }
//static double b_callback(double x, void *ab){ return x ? -x - *(double*)ab : 0; }
static double a_callback(double x, void *ab){ return log(x)- *(double*)ab; }
static double b_callback(double x, void *abv){ 
    double *ab = abv;
    return x/gsl_pow_2(ab[0]) - ab[1]; 
}

static void gamma_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  Nullcheck_p(p) 
    if (apop_settings_get_group(p, apop_rank))
       return gamma_rank_dlog_likelihood(d, gradient, p);
  double       	a    	= gsl_vector_get(p->parameters->vector, 0),
        		b    	= gsl_vector_get(p->parameters->vector, 1);
    //if (a <= 0 || b <= 0 || gsl_isnan(a) || gsl_isnan(b)) return GSL_POSINF;    
                        //a sign to the minimizer to look elsewhere.
  double psi_a_ln_b	= gsl_sf_psi(a) + log(b);
  double b_and_ab[2]    = {b, a/b};
    gsl_vector_set(gradient, 0, apop_map_sum(d, .fn_dp = a_callback, .param=&psi_a_ln_b));
    gsl_vector_set(gradient, 1, apop_map_sum(d, .fn_dp = b_callback, .param=&b_and_ab));
}

/* Just a wrapper for gsl_ran_gamma.

   cut & pasted from the GSL documentation:
\f$          p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx \f$

See the notes for \ref apop_exponential on a popular alternate form.
*/
static void gamma_rng( double *out, gsl_rng* r, apop_model *p){
    *out    = gsl_ran_gamma(r, gsl_vector_get(p->parameters->vector, 0), gsl_vector_get(p->parameters->vector, 1));
}

static double gamma_cdf(apop_data *d, apop_model *params){
  Nullcheck_m(params) Nullcheck_p(params) Nullcheck_d(d) 
  Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double alpha = gsl_vector_get(params->parameters->vector, 0);
    double beta = gsl_vector_get(params->parameters->vector, 1);
    return gsl_cdf_gamma_P(val, alpha, beta);
}

apop_model apop_gamma = {"Gamma distribution", 2,0,0, //estimate method is just the default MLE.
      .dsize=1, .log_likelihood = gamma_log_likelihood, 
     .score = gamma_dlog_likelihood, .constraint = beta_zero_and_one_greater_than_x_constraint, 
     .cdf = gamma_cdf, .draw = gamma_rng};
