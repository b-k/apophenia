/* \file apop_beta.c  The Beta distribution 
Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_beta The Beta distribution.

The beta distribution has two parameters and is restricted between zero and one. You
may also find \ref apop_beta_from_mean_var to be useful.

\adoc    Input_format  Any arrangement of scalar values. 
\adoc    Parameter_format   a vector, v[0]=\f$\alpha\f$; v[1]=\f$\beta\f$    
\adoc    RNG  Produces a scalar \f$\in[0,1]\f$. 
\adoc    settings None.  */

#include "apop_internal.h"

static double beta_log_likelihood(apop_data *d, apop_model *p);

/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static apop_model * beta_estimate(apop_data * data,  apop_model *est){
    Nullcheck_mpd(data, est, NULL);
    Get_vmsizes(data) //vsize, msize1,...
    double		mmean=0, mvar=0, vmean=0, vvar=0, alpha, beta;
    if (vsize){
        vmean = apop_mean(data->vector);
        vvar = apop_var(data->vector);
    }
    if (msize1)
        apop_matrix_mean_and_var(data->matrix, &mmean, &mvar);	
    double mean = mmean *(msize1*msize2/tsize) + vmean *(vsize/tsize);
    double var = mvar *(msize1*msize2/tsize) + vvar *(vsize/tsize);
    apop_data_add_names(est->parameters, 'r', "α", "β");
    alpha   = gsl_pow_2(mean) * ((1-mean)/var - 1/mean);
    beta    = alpha * (1-mean)/mean;
	gsl_vector_set(est->parameters->vector, 0, alpha);
	gsl_vector_set(est->parameters->vector, 1, beta);
    apop_data_add_named_elmt(est->info, "log likelihood", beta_log_likelihood(data, est));
    //apop_numerical_covariance_matrix(apop_beta, est, data);
	return est;
}

typedef struct{
    double alpha, beta; 
} ab_type;

static double betamap(double x, void *abin) {
    ab_type *ab = abin; 
    return (x < 0 || x > 1) ? 0
                : (ab->alpha-1) * log(x) + (ab->beta-1) *log(1-x); 
}

#define Get_ab(p) \
    ab_type ab = { .alpha = apop_data_get(p->parameters,0,-1), \
                   .beta  = apop_data_get(p->parameters,1,-1) };

static double beta_log_likelihood(apop_data *d, apop_model *p){
    Nullcheck_mpd(d, p, GSL_NAN); 
    Get_vmsizes(d) //tsize
    Get_ab(p) //ab
    Apop_stopif(isnan(ab.alpha+ab.beta), return GSL_NAN, 0, "NaN α or β input.");
	return apop_map_sum(d, .fn_dp = betamap, .param=&ab) - gsl_sf_lnbeta(ab.alpha, ab.beta) * tsize;
}

static double dbeta_callback(double x){ return log(1-x); }

static void beta_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    Nullcheck_mpd(d, m, )
    Get_vmsizes(d) //tsize
    Get_ab(m) //ab
    double lnsum = apop_map_sum(d, log);
    double ln_x_minus_1_sum = apop_map_sum(d, dbeta_callback);
	//Psi is the derivative of the log gamma function.
	gsl_vector_set(gradient, 0, lnsum  + (gsl_sf_psi(ab.alpha + ab.beta) - gsl_sf_psi(ab.alpha))*tsize);
	gsl_vector_set(gradient, 1, ln_x_minus_1_sum  + (gsl_sf_psi(ab.alpha + ab.beta) - gsl_sf_psi(ab.beta))*tsize);
}

static double beta_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1 and  0 < beta_2
    return apop_linear_constraint(v->parameters->vector, .margin= 1e-4);
}

static double beta_cdf(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN)
    Get_vmsizes(d)  //vsize
    Get_ab(params)
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    return gsl_cdf_beta_P(val, ab.alpha, ab.beta);
}

static void beta_rng(double *out, gsl_rng *r, apop_model* eps){
    Nullcheck_mp(eps, )
    Get_ab(eps)
    do {
    *out = gsl_ran_beta(r, ab.alpha, ab.beta);
    } while (*out <= 0 || *out >= 1);
}

apop_model apop_beta = {"Beta distribution", 2,0,0, .dsize=1, .estimate = beta_estimate, 
    .log_likelihood = beta_log_likelihood, .score = beta_dlog_likelihood, 
    .constraint = beta_constraint, .draw = beta_rng, .cdf = beta_cdf};
