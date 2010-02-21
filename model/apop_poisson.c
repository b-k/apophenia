/** \file apop_poisson.c

  The Poisson distribution.*/
/* Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

static double poisson_log_likelihood(apop_data *d, apop_model *);

static double data_mean(apop_data *d){
    Get_vmsizes(d)
    if (vsize && !msize1) return apop_vector_mean(d->vector);
    if (!vsize && msize1) return apop_matrix_mean(d->matrix);
    return apop_matrix_mean(d->matrix)*(msize1*msize2)/tsize +
               + apop_vector_mean(d->vector)*vsize/tsize;
}

static apop_model * poisson_estimate(apop_data * data,  apop_model *est){
  Nullcheck(data); Nullcheck_m(est);
  double		mean = data_mean(data);
    if (!est->parameters) 
        est->parameters   = apop_data_alloc(1,0,0);
	gsl_vector_set(est->parameters->vector, 0, mean);
    apop_data *info = apop_data_add_page(est->info, apop_data_alloc(1,0,0), "Info");
    apop_data_add_named_elmt(info, "log likelihood", poisson_log_likelihood(data, est));
    //to prevent an infinite loop, the jackknife needs to be flagged to
    //not run itself. We free-ride of the apop_lm_settings struct to signal.
    apop_lm_settings *dummy = apop_settings_get_group(est, "apop_lm");
    if (!dummy || dummy->want_cov=='y'){
        if (!dummy)
            Apop_model_add_group(est, apop_lm);
        Apop_settings_add(est, apop_lm, want_cov, 'n');
        apop_data_add_page(est->parameters, 
                apop_jackknife_cov(data, *est), 
                "Covariance");
    }
	return est;
}

static double beta_zero_greater_than_x_constraint(apop_data *returned_beta, apop_model *v){
    //constraint is 0 < beta_1
    return apop_linear_constraint(v->parameters->vector, .margin = 1e-4);
}


static double apply_me(double x, void *in){
  double    *ln_l = in;
    return x==0 ? 0 :  *ln_l *x - gsl_sf_lngamma(x+1);
}

static double poisson_log_likelihood(apop_data *d, apop_model * p){
  Get_vmsizes(d) //tsize
  Nullcheck(d); Nullcheck_m(p);
  double        lambda      = gsl_vector_get(p->parameters->vector, 0);
  double  ln_l 	= log(lambda);
  double  ll    = apop_map_sum(d, .fn_dp = apply_me, .param=&ln_l);
    return ll - tsize*lambda;
}

static void poisson_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
  Get_vmsizes(d) //tsize
  Nullcheck_v(d); Nullcheck_mv(p);
  double       	lambda  = gsl_vector_get(p->parameters->vector, 0);
  gsl_matrix      *data	= d->matrix;
  double           d_a  = apop_matrix_sum(data)/lambda;
    d_a -= tsize;
    gsl_vector_set(gradient,0, d_a);
}

/* Just a wrapper for gsl_ran_poisson.  */
static void poisson_rng(double *out, gsl_rng* r, apop_model *p){
    *out = gsl_ran_poisson(r, *p->parameters->vector->data);
}

apop_model apop_poisson = {"poisson", 1, 0,0, 
     .estimate = poisson_estimate, .log_likelihood = poisson_log_likelihood, 
     .score = poisson_dlog_likelihood, .constraint = beta_zero_greater_than_x_constraint, 
     .draw = poisson_rng};
