/* The Normal and Lognormal distributions.
 Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_normal The Normal (Gaussian) distribution

You know it, it's your attractor in the limit, it's the Gaussian distribution.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2)\f$

\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$

\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$

See also the \ref apop_multivariate_normal.

\adoc    Input_format     I use the elements of the matrix and vector, without regard to their order or position. 
\adoc    Settings   None.
\adoc    Parameter_format  
  As is custom, parameter zero (in the vector) is the mean, parmeter one is the standard deviation (i.e., the square root of the variance). 

\adoc    Predict  Returns the expected value. The <tt>->more</tt>
                 element holds a \ref apop_data set with the title <tt>\<Covariance\></tt>, whose 
                 matrix holds the covariance of the mean.
*/

#include "apop_internal.h"

static long double positive_sigma_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_2
    Staticdef(apop_data *, constraint, apop_data_falloc((1,1,2), 0, 0, 1));
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-5);
}

//This just takes the sum of (x-mu)^2. Using gsl_ran_gaussian_pdf
//would be to calculate log(exp((x-mu)^2)) == slow.
static double apply_me(double x, void *mu){ return x - *(double *)mu; }

static double apply_me2(double x, void *mu){ return gsl_pow_2(x - *(double *)mu); }

static long double normal_log_likelihood(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN);
    Get_vmsizes(d)
    double mu = gsl_vector_get(params->parameters->vector,0);
    double sd = gsl_vector_get(params->parameters->vector,1);
    long double ll  = -apop_map_sum(d, .fn_dp = apply_me2, .param = &mu)/(2*gsl_pow_2(sd));
    ll -= tsize*(M_LNPI+M_LN2+log(sd));
	return ll;
}

/*\adoc estimated_parameters Zeroth vector element is \f$\mu\f$, element 1 is \f$\sigma\f$.
 A page is added named <tt>\<Covariance\></tt> with the 2 \f$\times\f$ 2 covariance matrix for these two parameters
 \adoc estimated_info Reports the log likelihood.*/
static void normal_estimate(apop_data * data, apop_model *est){
    Nullcheck_mpd(data, est, );
    apop_prep(data, est);
    Get_vmsizes(data)
    double mmean=0, mvar=0, vmean=0, vvar=0;
    if (vsize){
        vmean = apop_mean(data->vector);
        vvar = apop_var(data->vector);
    }
    if (msize1) apop_matrix_mean_and_var(data->matrix, &mmean, &mvar);	
    double mean = mmean *(msize1*msize2/(tsize+0.0)) + vmean *(vsize/(tsize+0.0));
    double var = mvar *(msize1*msize2/(tsize+0.0)) + vvar *(vsize/(tsize+0.0));
    est->parameters->vector->data[0] = mean;
    est->parameters->vector->data[1] = sqrt(var);
	apop_name_add(est->parameters->names, "μ", 'r');
	apop_name_add(est->parameters->names, "σ",'r');

    apop_lm_settings *p = apop_settings_get_group(est, apop_lm);
    if (!p) p = Apop_model_add_group(est, apop_lm);
	if (!p || p->want_cov=='y'){
        apop_data *cov = apop_data_get_page(est->parameters, "<Covariance>");
        if (!cov) cov = apop_data_add_page(est->parameters, apop_data_calloc(2, 2), "<Covariance>");
        apop_data_set(cov, 0, 0, mean/tsize);
        apop_data_set(cov, 1, 1, 2*gsl_pow_2(var)/(tsize-1));
    }
    est->data = data;
    apop_data_add_named_elmt(est->info, "log likelihood", normal_log_likelihood(data, est));
}

static long double normal_cdf(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN)
    Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double mu = gsl_vector_get(params->parameters->vector, 0);
    double sd = gsl_vector_get(params->parameters->vector, 1);
    return gsl_cdf_gaussian_P(val-mu, sd);
}

static void normal_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *params){    
    Nullcheck_mpd(d, params, )
    Get_vmsizes(d)
    double mu = gsl_vector_get(params->parameters->vector,0),
           sd = gsl_vector_get(params->parameters->vector,1),
           dll, sll;
    dll = apop_map_sum(d, .fn_dp = apply_me, .param=&mu);
    sll = apop_map_sum(d, .fn_dp = apply_me2, .param=&mu);
    gsl_vector_set(gradient, 0, dll/gsl_pow_2(sd));
    gsl_vector_set(gradient, 1, sll/gsl_pow_3(sd)- tsize /sd);
}

/* \adoc predict Returns the mean, regardless of the input data you give (including
\c NULL). The second page is <tt>\<Covariance\></tt> of the mean.*/ 
apop_data * normal_predict(apop_data *dummy, apop_model *m){
    apop_data *out = apop_data_alloc(1,1);
    out->matrix->data[0] = m->parameters->vector->data[0];

    apop_data *cov = apop_data_get_page(out, "<Covariance>");
    if (!cov) cov = apop_data_add_page(out, apop_data_alloc(1,1), "<Covariance>");
    if (m->data){
           Get_vmsizes(m->data) //tsize
           cov->matrix->data[0] = m->parameters->vector->data[1]/ sqrt(tsize);
    } else cov->matrix->data[0] = 0;
    return out;
}

/*\adoc RNG An apophenia wrapper for the GSL's Normal RNG.

This one asks explicitly for a mean, and the GSL assumes zero and lets you add the mean yourself.  */
static int normal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0];
    return 0;
}

static void normal_prep(apop_data *data, apop_model *params){
    apop_score_vtable_add(normal_dlog_likelihood, apop_normal);
    apop_predict_vtable_add(normal_predict, apop_normal);
    apop_model_clear(data, params);
}

apop_model *apop_normal = &(apop_model){"Normal distribution", 2, 0, 0, .dsize=1,
 .estimate = normal_estimate, .log_likelihood = normal_log_likelihood, 
 .prep = normal_prep, .constraint = positive_sigma_constraint, 
 .draw = normal_rng, .cdf = normal_cdf};


/*\amodel apop_lognormal The Lognormal distribution

The log likelihood function for the lognormal distribution:

\f$f = exp(-(ln(x)-\mu)^2/(2\sigma^2))/ (x\sigma\sqrt{2\pi})\f$
\f$ln f = -(ln(x)-\mu)^2/(2\sigma^2) - ln(x) - ln(\sigma\sqrt{2\pi})\f$

\adoc    Input_format     I use the all elements of the matrix and vector, without regard to their order. 
\adoc    Parameter_format  Zeroth vector element is the mean (after logging); first is the std dev (after logging)    
\adoc    Estimate_results  Parameters are set. Log likelihood is calculated.    
\adoc    settings   None.    
*/

static double lnx_minus_mu_squared(double x, void *mu_in){
	return gsl_pow_2(log(x) - *(double *)mu_in);
}

static long double lognormal_log_likelihood(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN)
    Get_vmsizes(d) //tsize
    double mu = gsl_vector_get(params->parameters->vector, 0);
    double sd = gsl_vector_get(params->parameters->vector, 1);
    long double ll = -apop_map_sum(d, .fn_dp=lnx_minus_mu_squared, .param=&mu);
      ll /= (2*gsl_pow_2(sd));
      ll -= apop_map_sum(d, log, .part='m');
      ll -= tsize*(M_LNPI+M_LN2+log(sd));
	return ll;
}

/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static void lognormal_estimate(apop_data * data, apop_model *est){
    apop_data *cp = apop_data_copy(data);
    apop_prep(data, est);
    Apop_stopif(!cp->matrix && !cp->vector, est->error='d'; return, 
            0, "Neither matrix nor vector in the input data.");
    Get_vmsizes(cp); //vsize, msize1

    double mmean=0, mvar=0, vmean=0, vvar=0;
    if (vsize){
        apop_vector_log(cp->vector);
        vmean = apop_mean(cp->vector);
        vvar = apop_var(cp->vector);
    }
    if (msize2){
        for (int i=0; i< msize2; i++){
            Apop_col_v(cp, i, onecol);
            apop_vector_log(onecol);
        }
        apop_matrix_mean_and_var(cp->matrix, &mmean, &mvar);	
    }
    apop_data_free(cp);
    double mean = mmean *(msize1*msize2/tsize) + vmean *(vsize/tsize);
    double var = mvar *(msize1*msize2/tsize) + vvar *(vsize/tsize);
    est->parameters->vector->data[0] = mean;
    est->parameters->vector->data[1] = sqrt(var);

    apop_name_add(est->parameters->names, "μ", 'r');
    apop_name_add(est->parameters->names, "σ", 'r');
    apop_data_add_named_elmt(est->info, "log likelihood", lognormal_log_likelihood(data, est));
}

static long double lognormal_cdf(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN)
    Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double mu = gsl_vector_get(params->parameters->vector, 0);
    double sd = gsl_vector_get(params->parameters->vector, 1);
    return gsl_cdf_lognormal_P(val, mu, sd);
}

/* \adoc predict Returns the expeced value, \f$E(x) = e^(mu + sigma^2/2)\f$
  in the (0, 0)th element of the matrix, regardless of the input data you give (including \c NULL). */ 
apop_data * lognormal_predict(apop_data *dummy, apop_model *m){
    apop_data *out = apop_data_alloc(1,1);
    out->matrix->data[0] = exp(m->parameters->vector->data[0] 
                                + gsl_pow_2(m->parameters->vector->data[1])/2);
    return out;
}

double diff_sq(double x, void *mu){ return gsl_pow_2(log(x) - *(double*)mu); }

static void lognormal_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *params){    
    double mu = gsl_vector_get(params->parameters->vector,0),
           sd = gsl_vector_get(params->parameters->vector,1);
    Get_vmsizes(d); //tsize
    double dll = apop_map_sum(d, log) - mu*tsize;
    double sll = apop_map_sum(d, .fn_dp=diff_sq, .param=&mu);
    gsl_vector_set(gradient, 0, dll/gsl_pow_2(sd));
    gsl_vector_set(gradient, 1, sll/gsl_pow_3(sd)- tsize/sd);
}

/* \adoc RNG An Apophenia wrapper for the GSL's Normal RNG, exp'ed.  */
static int lognormal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = exp(gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0]);
    return 0;
}

static void lognormal_prep(apop_data *data, apop_model *params){
    apop_score_vtable_add(lognormal_dlog_likelihood, apop_lognormal);
    apop_model_clear(data, params);
}

apop_model *apop_lognormal = &(apop_model){"Lognormal distribution", 2, 0, 0, .dsize=1,
 .estimate = lognormal_estimate, .log_likelihood = lognormal_log_likelihood,
 .prep = lognormal_prep, .constraint = positive_sigma_constraint, 
 .draw = lognormal_rng, .cdf= lognormal_cdf};
