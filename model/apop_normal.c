/* The Normal and Lognormal distributions.
 Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_normal The Normal (Gaussian) distribution

You know it, it's your attractor in the limit, it's the Gaussian distribution.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2)\f$

\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$

\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$

\adoc    Input_format     I use the elements of the matrix, without regard to their order. 
\adoc    Parameter_format  
  As is custom, the first parameter (in the vector) is the mean, the second is the standard deviation (i.e., the square root of the variance). 

\adoc    Predict  Returns the expected value. The <tt>->more</tt>
                 element holds a \ref apop_data set with the title <tt>"Covariance"</tt>, whose matrix holds the covariance of the mean. If <tt>estimated_model->data</tt> then this is going to be \f$\mu/\sqrt{n}\f$; if the model's <tt>data == NULL</tt> then cov = 0. 
                 
                 Format subject to change. 
\adoc    settings   \ref apop_parts_wanted ; notably the \c want_cov element,
which refers to the covariance matrix for the mean and variance. You can set
that to \c 'n' if you are estimating millions of these and need to save time.
*/

#include "model.h"
#include "mapply.h"
#include "settings.h"
#include "internal.h"
#include "likelihoods.h"

static double beta_1_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint) {
        constraint = apop_data_calloc(1,1,2);
        apop_data_set(constraint, 0, 1, 1);
    }
    return apop_linear_constraint(v->parameters->vector, constraint, 1e-5);
}

//This just takes the sum of (x-mu)^2. Using gsl_ran_gaussian_pdf
//would be to calculate log(exp((x-mu)^2)) == slow.
static double apply_me(double x, void *mu){ return x - *(double *)mu; }
static double apply_me2(double x, void *mu){ return gsl_pow_2(x - *(double *)mu); }

/* The log likelihood function for the Normal.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-(x-\mu)^2 / 2\sigma^2)\f$
\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$
*/
static double normal_log_likelihood(apop_data *d, apop_model *params){
  assert(params->parameters);
  Get_vmsizes(d)
    double mu = gsl_vector_get(params->parameters->vector,0);
    double sd = gsl_vector_get(params->parameters->vector,1);
  long double   ll  = -apop_map_sum(d, .fn_dp = apply_me2, .param = &mu)/(2*gsl_pow_2(sd));
    ll    -=  tsize*(M_LNPI+M_LN2+log(sd));
	return ll;
}

/*\adoc estimated_parameters Zeroth vector element is \f$\mu\f$, element 1 is \f$\sigma\f$.
 A page is added named <tt>\<Covariance\></tt> with the 2 \f$\times\f$ 2 covariance matrix for these two parameters
 \adoc estimated_info Reports the log likelihood.*/
static apop_model * normal_estimate(apop_data * data, apop_model *est){
    Get_vmsizes(data)
  double		mmean=0, mvar=0, vmean=0, vvar=0;
  apop_lm_settings *p = apop_settings_get_group(est, apop_lm);
    if (!p) 
        p = Apop_model_add_group(est, apop_lm);
    if (vsize){
        vmean = apop_mean(data->vector);
        vvar = apop_var(data->vector);
    }
    if (msize1)
        apop_matrix_mean_and_var(data->matrix, &mmean, &mvar);	
    double mean = mmean *(msize1*msize2/tsize) + vmean *(vsize/tsize);
    double var = mvar *(msize1*msize2/tsize) + vvar *(vsize/tsize);
    est->parameters->vector->data[0] = mean;
    est->parameters->vector->data[1] = sqrt(var);
	apop_name_add(est->parameters->names, "mu", 'r');
	apop_name_add(est->parameters->names, "sigma",'r');
	if (!p || p->want_cov=='y'){
        apop_data *cov = apop_data_get_page(est->parameters, "<Covariance>");
        if (!cov)
            cov = apop_data_add_page(est->parameters, apop_data_calloc(0, 2, 2), "<Covariance>");
        apop_data_set(cov, 0, 0, mean/tsize);
        apop_data_set(cov, 1, 1, 2*gsl_pow_2(var)/(tsize-1));
    }
    est->data = data;
    est->info = apop_data_alloc();
    apop_data_add_named_elmt(est->info, "log likelihood", normal_log_likelihood(data, est));
	return est;
}

static double normal_cdf(apop_data *d, apop_model *params){
  Nullcheck_m(params) Nullcheck_p(params) Nullcheck_d(d) 
  Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double mu = gsl_vector_get(params->parameters->vector, 0);
    double sd = gsl_vector_get(params->parameters->vector, 1);
    return gsl_cdf_gaussian_P(val-mu, sd);
}

/** Gradient of the log likelihood function

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$
\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$
\f$d\ln N(\mu,\sigma)/d\sigma = ((x-\mu)^2 / \sigma^3) - 1/\sigma \f$
 */
static void normal_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *params){    
  Get_vmsizes(d)
  double  mu      = gsl_vector_get(params->parameters->vector,0),
          sd      = gsl_vector_get(params->parameters->vector,1),
          dll, sll;
    dll = apop_map_sum(d, .fn_dp = apply_me, .param=&mu);
    sll = apop_map_sum(d, .fn_dp = apply_me2, .param=&mu);
    gsl_vector_set(gradient, 0, dll/gsl_pow_2(sd));
    //gsl_vector_set(gradient, 1, sll/(2*gsl_pow_2(ss))- data->size1 * data->size2 * 0.5/ss);
    gsl_vector_set(gradient, 1, sll/gsl_pow_3(sd)- tsize /sd);
}

//Just the mean and the variance of the mean.
apop_data * normal_predict(apop_data *dummy, apop_model *m){
    apop_data *out = apop_data_alloc(0,1,1);
    out->matrix->data[0] = m->parameters->vector->data[0];

    apop_data *cov = apop_data_get_page(out, "<Covariance>");
    if (!cov)
        cov = apop_data_add_page(out, apop_data_alloc(0,1,1), "<Covariance>");
    if (m->data){
        Get_vmsizes(m->data) //tsize
        cov->matrix->data[0] = m->parameters->vector->data[1]/ sqrt(tsize);
    } else
        cov->matrix->data[0] = 0;
    return out;
}

/*\adoc RNG An apophenia wrapper for the GSL's Normal RNG.

Two differences: this one asks explicitly for a mean, and the GSL assumes zero and
makes you add the mean yourself; Apophenia tends to prefer the variance (\f$\sigma^2\f$)
wherever possible, while the GSL uses the standard deviation here (\f$\sigma\f$) */
static void normal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0];
}

apop_model apop_normal = {"Normal distribution", 2, 0, 0, .dsize=1,
 .estimate = normal_estimate, .log_likelihood = normal_log_likelihood, 
 .score = normal_dlog_likelihood, 
 .constraint = beta_1_greater_than_x_constraint, .draw = normal_rng, 
 .cdf = normal_cdf, .predict = normal_predict};


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

static double lognormal_log_likelihood(apop_data *d, apop_model *params){
    Get_vmsizes(d) //tsize
    Nullcheck(params)
    Nullcheck_p(params)
    double mu	   = gsl_vector_get(params->parameters->vector, 0);
    double sd      = gsl_vector_get(params->parameters->vector, 1);
    long double ll = -apop_map_sum(d, .fn_dp=lnx_minus_mu_squared, .param=&mu);
      ll /= (2*gsl_pow_2(sd));
      ll -= apop_map_sum(d, log, .part='m');
      ll -= tsize*(M_LNPI+M_LN2+log(sd));
	return ll;
}

/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static apop_model * lognormal_estimate(apop_data * data, apop_model *parameters){
  apop_model 	*est = apop_model_copy(*parameters);
  double   mean    = 0,
           var     = 0; 
  apop_lm_settings *p = apop_settings_get_group(est, apop_lm);
    if (!p) 
        p = Apop_model_add_group(est, apop_lm);
    apop_matrix_mean_and_var(data->matrix, &mean, &var);
    if (!est->parameters)
        est->parameters = apop_data_alloc(2);
    double sigsq   = log(1+ var/gsl_pow_2(mean));
	gsl_vector_set(est->parameters->vector, 0, log(mean)- sigsq/2);
	gsl_vector_set(est->parameters->vector, 1, sqrt(sigsq));
    est->data           = data;
    parameters->info = apop_data_alloc();
    apop_data_add_named_elmt(parameters->info, "log likelihood", lognormal_log_likelihood(data, parameters));
	return est;
}


/*\adoc CDF Yes. */
static double lognormal_cdf(apop_data *d, apop_model *params){
  Nullcheck_m(params) Nullcheck_p(params) Nullcheck_d(d) 
  Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double mu = gsl_vector_get(params->parameters->vector, 0);
    double sd = gsl_vector_get(params->parameters->vector, 1);
    return gsl_cdf_lognormal_P(val, mu, sd);
}

apop_data * lognormal_predict(apop_data *dummy, apop_model *m){
    //E(x) = e^(mu + sigma^2/2)
    apop_data *out = apop_data_alloc(1,1);
    out->matrix->data[0] = exp(m->parameters->vector->data[0] 
                                + gsl_pow_2(m->parameters->vector->data[1])/2);
    return out;
}

/* This is copied from the Normal. The first one who needs it gets to
 * fix it. 
static void lognormal_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *params){    
              mu      = gsl_vector_get(params->parameters->vector,0);
  double      sd      = gsl_vector_get(params->parameters->vector,1),
              dll     = 0,
              sll     = 0,
              x;
  int         i,j;
  gsl_matrix  *data   = d->matrix;
    for (i=0;i< data->size1; i++)
        for (j=0;j< data->size2; j++){
            x    = gsl_matrix_get(data, i, j);
            dll += (x - mu);
            sll += gsl_pow_2(x - mu);
        }
    gsl_vector_set(gradient, 0, dll/gsl_pow_2(sd));
    //gsl_vector_set(gradient, 1, sll/(2*gsl_pow_2(ss))- data->size1 * data->size2 * 0.5/ss);
    gsl_vector_set(gradient, 1, sll/gsl_pow_3(sd)- data->size1 * data->size2 /sd);
}
*/

/* \adoc RNG An Apophenia wrapper for the GSL's Normal RNG, logged.

\param r	a gsl_rng already allocated
\param *out	To where I will write the drawn number
\param *p   A pointer to the model.
 */
static void lognormal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = exp(gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0]);
}

apop_model apop_lognormal = {"Lognormal distribution", 2, 0, 0, .dsize=1,
 .estimate = lognormal_estimate, .log_likelihood = lognormal_log_likelihood, /*.score = lognormal_dlog_likelihood,*/ 
 .constraint = beta_1_greater_than_x_constraint, .draw = lognormal_rng,
  .cdf= lognormal_cdf};
