/** \file apop_normal.c

The Normal and Lognormal distributions.*/
/* Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "settings.h"
#include "internal.h"
#include "likelihoods.h"
static double normal_log_likelihood(apop_data *d, apop_model *params);

/** The normal estimate */
static apop_model * normal_estimate(apop_data * data, apop_model *parameters){
    Get_vmsizes(data)
  double		mmean=0, mvar=0, vmean=0, vvar=0;
  apop_model 	*est = apop_model_copy(*parameters);
  apop_ls_settings *p = apop_settings_get_group(est, "apop_ls");
    if (!p) {
        Apop_model_add_group(est, apop_ls);
        p = apop_settings_get_group(est, "apop_ls");
    }
    if (vsize){
        vmean = apop_mean(data->vector);
        vvar = apop_var(data->vector);
    }
    if (msize1)
        apop_matrix_mean_and_var(data->matrix, &mmean, &mvar);	
    double mean = mmean *(msize1*msize2/tsize) + vmean *(vsize/tsize);
    double var = mvar *(msize1*msize2/tsize) + vvar *(vsize/tsize);
    apop_model_clear(data, est);
	apop_data_add_named_elmt(est->parameters,"mu", mean);
	apop_data_add_named_elmt(est->parameters,"sigma", sqrt(var));
    est->llikelihood	= normal_log_likelihood(data, est);
	if (!p || p->want_cov=='y'){
        est->covariance   = apop_data_calloc(0, 2, 2);
        apop_data_set(est->covariance, 0, 0, mean/tsize);
        apop_data_set(est->covariance, 1, 1, 2*gsl_pow_2(var)/(tsize-1));
    }
    est->data = data;
	return est;
}

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

The log likelihood function and dlog likelihood don't care about your
rows of data; if you have an 8 x 7 data set, it will give you the log
likelihood of those 56 observations given the mean and variance (i.e.,
\f$\sigma^2\f$, not std deviation=\f$\sigma\f$) you provide.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-(x-\mu)^2 / 2\sigma^2)\f$
\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\param beta	beta[0]=the mean; beta[1]=the variance
\param d	the set of data points; see notes.
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


/** An apophenia wrapper for the GSL's Normal RNG.

Two differences: this one asks explicitly for a mean, and the GSL
assumes zero and makes you add the mean yourself; Apophenia tends to
prefer the variance (\f$\sigma^2\f$) wherever possible, while the GSL
uses the standard deviation here (\f$\sigma\f$)

\param r	a gsl_rng already allocated
\param *out	To where I will write the drawn number
\param *p   A pointer to the model.
 */
static void normal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0];
}

/** You know it, it's your attractor in the limit, it's the Gaussian distribution.

  As is custom, the first parameter is the mean, the second is the standard deviation (i.e., the square root of the variance).

The log likelihood function and dlog likelihood don't care about your
rows of data; if you have an 8 x 7 data set, it will give you the log
likelihood of those 56 observations given the mean and variance you provide.

\f$N(\mu,\sigma^2) = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2)\f$

\f$\ln N(\mu,\sigma^2) = (-(x-\mu)^2 / 2\sigma^2) - \ln (2 \pi \sigma^2)/2 \f$

\f$d\ln N(\mu,\sigma^2)/d\mu = (x-\mu) / \sigma^2 \f$

\f$d\ln N(\mu,\sigma^2)/d\sigma^2 = ((x-\mu)^2 / 2(\sigma^2)^2) - 1/2\sigma^2 \f$

The \c apop_ls settings group includes a \c want_cov element, which refers to the covariance matrix for the mean and variance. You can set that to \c 'n' if you are estimating millions of these and need to save time (i.e. <tt>Apop_model_add_group(your_model, apop_ls, .want_cov = 'n');</tt>).
\hideinitializer
\ingroup models
*/
apop_model apop_normal = {"Normal distribution", 2, 0, 0,
 .estimate = normal_estimate, .log_likelihood = normal_log_likelihood, .score = normal_dlog_likelihood, 
 .constraint = beta_1_greater_than_x_constraint, .draw = normal_rng};


//The Lognormal distribution


static double lognormal_log_likelihood(apop_data *d, apop_model *params);

static apop_model * lognormal_estimate(apop_data * data, apop_model *parameters){
  apop_model 	*est = apop_model_copy(*parameters);
  double   mean    = 0,
           var     = 0; 
  apop_ls_settings *p = apop_settings_get_group(est, "apop_ls");
    if (!p) {
        Apop_model_add_group(est, apop_ls);
        p = apop_settings_get_group(est, "apop_ls");
    }
    apop_matrix_mean_and_var(data->matrix, &mean, &var);
    if (!est->parameters)
        est->parameters = apop_data_alloc(2, 0, 0);
    double sigsq   = log(1+ var/gsl_pow_2(mean));
	gsl_vector_set(est->parameters->vector, 0, log(mean)- sigsq/2);
	gsl_vector_set(est->parameters->vector, 1, sqrt(sigsq));
    est->llikelihood	= lognormal_log_likelihood(data, est);
    est->data           = data;
	return est;
}

static double lnx_minus_mu_squared(double x, void *mu_in){
	return gsl_pow_2(log(x) - *(double *)mu_in);
}

/* The log likelihood function for the lognormal distribution.

\$f = exp(-(ln(x)-\mu)^2/(2\sigma^2))/ (x\sigma\sqrt{2\pi})\$
\$ln f = -(ln(x)-\mu)^2/(2\sigma^2) - ln(x) - ln(\sigma\sqrt{2\pi})\$

\param beta	beta[0]=the mean (after logging); beta[1]=the variance (after logging)
\param d	the set of data points; see notes.
*/
static double lognormal_log_likelihood(apop_data *d, apop_model *params){
    Get_vmsizes(d)
    Nullcheck(params)
    Nullcheck_p(params)
    double mu	   = gsl_vector_get(params->parameters->vector,0);
    double sd      = gsl_vector_get(params->parameters->vector,1);
    long double ll = -apop_map_sum(d, .fn_dp=lnx_minus_mu_squared, .param=&mu);
      ll /= (2*gsl_pow_2(sd));
      ll -= apop_map_sum(d, log, .part='m');
      ll -= tsize*(M_LNPI+M_LN2+log(sd));
	return ll;
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

/** An Apophenia wrapper for the GSL's Normal RNG, logged.

\param r	a gsl_rng already allocated
\param *out	To where I will write the drawn number
\param *p   A pointer to the model.
 */
static void lognormal_rng(double *out, gsl_rng *r, apop_model *p){
	*out = exp(gsl_ran_gaussian(r, p->parameters->vector->data[1]) + p->parameters->vector->data[0]);
}

/** The lognormal distribution. 

\hideinitializer
\ingroup models
*/
apop_model apop_lognormal = {"Lognormal distribution", 2, 0, 0,
 .estimate = lognormal_estimate, .log_likelihood = lognormal_log_likelihood, /*.score = lognormal_dlog_likelihood,*/ 
 .constraint = beta_1_greater_than_x_constraint, .draw = lognormal_rng};

