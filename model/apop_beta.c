/** \file apop_beta.c  The Beta distribution */
/*Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_beta The Beta distribution.

The beta distribution has two parameters and is restricted between zero and one. You may also find \ref
apop_beta_from_mean_var to be useful. 

\adoc    Input_format  Any arrangement of scalar values. 
\adoc    Parameter_format   a vector, v[0]=\f$\alpha\f$; v[1]=\f$\beta\f$    
\adoc    RNG  Produces a scalar \f$\in[0,1]\f$. 
\adoc    settings None. 
*/

#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

/** The Beta distribution is useful for modeling because it is bounded between zero and one, and can be either unimodal (if the variance is low) or bimodal (if the variance is high), and can have either a slant toward the bottom or top of the range (depending on the mean).

The distribution has two parameters, typically named \f$\alpha\f$ and \f$\beta\f$, which can be difficult to interpret. However, there is a one-to-one mapping between (alpha, beta) pairs and (mean, variance) pairs. Since we have good intuition about the meaning of means and variances, this function takes in a mean and variance, calculates alpha and beta behind the scenes, and returns a random draw from the appropriate Beta distribution.

\param m
The mean the Beta distribution should have. Notice that m
is in [0,1].

\param v
The variance which the Beta distribution should have. It is in (0, 1/12), where (1/12) is the variance of a Uniform(0,1) distribution. Funny things happen with variance near 1/12 and mean far from 1/2.

\return
Returns an \c apop_beta model with its parameters appropriately set.

*/
apop_model *apop_beta_from_mean_var(double m, double v){
    apop_assert_s(m<1&&m > 0, "You asked for a beta distribution "
                        "with mean %g, but the mean of the beta will always "
                        "be strictly between zero and one.", m);
    double k     = (m * (1- m)/ v) -1 ;
    double alpha = m*k;
    double beta  = k * (1-m);
    return apop_model_set_parameters(apop_beta, alpha, beta);
}

static double beta_log_likelihood(apop_data *d, apop_model *p);

/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static apop_model * beta_estimate(apop_data * data,  apop_model *est){
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

static double beta_log_likelihood(apop_data *d, apop_model *p){
    Get_vmsizes(d) //tsize
    Nullcheck_m(p); Nullcheck_p(p);
    ab_type ab = { .alpha = apop_data_get(p->parameters,0,-1),
                   .beta  = apop_data_get(p->parameters,1,-1)
    };
	return apop_map_sum(d, .fn_dp = betamap, .param=&ab) + gsl_sf_lnbeta(ab.alpha, ab.beta) * tsize;
}

static double dbeta_callback(double x){ return log(1-x); }

static void beta_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    Get_vmsizes(d) //tsize
    Nullcheck_m(m); Nullcheck_p(m);
  double bb	= gsl_vector_get(m->parameters->vector, 0);
  double a	= gsl_vector_get(m->parameters->vector, 1);
  double lnsum = apop_map_sum(d, log);
  double ln_x_minus_1_sum = apop_map_sum(d, dbeta_callback);
	//Psi is the derivative of the log gamma function.
	gsl_vector_set(gradient, 0, lnsum  + (-gsl_sf_psi(a) + gsl_sf_psi(a+bb))*tsize);
	gsl_vector_set(gradient, 1, ln_x_minus_1_sum  + (-gsl_sf_psi(bb) + gsl_sf_psi(a+bb))*tsize);
}

static double beta_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1 and  0 < beta_2
    return apop_linear_constraint(v->parameters->vector, .margin= 1e-4);
}

static double beta_cdf(apop_data *d, apop_model *params){
  Nullcheck_m(params) Nullcheck_p(params) Nullcheck_d(d) 
  Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double alpha = gsl_vector_get(params->parameters->vector, 0);
    double beta = gsl_vector_get(params->parameters->vector, 1);
    return gsl_cdf_beta_P(val, alpha, beta);
}

static void beta_rng(double *out, gsl_rng *r, apop_model* eps){
    *out = gsl_ran_beta(r, apop_data_get(eps->parameters,0,-1), apop_data_get(eps->parameters,1,-1));
}

static void beta_print(apop_model *m){
    fprintf(apop_opts.output_pipe,
            "Beta distribution with alpha = %g and beta = %g.\n", apop_data_get(m->parameters,0,-1)
                                                               , apop_data_get(m->parameters,1,-1));
}

apop_model apop_beta = {"Beta distribution", 2,0,0, .dsize=1, .estimate = beta_estimate, 
    .log_likelihood = beta_log_likelihood, .score = beta_dlog_likelihood, 
    .constraint = beta_constraint, .draw = beta_rng, .cdf = beta_cdf, .print=beta_print};
