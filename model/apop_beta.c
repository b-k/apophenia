/** \file apop_beta.c  The Beta distribution */
/*Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
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
    apop_assert(m<1&&m > 0, NULL, 0, 's', "You asked for a beta distribution "
                        "with mean %g, but the mean of the beta will always "
                        "be strictly between zero and one.", m);
    apop_assert(v <= 1./12, NULL, 0, 's', "You asked for a beta distribution "
                        "with variance %g, but the variance of the beta will always "
                        "be less than 1/12 (=0.08333).", v);
    double k     = (m * (1- m)/ v) -1 ;
    double alpha = m*k;
    double beta  = k * (1-m);
    return apop_model_set_parameters(apop_beta, alpha, beta);
}

static double beta_log_likelihood(apop_data *d, apop_model *p);

static apop_model * beta_estimate(apop_data * data,  apop_model *parameters){
  apop_model 	*est= parameters ? parameters : apop_model_copy(apop_beta);
  apop_model_clear(data, est);
  double		mean, var, alpha, beta;
	apop_matrix_mean_and_var(data->matrix, &mean, &var);	
    alpha   = gsl_pow_2(mean) * ((1-mean)/var - 1/mean);
    beta    = alpha * (1-mean)/mean;
	gsl_vector_set(est->parameters->vector, 0, alpha);
	gsl_vector_set(est->parameters->vector, 1, beta);
    est->llikelihood	= beta_log_likelihood(data, parameters);
    //apop_numerical_covariance_matrix(apop_beta, est, data);
	return est;
}

typedef struct{
    double alpha, beta; 
} ab_type;

static double betamap(double x, void *abin) {
    ab_type *ab = abin; 
    return (ab->alpha-1) * log(x) + (ab->beta-1) *log(1-x); 
}

static double beta_log_likelihood(apop_data *d, apop_model *p){
    apop_assert(p->parameters,  0, 0, 's', "You asked me to evaluate an un-parametrized model.");
    ab_type ab;
    size_t size = (d->vector ? d->vector->size : 0) + (d->matrix ? d->matrix->size1 + d->matrix->size2 : 0);
    ab.alpha       = apop_data_get(p->parameters,0,-1),
    ab.beta        = apop_data_get(p->parameters,1,-1);
	return apop_map_sum(d, .fn_dp = betamap, .param=&ab) + gsl_sf_lnbeta(ab.alpha, ab.beta) * size;
}

static double dbeta_callback(double x){ return log(1-x); }

static void beta_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
  size_t size = (d->vector ? d->vector->size : 0) + (d->matrix ? d->matrix->size1 + d->matrix->size2 : 0);
  apop_assert_void(m->parameters, 0,'s', "You asked me to evaluate an un-parametrized model.");
  double bb	= gsl_vector_get(m->parameters->vector, 0);
  double a	= gsl_vector_get(m->parameters->vector, 1);
  double lnsum = apop_map_sum(d, log);
  double ln_x_minus_1_sum = apop_map_sum(d, dbeta_callback);
	//Psi is the derivative of the log gamma function.
	gsl_vector_set(gradient, 0, lnsum  + (-gsl_sf_psi(a) + gsl_sf_psi(a+bb))*size);
	gsl_vector_set(gradient, 1, ln_x_minus_1_sum  + (-gsl_sf_psi(bb) + gsl_sf_psi(a+bb))*size);
}

static double beta_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1 and  0 < beta_2
    return apop_linear_constraint(v->parameters->vector, .margin= 1e-4);
}

static void beta_rng(double *out, gsl_rng *r, apop_model* eps){
    *out = gsl_ran_beta(r, apop_data_get(eps->parameters,0,-1), apop_data_get(eps->parameters,1,-1));
}

/** The beta distribution.
\hideinitializer
\ingroup models */
apop_model apop_beta = {"Beta distribution", 2,0,0, .estimate = beta_estimate, 
    .log_likelihood = beta_log_likelihood, .score = beta_dlog_likelihood, 
    .constraint = beta_constraint, .draw = beta_rng};
