/* The Exponential distribution.
 Copyright (c) 2005--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

 \amodel apop_exponential The Exponential distribution.

\f$Z(\mu,k) 	= \sum_k 1/\mu e^{-k/\mu} 			\f$ <br>
\f$ln Z(\mu,k) 	= \sum_k -\ln(\mu) - k/\mu			\f$ <br>
\f$dln Z(\mu,k)/d\mu 	= \sum_k -1/\mu + k/(\mu^2)			\f$ <br>

Some write the function as:
\f$Z(C,k) dx = \ln C C^{-k}. \f$
If you prefer this form, just convert your parameter via \f$\mu = {1\over \ln C}\f$
(and convert back from the parameters this function gives you via \f$C=\exp(1/\mu)\f$).

\adoc    Input_format  
Ignores the matrix structure of the input data, so send in a 1 x N, an N x 1, or an N x M.
\li See also \ref apop_data_rank_compress for means of dealing with one more input data format.
                    
\adoc    Parameter_format   \f$\mu\f$ is in the zeroth element of the vector.   
\adoc    CDF  Produces a single number.
\adoc    settings   None.  */

#include "apop_internal.h"

static long double beta_greater_than_x_constraint(apop_data *data, apop_model *v){
    //constraint is 0 < beta_1
    return apop_linear_constraint(v->parameters->vector, .margin = 1e-3);
}

static long double exponential_log_likelihood(apop_data *d, apop_model *p){
    Nullcheck_mpd(d, p, GSL_NAN);
    Get_vmsizes(d) //tsize
    double mu = gsl_vector_get(p->parameters->vector, 0);
    double llikelihood = -((d->matrix ? apop_matrix_sum(d->matrix):0) + (d->vector ? apop_sum(d->vector) : 0))/ mu;
	llikelihood	-= tsize * log(mu);
	return llikelihood;
}

static void exponential_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *p){
    Nullcheck_mpd(d, p, );
    Get_vmsizes(d) //tsize
    double mu = gsl_vector_get(p->parameters->vector, 0);
    double d_likelihood = (d->matrix ? apop_matrix_sum(d->matrix):0) + (d->vector ? apop_sum(d->vector) : 0);
	d_likelihood /= gsl_pow_2(mu);
	d_likelihood -= tsize /mu;
	gsl_vector_set(gradient,0, d_likelihood);
}

/* \adoc estimated_info   Reports <tt>log likelihood</tt>. */
static void exponential_estimate(apop_data * data,  apop_model *est){
    apop_score_vtable_add(exponential_dlog_likelihood, apop_exponential);
    apop_name_add(est->parameters->names, "μ", 'r');
    Get_vmsizes(data); //msize1, msize2, vsize, tsize
    double mu =  (vsize ? vsize * apop_vector_mean(data->vector):0
                + msize1 ? msize1*msize2 * apop_matrix_mean(data->matrix):0)/tsize;
	gsl_vector_set(est->parameters->vector, 0, mu);
    apop_data_add_named_elmt(est->info, "log likelihood", exponential_log_likelihood(data, est));
}

static long double expo_cdf(apop_data *d, apop_model *params){
    Nullcheck_mpd(d, params, GSL_NAN);
    Get_vmsizes(d)  //vsize
    double val = apop_data_get(d, 0, vsize ? -1 : 0);
    double lambda = gsl_vector_get(params->parameters->vector, 0);
    return gsl_cdf_exponential_P(val, lambda);
}

/* \adoc RNG Just a wrapper for \c gsl_ran_exponential.  */
static int exponential_rng(double *out, gsl_rng* r, apop_model *p){
	*out = gsl_ran_exponential(r, p->parameters->vector->data[0]);
    return 0;
}

static void exponential_prep(apop_data *data, apop_model *params){
    apop_score_vtable_add(exponential_dlog_likelihood, apop_exponential);
    apop_model_clear(data, params);
}

apop_model *apop_exponential = &(apop_model){"Exponential distribution", 1,0,0,.dsize=1,
	 .estimate = exponential_estimate, .log_likelihood = exponential_log_likelihood, 
     .prep = exponential_prep, .constraint = beta_greater_than_x_constraint, 
     .draw = exponential_rng, .cdf = expo_cdf};
