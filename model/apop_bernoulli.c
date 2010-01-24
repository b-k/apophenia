/** \file apop_bernoulli.c 
 
  The Bernoulli distribution as an \ref apop_model.*/
/*Copyright (c) 2007--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "likelihoods.h"

static double bernie_ll(double x, void * pin){ 
    double *p = pin; 
    return x ? log(*p) : log(1-*p); 
}

static double bernoulli_log_likelihood(apop_data *d, apop_model *params){
    apop_assert(params->parameters, 0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    double p   = apop_data_get(params->parameters,0,-1);
	return apop_map_sum(d, .fn_dp = bernie_ll, .param=&p);
}

static double nonzero (double in) { return in !=0; }

static apop_model * bernoulli_estimate(apop_data * data,  apop_model *est){
  double		p       = 0;
  double        n       = (data->vector ? data->vector->size : 0)
                        + (data->matrix ? data->matrix->size1*data->matrix->size2 : 0);
    p   = apop_map_sum(data, nonzero)/n;
	gsl_vector_set(est->parameters->vector, 0, p);
    est->llikelihood	= bernoulli_log_likelihood(data, est);
    apop_data *cov = apop_data_alloc(0,1,1);
    apop_data_set(cov, 0,0, p*(1-p));
    apop_data_add_page(est->parameters, cov, "Covariance");
	return est;
}

static double bernoulli_constraint(apop_data *data, apop_model *inmodel){
    //constraint is 0 < b and  1 > b
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_calloc(2,2,1);
        apop_data_fill(constraint, 0., 1.,
                                -1., -1.);
    }
    return apop_linear_constraint(inmodel->parameters->vector, constraint, 1e-3);
}

static void bernoulli_rng(double *out, gsl_rng *r, apop_model* eps){
    *out = gsl_rng_uniform (r) < eps->parameters->vector->data[0]; 
}

apop_model apop_bernoulli = {"Bernoulli distribution", 1,0,0,
	.estimate = bernoulli_estimate, .log_likelihood = bernoulli_log_likelihood, 
   .constraint =  bernoulli_constraint, .draw = bernoulli_rng};
