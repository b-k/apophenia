/** \file apop_bernoulli.c 
 
  The Bernoulli distribution as an \ref apop_model.*/
/*Copyright (c) 2007--2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "likelihoods.h"

static double bernoulli_log_likelihood(apop_data*, apop_model*);

static apop_model * bernoulli_estimate(apop_data * data,  apop_model *parameters){
  apop_model 	*est= parameters ? parameters : apop_model_copy(apop_bernoulli);
  apop_model_clear(data, est);
  double		p, sum  = 0;
  int           i, j,
                n       = (data->matrix->size1*data->matrix->size2);
    for (i=0; i< data->matrix->size1; i++)
        for (j=0; j< data->matrix->size2; j++)
            if (gsl_matrix_get(data->matrix, i, j))
                sum++;
    p   = (sum+0.0)/n;
	gsl_vector_set(est->parameters->vector, 0, p);
    est->llikelihood	= bernoulli_log_likelihood(data, est);
    est->covariance = apop_data_alloc(0,1,1);
    apop_data_set(est->covariance, 0,0, p*(1-p));
	return est;
}

static double bernie_ll(double x, void * pin){ double *p = pin; return x ? log(*p) : log(1-*p); }

static double bernoulli_log_likelihood(apop_data *d, apop_model *params){
    apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    double p   = apop_data_get(params->parameters,0,-1);
	return apop_map_sum(d, .fn_dp = bernie_ll, .param=&p);
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

/** The Bernoulli model.

  Data format: the matrix can have any size, and I just count up zeros
  and non-zeros. The bernoulli paramter $p$ is the percentage of non-zero
  values in the matrix. Its variance is $p(1-p)$.

\ingroup models
*/
apop_model apop_bernoulli = {"Bernoulli distribution", 1,0,0,
	.estimate = bernoulli_estimate, .log_likelihood = bernoulli_log_likelihood, 
   .constraint =  bernoulli_constraint, .draw = bernoulli_rng};
