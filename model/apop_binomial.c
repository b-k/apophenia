/** \file apop_binomial.c 
 
  The binomial distribution as an \c apop_model.

Copyright (c) 2007 by Ben Klemens. Licensed under the GNU GPL version 2.
*/
#include "model.h"

//The default list. You probably don't need them all.
#include "types.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>

static double binomial_log_likelihood(const apop_data*, apop_data*, apop_params*);

/*

\todo Look up  the covariance matrix for the parameters of the Binomial */
static apop_params * binomial_estimate(apop_data * data,  apop_params *parameters){
apop_params *est	    = apop_params_alloc(data, &apop_binomial, parameters, NULL);
double		mean, var;
	apop_matrix_mean_and_var(data->matrix, &mean, &var);	
	gsl_vector_set(est->parameters->vector, 0, gsl_pow_2(mean)/(mean-var));      //n
	gsl_vector_set(est->parameters->vector, 1, 1- var/mean); //p
	if (est->uses.log_likelihood)
		est->log_likelihood	= binomial_log_likelihood(est->parameters, data, parameters);
	if (est->uses.covariance)
		apop_numerical_covariance_matrix(apop_binomial, est, data);
	return est;
}

static double binomial_log_likelihood(const apop_data *in, apop_data *d, apop_params *params){
int		    i,j;
double	    x, 
            loglike = 0,
            n       = apop_data_get(in ? in: params->parameters,0,-1),
            p       = apop_data_get(in ? in: params->parameters,1,-1);
gsl_matrix 	*data 	= d->matrix;
	for(i=0;i< data->size1; i++)
        for(j=0; j< data->size2; j++){
            x       = gsl_matrix_get(data, i, j);
		    loglike += log(gsl_ran_binomial_pdf(x, p, n));
        }
	return loglike;
}

static double binomial_p(const apop_data *binomial, apop_data *d, apop_params *p){
    return exp(binomial_log_likelihood(binomial, d, p));
}

static double binomial_constraint(const apop_data *binomial, apop_data *returned_binomial, apop_params *v){
    //constraint is 0 < binomial_1 and  0 < binomial_2
  static apop_data *constraint = NULL;
    if (!constraint)constraint= apop_data_calloc(2,2,1);
    apop_data_set(constraint, 0, 0, 1);
    apop_data_set(constraint, 1, 0, 1);
    return apop_linear_constraint(binomial->vector, constraint, 1e-3, returned_binomial->vector);
}

static void binomial_rng(double *out, gsl_rng *r, apop_params* eps){
    *out =   gsl_ran_binomial(r, eps->parameters->vector->data[1],eps->parameters->vector->data[0]); }

/** The binomial model.

The parameters are kept in the vector element of the \c apop_params parameters element. \c parameters->vector->data[0]==n;
\c parameters->vector->data[1]==p.

\ingroup models
*/
apop_model apop_binomial = {"Binomial distribution", 2,0,0,
	binomial_estimate, binomial_p, binomial_log_likelihood, NULL, binomial_constraint, binomial_rng};
