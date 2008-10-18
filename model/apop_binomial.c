/** \file apop_binomial.c 
 
  The binomial distribution as an \c apop_model.

Copyright (c) 2006--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

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

static double binomial_p(apop_data *d, apop_model *p);
static double binomial_log_likelihood(apop_data*, apop_model*);

static void get_hits_and_misses(const apop_data *data, char method, double *hitcount, double *misscount){
    if (method == 't'){
        size_t        i, j;
        *hitcount=0, *misscount=0;
        for(i=0; i < data->matrix->size1; i ++)
            for(j=0; j < data->matrix->size2; j ++)
                if (gsl_matrix_get(data->matrix, i, j))
                    (*hitcount)    ++;
                else
                    (*misscount)   ++;
    } else {
        APOP_COL(data, 0, misses);
        APOP_COL(data, 1, hits);
        *hitcount = apop_vector_sum(hits);
        *misscount = apop_vector_sum(misses);
    }
}

/*

\todo Look up  the covariance matrix for the parameters of the Binomial */
static apop_model * binomial_estimate(apop_data * data,  apop_model *parameters){
  apop_assert(data,  0, 0,'s', "You asked me to estimate the parameters of a model but sent NULL data.");
  apop_model 	*est= parameters ? parameters : apop_model_copy(apop_binomial);
  apop_model_clear(data, est);
  double hitcount, misscount;
  char method = apop_settings_get_group(parameters, "apop_rank") ? 'b' : 't';
    get_hits_and_misses(data, method, &hitcount, &misscount);   
    apop_data_add_named_elmt(est->parameters, "n", hitcount+misscount);
    apop_data_add_named_elmt(est->parameters, "p", hitcount/(hitcount + misscount));
    est->llikelihood	= binomial_log_likelihood(data, parameters);
    //apop_numerical_covariance_matrix(apop_binomial, est, data);
    return est;
}

static double binomial_p(apop_data *d, apop_model *params){
  apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
    return exp(binomial_log_likelihood(d, params));
}

static double binomial_log_likelihood(apop_data *d, apop_model *params){
  apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  double	  n       = apop_data_get(params->parameters, 0, -1),
              p       = apop_data_get(params->parameters, 1, -1);
  double hitcount, misscount, ll = 0;
  int    i;
  char method = apop_settings_get_group(params, "apop_rank") ? 'b' : 't';
    if (method == 't'){
        get_hits_and_misses(d, method, &hitcount, &misscount);
        return log(gsl_ran_binomial_pdf(hitcount, p, n));
    } else {
        for (i=0; i< d->matrix->size1; i++){
            hitcount = gsl_matrix_get(d->matrix, i, 1);
            ll += log(gsl_ran_binomial_pdf(hitcount, p, n));
        }
        return ll;
    }
}

static double binomial_constraint(apop_data *data, apop_model *b){
    //constraint is 0 < binomial_1 and  0 < binomial_2
  static apop_data *constraint = NULL;
    if (!constraint){
        constraint= apop_data_alloc(2,2,1);
        apop_data_fill(constraint, 0., 1., 0.,
                                   0., 0., 1.);
    }
    return apop_linear_constraint(b->parameters->vector, constraint, 1e-3);
}

static void binomial_rng(double *out, gsl_rng *r, apop_model* eps){
    *out =   gsl_ran_binomial(r, eps->parameters->vector->data[1],eps->parameters->vector->data[0]); }

/** The binomial model.

The parameters are kept in the vector element of the \c apop_model parameters element. \c parameters->vector->data[0]==n;
\c parameters->vector->data[1]==p.

Input data can take two forms:

The default is to take the data to have a binary form, meaning that
the system counts zeros as failures and non-zeros as successes. \f$N\f$
is the size of the matrix.

In rank-type format, the data is taken to be the two-column miss-hit
format: column zero of the matrix represents failures and column one
represents successes. Set this using, e.g.,
\code
apop_model *estimate_me = apop_model_copy(apop_binomial);
Apop_settings_alloc(estimate_me, apop_rank, NULL);
apop_estimate(your_data, estimate_me);
\endcode

\ingroup models
*/
apop_model apop_binomial = {"Binomial distribution", 2,0,0,
	.estimate = binomial_estimate, .p = binomial_p, .log_likelihood = binomial_log_likelihood, 
   .constraint = binomial_constraint, .draw = binomial_rng};
