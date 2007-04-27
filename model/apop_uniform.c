/** \file apop_uniform.c 

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

apop_model apop_uniform;

static apop_params * uniform_estimate(apop_data * data,  apop_params *parameters){
  apop_params 	*est	    = apop_params_alloc(data, &apop_uniform, parameters, NULL);
    est->parameters->vector->data[0]    = gsl_matrix_min(data->matrix);
    est->parameters->vector->data[1]    = gsl_matrix_max(data->matrix);
    return est;
}


static double unif_ll(const apop_data *params, apop_data *d, apop_params *v){
    if (gsl_matrix_min(d->matrix)> params->vector->data[0] && gsl_matrix_max(d->matrix)< params->vector->data[1])
        return log(params->vector->data[1] - params->vector->data[0]) *  d->matrix->size1 * d->matrix->size2;
    else
        return GSL_NEGINF;
}

static double unif_p(const apop_data *params, apop_data *d, apop_params *v){
    if (gsl_matrix_min(d->matrix)> params->vector->data[0] && gsl_matrix_max(d->matrix)< params->vector->data[1])
        return pow(params->vector->data[1] - params->vector->data[0],  d->matrix->size1 * d->matrix->size2);
    else
        return 0;
}


static void uniform_rng( double *out, gsl_rng *r, apop_params* eps){
    *out =  gsl_rng_uniform(r) *(eps->parameters->vector->data[1]- eps->parameters->vector->data[0])+ eps->parameters->vector->data[0];
}

/** The uniform model.
This is the two-parameter version of the uniform, expressing a uniform distribution over [a, b].

The MLE of this distribution is simply a = min(your data); b = max(your data).
Primarily useful for the RNG, such as when you have a Uniform prior model.

\ingroup models
*/
apop_model apop_uniform = {"Uniform distribution", 2, 0, 0,  
    .estimate = uniform_estimate,  .p = unif_p,.log_likelihood = unif_ll,  .draw = uniform_rng};
