/** \file apop_uniform.c  */
/* Copyright (c) 2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "internal.h"
#include "likelihoods.h"

apop_model apop_uniform;

static apop_model * uniform_estimate(apop_data * data,  apop_model *parameters){
  apop_model 	*est= parameters ? parameters : apop_model_copy(apop_uniform);
  apop_model_clear(data, est);
    est->parameters->vector->data[0]    = gsl_matrix_min(data->matrix);
    est->parameters->vector->data[1]    = gsl_matrix_max(data->matrix);
    return est;
}

static double unif_ll(apop_data *d, apop_model *m){
  Get_vmsizes(d) //tsize
  Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
    double min = GSL_MIN(msize1 ? gsl_matrix_min(d->matrix) : GSL_POSINF,
                          vsize ? gsl_vector_min(d->vector) : GSL_POSINF);
    double max = GSL_MAX(msize1 ? gsl_matrix_max(d->matrix) : GSL_NEGINF,
                          vsize ? gsl_vector_max(d->vector) : GSL_NEGINF);
    if (min> m->parameters->vector->data[0] && max < m->parameters->vector->data[1])
        return -log(m->parameters->vector->data[1] - m->parameters->vector->data[0]) * tsize;
    return GSL_NEGINF;
}

static double unif_p(apop_data *d, apop_model *m){
  Get_vmsizes(d) //tsize
  Nullcheck(d); Nullcheck_m(m); Nullcheck_p(m);
    double min = GSL_MIN(msize1 ? gsl_matrix_min(d->matrix) : GSL_POSINF,
                          vsize ? gsl_vector_min(d->vector) : GSL_POSINF);
    double max = GSL_MAX(msize1 ? gsl_matrix_max(d->matrix) : GSL_NEGINF,
                          vsize ? gsl_vector_max(d->vector) : GSL_NEGINF);
    if (min> m->parameters->vector->data[0] && max< m->parameters->vector->data[1])
        return pow(m->parameters->vector->data[1] - m->parameters->vector->data[0], -tsize);
    return 0;
}

static void uniform_rng(double *out, gsl_rng *r, apop_model* eps){
    *out =  gsl_rng_uniform(r) *(eps->parameters->vector->data[1]- eps->parameters->vector->data[0])+ eps->parameters->vector->data[0];
}

apop_model apop_uniform = {"Uniform distribution", 2, 0, 0,  
    .estimate = uniform_estimate,  .p = unif_p,.log_likelihood = unif_ll,  .draw = uniform_rng};



static apop_model * improper_uniform_estimate(apop_data * data,  apop_model *parameters){
    return parameters; }

static double improper_unif_ll(apop_data *d, apop_model *m){ return 0; }
static double improper_unif_p (apop_data *d, apop_model *m){ return 1; }

static void improper_uniform_rng(double *out, gsl_rng *r, apop_model* eps){
    apop_assert_void(0, 0, 's', "It doesn't make sense to make random draws from an improper Uniform.");
}

apop_model apop_improper_uniform = {"Improper uniform distribution", 2, 0, 0,  
    .estimate = improper_uniform_estimate,  .p = improper_unif_p,
    .log_likelihood = improper_unif_ll,  .draw = improper_uniform_rng};
