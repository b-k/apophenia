/** \file 
        The Dirichlet distribution */
/*Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "mapply.h"
#include "internal.h"
#include "likelihoods.h"

static double dirichletlnmap(gsl_vector *v, void *pin) {
    gsl_vector *params = pin;
    if (v->stride ==1)
        return gsl_ran_dirichlet_lnpdf (params->size, params->data, v->data);
    //else:
    double t[v->size]; 
    for(size_t i=0; i < v->size; i ++)
        t[i] = v->data[i];
    return gsl_ran_dirichlet_lnpdf (params->size, params->data, t);
}

static double dirichlet_log_likelihood(apop_data *d, apop_model *p){
    Nullcheck(d); Nullcheck_m(p); Nullcheck_p(p);
	return apop_map_sum(d, .fn_vp = dirichletlnmap, .param=p->parameters->vector, .part='r');
}

static void dirichlet_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    Nullcheck_mv(m); Nullcheck_pv(m);
    double param_sum = apop_sum(m->parameters->vector);
    int n = d->matrix->size1;
    for(size_t i=0; i < m->parameters->vector->size; i ++){
        double thisparam = gsl_vector_get(m->parameters->vector, i);
        Apop_col(d, i, onecol);
        gsl_vector_set(gradient, i,  //Psi is the derivative of the log gamma function.
                apop_vector_map_sum(onecol, log) + n*gsl_sf_psi(param_sum) - n*gsl_sf_psi(thisparam));
    }
}

static double dirichlet_constraint(apop_data *data, apop_model *v){
    //all elements are > 0.
    return apop_linear_constraint(v->parameters->vector, .margin= 1e-4);
}

static void dirichlet_rng(double *out, gsl_rng *r, apop_model* eps){
    gsl_ran_dirichlet(r, eps->parameters->vector->size, eps->parameters->vector->data, out);
}

apop_model apop_dirichlet = {"Dirichlet distribution", -1,0,0, .dsize=-1,
    .log_likelihood = dirichlet_log_likelihood, .score = dirichlet_dlog_likelihood,
    .constraint = dirichlet_constraint, .draw = dirichlet_rng};
