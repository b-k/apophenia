/* The Dirichlet distribution 
Copyright (c) 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  

\amodel apop_dirichlet A multivariate generalization of the \ref apop_beta "Beta distribution".

\adoc    Input_format      Each row of your data matrix is a single observation.  
\adoc    Parameter_format   The estimated parameters are in the output model's <tt>parameters->vector</tt>. The size of the model is determined by the width of your input data set, so later RNG draws, \&c will match in size.
\adoc    settings   MLE-type: \ref apop_mle_settings, \ref apop_parts_wanted_settings   
*/

#include "apop_internal.h"

static double dirichletlnmap(gsl_vector *v, void *pin) {
    //we used gsl_matrix_row to get here==>guaranteed that v->stride==1.
    gsl_vector *params = pin;
    return gsl_ran_dirichlet_lnpdf (params->size, params->data, v->data);
}

static long double dirichlet_log_likelihood(apop_data *d, apop_model *p){
    Nullcheck_mpd(d, p, GSL_NAN);
    Apop_stopif(!p->parameters->vector, return GSL_NAN, 0, "parameters should be in inmodel->parameters->vector.");
    double paramsum = apop_sum(p->parameters->vector);
    Apop_stopif(fabs(paramsum)<1e-5, return GSL_NAN, 0, "Parameter total is too close to zero.");
    Apop_stopif(isnan(paramsum), return GSL_NAN, 0, "NaN parameter.");
	return apop_map_sum(d, .fn_vp = dirichletlnmap, .param=p->parameters->vector, .part='r');
}

static void dirichlet_dlog_likelihood(apop_data *d, gsl_vector *gradient, apop_model *m){
    Nullcheck_mpd(d, m, );
    double param_sum = apop_sum(m->parameters->vector);
    int n = d->matrix->size1;
    for(size_t i=0; i < m->parameters->vector->size; i ++){
        double thisparam = gsl_vector_get(m->parameters->vector, i);
        Apop_col_v(d, i, onecol);
        gsl_vector_set(gradient, i,  //Psi is the derivative of the log gamma function.
                apop_vector_map_sum(onecol, log) + n*gsl_sf_psi(param_sum) - n*gsl_sf_psi(thisparam));
    }
}

static long double dirichlet_constraint(apop_data *data, apop_model *v){
    //all elements are > 0.
    return apop_linear_constraint(v->parameters->vector, .margin= 1e-4);
}

/*\adoc    RNG  A call to \c gsl_ran_dirichlet.*/
static int dirichlet_rng(double *out, gsl_rng *r, apop_model* eps){
    gsl_ran_dirichlet(r, eps->parameters->vector->size, eps->parameters->vector->data, out);
    return 0;
}

static void dirichlet_prep(apop_data *data, apop_model *params){
    apop_score_vtable_add(dirichlet_dlog_likelihood, apop_dirichlet);
    apop_model_clear(data, params);
}

apop_model *apop_dirichlet = &(apop_model){"Dirichlet distribution", -1,0,0, .dsize=-1,
    .log_likelihood = dirichlet_log_likelihood, .prep = dirichlet_prep,
    .constraint = dirichlet_constraint, .draw = dirichlet_rng};
