/** \file apop_model_fix_params.c 
 Set some of the parameters of a model to fixed values.*/

/* There's only one function here. Its header is in likelihoods.h
 
Copyright (c) 2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "asst.h"
#include "model.h"
#include "mapply.h"
#include "output.h"
#include "settings.h"
#include "likelihoods.h"

static double find_nans(double in){ return !gsl_isnan(in); }

typedef struct {
    apop_model *base_model;
    size_t      *row, *col, *page;
    int         ct;
} apop_model_fixed_params_settings;

static void addin(apop_model_fixed_params_settings *m, size_t i, size_t j, size_t page){
    m->row  = realloc(m->row, ++(m->ct) * sizeof(size_t));
    m->col  = realloc(m->col, m->ct * sizeof(size_t));
    m->page  = realloc(m->page, m->ct * sizeof(size_t));
    m->row[m->ct-1]    = i;
    m->col[m->ct-1]    = j;
    m->page[m->ct-1]   = page;
}

static int find_missing(const apop_data *mask, apop_model *mc, size_t page){
    //generate a list of fixed-parameter positions, and their paramvals.
  apop_model_fixed_params_settings  *mset = apop_settings_get_group(mc, apop_model_fixed_params);
    if (page==0)
        mset->ct = 0;
    //find out where the NaNs are
    for (size_t i=0; mask->vector && i< mask->vector->size; i++)
            if (!apop_data_get(mask, i, -1))
                addin(mset, i, -1, page);
    for (size_t i=0; mask->matrix && i< mask->matrix->size1; i++)
        for (int j=0; j <mask->matrix->size2; j++)
            if (!apop_data_get(mask, i, j))
                addin(mset, i, j, page);
    if (page == 0)
        apop_assert(mset->ct, 0, 0,'s',"You're asking me to estimate a model where every single parameter is fixed.");
    if (mask->more)
        find_missing(mask->more, mc, page+1);
    mc->vbase = mset->ct;
    return mset->ct;
}

static void unpack(const apop_data *v, apop_model *m){
    apop_model_fixed_params_settings * mset = Apop_settings_get_group(m, apop_model_fixed_params);
    apop_data *page =mset->base_model->parameters;
    for (int i=0; i< mset->ct; i++){
        apop_data_set(page, mset->row[i], mset->col[i], gsl_vector_get(v->vector,i));
        if (i< mset->ct-1 && mset->page[i+1] != mset->page[i])
            page = page->more;
    }
}

static void  pack(apop_data *out,const  apop_data  *in, apop_model *m){
   apop_model_fixed_params_settings *mset = Apop_settings_get_group(m, apop_model_fixed_params);
    for(int i =0; i< mset->ct; i++){
        apop_data_set(out, i, -1, apop_data_get(in, mset->row[i], mset->col[i]));
        if (i< mset->ct-1 && mset->page[i+1] != mset->page[i])
            out = out->more;
    }
}

static void *apop_model_fixed_params_settings_copy (apop_model_fixed_params_settings *in ){ 
    apop_model_fixed_params_settings *out = malloc(sizeof(apop_model_fixed_params_settings));
    *out = *in;
    return out;
} 

static void apop_model_fixed_params_settings_free (apop_model_fixed_params_settings *in ){ 
    free(in); }

static apop_model_fixed_params_settings *apop_model_fixed_params_settings_init (apop_model_fixed_params_settings in){
    apop_model_fixed_params_settings *out = malloc(sizeof(apop_model_fixed_params_settings));
    *out = (apop_model_fixed_params_settings){ };
    apop_varad_setting(in, out, base_model, NULL);
    return out;
}

static double i_ll(apop_data *d, apop_model *fixed_model){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(fixed_model, apop_model_fixed_params);
    unpack(fixed_model->parameters, fixed_model);
    return apop_log_likelihood(d, p->base_model);
}

static double i_p(apop_data *d, apop_model *fixed_model){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(fixed_model, apop_model_fixed_params);
    unpack(fixed_model->parameters, fixed_model);
    return apop_p(d, p->base_model);
}

static double  i_constraint(apop_data *data, apop_model *fixed_model){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(fixed_model, apop_model_fixed_params);
    unpack(fixed_model->parameters, fixed_model);
  double out = p->base_model->constraint(data, p->base_model);
    if (out) 
        pack(p->base_model->parameters, fixed_model->parameters, fixed_model);
    return out;
}

static void i_draw(double *out, gsl_rng* r, apop_model *eps){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(eps, apop_model_fixed_params);
  apop_data             *tmp    = p->base_model->parameters;
    unpack(eps->parameters, eps);
    p->base_model->draw(out, r, p->base_model);
    p->base_model->parameters   = tmp;
}

static apop_model *fixed_est(apop_data * data, apop_model *params){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(params, apop_model_fixed_params);
    if (!data)
        data    = params->data;
    apop_model *e = apop_maximum_likelihood(data, params);
    unpack(e->parameters, params);
    apop_data_free(e->parameters);
    e->parameters   = apop_data_copy(p->base_model->parameters);
    return e;
}

static apop_model fixed_param_model = {"Fill me", .estimate=fixed_est, .p = i_p, .log_likelihood=i_ll, 
                                    .constraint= i_constraint, .draw=i_draw};

/** Produce a model based on another model, but with some of the 
  parameters fixed at a given value. 
  
You will send me the model whose parameters you want fixed, with the \c parameters element
set as follows. For the fixed parameters, simply give the values to which they will
be fixed. Set the free parameters to \c NaN.

The output is an \c apop_model that can be estimated, Bayesian updated, et cetera.

\li The \c estimate method always uses an MLE, and it never calls the base model's \c estimate method.

\li If the input model has MLE-style settings attached, I'll use them for the \c estimate method. Otherwise, I'll set my own.

\li If the parameter input has non-NaN values at the free parameters, then I'll use those as the starting point for any search; else I'll start from <b>1</b> as usual.

\li I do check the \c more pointer of the \c parameters for additional pages and <tt>NaN</tt>s on those pages.

Here is a sample program. It produces a few thousand draws from a Multivariate Normal distribution,
and then tries to recover the means given a var/covar matrix fixed at the correct variance.

\include fix_params.c
  
  \param model_in   The base model
 \return a model that can be used like any other, with the given params fixed or free.
  */
apop_model * apop_model_fix_params(apop_model *model_in){
   apop_assert(model_in, NULL, 0, 's', "You sent me a NULL model.");
    apop_data *paramvals = model_in->parameters; //just an alias
   apop_assert(model_in->parameters, NULL, 0, 's', "I need parameters passed in either via the model->parameters or as an argument to the function");
   apop_data * mask = apop_map(paramvals, find_nans, .all_pages='y');
    apop_model *model_out  = apop_model_copy(fixed_param_model);
    apop_model *base = apop_model_copy(*model_in);
    Apop_model_add_group(model_out, apop_model_fixed_params, .base_model = base);
    find_missing(mask, model_out, 0);
    if (!Apop_settings_get_group(model_out, apop_mle))
        Apop_model_add_group(model_out, apop_mle, .parent= model_out, .method=APOP_CG_PR,
                                     .want_cov='n', .step_size=1, .tolerance=0.2);
    if (!model_in->p) model_out->p = NULL;
    if (!model_in->log_likelihood) model_out->log_likelihood = NULL;
    if (!model_in->score) model_out->score = NULL;
    if (!model_in->constraint) model_out->constraint = NULL;
    if (!model_in->draw) model_out->draw = NULL;
    snprintf(model_out->name, 100, "%s, with some params fixed", model_in->name);
    return model_out;
}
