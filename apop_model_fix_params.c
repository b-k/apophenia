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

typedef struct {
    apop_model *base_model;
    size_t      *row, *col;
    int         ct;
} apop_model_fixed_params_settings;

static void addin(apop_model_fixed_params_settings *m, size_t i, size_t j){
    m->row  = realloc(m->row, ++(m->ct) * sizeof(size_t));
    m->col  = realloc(m->col, m->ct * sizeof(size_t));
    m->row[m->ct-1]    = i;
    m->col[m->ct-1]    = j;
}

static int find_missing(const apop_data *mask, apop_model *mc){
    //generate a list of fixed-parameter positions, and their paramvals.
  apop_model_fixed_params_settings  *mset = apop_settings_get_group(mc, "apop_model_fixed_params");
  mset->ct = 0;
    //find out where the NaNs are
    for (size_t i=0; mask->vector && i< mask->vector->size; i++)
            if (!apop_data_get(mask, i, -1))
                addin(mset, i, -1);
    for (size_t i=0; mask->matrix && i< mask->matrix->size1; i++)
        for (int j=0; j <mask->matrix->size2; j++)
            if (!apop_data_get(mask, i, j))
                addin(mset, i, j);
    apop_assert(mset->ct, 0, 0,'s',"You're asking me to estimate a model where every single parameter is fixed.");
    mc->vbase = mset->ct;
    return mset->ct;
}

static void unpack(const apop_data *v, apop_model *m){
apop_model_fixed_params_settings * mset = Apop_settings_get_group(m, apop_model_fixed_params);
    for (int i=0; i< mset->ct; i++)
        apop_data_set(mset->base_model->parameters, mset->row[i], mset->col[i], gsl_vector_get(v->vector,i));
}

static void  pack(apop_data *out,const  apop_data  *in, apop_model *m){
   apop_model_fixed_params_settings *mset = Apop_settings_get_group(m, apop_model_fixed_params);
    for(int i =0; i< mset->ct; i++)
        apop_data_set(out, i, -1, apop_data_get(in, mset->row[i], mset->col[i]));

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
  apop_model_fixed_params_settings *p    = apop_settings_get_group(fixed_model, "apop_model_fixed_params");
    unpack(fixed_model->parameters, fixed_model);
    return apop_log_likelihood(d, p->base_model);
}

static double i_p(apop_data *d, apop_model *fixed_model){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(fixed_model, "apop_model_fixed_params");
    unpack(fixed_model->parameters, fixed_model);
    return apop_p(d, p->base_model);
}

static double  i_constraint(apop_data *data, apop_model *fixed_model){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(fixed_model, "apop_model_fixed_params");
    unpack(fixed_model->parameters, fixed_model);
  double out = p->base_model->constraint(data, p->base_model);
    if (out) 
        pack(p->base_model->parameters, fixed_model->parameters, fixed_model);
    return out;
}

static void i_draw(double *out, gsl_rng* r, apop_model *eps){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(eps, "apop_model_fixed_params");
  apop_data             *tmp    = p->base_model->parameters;
    unpack(eps->parameters, eps);
    p->base_model->draw(out, r, p->base_model);
    p->base_model->parameters   = tmp;
}

static apop_model *fixed_est(apop_data * data, apop_model *params){
  apop_model_fixed_params_settings *p    = apop_settings_get_group(params, "apop_model_fixed_params");
    if (!data)
        data    = params->data;
    apop_model *e = apop_maximum_likelihood(data, params);
    unpack(e->parameters, params);
    apop_data_free(e->parameters);
    e->parameters   = apop_data_copy(p->base_model->parameters);
    return e;
}

static double find_nans(double in){ return !gsl_isnan(in); }

static apop_model fixed_param_model = {"Fill me", .estimate=fixed_est, .p = i_p, .log_likelihood=i_ll, 
                                    .constraint= i_constraint, .draw=i_draw};

/** Produce a model based on another model, but with some of the 
  parameters fixed at a given value. 
  
As well as a pointer to the model whose parameters are to be fixed, I need two sets of data.

For the fixed parameters, I need to know the values at which they'll be fixed. Send me an \ref apop_data set that has the same shape as the parameters of your model; at the positions of the fixed parameters, give the values to which they will be fixed. For the free parameters, I (mostly) don't care what value they have. This set of parameters can be either set as the <tt>model_in->parameters</tt> element, or as an argument to the model. [If you give me both, I will use the one explicitly sent in rather than the one attached to the input model.]

I also need to know which parameters to fix, which requires a mask that I can hold over the parameter set. Again, the mask is an \ref apop_data set of the same size and shape as your data. Where there is a nonzero marker, I will fix the parameter.

You again have two options for giving me this information. You can use the parameter matrix as the mask: just set parameters to be left free to a nonzero value (including \c GSL_NAN). Or, you can explicitly send in a mask, with ones at params to be fixed and zero elsewhere. Again, I will try the explicit mask first.

  
You also need to input the base model, which I will copy (along with its settings groups) to form the new model.

The output is an \c apop_model that can be estimated, Bayesian updated, et cetera.

\li The \c estimate method always uses an MLE, and it never calls the base model's \c estimate method.

\li If the input model has MLE-style settings attached, I'll use them for the \c estimate method. Otherwise, I'll set my own.

\li If the parameter input has non-NaN values at the free parameters, then I'll use those as the starting point for any search; else I'll start from <b>1</b> as usual.

Here is a sample program. It produces a few thousand draws from a multivariate normal distribution,
and then tries to recover the means given a var/covar matrix fixed at the correct variance.

\include fix_params.c
  
  \param model_in   The base model
 \return a model that can be used like any other, with the given params fixed or free.
  */
apop_model * apop_model_fix_params(apop_model *model_in){
   apop_assert(model_in, NULL, 0, 's', "You sent me a NULL model.");
    apop_data *paramvals = model_in->parameters; //just an alias
   apop_assert(model_in->parameters, NULL, 0, 's', "I need parameters passed in either via the model->parameters or as an argument to the function");
   apop_data * mask = apop_map(paramvals, find_nans);
    apop_model *model_out  = apop_model_copy(fixed_param_model);
    apop_model *base = apop_model_copy(*model_in);
    Apop_model_add_group(model_out, apop_model_fixed_params, .base_model = base);
    find_missing(mask, model_out);
    if (!Apop_settings_get_group(model_out, apop_mle))
        Apop_model_add_group(model_out, apop_mle, .parent= model_out, .method=APOP_CG_PR,
                                     .step_size=1, .tolerance=0.2);
    if (!model_in->p) model_out->p = NULL;
    if (!model_in->log_likelihood) model_out->log_likelihood = NULL;
    if (!model_in->score) model_out->score = NULL;
    if (!model_in->constraint) model_out->constraint = NULL;
    if (!model_in->draw) model_out->draw = NULL;
    snprintf(model_out->name, 100, "%s, with some params fixed", model_in->name);
    return model_out;
}
