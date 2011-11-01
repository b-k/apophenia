/** \file 
 Set some of the parameters of a model to fixed values.*/

/* There's only one function here. Its header is in likelihoods.h
 
Copyright (c) 2007, 2009, 2011 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
void apop_data_predict_fill(apop_data *data, apop_data *predict);
apop_data *apop_predict_table_prep(apop_data *in, char fill_with_nans);

typedef struct {
    apop_model *base_model;
    apop_data   *predict;
    int         ct;
} apop_model_fixed_params_settings;

static void unpack(apop_data *out, apop_model *m){
    //real param set --> predict table 
   apop_model_fixed_params_settings *mset = Apop_settings_get_group(m, apop_model_fixed_params);
   Apop_col_t(mset->predict, "predict", p_in_tab);
   gsl_vector_memcpy(p_in_tab, m->parameters->vector);
   apop_data_predict_fill(out, mset->predict);
}

static void pack(apop_data *in, apop_model *m){
    //predict table --> real param set 
   apop_model_fixed_params_settings *mset = Apop_settings_get_group(m, apop_model_fixed_params);
   apop_data *predict = mset->predict;
    for(int i =0; i< predict->matrix->size1; i++){
        apop_data_set(predict, .row =i, .colname="predict", .val=apop_data_get(in, 
                                                        apop_data_get(predict, .row=i, .colname="row"),
                                                        apop_data_get(predict, .row=i, .colname="col")));
        if (i< mset->ct-1 && apop_data_get(predict, .row= i+1, .colname="page") 
                                != apop_data_get(predict, .row= i, .colname="page"))
            in = in->more;
    }
   Apop_col_t(mset->predict, "predict", p_in_tab);
   gsl_vector_memcpy(m->parameters->vector, p_in_tab);
}

//The macros generating the fixed_param_settings group's init/copy/free functions:
Apop_settings_init(apop_model_fixed_params, 
    Apop_assert(in.base_model, "I can't fix a NULL model's parameters.");
)
Apop_settings_copy(apop_model_fixed_params, )
Apop_settings_free(apop_model_fixed_params, )

static double fix_params_ll(apop_data *d, apop_model *fixed_model){
    apop_model *base_model = Apop_settings_get(fixed_model, apop_model_fixed_params, base_model);
    unpack(base_model->parameters, fixed_model);
    return apop_log_likelihood(d, base_model);
}

static double fix_params_p(apop_data *d, apop_model *fixed_model){
    apop_model *base_model = Apop_settings_get(fixed_model, apop_model_fixed_params, base_model);
    unpack(base_model->parameters, fixed_model);
    return apop_p(d, base_model);
}

static double  fix_params_constraint(apop_data *data, apop_model *fixed_model){
    apop_model *base_model = Apop_settings_get(fixed_model, apop_model_fixed_params, base_model);
    unpack(base_model->parameters, fixed_model);
    double out = base_model->constraint(data, base_model);
    if (out) 
        pack(base_model->parameters, fixed_model);
    return out;
}

static void fix_params_draw(double *out, gsl_rng* r, apop_model *eps){
    apop_model *base_model = Apop_settings_get(eps, apop_model_fixed_params, base_model);
    apop_data *tmp    = base_model->parameters;
    unpack(base_model->parameters, eps);
    base_model->draw(out, r, base_model);
    base_model->parameters   = tmp;
}

static apop_model *fixed_est(apop_data * data, apop_model *params){
    apop_model *base_model = Apop_settings_get(params, apop_model_fixed_params, base_model);
    if (!data) data = params->data;
    apop_model *e = apop_maximum_likelihood(data, params);
    unpack(base_model->parameters, e);
    apop_data_free(e->parameters);
    e->parameters   = apop_data_copy(base_model->parameters);
    return e;
}

static void fixed_param_show(apop_model *m){
    apop_model_fixed_params_settings *mset = Apop_settings_get_group(m, apop_model_fixed_params);
    printf("The fill-in table:\n");
    apop_data_show(mset->predict);
    printf("The base model, after unpacking:\n");
    unpack(mset->base_model->parameters, m);
    apop_model_print(mset->base_model);
}

static apop_model fixed_param_model = {"Fill me", .estimate=fixed_est, .p = fix_params_p, 
            .log_likelihood=fix_params_ll, .constraint= fix_params_constraint, 
            .draw=fix_params_draw, .print=fixed_param_show};

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

\li If you need to do anything with the base model

Here is a sample program. It produces a few thousand draws from a Multivariate Normal distribution,
and then tries to recover the means given a var/covar matrix fixed at the correct variance.

\include fix_params.c
  
  \param model_in   The base model
 \return a model that can be used like any other, with the given params fixed or free.
  */
apop_model * apop_model_fix_params(apop_model *model_in){
    Nullcheck_mp(model_in)
    apop_model *model_out  = apop_model_copy(fixed_param_model);
    apop_model *base = apop_model_copy(*model_in);
    Apop_model_add_group(model_out, apop_model_fixed_params, .base_model = base);

    apop_data *predict_tab; //Keep the predict tab on the data set and in the settings struct
    if (!(predict_tab = apop_data_get_page(model_in->parameters, "<predict>")))
        predict_tab = apop_predict_table_prep(model_in->parameters, 'y');
    if (!predict_tab || !predict_tab->matrix|| !predict_tab->matrix->size2){
        apop_data_free(predict_tab);
        apop_model_free(model_out);
        Apop_assert_c(0, base, 1, "No free parameters (which would be marked with"
                " an NaN). Returning a copy of the input model.");
    }
    apop_settings_set(model_out, apop_model_fixed_params, predict, predict_tab);

    if (Apop_settings_get_group(base, apop_mle))
        apop_settings_copy_group(model_out, base, "apop_mle");
    else Apop_model_add_group(model_out, apop_mle, .parent= model_out, .method=APOP_CG_PR,
                                     .want_cov='n', .step_size=1, .tolerance=0.2);

#define cut_if_missing(method) if (!model_in->method) model_out->method = NULL;
    cut_if_missing(p);
    cut_if_missing(draw);
    cut_if_missing(score);
    cut_if_missing(constraint);
    cut_if_missing(log_likelihood);
    model_out->vbase = predict_tab->matrix->size1;
    snprintf(model_out->name, 100, "%s, with some params fixed", model_in->name);
    return model_out;
}
