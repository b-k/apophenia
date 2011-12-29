/** \file 
 Set some of the parameters of a model to fixed values.*/

/* There's only one function here. Its header is in likelihoods.h
 
Copyright (c) 2007, 2009, 2011 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */


/*

   The transform model is a shell model for your actual models. Here are some examples:

   You want to find the most likely values for missing data. So: find the empty spots, put them in a table,
   and use that table as the parameter set to be searched by an MLE. 
        [this part's a mess]
       setting up the model: put the missing table in the parameters slot
       [now do ML step on the original model]

       There's only one case where we go in the other direction: if there is a binding constraint, the rule is that the parameters are changed in place; that has to be communicated back to the parent function. Thus:
       taking down the model: copy the parameters back to the parameters slot.

    You want to constrain the search to just a few of the parameters; other parameters are fixed at given values
    an ML search offers parameters, then:
        setting up: the parameters are copied to the not-fixed slots in the base model
        base_model.log_likelihood is called

        taking down, in case of a constraint: copy the appropriate slots from the base model to the parameter list


    One-parameter Normal distribution: variance is always 0.1*Mean.
        setting up: given one parameter (the mean), write the mean and variance.
        
        take down: you've gotta come up with a rule, I suppose.

    coordinate transform: space one is in imperial inches, space two is in metric centimeters
        setting up: divide by GSL_CONST_MKSA_INCH

        take down: multiply by GSL_CONST_MKSA_INCH

 */

#include "apop_internal.h"
void apop_data_predict_fill(apop_data *data, apop_data *predict);
apop_data *apop_predict_table_prep(apop_data *in, char fill_with_nans);

typedef void (*model_transform_fn)(apop_model *, apop_model *);

typedef struct {
    apop_model *base_model;
    model_transform_fn base_to_shell, shell_to_base;
} apop_model_transform_settings;

/*
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
*/

//The macros generating the fixed_param_settings group's init/copy/free functions:
Apop_settings_init(apop_model_transform, 
    Apop_assert(in.base_model, "I can't fix a NULL model's parameters.");
)
Apop_settings_copy(apop_model_transform, )
Apop_settings_free(apop_model_transform, )

static double transform_ll(apop_data *d, apop_model *shell_model){
    apop_model *base_model = Apop_settings_get(shell_model, apop_model_transform, base_model);
    Apop_settings_get(shell_model, apop_model_transform, shell_to_base)(shell_model, base_model);
    //unpack(base_model->parameters, fixed_model);
    return apop_log_likelihood(d, base_model);
}

static double transform_p(apop_data *d, apop_model *shell_model){
    apop_model *base_model = Apop_settings_get(shell_model, apop_model_transform, base_model);
    Apop_settings_get(shell_model, apop_model_transform, shell_to_base)(shell_model, base_model);
    //unpack(base_model->parameters, fixed_model);
    return apop_p(d, base_model);
}

static double  transform_constraint(apop_data *data, apop_model *shell_model){
    apop_model *base_model = Apop_settings_get(shell_model, apop_model_transform, base_model);
    Apop_settings_get(shell_model, apop_model_transform, shell_to_base)(shell_model, base_model);
    //unpack(base_model->parameters, fixed_model);
    double out = base_model->constraint(data, base_model);
    if (out) 
        Apop_settings_get(shell_model, apop_model_transform, base_to_shell)(base_model, shell_model);
        //pack(base_model->parameters, shell_model);
    return out;
}

static void transform_draw(double *out, gsl_rng* r, apop_model *shell_model){
    apop_model *base_model = Apop_settings_get(shell_model, apop_model_transform, base_model);
    Apop_settings_get(shell_model, apop_model_transform, shell_to_base)(shell_model, base_model);
    //apop_data *tmp    = base_model->parameters;           //Why do I save the prior state?
    //unpack(base_model->parameters, eps);
    base_model->draw(out, r, base_model);
                                                //Do I want to transform the data here?
    //base_model->parameters   = tmp;
}

static apop_model *transform_est(apop_data * data, apop_model *shell_model){
    apop_model *base_model = Apop_settings_get(shell_model, apop_model_transform, base_model);
    //if (!data) data = params->data;
    apop_model *e = apop_maximum_likelihood(data, shell_model);
    Apop_settings_get(shell_model, apop_model_transform, shell_to_base)(shell_model, base_model);
    //unpack(base_model->parameters, e);
    /*apop_data_free(e->parameters);                        //Again, why?
    e->parameters   = apop_data_copy(base_model->parameters);
    */
    return e;
}

static void transform_show(apop_model *m){
    apop_model_transform_settings *mset = Apop_settings_get_group(m, apop_model_transform);
    /*printf("The fill-in table:\n");
    apop_data_show(mset->predict);
    printf("The base model, after unpacking:\n");
    unpack(mset->base_model->parameters, m);*/
    apop_model_print(mset->base_model);
}

static apop_model transformed_model = {"Fill me", .estimate=transform_est, .p = transform_p, 
            .log_likelihood=transform_ll, .constraint= transform_constraint, 
            .draw=transform_draw, .print=transform_show};


apop_model * apop_model_transform(apop_model model_in, model_transform_fn shell_to_base, model_transform_fn base_to_shell){
    apop_model *model_out  = apop_model_copy(transformed_model);
    apop_model *base = apop_model_copy(model_in);
    Apop_model_add_group(model_out, apop_model_transform, .base_model = base, 
                                        .shell_to_base=shell_to_base, .base_to_shell=base_to_shell);


    #define cut_if_missing(method) if (!model_in.method) model_out->method = NULL;
    cut_if_missing(p);
    cut_if_missing(draw);
    cut_if_missing(score);
    cut_if_missing(constraint);
    cut_if_missing(log_likelihood);
    //model_out->vbase = predict_tab->matrix->size1;
    snprintf(model_out->name, 100, "%s, transformed", model_in.name);
    return model_out;
}

#include <gsl/gsl_const_mksa.h>
void to_cm(apop_model *i, apop_model *m){
    if (!i || !i->parameters) return;
    if (!m->parameters) m->parameters = apop_data_copy(i->parameters);
    else apop_data_memcpy(m->parameters, i->parameters);
    if (m->parameters->vector) gsl_vector_scale(m->parameters->vector, 1*GSL_CONST_MKSA_INCH);
    if (m->parameters->matrix) gsl_matrix_scale(m->parameters->matrix, 1*GSL_CONST_MKSA_INCH);
    //data
    if (i->data && !m->data) m->data = i->data;
}

void from_cm(apop_model *i, apop_model *m){
    if (!i || !i->parameters) return;
    if (!m->parameters) m->parameters = apop_data_copy(i->parameters);
    else apop_data_memcpy(m->parameters, i->parameters);
    if (m->parameters->vector) gsl_vector_scale(m->parameters->vector, 1*GSL_CONST_MKSA_INCH);
    if (m->parameters->matrix) gsl_matrix_scale(m->parameters->matrix, 1*GSL_CONST_MKSA_INCH);
}

//If the base model has no data but the parent model does, then set the base model's
//pointer.
void test_transform(){
    apop_model *mm = apop_model_transform(apop_normal, to_cm, from_cm);
    apop_data *test_data=apop_data_fill(apop_data_alloc(5,2),
                                36, 12,
                                32, 16,
                                41, 10,
                                34, 14,
                                30, 10);
    apop_estimate(test_data, *mm);
}


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
/*
apop_model * apop_model_fix_params(apop_model *model_in){
    apop_model *out = apop_model_transform(in, unpack, pack);

    apop_data *predict_tab; //Keep the predict tab on the data set and in the settings struct
    if (!(predict_tab = apop_data_get_page(model_in->parameters, "<predict>")))
        predict_tab = apop_predict_table_prep(model_in->parameters, 'y');
    if (!predict_tab || !predict_tab->matrix|| !predict_tab->matrix->size2){
        apop_data_free(predict_tab);
        apop_model_free(model_out);
        Apop_assert_c(0, base, 1, "No free parameters (which would be marked with"
                " an NaN). Returning a copy of the input model.");
    }
    apop_settings_set(model_out, apop_model_transform, predict, predict_tab);

    if (Apop_settings_get_group(base, apop_mle))
        apop_settings_copy_group(model_out, base, "apop_mle");
    else Apop_model_add_group(model_out, apop_mle, .parent= model_out, .method=APOP_CG_PR,
                                     .want_cov='n', .step_size=1, .tolerance=0.2);

}
*/
