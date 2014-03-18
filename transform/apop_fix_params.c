/** \file 
 Set some of the parameters of a model to fixed values.*/

/* There's only one public function here. Its header is in likelihoods.h
 
Copyright (c) 2007, 2009, 2011 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "apop_internal.h"
static apop_model *fixed_param_model;

//The model keeps a table of what the blanks should be filled in with.
//This first section does the work for that part.
static double find_nans(double in){ return isnan(in); }

static void addin(apop_data *predict, size_t i, int j, size_t page){
    int len;
    if (!predict->matrix){
        predict->matrix = gsl_matrix_alloc(1,5); 
        len = 0;
    } else 
        len = predict->matrix->size1;
    apop_matrix_realloc(predict->matrix, len + 1, predict->matrix->size2);
    apop_data_set(predict, .row=len, .colname="row", .val=i);
    apop_data_set(predict, .row=len, .colname="col", .val=j);
    apop_data_set(predict, .row=len, .colname="page", .val=page);
}

static int find_missing(const apop_data *data, apop_data *predict, size_t page, int ct){
    //generate a list of fixed-parameter positions, and their paramvals.
   apop_data * mask = apop_map((apop_data*)data, find_nans, .all_pages='y');
    //find out where the NaNs are
    for (size_t i=0; mask->vector && i< mask->vector->size; i++)
            if (apop_data_get(mask, i, -1))
                addin(predict, i, -1, page);
    for (size_t i=0; mask->matrix && i< mask->matrix->size1; i++)
        for (int j=0; j <mask->matrix->size2; j++)
            if (apop_data_get(mask, i, j))
                addin(predict, i, j, page);
    if (mask->more)
        ct += mask->vector ? apop_sum(mask->vector) : 0
              + mask->matrix ? apop_matrix_sum(mask->matrix) : 0
              + find_missing(mask->more, predict, page+1, ct);
    apop_data_free(mask);
    return ct;
}

static apop_data *apop_predict_table_prep(apop_data *in, char fill_with_nans){
    apop_data *out = apop_data_alloc( );
    apop_name_add(out->names, "<fillins>", 'h');
    apop_name_add(out->names, "row", 'c');
    apop_name_add(out->names, "col", 'c');
    apop_name_add(out->names, "page", 'c');
    apop_name_add(out->names, "value", 'c');
    if (fill_with_nans == 'y')
        find_missing(in, out, 0, 0);
    return out;
}

/* Take a \c predict table and set the entries in the data set to the given predicted
  value. Functions for prediction and imputation use this internally, and append to your
  data a \c predict table of the right form.  For example, \c apop_ml_impute uses
  this internally.
  
  I assume that the ordering of elements in the \c predict table include everything on the
  first page, then everything on the second, et cetera. 

\param data The data set to be filled in. 
\param predict The set of fillins.
*/
static void apop_data_predict_fill(apop_data *data, apop_data *predict){
    if (!predict) return;
    int this_page_ct = 0;
    apop_data *this_page = data;
    for (int i=0; i < predict->matrix->size1; i++){
        int p = apop_data_get(predict, .row=i, .colname="page");
        if (p != this_page_ct){
            this_page_ct = p; 
            this_page = this_page->more;
        }
        apop_data_set(this_page, .row= apop_data_get(predict, .row=i, .colname="row"),
                                 .col= apop_data_get(predict, .row=i, .colname="col"),
                                 .val= apop_data_get(predict, .row=i, .colname="value"));
    }
}

/////////End predict table machinery.

typedef struct {
    apop_model *base_model;
    apop_data *predict;
    int ct;
} apop_fix_params_settings;

static void unpack(apop_data *out, apop_model *m){
    //real param set --> predict table 
   apop_fix_params_settings *mset = Apop_settings_get_group(m, apop_fix_params);
   Apop_col_tv(mset->predict, "value", p_in_tab);
   gsl_vector_memcpy(p_in_tab, m->parameters->vector);
   apop_data_predict_fill(out, mset->predict);
}

static void pack(apop_data *in, apop_model *m){
    //predict table --> real param set 
   apop_fix_params_settings *mset = Apop_settings_get_group(m, apop_fix_params);
   apop_data *predict = mset->predict;
    for(int i =0; i< predict->matrix->size1; i++){
        apop_data_set(predict, .row =i, .colname="value", .val=apop_data_get(in, 
                                                        apop_data_get(predict, .row=i, .colname="row"),
                                                        apop_data_get(predict, .row=i, .colname="col")));
        if (i< mset->ct-1 && apop_data_get(predict, .row= i+1, .colname="page") 
                                != apop_data_get(predict, .row= i, .colname="page"))
            in = in->more;
    }
   Apop_col_tv(mset->predict, "value", p_in_tab);
   gsl_vector_memcpy(m->parameters->vector, p_in_tab);
}

//The macros generating the fixed_param_settings group's init/copy/free functions:
Apop_settings_init(apop_fix_params, 
    Apop_assert(in.base_model, "I can't fix a NULL model's parameters.");
)
Apop_settings_copy(apop_fix_params, )
Apop_settings_free(apop_fix_params, )

static long double fix_params_ll(apop_data *d, apop_model *fixed_model){
    apop_model *base_model = Apop_settings_get(fixed_model, apop_fix_params, base_model);
    unpack(base_model->parameters, fixed_model);
    return apop_log_likelihood(d, base_model);
}

static long double fix_params_p(apop_data *d, apop_model *fixed_model){
    apop_model *base_model = Apop_settings_get(fixed_model, apop_fix_params, base_model);
    unpack(base_model->parameters, fixed_model);
    return apop_p(d, base_model);
}

static long double fix_params_constraint(apop_data *data, apop_model *fixed_model){
    apop_model *base_model = Apop_settings_get(fixed_model, apop_fix_params, base_model);
    unpack(base_model->parameters, fixed_model);
    long double out = base_model->constraint(data, base_model);
    if (out) pack(base_model->parameters, fixed_model);
    return out;
}

static int fix_params_draw(double *out, gsl_rng* r, apop_model *eps){
    apop_model *base_model = Apop_settings_get(eps, apop_fix_params, base_model);
    unpack(base_model->parameters, eps);
    return base_model->draw(out, r, base_model);
}

static void fixed_est(apop_data * data, apop_model *params){
    if (!data) data = params->data;
    apop_maximum_likelihood(data, params);
    apop_model *base_model = Apop_settings_get(params, apop_fix_params, base_model);
    unpack(base_model->parameters, params);
}

static void fixed_param_show(apop_model *m, FILE *out){
    apop_fix_params_settings *mset = Apop_settings_get_group(m, apop_fix_params);
    fprintf(out, "The fill-in table:\n");
    apop_data_print(mset->predict, .output_pipe=out);
    if (!m->parameters) printf("This copy of the model has not yet been estimated.\n");
    else {
        fprintf(out, "The base model, after unpacking:\n");
        unpack(mset->base_model->parameters, m);
    }
    apop_model_print(mset->base_model, out);
}

static void fixed_param_prep(apop_data *data, apop_model *params){
    apop_model_print_vtable_add(fixed_param_show, fixed_param_model);
    apop_model_clear(data, params);
    //apop_model *base_model = Apop_settings_get(params, apop_fix_params, base_model);
    //apop_prep(data, base_model);
}

static apop_model *fixed_param_model = &(apop_model){"Fill me", .estimate=fixed_est, .p = fix_params_p, 
            .log_likelihood=fix_params_ll, .constraint= fix_params_constraint, 
            .draw=fix_params_draw, .prep=fixed_param_prep};

/** Produce a model based on another model, but with some of the parameters fixed at a given value. 
  
You will send me the model whose parameters you want fixed, with the \c parameters element
set as follows. For the fixed parameters, simply give the values to which they will
be fixed. Set the free parameters to \c NaN.

For example, here is a Binomial distribution with a fixed \f$n=30\f$ but \f$p_1\f$ allowed to float freely:

\code
apop_model *bi30 = apop_model_fix_params(apop_model_set_parameters(apop_binomial, 30, GSL_NAN));
Apop_model_add_group(bi30, apop_mle, .starting_pt=(double[]){.5}); // The Binomial doesn't like the default 
                                                                   // starting point of 1.
apop_model *out = apop_estimate(your_data, bi30);
\endcode

The output is an \c apop_model that can be estimated, Bayesian updated, et cetera.

\li Rather than using this model, you may simply want a now-filled-in copy of the original model. Use \ref  apop_model_fix_params_get_base to retrieve the original model's parameters.

\li The \c estimate method always uses an MLE, and it never calls the base model's \c estimate method.

\li If the input model has MLE-style settings attached, I'll use them for the \c estimate method. Otherwise, I'll set my own.

\li If the parameter input has non-NaN values at the free parameters, then I'll use those as the starting point for any MLE search; the defaults for the variables without fixed values starts from <b>1</b> as usual.

\li I do check the \c more pointer of the \c parameters for additional pages and <tt>NaN</tt>s on those pages.

Here is a sample program. It produces a few thousand draws from a Multivariate Normal distribution,
and then tries to recover the means given a var/covar matrix fixed at the correct variance.

\include fix_params.c
  
\param model_in   The base model
\return a model that can be used like any other, with the given params fixed or free.
*/
apop_model * apop_model_fix_params(apop_model *model_in){
    Nullcheck_mp(model_in, NULL)
    apop_model *model_out  = Apop_model_copy_set(fixed_param_model,
                                apop_fix_params, .base_model = model_in);

    apop_data *predict_tab; //Keep the predict tab on the data set and in the settings struct
    predict_tab = apop_predict_table_prep(model_in->parameters, 'y');
    Apop_stopif (!predict_tab || !predict_tab->matrix|| !predict_tab->matrix->size2,
        apop_data_free(predict_tab);
        apop_model_free(model_out);
        return apop_model_copy(model_in);
        , 1, "No free parameters (which would be marked with a NaN). "
                "Returning a copy of the input model."
    );
    apop_settings_set(model_out, apop_fix_params, predict, predict_tab);

    if (Apop_settings_get_group(model_in, apop_mle))
        apop_settings_copy_group(model_out, model_in, "apop_mle");
    else Apop_model_add_group(model_out, apop_mle, .method="PR cg",
                                     .step_size=1, .tolerance=0.2);

    #define cut_if_missing(method) if (!model_in->method) model_out->method = NULL;
    cut_if_missing(p);
    cut_if_missing(draw);
    cut_if_missing(constraint);
    cut_if_missing(log_likelihood);
    model_out->vsize = predict_tab->matrix->size1;
    model_out->dsize = model_in->dsize;
    snprintf(model_out->name, 100, "%s, with some params fixed", model_in->name);
    return model_out;
}

/** The \ref apop_model_fix_params function produces a model that has only the non-fixed
  parameters of the model. After estimation of the fixed-parameter model, this function
  fills the \c parameters element of the base model and returns a pointer to the
  base model.
*/
apop_model * apop_model_fix_params_get_base(apop_model *fixed_model){
    apop_model *base_model = Apop_settings_get(fixed_model, apop_fix_params, base_model);
    unpack(base_model->parameters, fixed_model);
    return base_model;
}
