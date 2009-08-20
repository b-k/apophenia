/** \file apop_missing_data.c Some missing data handlers. */
/* Copyright (c) 2007, 2009 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "internal.h"
#include "likelihoods.h"

/** If there is an NaN anywhere in the row of data (including the matrix and the vector) then delete the row from the data set.

The function returns a new data set with the NaNs removed, so the original data set is left unmolested. You may want to \c apop_data_free the original immediately after this function.

If every row has an NaN, then this returns \c NULL.

    \param d    The data, with NaNs
    \return     A (potentially shorter) copy of the data set, without NaNs.

    \todo  Doesn't handle text.
*/
apop_data * apop_data_listwise_delete(apop_data *d){
    Get_vmsizes(d) //defines vsize, msize1, msize2.
    int i, j, to_rm;
    int max = d->matrix ? d->matrix->size2: 0;
    int min = d->vector ? -1 : 0;
    apop_assert(msize1 || vsize, NULL, 0, 'c', 
            "You sent to apop_data_listwise_delete a data set with void matrix and vector. Confused, it is returning NULL.\n");
    //find out where the NaNs are
  gsl_vector *marked = gsl_vector_calloc(msize1);
    for (i=0; i< d->matrix->size1; i++)
        for (j=min; j <max; j++)
            if (gsl_isnan(apop_data_get(d, i, j))){
                    gsl_vector_set(marked, i, 1);
                    break;
            }
    to_rm   = apop_sum(marked);
    //copy the good data.
    if (to_rm  == msize1)
        return NULL;
  apop_data *out = apop_data_alloc(0, msize1-to_rm, msize1 ? max : -1);
    out->names  = apop_name_copy(d->names); 
    if (vsize && msize1)
        out->vector = gsl_vector_alloc(msize1 - to_rm);
    j   = 0;
    for (i=0; i< msize1; i++){
        if (!gsl_vector_get(marked, i)){
            if (vsize)
                gsl_vector_set(out->vector, j, gsl_vector_get(d->vector, i));
            if (msize1){
                APOP_ROW(d, i, v);
                gsl_matrix_set_row(out->matrix, j, v);
                if (d->names->row && d->names->rowct > i)
                    apop_name_add(out->names, d->names->row[i], 'r');
            }
            j++;
        }
    }
    gsl_vector_free(marked);
    return out;
}


//ML imputation

typedef struct {
    size_t      *row, *col;
    int         ct;
    apop_model  *local_mvn;
} apop_ml_imputation_settings;

apop_ml_imputation_settings * apop_ml_imputation_settings_alloc(apop_model *local){
    apop_ml_imputation_settings *out = malloc(sizeof(apop_ml_imputation_settings));
    out->local_mvn = local;
    return out;
}

void  apop_ml_imputation_settings_free(apop_ml_imputation_settings *in){
    free(in->row);  free(in->col); free(in);}

void * apop_ml_imputation_settings_copy(apop_ml_imputation_settings *in){
    apop_ml_imputation_settings *out = malloc(sizeof(apop_ml_imputation_settings));
    return (out = in);
}

static void addin(apop_ml_imputation_settings *m, size_t i, size_t j, double** starting_pt, gsl_vector *imean){
    m->row  = realloc(m->row, ++(m->ct) * sizeof(size_t));
    m->col  = realloc(m->col, m->ct * sizeof(size_t));
    *starting_pt = realloc(*starting_pt, m->ct*sizeof(double));
    m->row[m->ct-1]    = i;
    m->col[m->ct-1]    = j;
    (*starting_pt)[m->ct-1] = 1;
}

static int  find_missing(apop_data *d, apop_model *mc, gsl_vector *initialmean){
  apop_ml_imputation_settings  *mask = apop_settings_get_group(mc, "apop_ml_imputation");
  double ** starting_pt = &(Apop_settings_get(mc, apop_mle, starting_pt));
  int i, j, min = 0, max = 0, ct = 0;
    //get to know the input.
    if (d->matrix)
        max = d->matrix->size2;
    if (d->vector)
        min = -1;
    mask->row   = 
    mask->col   = NULL;
    mask->ct    = 0;
    //find out where the NaNs are
    for (i=0; i< d->matrix->size1; i++)
        for (j=min; j <max; j++)
            if (gsl_isnan(apop_data_get(d, i, j))){
                ct ++;
                addin(mask, i, j, starting_pt, initialmean);
            }
    return ct;
}

static void unpack(const apop_data *v, apop_data *x, apop_ml_imputation_settings * m){
    for (int i=0; i< m->ct; i++)
        apop_data_set(x, m->row[i], m->col[i], gsl_vector_get(v->vector,i));
}

//The model to send to the optimization

static double ll(apop_data *d, apop_model * ep){
  apop_ml_imputation_settings  *m = apop_settings_get_group(ep, "apop_ml_imputation");
    unpack(ep->parameters, d, m);
    return apop_log_likelihood(d, m->local_mvn);
}

static apop_model apop_ml_imputation_model= {"Impute missing data via maximum likelihood", .log_likelihood= ll};

/**
    Impute the most likely data points to replace NaNs in the data, and
    insert them into the given data. That is, the data set is modified
    in place.

\param  d       The data set. It comes in with NaNs and leaves entirely filled in.
\param  mvn A parametrized \c apop_model from which you expect the data was derived.
if \c NULL, then I'll use the Multivariate Normal that best fits the data after listwise deletion.

\return An estimated \ref apop_ml_imputation_model . Also, the data input will be filled in and ready to use.
*/
apop_model * apop_ml_imputation(apop_data *d,  apop_model* mvn){
    if (!mvn){
        apop_data *list_d = apop_data_listwise_delete(d);
        apop_assert(list_d, NULL, 0, 's', "Listwise deletion returned no whole rows, "
                            "so I couldn't fit a Multivariate Normal to your data. "
                            "Please provide a pre-estimated initial model.");
        mvn = apop_estimate(list_d, apop_multivariate_normal);
        apop_data_free(list_d);
    }
  apop_model *mc       = apop_model_copy(apop_ml_imputation_model);
    Apop_settings_add_group(mc, apop_ml_imputation, mvn);
    if (Apop_settings_get_group(mvn, apop_mle))
        apop_settings_copy_group(mc, mvn, "apop_mle");
    else 
        Apop_model_add_group(mc, apop_mle, .parent=mc, .method = APOP_CG_PR,
                                    .step_size=2, .tolerance=0.2);
    int missing_ct = find_missing(d, mc, NULL);
    apop_assert(missing_ct, NULL, 1, 'c', "You sent apop_ml_imputation a data set with no NANs");
    mc->vbase          = Apop_settings_get(mc, apop_ml_imputation, ct);
    return  apop_maximum_likelihood(d, *mc);//already unpacked.
}
