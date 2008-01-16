/** \file apop_missing_data.c
 
  Some missing data handlers.

Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "likelihoods.h"

/** If there is an NaN anywhere in the row of data (including the matrix
    and the vector) then delete the row from the data set.

    The function returns a new data set with the NaNs removed, so
    the original data set is left unmolested. You may want to \c
    apop_data_free the original immediately after this function.

    If every row has an NaN, then this returns NULL; you may want to
    check for this after the function returns.

    \param d    The data, with NaNs
    \return     A (potentially shorter) copy of the data set, without NaNs.

    \todo  Doesn't handle names or text.
*/
apop_data * apop_data_listwise_delete(apop_data *d){
  int i, j, min = 0, max = 0, height=0, has_vector=0, has_matrix=0, to_rm;
    //get to know the input.
    if (d->matrix){
        height      = d->matrix->size1;
        max         = d->matrix->size2;
        has_matrix  ++;
    } 
    if (d->vector){
        height      = height ?  height : d->vector->size;
        min         = -1;
        has_vector  ++;
    } 
    if (!has_matrix && !has_vector) {
        apop_error(0, 'c', "You sent to apop_data_listwise_delete a data set with void matrix and vector. Confused, it is returning NULL.\n");
        return NULL;
        }
    //find out where the NaNs are
  gsl_vector *marked = gsl_vector_calloc(height);
    for (i=0; i< d->matrix->size1; i++)
        for (j=min; j <max; j++)
            if (gsl_isnan(apop_data_get(d, i, j))){
                    gsl_vector_set(marked, i, 1);
                    break;
            }
    to_rm   = apop_sum(marked);
    //copy the good data.
    if (to_rm  == height)
        return NULL;
  apop_data *out = apop_data_alloc(0,height-to_rm, has_matrix ? max : -1);
    out->names  = apop_name_copy(d->names);                           //You loser!!! Fix this. And add text!!!!
    if (has_vector && has_matrix)
        out->vector = gsl_vector_alloc(height - to_rm);
    j   = 0;
    for (i=0; i< height; i++){
        if (!gsl_vector_get(marked, i)){
            if (has_vector)
                gsl_vector_set(out->vector, j, gsl_vector_get(d->vector, i));
            if (has_matrix){
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

static apop_model apop_ml_imputation_model;

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

void  apop_ml_imputation_settings_free(apop_ml_imputation_settings *in){free(in);}

void * apop_ml_imputation_settings_copy(apop_ml_imputation_settings *in){
    apop_ml_imputation_settings *out = malloc(sizeof(apop_ml_imputation_settings));
    out->local_mvn = in->local_mvn;
    out->ct = in->ct;
    out->row = in->row;
    out->col = in->col;
    return out;
}

static void addin(apop_ml_imputation_settings *m, size_t i, size_t j, double** starting_pt, apop_data *meanvar){
    m->row  = realloc(m->row, ++(m->ct) * sizeof(size_t));
    m->col  = realloc(m->col, m->ct * sizeof(size_t));
    *starting_pt = realloc(*starting_pt, m->ct*sizeof(double));
    m->row[m->ct-1]    = i;
    m->col[m->ct-1]    = j;
    (*starting_pt)[m->ct-1] = gsl_vector_get(meanvar->vector, j);
}

static void  find_missing(apop_data *d, apop_model *mc, apop_data *meanvar){
  apop_ml_imputation_settings  *mask = apop_settings_get_group(mc, "apop_ml_imputation");
  double ** starting_pt = &(Apop_settings_get(mc, apop_mle, starting_pt));
  int i, j, min = 0, max = 0;
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
            if (gsl_isnan(apop_data_get(d, i, j)))
                addin(mask, i, j, starting_pt, meanvar);
}

static void unpack(const apop_data *v, apop_data *x, apop_ml_imputation_settings * m){
  int i;
    for (i=0; i< m->ct; i++){
        apop_data_set(x, m->row[i], m->col[i], gsl_vector_get(v->vector,i));
    }
}

//The model to send to the optimization

static double ll(apop_data *d, apop_model * ep){
  apop_ml_imputation_settings  *m = apop_settings_get_group(ep, "apop_ml_imputation");
    unpack(ep->parameters, d, m);
    return apop_multivariate_normal.log_likelihood(d, m->local_mvn);
}


static apop_model apop_ml_imputation_model= {"Impute missing data via maximum likelihood", 0,0,0, .log_likelihood= ll};


/**
    Impute the most likely data points to replace NaNs in the data, and
    insert them into the given data. That is, the data set is modified
    in place.


\param  d       The data set. It comes in with NaNs and leaves entirely filled in.
\param  meanvar An \c apop_data set where the vector is the mean of each column and the matrix is the covariance matrix
\return A model ready for estimation. It includes a set of \ci apop_mle settings, which you will almost certainly want to tune.

*/
apop_model * apop_ml_imputation(apop_data *d,  apop_data* meanvar){
  apop_model *mc       = apop_model_copy(apop_ml_imputation_model);
  Apop_settings_add_group(mc, apop_mle, mc);
    apop_model *local_mvn             = apop_model_copy(apop_multivariate_normal);
    local_mvn->parameters = meanvar;
    Apop_settings_add_group(mc, apop_ml_imputation, local_mvn);
    Apop_settings_add(mc, apop_mle, method, APOP_SIMPLEX_NM);
    Apop_settings_add(mc, apop_mle, step_size, 2);
    Apop_settings_add(mc, apop_mle, tolerance, 0.2);
    find_missing(d, mc, meanvar);
    mc->vbase          = Apop_settings_get(mc, apop_ml_imputation, ct);
    return mc;
}
