/** \file apop_missing_data.c
 
  Some missing data handlers.

Copyright (c) 2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "mapply.h"
#include "stats.h"
#include "output.h"
#include "model.h"
#include "types.h"
#include "settings.h"
#include "regression.h"
#include "conversions.h"
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

    \todo  Doesn't handle text; doesn't delete row names as necessary.
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
    apop_assert(has_matrix || has_vector, NULL, 0, 'c', 
            "You sent to apop_data_listwise_delete a data set with void matrix and vector. Confused, it is returning NULL.\n");
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
    out->names  = apop_name_copy(d->names); 
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

void  apop_ml_imputation_settings_free(apop_ml_imputation_settings *in){
    free(in->row);  free(in->col); free(in);}

void * apop_ml_imputation_settings_copy(apop_ml_imputation_settings *in){
    apop_ml_imputation_settings *out = malloc(sizeof(apop_ml_imputation_settings));
    out->local_mvn = in->local_mvn;
    out->ct = in->ct;
    out->row = in->row;
    out->col = in->col;
    return out;
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
  int i;
    for (i=0; i< m->ct; i++){
        apop_data_set(x, m->row[i], m->col[i], gsl_vector_get(v->vector,i));
    }
}

//The model to send to the optimization

static double ll(apop_data *d, apop_model * ep){
  apop_ml_imputation_settings  *m = apop_settings_get_group(ep, "apop_ml_imputation");
    unpack(ep->parameters, d, m);
    return apop_log_likelihood(d, m->local_mvn);
}


static apop_model apop_ml_imputation_model= {"Impute missing data via maximum likelihood", 0,0,0, .log_likelihood= ll};


/*
static double no_nan_val(const double in){ return isnan(in)? 0 : in;}
static double no_nan_col(const double in){ return isnan(in);}

static double apop_mean_no_nans(gsl_vector *in){
    return apop_vector_map_sum(in, no_nan_val)/apop_vector_map_sum(in, no_nan_col);
}*/

/**
    Impute the most likely data points to replace NaNs in the data, and
    insert them into the given data. That is, the data set is modified
    in place.


\param  d       The data set. It comes in with NaNs and leaves entirely filled in.
\param  mvn A parametrized \c apop_model from which you expect the data was derived.
This is traditionally a Multivariate Normal.

\return A model ready for estimation. It includes a set of \c apop_mle settings, which you will almost certainly want to tune.

*/
apop_model * apop_ml_imputation(apop_data *d,  apop_model* mvn){
  apop_model *mc       = apop_model_copy(apop_ml_imputation_model);
    Apop_settings_add_group(mc, apop_ml_imputation, mvn);
    if (Apop_settings_get_group(mvn, apop_mle))
        apop_settings_copy_group(mc, mvn, "apop_mle");
    else {
        Apop_settings_add_group(mc, apop_mle, mc);
        Apop_settings_add(mc, apop_mle, method, APOP_CG_PR);
        Apop_settings_add(mc, apop_mle, step_size, 2);
        Apop_settings_add(mc, apop_mle, tolerance, 0.2);
    }
    int missing_ct = find_missing(d, mc, NULL);
    apop_assert(missing_ct, NULL, 1, 'c', "You sent apop_ml_imputation a data set with no NANs");
    mc->vbase          = Apop_settings_get(mc, apop_ml_imputation, ct);
    apop_model *out = apop_maximum_likelihood(d, *mc);
    //unpack(out->parameters, d, apop_settings_get_group(mc, "apop_ml_imputation"));//already unpacked.
    return out;
}
