/** \file apop_WLS.c

  Weighted Least Squares. Much of the real work is done in regression.c.*/
/* Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "asst.h"
#include "model.h"
#include "stats.h"
#include "settings.h"
#include "regression.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>

static double wls_log_likelihood (apop_data *, apop_model *);

/* The procedure here is to simply modify the input data, run OLS on
the modified data, and then claim that the output was from WLS.  */
static apop_model * wls_estimate(apop_data *inset, apop_model *epin){
  apop_model    *epcopy = malloc(sizeof(*epcopy)),
                            *ep     = epin;
  apop_data                 *set;
  int                       i;
  gsl_vector                v;
  apop_ls_settings *insettings = apop_settings_get_group(epin, "apop_ls");
    if (!inset->weights){
        printf("You need to specify weights in the data set to use apop_wls.\n");
        return NULL;
    }
    set = (insettings && insettings->destroy_data)
        ? inset
        : apop_data_copy(inset);
    for (i=0; i< inset->matrix->size2; i++){
        v   = gsl_matrix_column(set->matrix, i).vector;
        gsl_vector_mul(&v, inset->weights);
    }
    memcpy(epcopy, ep, sizeof(*epcopy));
    apop_model *od = apop_model_copy(apop_ols);
    Apop_settings_alloc_add(od, apop_ls, destroy_data, 1, inset);
    apop_model *out     = apop_estimate(set, *od);
    out->estimate       = wls_estimate;
    out->log_likelihood = wls_log_likelihood;
    apop_model_free(od);
    return out;
}

/** The assumption that makes a log likelihood possible is that the
errors are normally distributed.

This function is a bit inefficient, in that it calculates the error terms,
which you may have already done in the OLS estimation.

 */
static double wls_log_likelihood (apop_data *d, apop_model *params){ 
  apop_assert(params->parameters,  0, 0,'s', "You asked me to evaluate an un-parametrized model.");
  int           i; 
  long double   total_prob  = 0; 
  double        sigma, expected, actual, weight;
  gsl_matrix	*data		= d->matrix;
  gsl_vector    v;
  gsl_vector    *errors     = gsl_vector_alloc(data->size1);
    if (!d->weights){
        printf("You need to specify weights to use apop_wls.\n");
        return 0;
    }
	for(i=0;i< data->size1; i++){
        v            = gsl_matrix_row(data, i).vector;
        gsl_blas_ddot(params->parameters->vector, &v, &expected);
        actual       = gsl_matrix_get(data,i, 0);
        expected    += gsl_vector_get(params->parameters->vector,0) * (1 - actual); //data isn't affine.
        weight       = gsl_vector_get(d->weights, i); //This is the only change from ols_ll.
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
	for(i=0;i< data->size1; i++){
        total_prob  += log(gsl_ran_gaussian_pdf(gsl_vector_get(errors, i), sigma));
	} 
    gsl_vector_free(errors);
    return total_prob;
}

/** The WLS model

  You will need to provide the weights in yourdata->weights. Otherwise, this model is just like \ref apop_ols.
\ingroup models
*/
apop_model apop_wls = {"Weighted Least Squares", -1,0,0, .estimate = wls_estimate, .log_likelihood = wls_log_likelihood};
