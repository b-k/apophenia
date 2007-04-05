/** \file apop_OLS.c

  OLS models. Much of the real work is done in regression.c.

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/

#include "model.h"
#include "regression.h"
#include "stats.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>


/** The procedure here is to simply modify the input data, run OLS on
the modified data, and then claim that the output was from WLS.
*/
apop_params * wls_estimate(apop_data *inset, apop_params *epin){
  apop_params    *epcopy = malloc(sizeof(*epcopy)),
                            *ep     = epin;
  apop_params             *out;
  apop_data                 *set;
  int                       i;
  gsl_vector                v;
  apop_OLS_params           *op = malloc(sizeof(*op));
    if (!inset->weights){
        printf("You need to specify weights in the data set to use apop_WLS.\n");
        return NULL;
    }
    if (!ep || strcmp(ep->method_name, "OLS") || ((apop_OLS_params*)ep->method_params)->destroy_data==0)
        set = apop_data_copy(inset);
    else
        set = inset;
    for (i=0; i< inset->matrix->size2; i++){
        v   = gsl_matrix_column(set->matrix, i).vector;
        gsl_vector_mul(&v, inset->weights);
    }
    memcpy(epcopy, ep, sizeof(*epcopy));
    op->destroy_data    = 1;
    epcopy->method_params = op;
    op->ep              = ep;
    out                 = apop_OLS.estimate(set, epcopy);
    free(out->model);
    out->model          = apop_model_copy(apop_WLS);
    return out;
}

/** The assumption that makes a log likelihood possible is that the
errors are normally distributed.

This function is a bit inefficient, in that it calculates the error terms,
which you may have already done in the OLS estimation.

 */
static double wls_log_likelihood (const apop_data *beta, apop_data *d, apop_params *params){ 
  int           i; 
  long double   total_prob  = 0; 
  double        sigma, expected, actual, weight;
  gsl_matrix	*data		= d->matrix;
  gsl_vector    v;
  gsl_vector    *errors     = gsl_vector_alloc(data->size1);
    if (!d->weights){
        printf("You need to specify weights to use apop_WLS.\n");
        return 0;
    }
	for(i=0;i< data->size1; i++){
        v            = gsl_matrix_row(data, i).vector;
        gsl_blas_ddot(beta->vector, &v, &expected);
        actual       = gsl_matrix_get(data,i, 0);
        expected    += gsl_vector_get(beta->vector,0) * (1 - actual); //data isn't affine.
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


static double wls_p (const apop_data *beta, apop_data *d, apop_params *p){ 
    return exp(wls_log_likelihood(beta, d, p));
            }

/** The WLS model

  You will need to provide the weights in data->weights. Otherwise, this model is just like \ref apop_OLS.
\ingroup models
*/
apop_model apop_WLS = {"WLS", -1,0,0, wls_estimate, wls_p, wls_log_likelihood};
