/** \file apop_OLS.c

  OLS models. Much of the real work is done in regression.c.

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "regression.h"
#include "stats.h"
#include <math.h>
#include <assert.h>
#include <apophenia/asst.h>
#include <gsl/gsl_blas.h>

/** The assumption that makes a log likelihood possible is that the
errors are normally distributed.

This function is a bit inefficient, in that it calculates the error terms,
which you may have already done in the OLS estimation.

 */
static double ols_log_likelihood (apop_data *d, apop_model *p){ 
  if (!p->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  int         i; 
  long double	total_prob  = 0; 
  double      sigma, expected, actual;
  gsl_matrix	*data		    = d->matrix;
  gsl_vector  v;
  gsl_vector  *errors         = gsl_vector_alloc(data->size1);
	for(i=0;i< data->size1; i++){
        v            = gsl_matrix_row(data, i).vector;
        gsl_blas_ddot(p->parameters->vector, &v, &expected);
        actual       = gsl_matrix_get(data,i, 0);
        expected    += gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
	for(i=0;i< data->size1; i++){
        total_prob  += log(gsl_ran_gaussian_pdf(gsl_vector_get(errors, i), sigma));
	} 
    gsl_vector_free(errors);
    return total_prob;
}

static double ols_p (apop_data *d, apop_model *p){ 
    return exp(ols_log_likelihood(d, p)); }

/** The IV model

  This is basically a wrapper for the IV regression function, \ref apop_estimate_IV
\ingroup models
*/
apop_model apop_IV = {.name="instrumental variables", .vbase = -1, .estimate =apop_estimate_IV, .p=ols_p,
                            .log_likelihood = ols_log_likelihood};
