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

static void ols_prep(apop_data *d, apop_model *m){
    if (!d->vector){
        APOP_COL(d, 0, independent);
        d->vector = apop_vector_copy(independent);
        gsl_vector_set_all(independent, 1);
        if (d->names->colct > 0) {		
            apop_name_add(d->names, d->names->column[0], 'v');
            sprintf(d->names->column[0], "1");
        }
    }
    void *mpt = m->prep; //also use the defaults.
    m->prep = NULL;
    apop_model_prep(d, m);
    m->prep = mpt;
}

/** The assumption that makes a log likelihood possible is that the
errors are normally distributed.

This function is a bit inefficient, in that it calculates the error terms,
which you may have already done in the OLS estimation.

 */
static double ols_log_likelihood (const apop_data *d, apop_model *p){ 
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

static double ols_p (const apop_data *d, apop_model *p){ 
    return exp(ols_log_likelihood(d, p)); }

/** The OLS model

  This is basically a wrapper for the OLS regression function, \ref apop_params_OLS.
\ingroup models
*/
apop_model apop_OLS = {.name="OLS", .vbase = -1, .estimate =apop_estimate_OLS, .p=ols_p,
                            .log_likelihood = ols_log_likelihood, .prep = ols_prep};

/** The GLS model

  This is basically a wrapper for the GLS regression function, \ref apop_params_GLS.
\ingroup models
*/
//apop_model apop_GLS = {"GLS", -1, apop_params_GLS, NULL, NULL, NULL, NULL, NULL};
