/** \file apop_OLS.c

  OLS models. Much of the real work is done in regression.c.

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */

#include "model.h"
#include "regression.h"
#include "stats.h"
#include "asst.h"

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
static double ols_log_likelihood (apop_data *d, apop_model *p){ 
  apop_assert(p->parameters, 0, 0,'s', "You asked me to evaluate an un-parametrized model. Returning zero.");
  int         i; 
  long double	ll  = 0; 
  double      sigma, expected, actual;
  gsl_matrix	*data		    = d->matrix;
  gsl_vector  *errors         = gsl_vector_alloc(data->size1);
	for(i=0;i< data->size1; i++){
        APOP_ROW(d, i, datarow);
        gsl_blas_ddot(p->parameters->vector, datarow, &expected);
        if (d->vector){ //then this has been prepped
            actual       = apop_data_get(d,i, -1);
        } else {
            actual       = gsl_matrix_get(data,i, 0);
            expected    += gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        }
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
	for(i=0;i< data->size1; i++){
        ll  += log(gsl_ran_gaussian_pdf(gsl_vector_get(errors, i), sigma));
	} 
    gsl_vector_free(errors);
    return ll;
}

/* $\partial {\cal N}(x\beta - y)/\partial \beta_i = \sum{x_i} \partial {\cal N}(K)/\partial K$ (at $K=x\beta -y$) */
static void ols_score(apop_data *d, gsl_vector *gradient, apop_model *p){ 
  apop_assert_void(p->parameters, 0,'s', "You asked me to evaluate an un-parametrized model. Not changing the gradient");
  size_t         i, j; 
  double      sigma, expected, actual;
  gsl_matrix	*data		    = d->matrix;
  gsl_vector  *errors         = gsl_vector_alloc(data->size1);
  gsl_vector  *normscore      = gsl_vector_alloc(2);
  apop_data  *subdata      = apop_data_alloc(0,1,1);
	for(i=0;i< data->size1; i++){
        APOP_ROW(d, i, datarow);
        gsl_blas_ddot(p->parameters->vector, datarow, &expected);
        if (d->vector){ //then this has been prepped
            actual       = apop_data_get(d,i, -1);
        } else {
            actual       = gsl_matrix_get(data,i, 0);
            expected    += gsl_vector_get(p->parameters->vector,0) * (1 - actual); //data isn't affine.
        }
        gsl_vector_set(errors, i, expected-actual);
    }
    sigma   = sqrt(apop_vector_var(errors));
    apop_model *norm = apop_model_set_parameters(apop_normal, 0.0, sigma);
    gsl_vector_set_all(gradient, 0);
	for(i=0;i< data->size1; i++){
        apop_data_set(subdata, 0, 0, gsl_vector_get(errors, i));
        apop_score(subdata, normscore, norm);
        for(j=0; j< data->size2; j++)
            apop_vector_increment(gradient, j, apop_data_get(d, i, j) * gsl_vector_get(normscore, 0));
	} 
    gsl_vector_free(errors);
    apop_model_free(norm);
}

/** The OLS model
  This one needs no introduction.

\ingroup models
*/
apop_model apop_OLS = {.name="OLS", .vbase = -1, .estimate =apop_estimate_OLS, 
                            .log_likelihood = ols_log_likelihood, .score=ols_score, .prep = ols_prep};

/** The GLS model

  This is basically a wrapper for the GLS regression function, \ref apop_params_GLS.
\ingroup models
*/
//apop_model apop_GLS = {"GLS", -1, apop_params_GLS, NULL, NULL, NULL, NULL, NULL};
