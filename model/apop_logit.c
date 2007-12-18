/** \file apop_logit.c

Copyright (c) 2005--2007 by Ben Klemens.  Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  */
#include "model.h"



//The default list. You probably don't need them all.
#include "types.h"
#include "conversions.h"
#include "likelihoods.h"
#include "model.h"
#include "linear_algebra.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdio.h>
#include <assert.h>

static void modify_in_data(apop_data *d){
    if (!d->vector){
        APOP_COL(d, 0, independent);
        d->vector = apop_vector_copy(independent);
        gsl_vector_set_all(independent, 1);
    }
}

static apop_model * logit_estimate(apop_data * data,  apop_model *parameters){
    modify_in_data(data);
	return apop_maximum_likelihood(data,  *parameters);
}

/*For the sake of the fdf function, we keep xdotbeta global.
  */

gsl_vector  *xdotbeta;
int         calculate_xdotbeta  = 1;
int         keep_xdotbeta       = 0;

static double logit_log_likelihood(apop_data *d, apop_model *p){
  if (!p->parameters)
      apop_error(0,'s', "%s: You asked me to evaluate an un-parametrized model.", __func__);
  size_t	    i;
  double	    loglike 	= 0,
                xb, exb;
  gsl_matrix 	*data 		= d->matrix;
    if (calculate_xdotbeta){
        xdotbeta    = gsl_vector_calloc(data->size1);
        gsl_blas_dgemv(CblasNoTrans, 1.0, data, p->parameters->vector, 0, xdotbeta);
    }
	for (i=0;i< data->size1; i++) {
        xb    = gsl_vector_get(xdotbeta, i);
        exb   = exp(gsl_vector_get(xdotbeta, i));
		if (gsl_vector_get(d->vector, i))
            loglike   += xb - log(1 + exb);
        else
            loglike   +=  - log(1 + exb);
	}
    if (!keep_xdotbeta)
        gsl_vector_free(xdotbeta);
	return loglike;
}

static double logit_p(apop_data *d, apop_model *p){
    return exp(logit_log_likelihood(d, p));
}

/** 
  Simple, but some trickery to keep xdotbeta. Notice that the two switches
  leave the function with the same values with which they came in.

static void logit_fdf(gsl_vector *beta, apop_data *d, double *f, gsl_vector *df, apop_model *p){
    keep_xdotbeta       = 1;
	*f	= logit_log_likelihood(beta, d, NULL);
    calculate_xdotbeta  = 0;
    keep_xdotbeta       = 0;
	logit_dlog_likelihood(beta, d, df, NULL);
    calculate_xdotbeta  = 1;
}
	*/



/** The binary logit model.
 The first column of the data matrix this model expects is ones and zeros;
 the remaining columns are values of the independent variables. As with other regression-type models, the 
 system will move the first column to the vector and insert a column of ones as the first variable.

\ingroup models
*/
apop_model apop_logit = {"Logit",-1,0,0,  
    .estimate = logit_estimate, .p = logit_p, .log_likelihood = logit_log_likelihood};
    //.score = logit_dlog_likelihood};
