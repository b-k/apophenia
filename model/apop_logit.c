/** \file apop_logit.c

Copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
*/
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


/* 
This is an MLE, so this is a one-liner.
*/
static apop_params * logit_estimate(apop_data * data,  apop_params *parameters){
	return apop_maximum_likelihood(data,  apop_logit, parameters);
}


/*For the sake of the fdf function, we keep xdotbeta global.
  */

gsl_vector  *xdotbeta;
int         calculate_xdotbeta  = 1;
int         keep_xdotbeta       = 0;

/*



*/
static double logit_log_likelihood(const apop_data *beta, apop_data *d, apop_params *p){
size_t	    i;
double	    loglike 	= 0,
            xb;
gsl_matrix 	*data 		= d->matrix;
gsl_matrix  independent = gsl_matrix_submatrix(data, 0,1, 
                                    data->size1, data->size2 -1).matrix;
    if (calculate_xdotbeta){
        xdotbeta    = gsl_vector_calloc(data->size1);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &independent, beta->vector, 0, xdotbeta);
        }
	for(i=0;i< data->size1; i++){
        xb    = gsl_vector_get(xdotbeta, i);
		if (gsl_matrix_get(data, i, 0))
            loglike   += xb/(1.+xb);
        else
            loglike   += 1./(1.+xb);
	}
    if (!keep_xdotbeta)
        gsl_vector_free(xdotbeta);
	return loglike;
}

static double logit_p(const apop_data *beta, apop_data *d, apop_params *p){
    return exp(logit_log_likelihood(beta, d, p));
}

/* The derivative of the logit distribution, for use in likelihood
  minimization. 
  The format is often the same as above: go line by line through a gsl_matrix.
  The sample is a three-dimensional parameter vector.
 */
static void logit_dlog_likelihood(const apop_data *beta, apop_data *d, gsl_vector *gradient, apop_params *p){
int		    i,j;
double	    dtotal[3];
gsl_matrix 	*data 	= d->matrix;
    dtotal[0]  = 0,
    dtotal[1]  = 0,
    dtotal[2]  = 0;
	for(i=0; i< data->size1; i++){
		dtotal[0]  += 0; //PLACE MATH HERE
		dtotal[1]  += 0; //PLACE MATH HERE
		dtotal[2]  += 0; //PLACE MATH HERE
	}
	for(j=0; j< beta->vector->size; j++){
	    gsl_vector_set(gradient,j,dtotal[j]);
    }
}


/** 
  Simple, but some trickery to keep xdotbeta. Notice that the two switches
  leave the function with the same values with which they came in.

static void logit_fdf( const gsl_vector *beta, apop_data *d, double *f, gsl_vector *df, apop_params *p){
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
 the remaining columns are values of the independent variables. Thus,
 the model will return (data columns)-1 parameters.

\ingroup models
*/
apop_model apop_logit = {"logit",-1,0,0,  
    .estimate = logit_estimate, .p = logit_p, .log_likelihood = logit_log_likelihood, 
    .score = logit_dlog_likelihood};
