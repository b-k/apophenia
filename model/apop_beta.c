/** apop_beta.c 
 
 (This file is not handled by Doxygen)

Base file copyright (c) 2005 by Ben Klemens. Licensed under the GNU GPL version 2.
Feel free to augment this with your own copyright: modifications (c) you, today.
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

static double beta_log_likelihood(const apop_data *beta, apop_data *d, void *p);

/* Every model should have an estimate function. If you are doing an
 MLE, then it should be the one line function here, and then you will need
 to fill in the beta_log_likelihood function below. If you are not
 doing an MLE, then you won't need the beta_log_likelihood function,
 but will probably instead do some substantial math here in this function.

 a
*/
static apop_estimate * beta_estimate(apop_data * data,  void *parameters){
apop_estimate 	*est	    = apop_estimate_alloc(data, apop_beta, parameters);
double		mean, var, alpha, beta;
	apop_matrix_mean_and_var(data->matrix, &mean, &var);	
    alpha   = gsl_pow_2(mean) * ((1-mean)/var - 1/mean);
    beta    = alpha * (1-mean)/mean;
	gsl_vector_set(est->parameters->vector, 0, alpha);
	gsl_vector_set(est->parameters->vector, 1, beta);
	if (est->ep.uses.log_likelihood)
		est->log_likelihood	= beta_log_likelihood(est->parameters, data, parameters);
	if (est->ep.uses.covariance)
		apop_numerical_covariance_matrix(apop_beta, est, data);
	return est;
    
	return NULL;
}


/*
Often, the input data is a gsl_matrix, and you need to go line by line
through the matrix, calculating the log likelihood. That's the sample
code here.  
*/
static double beta_log_likelihood(const apop_data *in, apop_data *d, void *p){
int		    i,j;
double	    x, 
            loglike    	= 0,
            alpha       = apop_data_get(in,0,-1),
            beta        = apop_data_get(in,1,-1);
gsl_matrix 	*data 		= d->matrix;
	for(i=0;i< data->size1; i++)
        for(j=0; j< data->size2; j++){
            x       = gsl_matrix_get(data, i, j);
		    loglike += alpha * log(x) + beta *log(1-x);
        }
	return loglike  + gsl_sf_lnbeta (alpha, beta)*data->size1*data->size2;
}

static double beta_p(const apop_data *beta, apop_data *d, void *p){
    return exp(beta_log_likelihood(beta, d, p));
}

/* The derivative of the beta distribution, for use in likelihood
  minimization. 
  The format is often the same as above: go line by line through a gsl_matrix.
  The sample is a three-dimensional parameter vector.

You can delete this function entirely if so inclined. If so, remember
to replace this function with NULL in the model definition below.
static void beta_dlog_likelihood(const gsl_vector *beta, apop_data *d, gsl_vector *gradient){
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
	for(j=0; j< beta->size; j++){
	    gsl_vector_set(gradient,j,dtotal[j]);
    }
}
 */


/* For constrained optimizations, you will need a constraint function.
 You can use \c apop_linear_constraint if your constraint can be expressed as a set of hyperplanes.

You can delete this function entirely if so inclined. If so, remember
to replace this function with NULL in the model definition below.
 */
static double beta_constraint(const apop_data *beta, void * d, apop_data *returned_beta, void *v){
    //constraint is 0 < beta_1 and  0 < beta_2
  static apop_data *constraint = NULL;
    if (!constraint)constraint= apop_data_calloc(2,2,1);
    apop_data_set(constraint, 0, 0, 1);
    apop_data_set(constraint, 1, 0, 1);
    return apop_linear_constraint(beta->vector, constraint, 1e-3, returned_beta->vector);
}


static void beta_rng(double *out, apop_data * a, apop_ep* eps, gsl_rng *r){
    *out = gsl_ran_beta(r, apop_data_get(a,0,-1), apop_data_get(a,1,-1));
}

/** The beta model.
You should describe the format of the input data here.

--The second parameter is either the number of parameters your model has, or a negative number indicating a special case; see the manual.
--If you deleted any of the functions above, replace their names with NULL here.

\ingroup models
*/
apop_model apop_beta = {"Beta distribution", -1, 0,0,
	beta_estimate, beta_p, beta_log_likelihood, NULL, beta_constraint, beta_rng};
